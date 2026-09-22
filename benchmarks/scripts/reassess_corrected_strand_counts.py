#!/usr/bin/env python3
"""Reassess the fixed 96-run candidate after the research extraction correction.

This is a new, post-hoc implementation reassessment. It cannot tune a candidate
or supply untouched evidence. Historical computations are never overwritten.
"""

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_direct_bias_stability as base
import extend_direct_bias_stability as extension
from audit_chrombpnet_comparison import canonical_id, checked_path, data_path, file_sha256
from audit_direct_bias_safety import record, write_immutable

SOURCE = ROOT / "benchmarks/results/direct_bias_stability_96_20260908_v1"
SOURCE_HASHES = {
    "configuration.freeze.json": "9e870f9a8b66bf3e37f7c330e097ccbecef8c83eb36f4b0985b97e9cb7f1d599",
    "cohorts.freeze.json": "3f9ce350545bde4931a923f74ba92039b80b50c81df7f618d033ec3ab3216ea8",
    "stability_manifest.json": "52fc25ca38d8a83b301ace2980a4947583a5605743cf729a66c6d9c8051da798",
}
SCHEMA = "fp-tools-corrected-count-reassessment-v1"


def source_documents():
    return {name: json.loads(checked_path({"path": str(SOURCE / name), "sha256": digest}, ROOT).read_text())
            for name, digest in SOURCE_HASHES.items()}


def freeze(root, output):
    if output.exists():
        raise ValueError("freeze requires a new output directory")
    docs = source_documents()
    old = docs["configuration.freeze.json"]
    matrix = pd.DataFrame(old["matrix"])
    extension.validate_matrix(matrix)
    inputs = {}

    def use(entry):
        path = checked_path(entry, root)
        inputs[str(path.resolve())] = record(path)
        return path

    for name, digest in SOURCE_HASHES.items():
        use({"path": str(SOURCE / name), "sha256": digest})
    source = json.loads(use(old["settings"]["source_configuration"]).read_text())
    settings = old["scientific_settings"]
    if settings != source["settings"]:
        raise ValueError("source scientific settings changed")
    matched_manifest = json.loads(use(settings["matched_cohorts"]).read_text())
    matched = pd.read_csv(use(matched_manifest["outputs"]["cohorts"]), sep="\t")
    output.mkdir(parents=True)
    contexts = {}
    for cell in ("HepG2", "K562"):
        ref = source["references"][cell]
        doc = ref["profile_manifest"]
        for name, digest in (("sites", "sites_sha256"), ("profiles_npz", "profiles_sha256")):
            use({"path": doc[name], "sha256": doc[digest]})
        use(ref["inputs"]["scores"])
        use(ref["inputs"]["tracks"]["dwm"])
        use({"path": doc["metadata"]["bias_model"], "sha256": settings["bias_model_sha256"]})
        sites, arrays, _, _ = base.load_reference(cell, ref, matched, root)
        support = json.loads(use(docs["cohorts.freeze.json"]["supports"][cell]).read_text())
        coverage = pd.read_csv(use(support["outputs"]["coverage"]), sep="\t")
        order = base.align_sites(coverage, sites)
        if len(order) != len(coverage):
            raise ValueError("coverage site set changed")
        strict = matched.loc[matched.cell.eq(cell) & matched.policy.eq("motif_accessibility_observed_total")]
        selected, cohorts = base.qualify_support(sites, strict, coverage.iloc[order].all_run_dwm_present.to_numpy())
        selected_sites = sites.iloc[selected].reset_index(drop=True)
        if extension.cohort_identity(selected_sites, cohorts) != docs["cohorts.freeze.json"]["memberships"][cell]["site_and_cohort_identity"]:
            raise ValueError("reassessment must preserve exact site/label order and cohort memberships")
        track_path = use(support["outputs"]["expected_tracks"])
        with np.load(track_path, allow_pickle=False) as tracks:
            if not np.array_equal(tracks["site_hash"], base.site_hashes(sites)):
                raise ValueError("expected track site order changed")
        directory = output / cell
        directory.mkdir()
        sites_path = directory / "sites.tsv.gz"
        selected_sites.to_csv(sites_path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        context_path = directory / "context.npz"
        np.savez_compressed(context_path, selected=selected, site_hash=base.site_hashes(selected_sites),
                            plus_log_bias=arrays["plus_log_bias"][selected], minus_log_bias=arrays["minus_log_bias"][selected])
        contexts[cell] = {"sites": record(sites_path), "arrays": record(context_path),
                          "tracks": record(track_path), "cohorts": {n: v.tolist() for n, v in cohorts.items()},
                          "runs": support["runs"]}
    # Freeze both cells and every run identity before any corrected extraction.
    for row in matrix.itertuples(index=False):
        for name in ("bam", "bai", "dwm_expected"):
            use({"path": getattr(row, name), "sha256": getattr(row, name + "_sha256")})
        print(f"Verified freeze inputs: {row.sample} {row.depth} {row.seed}", flush=True)
    config = {"schema": SCHEMA, "evidence": "post_hoc_fixed_candidate_corrected_count_reassessment",
              "matrix": old["matrix"], "scientific_settings": settings,
              "contexts": contexts, "inputs": list(inputs.values()),
              "code": [record(Path(entry["path"])) for entry in old["code"]] + [record(Path(__file__))],
              "runtime": extension.runtime(), "source_result_sha256": SOURCE_HASHES["stability_manifest.json"],
              "candidate_retuned": False, "new_holdout_labels_read": False,
              "final_untouched_control_read": False, "promotion_eligible": False,
              "archived_observed_profiles_reused": False}
    config["configuration_id"] = canonical_id(config)
    write_immutable(output / "configuration.freeze.json", config)
    print(f"Frozen corrected-count reassessment: {config['configuration_id']}", flush=True)


def validate_configuration(config):
    if (config.get("schema") != SCHEMA or config.get("evidence") != "post_hoc_fixed_candidate_corrected_count_reassessment"
            or config.get("source_result_sha256") != SOURCE_HASHES["stability_manifest.json"]
            or any(config.get(k) is not False for k in ("candidate_retuned", "new_holdout_labels_read",
                "final_untouched_control_read", "promotion_eligible", "archived_observed_profiles_reused"))):
        raise ValueError("invalid corrected-count scope")
    docs = source_documents()
    if config["scientific_settings"] != docs["configuration.freeze.json"]["scientific_settings"] or config["matrix"] != docs["configuration.freeze.json"]["matrix"]:
        raise ValueError("candidate settings or matrix changed")


def corrected_arrays(sites, bam, arrays):
    plus, minus, valid = base.extract_strand_cut_profiles(
        sites, bam, flank=100, read_shift=(4, -5), minimum_mapq=30, keep_duplicates=False)
    if not valid.all():
        raise ValueError("corrected BAM lacks a frozen window")
    plus, minus, _ = base.orient_strand_log_bias(plus, minus, sites.TFBS_strand.tolist())
    expected = base.expected_for_observed(plus, minus, arrays["plus_log_bias"], arrays["minus_log_bias"])
    return plus + minus, expected


def execute(root, output, verify_only=False):
    config = extension.read_document(output / "configuration.freeze.json", "configuration_id")
    validate_configuration(config)
    if config["runtime"] != extension.runtime():
        raise ValueError("frozen reassessment runtime changed")
    for entry in config["code"]:
        checked_path(entry, root)
    if verify_only:
        for entry in config["inputs"]:
            checked_path(entry, root)
    # Small metadata and frozen scientific input contexts are always rechecked.
    completed = []
    for cell, context in config["contexts"].items():
        sites = pd.read_csv(checked_path(context["sites"], root), sep="\t")
        if set(sites.cell) != {cell} or not set(sites.TFBS_chr).issubset(base.TEST_CHROMOSOMES):
            raise ValueError("unexpected cell or chromosome in frozen reassessment")
        with np.load(checked_path(context["arrays"], root), allow_pickle=False) as data:
            arrays = {name: data[name] for name in data.files}
        if not np.array_equal(arrays["site_hash"], base.site_hashes(sites)):
            raise ValueError("context/site order mismatch")
        cohorts = {n: np.asarray(v, dtype=int) for n, v in context["cohorts"].items()}
        prior = source_documents()["cohorts.freeze.json"]["memberships"][cell]
        if extension.cohort_identity(sites, cohorts) != prior["site_and_cohort_identity"]:
            raise ValueError("cohort membership changed")
        with np.load(checked_path(context["tracks"], root), allow_pickle=False) as data:
            tracks = {n: data[n][arrays["selected"]] for n in data.files if n != "site_hash"}
        frame = pd.DataFrame(config["matrix"]).query("cell == @cell").sort_values(base.RUN_KEYS).reset_index(drop=True)
        if frame[base.RUN_KEYS].to_dict("records") != context["runs"]:
            raise ValueError("track/run order changed")
        for index, row in enumerate(frame.itertuples(index=False)):
            for name in ("bam", "bai", "dwm_expected"):
                checked_path({"path": getattr(row, name), "sha256": getattr(row, name + "_sha256")}, root)
            directory = output / row.sample / row.depth / f"seed_{row.seed}"
            run = {key: getattr(row, key) for key in base.RUN_KEYS}
            extraction = directory / "extraction.json"
            if extraction.exists():
                prior_stage = base.verify_stage(extraction, config["configuration_id"], root, expected_run=run)
                if prior_stage.get("archived_counts_reused") is not False:
                    raise ValueError("resumed extraction does not certify fresh corrected counts")
                with np.load(checked_path(prior_stage["outputs"]["profiles"], root), allow_pickle=False) as data:
                    if not np.array_equal(data["site_hash"], base.site_hashes(sites)):
                        raise ValueError("resumed extraction site order changed")
                    observed, expected = data["observed"], data["expected"]
            elif verify_only:
                raise ValueError("missing corrected extraction")
            else:
                directory.mkdir(parents=True, exist_ok=True)
                profile_path = directory / "profiles.npz"
                if profile_path.exists():
                    raise ValueError("incomplete extraction retained; do not overwrite")
                observed, expected = corrected_arrays(sites, data_path(row.bam, root), arrays)
                np.savez_compressed(profile_path, observed=observed, expected=expected, site_hash=base.site_hashes(sites))
                stage = {"configuration_id": config["configuration_id"], "run": run,
                         "outputs": {"profiles": record(profile_path)}, "archived_counts_reused": False}
                stage["stage_id"] = canonical_id(stage)
                write_immutable(extraction, stage)
            evaluation = directory / "evaluation.json"
            if evaluation.exists():
                base.verify_stage(evaluation, config["configuration_id"], root, expected_run=run)
            elif verify_only:
                raise ValueError("missing corrected evaluation")
            else:
                if list(directory.glob("*.tsv.gz")):
                    raise ValueError("incomplete evaluation retained; do not overwrite")
                methods = base.score_panel(observed, expected, tracks["fixed_dwm"], tracks[f"run_{index}"],
                                          config["scientific_settings"]["dispersion"][cell])
                base.evaluate_run(row, sites, cohorts, methods, observed, config["scientific_settings"],
                                  directory, config["configuration_id"], root)
            completed.append({"extraction": record(extraction), "evaluation": record(evaluation)})
    if verify_only:
        result = base.verify_stage(output / "reassessment_manifest.json", config["configuration_id"], root)
        if (result["stages"] != completed or result.get("schema") != SCHEMA
                or any(result.get(k) is not False for k in ("promotion_eligible", "candidate_retuned",
                    "new_holdout_labels_read", "final_untouched_control_read"))):
            raise ValueError("result stage matrix differs")
        print("Verified all 96 corrected-count reassessments", flush=True)
        return
    outputs = {}
    for name in ("metrics", "bootstrap", "profiles", "shapes", "coverage"):
        tables = []
        for entry in completed:
            stage = json.loads(Path(entry["evaluation"]["path"]).read_text())
            table = pd.read_csv(checked_path(stage["outputs"][name], root), sep="\t")
            if not table.empty:
                tables.append(table)
        path = output / f"stability_{name}.tsv.gz"
        if path.exists():
            raise ValueError("incomplete summary retained; do not overwrite")
        pd.concat(tables, ignore_index=True).to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(path)
    result = {"configuration_id": config["configuration_id"], "schema": SCHEMA,
              "stages": completed, "outputs": outputs, "promotion_eligible": False,
              "candidate_retuned": False, "new_holdout_labels_read": False,
              "final_untouched_control_read": False}
    result["stage_id"] = canonical_id(result)
    write_immutable(output / "reassessment_manifest.json", result)
    print("Completed all 96 fixed-candidate corrected-count reassessments", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--freeze-only", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    output = args.outdir.resolve()
    if not output.is_relative_to(ROOT / "benchmarks/results") or output == SOURCE:
        parser.error("use a new research result directory")
    if any(os.environ.get(k) != v for k, v in extension.THREADS.items()):
        parser.error("set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1")
    if args.freeze_only:
        freeze(args.data_root.resolve(), output)
    else:
        with extension.exclusive_run(output):
            execute(args.data_root.resolve(), output,
                    args.verify_only or (output / "reassessment_manifest.json").exists())
