#!/usr/bin/env python3
"""Test the recorded direct bias operator on existing naked-DNA replicate 2.

Freeze code and input identities, calibrate only on validation negatives, then
apply immutable thresholds to an already examined candidate-screen library.
No candidate fitting, new holdout labels or final untouched control is allowed.
"""

from __future__ import annotations

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from audit_chrombpnet_comparison import (
    PINNED_COMMIT, SITE_KEYS, canonical_id, check_identity, checked_path, data_path, file_sha256,
)
from build_strand_functional_profiles import site_hashes
from evaluate_chrombpnet_reference import extract_prediction_profiles, scale_to_observed_total
from evaluate_factorization_residual_safety import conformal_upper_threshold, wilson_interval
from evaluate_parametric_factorization import residual_score
from fp_tools.tools.frozen_bias_evaluation import motif_residual_effect

OPERATOR = "total_matched_direct_expected_signed_deviance_conventional_geometry"
SCHEMA = "fp-tools-direct-bias-candidate-safety-v1"
METHODS = ("parametric_direct", "DWM_expected")


def direct_scores(observed: np.ndarray, expected: np.ndarray, dispersion: float):
    if (observed.shape != expected.shape or observed.ndim != 2 or observed.shape[1] != 201
            or not np.isfinite(observed).all() or not np.isfinite(expected).all()
            or np.any(observed < 0) or np.any(expected < 0)
            or not np.isfinite(dispersion) or dispersion < 0):
        raise ValueError("direct scores require finite nonnegative 201-position profiles")
    scaled = scale_to_observed_total(observed, expected)
    score, residual = residual_score(observed, scaled, np.arange(-100, 101), "deviance", dispersion)
    return score, scaled, residual


def calibrate(sites, scores, valid, totals, *, target=0.05, minimum_negatives=200):
    rows = []
    for tf, task in sites.groupby("tf", sort=True):
        indexes = task.index.to_numpy(dtype=int)
        member = (task.TFBS_chr.isin(["chr16", "chr17", "chr18"]).to_numpy()
                  & task.chip_label.eq(0).to_numpy() & valid[indexes]
                  & np.isfinite(scores[indexes]) & (totals[indexes] > 0))
        values = scores[indexes][member]
        eligible = len(values) >= minimum_negatives
        rows.append({"tf": tf, "validation_negatives": len(values),
                     "threshold": conformal_upper_threshold(values, target) if eligible else np.nan,
                     "threshold_eligible": eligible})
    return pd.DataFrame(rows)


def safety_row(candidate, baseline, totals, *, candidate_threshold, baseline_threshold, threshold_eligible):
    candidate, baseline, totals = [np.asarray(v, dtype=float) for v in (candidate, baseline, totals)]
    if candidate.ndim != 1 or candidate.shape != baseline.shape or candidate.shape != totals.shape:
        raise ValueError("safety scores and totals must have identical one-dimensional support")
    finite = np.isfinite(candidate) & np.isfinite(baseline) & np.isfinite(totals)
    informative = finite & (totals > 0)
    candidate_calls = informative & (candidate > candidate_threshold)
    baseline_calls = informative & (baseline > baseline_threshold)
    row = {"false_positive_calls": int(candidate_calls.sum()), "dwm_calls": int(baseline_calls.sum()),
           "threshold_eligible": bool(threshold_eligible), "final_untouched_safety": False}
    passes = bool(threshold_eligible and np.isfinite(candidate_threshold) and np.isfinite(baseline_threshold))
    for prefix, mask in (("", finite), ("informative_", informative)):
        n = int(mask.sum())
        rate = candidate_calls.sum() / n if n else np.nan
        dwm_rate = baseline_calls.sum() / n if n else np.nan
        _, upper = wilson_interval(int(candidate_calls.sum()), n)
        row.update({prefix + "support": n, prefix + "false_positive_rate": rate,
                    prefix + "false_positive_rate_upper_95": upper,
                    prefix + "false_positive_rate_increase_over_dwm": rate - dwm_rate})
        passes &= bool(n > 0 and rate <= 0.05 and upper <= 0.05 and rate - dwm_rate <= 0.01)
    row["passes_candidate_screen"] = passes
    return row


def write_immutable(path: Path, document: dict) -> None:
    rendered = json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if path.exists() and path.read_text() != rendered:
        raise ValueError(f"immutable safety artifact differs: {path}")
    if not path.exists():
        path.write_text(rendered)


def record(path: Path) -> dict:
    return {"path": str(path.resolve()), "sha256": file_sha256(path)}


def one_record(records, predicate) -> dict:
    matched = [item for item in records if predicate(item)]
    if len(matched) != 1:
        raise ValueError("expected exactly one registered input artifact")
    return matched[0]


def load_artifact(item: dict, data_root: Path, cell: str, *, control: bool):
    document = json.loads(checked_path(item, data_root).read_text())
    metadata = document["metadata"]
    if metadata.get("labels_used") is not False or metadata.get("flank") != 100:
        raise ValueError("profile construction must be label-free with flank 100")
    if control and "NakedDNA_rep2" not in str(metadata.get("bam", metadata.get("observed", ""))):
        raise ValueError("only the recorded naked-DNA replicate 2 is allowed")
    sites = pd.read_csv(checked_path({"path": document["sites"], "sha256": document["sites_sha256"]}, data_root), sep="\t")
    if set(sites.cell) != {cell} or sites.duplicated(SITE_KEYS).any():
        raise ValueError("profile sites have unexpected cells or duplicate loci")
    if control:
        if any("chip" in column.lower() or "label" in column.lower() for column in sites.columns):
            raise ValueError("naked-DNA site table contains labels")
    elif (not sites.chip_label.isin([0, 1]).all()
          or not set(sites.TFBS_chr).issubset({f"chr{i}" for i in range(1, 19)})):
        raise ValueError("calibration source contains non-development chromosomes or labels")
    path = checked_path({"path": document["profiles_npz"], "sha256": document["profiles_sha256"]}, data_root)
    with np.load(path, allow_pickle=False) as loaded:
        arrays = {name: np.asarray(loaded[name]) for name in loaded.files}
    if not np.array_equal(arrays["site_hash"], site_hashes(sites)):
        raise ValueError("profile site order mismatch")
    if document["schema"] == "fp-tools-strand-functional-profiles-v1":
        if metadata.get("read_shift") != [4, -5]:
            raise ValueError("direct bias cut shift differs from the recorded comparison")
        observed = arrays["plus_observed"] + arrays["minus_observed"]
        expected = arrays["plus_expected"] + arrays["minus_expected"]
    elif document["schema"] == "fp-tools-combined-functional-profiles-v1" and control:
        observed, expected = arrays["observed"], arrays["expected"]
    else:
        raise ValueError("unsupported safety profile schema")
    return sites, observed, expected, arrays["valid"].astype(bool), metadata


def resolve_plan(plan: dict, root: Path):
    if (plan.get("schema") != "fp-tools-direct-bias-safety-plan-v1"
            or plan.get("operator") != OPERATOR or plan.get("naked_replicate") != "NakedDNA_rep2"
            or plan.get("new_holdout_labels_allowed") is not False
            or plan.get("final_untouched_control_allowed") is not False):
        raise ValueError("unsupported direct-bias candidate-screen plan")
    if (not plan.get("comparisons") or any(item.get("cell") not in {"K562", "HepG2"} for item in plan["comparisons"])
            or not 0 < plan["target_validation_negative_call_rate"] < 1
            or plan["minimum_validation_negatives"] < 200
            or plan["residual_bootstraps"] < 1000):
        raise ValueError("invalid development-only safety plan")
    validation = json.loads(checked_path(plan["validation_input_registry"], root).read_text())
    naked = json.loads(checked_path(plan["naked_input_registry"], root).read_text())
    if naked.get("replicate") != "NakedDNA_rep2" or naked.get("naked_dna_labels_used") is not False:
        raise ValueError("naked registry is not the candidate-screen library")
    cells = {}
    for item in plan["comparisons"]:
        cell = item["cell"]
        if cell not in {"K562", "HepG2"} or cell in cells:
            raise ValueError("safety audit accepts distinct development cells only")
        path = checked_path(item, root)
        comparison = json.loads(path.read_text())
        freeze = json.loads(path.with_name("chrombpnet_comparison_inputs.freeze.json").read_text())
        check_identity(comparison, "result_id")
        check_identity(freeze, "input_id")
        if (comparison["input_id"] != freeze["input_id"] or freeze["cells"] != [cell]
                or freeze.get("pinned_chrombpnet_commit") != PINNED_COMMIT):
            raise ValueError("comparison provenance differs from safety plan")
        artifact = one_record(validation["validation_inputs"], lambda r: Path(r["path"]).name == f"{cell}.json")
        test_artifact = one_record(freeze["inputs"], lambda r: r["purpose"] == f"{cell}-test-manifest")
        test_metadata = json.loads(checked_path(test_artifact, root).read_text())["metadata"]
        candidate = one_record(naked["inputs"], lambda r: r["purpose"] == f"{cell}-naked-profile-manifest")
        dwm = one_record(naked["inputs"], lambda r: r["purpose"] == f"{cell}-naked-DWM" and r["path"].endswith(".json"))
        track = one_record(freeze["inputs"], lambda r: r["purpose"] == f"{cell}-dwm")
        checked_path(track, root)
        for source in (artifact, candidate):
            metadata = json.loads(checked_path(source, root).read_text())["metadata"]
            for key in ("bias_model_sha256", "read_shift", "genome_sha256", "flank", "minimum_mapq", "keep_duplicates"):
                if metadata[key] != test_metadata[key]:
                    raise ValueError(f"safety artifact differs from comparison: {cell} {key}")
            if source is artifact and metadata["bam_sha256"] != test_metadata["bam_sha256"]:
                raise ValueError("validation and comparison must use the same observed library")
        checked_path(dwm, root)
        cells[cell] = {"validation": artifact, "naked": candidate, "naked_dwm": dwm,
                       "dwm_track": track, "dispersion": float(freeze["dispersion"]),
                       "comparison_result_id": comparison["result_id"],
                       "bias_model_sha256": test_metadata["bias_model_sha256"]}
    return cells


def read_thresholds(path: Path, configuration_id: str, root: Path):
    frozen = json.loads(path.read_text())
    check_identity(frozen, "threshold_id")
    if (frozen["configuration_id"] != configuration_id
            or frozen["calibration_chromosomes"] != ["chr16", "chr17", "chr18"]
            or frozen["naked_profiles_scored"] is not False):
        raise ValueError("threshold configuration mismatch")
    table = pd.read_csv(checked_path(frozen["thresholds"], root), sep="\t")
    if table.duplicated(["cell", "tf", "method"]).any():
        raise ValueError("duplicate safety threshold")
    return frozen, table


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    plan = json.loads(args.plan.read_text())
    cells = resolve_plan(plan, args.data_root)
    sources = [Path(__file__), SCRIPT_DIR / "audit_chrombpnet_comparison.py",
               SCRIPT_DIR / "evaluate_chrombpnet_reference.py", SCRIPT_DIR / "evaluate_parametric_factorization.py",
               SCRIPT_DIR / "evaluate_factorization_residual_safety.py", SCRIPT_DIR / "build_strand_functional_profiles.py",
               REPOSITORY / "src/fp_tools/tools/parametric_bias.py", REPOSITORY / "src/fp_tools/tools/frozen_bias_evaluation.py",
               REPOSITORY / "src/fp_tools/tools/parametric_factorization.py",
               REPOSITORY / "src/fp_tools/tools/functional_footprints.py", REPOSITORY / "src/fp_tools/utils/bigwig.py"]
    configuration = {"schema": "fp-tools-direct-bias-safety-configuration-v1", "plan": record(args.plan),
                     "operator": OPERATOR, "cells": cells, "code": [record(p) for p in sources],
                     "runtime": {"python": sys.version, **{name: version(name) for name in ("numpy", "pandas", "scipy", "pyBigWig")}}}
    configuration["configuration_id"] = canonical_id(configuration)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_immutable(args.outdir / "configuration.freeze.json", configuration)
    result_path = args.outdir / "safety_manifest.json"
    threshold_manifest = args.outdir / "thresholds.freeze.json"
    if threshold_manifest.exists():
        frozen, thresholds = read_thresholds(threshold_manifest, configuration["configuration_id"], args.data_root)
    if result_path.exists():
        if not threshold_manifest.exists():
            raise ValueError("completed safety result lacks frozen thresholds")
        result = json.loads(result_path.read_text())
        check_identity(result, "safety_id")
        if (result["configuration_id"] != configuration["configuration_id"]
                or result["threshold_id"] != frozen["threshold_id"]):
            raise ValueError("completed safety configuration mismatch")
        for output in result["outputs"].values():
            checked_path(output, args.data_root)
        print("Verified completed direct-bias safety audit; no stage rerun.")
        return 0
    if not threshold_manifest.exists():
        frames = []
        for cell, inputs in sorted(cells.items()):
            sites, observed, expected, valid, _ = load_artifact(inputs["validation"], args.data_root, cell, control=False)
            dwm, dwm_valid = extract_prediction_profiles(sites, checked_path(inputs["dwm_track"], args.data_root), 100, require_dense=False)
            valid &= dwm_valid
            # Invalid genomic windows remain excluded from both calibration arms.
            dwm = np.nan_to_num(dwm)
            for name, values in ((METHODS[0], expected), (METHODS[1], dwm)):
                scores = direct_scores(observed, values, inputs["dispersion"])[0]
                frame = calibrate(sites, scores, valid, observed.sum(axis=1),
                                  target=plan["target_validation_negative_call_rate"],
                                  minimum_negatives=plan["minimum_validation_negatives"])
                frame["cell"], frame["method"] = cell, name
                frames.append(frame)
        thresholds = pd.concat(frames, ignore_index=True)
        threshold_path = args.outdir / "validation_thresholds.tsv"
        thresholds.to_csv(threshold_path, sep="\t", index=False)
        frozen = {"configuration_id": configuration["configuration_id"], "thresholds": record(threshold_path),
                  "calibration_chromosomes": ["chr16", "chr17", "chr18"], "naked_profiles_scored": False}
        frozen["threshold_id"] = canonical_id(frozen)
        write_immutable(threshold_manifest, frozen)
    print("Validation thresholds frozen; evaluating existing naked-DNA replicate 2.", flush=True)
    lookup = thresholds.set_index(["cell", "tf", "method"])
    rate_rows, residual_rows, curve_frames, site_frames = [], [], [], []
    for cell, inputs in sorted(cells.items()):
        sites, observed, expected, valid, _ = load_artifact(inputs["naked"], args.data_root, cell, control=True)
        dwm_sites, dwm_observed, dwm_expected, dwm_valid, _ = load_artifact(inputs["naked_dwm"], args.data_root, cell, control=True)
        order = pd.Index(site_hashes(dwm_sites)).get_indexer(site_hashes(sites))
        if np.any(order < 0):
            raise ValueError("DWM control does not contain the same motif loci")
        dwm_expected, dwm_observed = dwm_expected[order], dwm_observed[order]
        valid &= dwm_valid[order]
        methods = {name: direct_scores(observed, values, inputs["dispersion"])
                   for name, values in ((METHODS[0], expected), (METHODS[1], dwm_expected))}
        totals = observed.sum(axis=1)
        for tf, task in sites.groupby("tf", sort=True):
            indexes = task.index.to_numpy(dtype=int)
            indexes = indexes[valid[indexes]]
            cutoffs = [lookup.loc[(cell, tf, name)] for name in METHODS]
            row = safety_row(methods[METHODS[0]][0][indexes], methods[METHODS[1]][0][indexes], totals[indexes],
                             candidate_threshold=cutoffs[0].threshold, baseline_threshold=cutoffs[1].threshold,
                             threshold_eligible=all(bool(r.threshold_eligible) for r in cutoffs))
            row.update(cell=cell, tf=tf, dispersion=inputs["dispersion"],
                       validation_negatives=int(cutoffs[0].validation_negatives),
                       candidate_threshold=cutoffs[0].threshold, baseline_threshold=cutoffs[1].threshold,
                       zero_dwm_expected_with_cuts=int(((dwm_expected[indexes].sum(axis=1) == 0) & (totals[indexes] > 0)).sum()),
                       dwm_observed_count_disagreements=int(np.any(observed[indexes] != dwm_observed[indexes], axis=1).sum()))
            rate_rows.append(row)
            frame = sites.loc[indexes, SITE_KEYS].copy()
            frame["total_cuts"] = totals[indexes]
            for name, (scores, scaled, residuals) in methods.items():
                frame[name + "_score"] = scores[indexes]
                informative = indexes[totals[indexes] > 0]
                if len(informative) >= 2:
                    residual_rows.append({"cell": cell, "tf": tf, "method": name,
                                          **motif_residual_effect(observed[informative], scaled[informative], np.arange(-100, 101),
                                              bootstraps=plan["residual_bootstraps"], seed=plan["seed"],
                                              threshold=plan["motif_residual_absolute_center_flank_limit"])})
                    curve_frames.append(pd.DataFrame({"cell": cell, "tf": tf, "method": name,
                        "position": np.arange(-100, 101), "mean_observed": observed[informative].mean(axis=0),
                        "mean_expected": scaled[informative].mean(axis=0), "mean_deviance": residuals[informative].mean(axis=0)}))
            site_frames.append(frame)
    tables = {"false_positive_rates": pd.DataFrame(rate_rows), "motif_residuals": pd.DataFrame(residual_rows),
              "naked_profiles": pd.concat(curve_frames, ignore_index=True), "site_scores": pd.concat(site_frames, ignore_index=True)}
    outputs = {}
    for name, table in tables.items():
        path = args.outdir / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        outputs[name] = record(path)
    result = {"schema": SCHEMA, "configuration_id": configuration["configuration_id"],
              "configuration": record(args.outdir / "configuration.freeze.json"),
              "thresholds": record(threshold_manifest),
              "threshold_id": frozen["threshold_id"], "operator": OPERATOR,
              "evidence": "retrospective_development_candidate_screen", "replicate": "NakedDNA_rep2",
              "new_holdout_labels_read": False, "final_untouched_safety": False,
              "promotion_eligible": False, "outputs": outputs}
    result["safety_id"] = canonical_id(result)
    write_immutable(result_path, result)
    print(tables["false_positive_rates"][["cell", "tf", "support", "informative_support", "false_positive_calls", "passes_candidate_screen"]].to_string(index=False))
    print(result_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
