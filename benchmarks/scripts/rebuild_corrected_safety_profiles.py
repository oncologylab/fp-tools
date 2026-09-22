#!/usr/bin/env python3
"""Reextract frozen development/control-2 sites without fitting or scoring."""

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
from audit_chrombpnet_comparison import canonical_id, checked_path
from audit_direct_bias_safety import record, write_immutable

PROTOCOL = ROOT / "benchmarks/manifests/corrected_count_safety_protocol_v1.json"
PROTOCOL_SHA = "be8ccb69af72b664264532db8f1c8a06011e27581a7fd2404efb818663c8912b"


def validate_sites(sites, arrays, cell, role):
    if cell not in ("HepG2", "K562") or role not in ("validation", "naked"):
        raise ValueError("only development cells and control replicate 2 are allowed")
    if sites.empty or set(sites.cell) != {cell} or sites.duplicated(base.SITE_KEYS).any():
        raise ValueError("invalid or duplicate site support")
    if role == "validation":
        if not set(sites.TFBS_chr).issubset({f"chr{i}" for i in range(1, 19)}) or not sites.chip_label.isin([0, 1]).all():
            raise ValueError("calibration inputs contain non-development labels")
    elif any("chip" in name.lower() or "label" in name.lower() for name in sites.columns):
        raise ValueError("control sites must not contain labels")
    if not np.array_equal(arrays["site_hash"], base.site_hashes(sites)):
        raise ValueError("frozen site order changed")
    for name in ("plus_log_bias", "minus_log_bias"):
        if arrays[name].shape != (len(sites), 201):
            raise ValueError("frozen sequence arrays must match sites and width")
    if arrays["valid"].shape != (len(sites),):
        raise ValueError("invalid archived validity mask")


def rebuild(sites, arrays, bam):
    plus, minus, bam_valid = base.extract_strand_cut_profiles(
        sites, bam, flank=100, read_shift=(4, -5), minimum_mapq=30, keep_duplicates=False)
    plus, minus, _ = base.orient_strand_log_bias(plus, minus, sites.TFBS_strand.tolist())
    first, first_valid = base.expected_from_cached_log_bias(plus, arrays["plus_log_bias"])
    second, second_valid = base.expected_from_cached_log_bias(minus, arrays["minus_log_bias"])
    expected = first + second
    # Original sequence validity is retained separately; no site is dropped.
    sequence_valid = first_valid & second_valid
    return dict(observed=plus + minus, expected=expected, plus_observed=plus,
                minus_observed=minus, valid=bam_valid & sequence_valid,
                archived_valid=arrays["valid"].astype(bool), site_hash=base.site_hashes(sites))


def freeze(data_root, output):
    protocol = json.loads(checked_path({"path": str(PROTOCOL), "sha256": PROTOCOL_SHA}, ROOT).read_text())
    source = json.loads(checked_path(protocol["source_configuration"], ROOT).read_text())
    for key in ("source_thresholds", "source_result"):
        checked_path(protocol[key], ROOT)
    if output.exists():
        raise ValueError("freeze requires a fresh output directory")
    items = []
    for cell in protocol["cells"]:
        for role in ("validation", "naked"):
            manifest = checked_path(source["cells"][cell][role], data_root)
            doc = json.loads(manifest.read_text())
            meta = doc["metadata"]
            if (doc["schema"] != "fp-tools-strand-functional-profiles-v1" or meta.get("labels_used") is not False
                    or any(meta[k] != protocol[k] for k in ("bias_model_sha256", "flank", "read_shift", "minimum_mapq", "keep_duplicates"))):
                raise ValueError("source profile construction differs from protocol")
            if role == "naked" and "NakedDNA_rep2" not in meta["bam"]:
                raise ValueError("only archived naked-DNA replicate 2 is allowed")
            paths = {name: checked_path({"path": doc[key], "sha256": doc[digest]}, data_root)
                     for name, key, digest in (("sites", "sites", "sites_sha256"), ("arrays", "profiles_npz", "profiles_sha256"))}
            sites = pd.read_csv(paths["sites"], sep="\t")
            with np.load(paths["arrays"], allow_pickle=False) as arrays:
                validate_sites(sites, arrays, cell, role)
            bam = checked_path({"path": meta["bam"], "sha256": meta["bam_sha256"]}, data_root)
            bai = Path(str(bam) + ".bai")
            if not bai.is_file():
                raise ValueError("matching BAM index is required")
            items.append(dict(cell=cell, role=role, manifest=record(manifest), bam=record(bam), bai=record(bai),
                              **{k: record(p) for k, p in paths.items()}))
    code_paths = [Path(__file__), Path(base.__file__), Path(extension.__file__),
                  ROOT / "benchmarks/scripts/build_strand_functional_profiles.py",
                  ROOT / "benchmarks/scripts/audit_chrombpnet_comparison.py",
                  ROOT / "benchmarks/scripts/audit_direct_bias_safety.py",
                  ROOT / "src/fp_tools/tools/parametric_bias.py",
                  ROOT / "benchmarks/scripts/evaluate_frozen_functional_depth_matrix.py"]
    code_paths = sorted(set(code_paths + [Path(e["path"]) for e in source["code"]]))
    config = dict(schema="fp-tools-corrected-safety-profiles-v1", protocol=record(PROTOCOL),
                  source_configuration=protocol["source_configuration"], items=items,
                  code=[record(p) for p in code_paths], runtime=extension.runtime(),
                  profiles_scored=False, candidate_fitted=False, final_control_read=False)
    config["configuration_id"] = canonical_id(config)
    output.mkdir(parents=True)
    write_immutable(output / "configuration.freeze.json", config)
    print("Frozen four corrected safety-profile inputs", flush=True)


def execute(data_root, output, verify_only=False):
    config = extension.read_document(output / "configuration.freeze.json", "configuration_id")
    if (config["schema"] != "fp-tools-corrected-safety-profiles-v1" or config["runtime"] != extension.runtime()
            or config["protocol"]["sha256"] != PROTOCOL_SHA
            or any(config[k] is not False for k in ("profiles_scored", "candidate_fitted", "final_control_read"))
            or [(i["cell"], i["role"]) for i in config["items"]] != [(c, r) for c in ("HepG2", "K562") for r in ("validation", "naked")]):
        raise ValueError("frozen execution scope or runtime changed")
    checked_path(config["protocol"], ROOT)
    checked_path(config["source_configuration"], ROOT)
    for entry in config["code"]:
        checked_path(entry, ROOT)
    stages = []
    for item in config["items"]:
        paths = {k: checked_path(item[k], data_root) for k in ("manifest", "bam", "bai", "sites", "arrays")}
        sites = pd.read_csv(paths["sites"], sep="\t")
        with np.load(paths["arrays"], allow_pickle=False) as loaded:
            validate_sites(sites, loaded, item["cell"], item["role"])
            arrays = {k: loaded[k] for k in ("plus_log_bias", "minus_log_bias", "valid")}
        prefix = output / (item["cell"] + "." + item["role"])
        receipt, profiles = prefix.with_suffix(prefix.suffix + ".json"), prefix.with_suffix(prefix.suffix + ".npz")
        if receipt.exists():
            done = base.verify_stage(receipt, config["configuration_id"], ROOT)
            if done["cell"] != item["cell"] or done["role"] != item["role"]:
                raise ValueError("resumed stage scope mismatch")
        elif verify_only:
            raise ValueError("missing corrected profile stage")
        else:
            if profiles.exists():
                raise ValueError("partial output preserved; do not overwrite")
            np.savez_compressed(profiles, **rebuild(sites, arrays, paths["bam"]))
            done = dict(configuration_id=config["configuration_id"], cell=item["cell"], role=item["role"],
                        outputs={"profiles": record(profiles), "sites": item["sites"]})
            done["stage_id"] = canonical_id(done)
            write_immutable(receipt, done)
        stages.append(record(receipt))
        print("Verified corrected profiles:", item["cell"], item["role"], flush=True)
    result = dict(configuration_id=config["configuration_id"], stages=stages, outputs={},
                  profiles_scored=False, candidate_fitted=False, final_control_read=False)
    result["stage_id"] = canonical_id(result)
    path = output / "profiles_manifest.json"
    if verify_only:
        if base.verify_stage(path, config["configuration_id"], ROOT) != result:
            raise ValueError("completed profile matrix differs")
    else:
        write_immutable(path, result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--freeze-only", action="store_true")
    mode.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    output = args.outdir.resolve()
    if any(os.environ.get(k) != v for k, v in extension.THREADS.items()):
        parser.error("set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1")
    if not output.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("output must stay in research benchmark results")
    if args.freeze_only:
        freeze(args.data_root.resolve(), output)
    else:
        with extension.exclusive_run(output):
            execute(args.data_root.resolve(), output, args.verify_only)
