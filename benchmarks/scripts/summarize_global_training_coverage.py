#!/usr/bin/env python3
"""Audit global-operator training coverage without labels or score evaluation."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_global_parametric_training_profiles as builder
import global_parametric_scoring as scoring
from audit_chrombpnet_comparison import canonical_id, checked_path
from audit_direct_bias_safety import record, write_immutable

ROOT = builder.ROOT
SOURCE = ROOT / "benchmarks/results/global_parametric_training_profiles_20260922_v1"
SOURCE_SHA = "d903eec175130993afe35e6f3c64701e25399d1e29f6d2b98f8750825849b328"


def neighborhood_sum(values):
    values = np.asarray(values)
    if values.ndim != 2 or values.shape[1] != 401:
        raise ValueError("401-position arrays required")
    prefix = np.pad(np.cumsum(values, axis=1), ((0, 0), (1, 0)))
    return prefix[:, 201:] - prefix[:, :-201]


def annular_sum(values):
    values = np.asarray(values)
    if values.ndim != 2 or values.shape[1] != 401:
        raise ValueError("401-position arrays required")
    prefix = np.pad(np.cumsum(values, axis=1), ((0, 0), (1, 0)))
    return (prefix[:, 51:252]-prefix[:, :201]) + (prefix[:, 201:402]-prefix[:, 150:351])


def summarize(sites, arrays):
    if any("label" in c.lower() or "chip" in c.lower() for c in sites.columns):
        raise ValueError("coverage audit must not read labels")
    missing_sequence = ((neighborhood_sum(np.isfinite(arrays["plus_log_bias"])) != 201)
                        | (neighborhood_sum(np.isfinite(arrays["minus_log_bias"])) != 201))
    empty_plus = annular_sum(arrays["plus_observed"]) == 0
    empty_minus = annular_sum(arrays["minus_observed"]) == 0
    valid = arrays["bam_valid"][:, None] & ~missing_sequence & ~empty_plus & ~empty_minus
    np.testing.assert_array_equal(valid, arrays["position_valid"])
    np.testing.assert_array_equal(valid.all(axis=1), arrays["valid"])
    geometry = scoring.conventional_weights() != 0
    rows, positions = [], []
    for (cell, tf), task in sites.groupby(["cell", "tf"]):
        ix = task.index.to_numpy()
        rows.append(dict(cell=cell, tf=tf, sites=len(ix), bam_valid=int(arrays["bam_valid"][ix].sum()),
            positive_observed_sites=int((arrays["observed"][ix].sum(axis=1) > 0).sum()),
            central_expectation_sites=int(valid[ix, 100].sum()), full_profile_sites=int(valid[ix].all(axis=1).sum()),
            conventional_geometry_sites=int(valid[ix, 60:141][:, geometry].all(axis=1).sum()),
            central_missing_sequence=int(missing_sequence[ix, 100].sum()),
            central_empty_plus_background=int(empty_plus[ix, 100].sum()),
            central_empty_minus_background=int(empty_minus[ix, 100].sum()),
            central_both_backgrounds_empty=int((empty_plus[ix, 100] & empty_minus[ix, 100]).sum())))
        for offset in range(201):
            positions.append(dict(cell=cell, tf=tf, position=offset-100, sites=len(ix),
                supported=int(valid[ix, offset].sum()), missing_sequence=int(missing_sequence[ix, offset].sum()),
                empty_plus_background=int(empty_plus[ix, offset].sum()), empty_minus_background=int(empty_minus[ix, offset].sum())))
    return pd.DataFrame(rows), pd.DataFrame(positions)


def freeze(root, out):
    checked_path(dict(path=str(SOURCE / "configuration.freeze.json"), sha256=SOURCE_SHA), ROOT)
    config = builder.parent.extension.read_document(SOURCE / "configuration.freeze.json", "configuration_id")
    result = builder.parent.base.verify_stage(SOURCE / "profiles_manifest.json", config["configuration_id"], root)
    if result.get("profiles_scored") is not False or result.get("labels_read") is not False:
        raise ValueError("unexpected source scope")
    if out.exists():
        raise ValueError("use a fresh output directory")
    protocol_path = checked_path(dict(path=str(builder.PROTOCOL), sha256=builder.PROTOCOL_SHA), ROOT)
    inputs = [record(SOURCE / "configuration.freeze.json"), record(SOURCE / "profiles_manifest.json"), record(protocol_path)]
    for entry in result["stages"]:
        stage = builder.parent.base.verify_stage(checked_path(entry, root), config["configuration_id"], root)
        inputs.append(entry)
        inputs.extend(stage["outputs"].values())
    frozen = dict(inputs=inputs, stages=result["stages"], source_configuration_id=config["configuration_id"],
        chromosomes=json.loads(protocol_path.read_text())["initial_evaluation_chromosomes"],
        code=[record(Path(__file__)), record(Path(scoring.__file__))], runtime=builder.parent.extension.runtime(),
        labels_read=False, scores_evaluated=False, geometry_selected=False, promotion_eligible=False,
        reason_counts="nonexclusive; sequence and empty-strand reasons can overlap")
    frozen["configuration_id"] = canonical_id(frozen)
    out.mkdir(parents=True)
    write_immutable(out / "configuration.freeze.json", frozen)


def execute(root, out, verify_only=False):
    config = builder.parent.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    if (config["runtime"] != builder.parent.extension.runtime()
            or any(config[k] is not False for k in ("labels_read", "scores_evaluated", "geometry_selected", "promotion_eligible"))):
        raise ValueError("scope or runtime changed")
    for entry in config["inputs"] + config["code"]:
        checked_path(entry, root)
    path = out / "coverage_manifest.json"
    if path.exists():
        builder.parent.base.verify_stage(path, config["configuration_id"], root)
        print("Verified label-free training coverage summary", flush=True)
        return
    if verify_only or list(out.glob("*.tsv")):
        raise ValueError("missing result or partial outputs preserved")
    summaries, positions = [], []
    for entry in config["stages"]:
        stage = builder.parent.base.verify_stage(checked_path(entry, root), config["source_configuration_id"], root)
        sites = pd.read_csv(checked_path(stage["outputs"]["sites"], root), sep="\t")
        if not sites.TFBS_chr.isin(config["chromosomes"]).all():
            raise ValueError("non-training chromosomes")
        with np.load(checked_path(stage["outputs"]["profiles"], root), allow_pickle=False) as data:
            np.testing.assert_array_equal(data["site_hash"], builder.strands.site_hashes(sites))
            summary, position = summarize(sites, data)
        summaries.append(summary); positions.append(position)
    summary = pd.concat(summaries, ignore_index=True)
    if len(summary) != 21 or summary.duplicated(["cell", "tf"]).any():
        raise ValueError("task denominator differs")
    outputs = {}
    for name, table in (("tasks", summary), ("positions", pd.concat(positions, ignore_index=True))):
        target = out / f"{name}.tsv"
        table.to_csv(target, sep="\t", index=False)
        outputs[name] = record(target)
    result = dict(configuration_id=config["configuration_id"], outputs=outputs,
        labels_read=False, scores_evaluated=False, promotion_eligible=False)
    result["stage_id"] = canonical_id(result)
    write_immutable(path, result)
    print("Completed all 21 training-task coverage rows without scores or labels", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--freeze-only", action="store_true")
    modes.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(); out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("use research results")
    if args.freeze_only:
        freeze(args.data_root.resolve(), out)
    else:
        with builder.parent.extension.exclusive_run(out):
            execute(args.data_root.resolve(), out, args.verify_only)
