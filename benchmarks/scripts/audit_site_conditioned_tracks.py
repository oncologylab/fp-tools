#!/usr/bin/env python3
"""Audit whether site-conditioned training profiles admit one genomic track."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reassess_corrected_safety as source
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = source.ROOT
PROFILES = ROOT / "benchmarks/results/corrected_safety_profiles_20260922_v1"
PROFILE_SHA = "dc053ee9d0e0d6154cabe0e3915a6a78ba39887ffe28e2850a5059db4aa3a9e9"
TRAIN = [f"chr{i}" for i in range(1, 16)]


def adjacent_overlaps(sites, observed, expected, valid, flank=100):
    """Deterministic adjacent distinct centers; restore genomic orientation.

    This subset can disprove a single-track representation, but absence of a
    conflict cannot prove compatibility across every overlapping window.
    """
    n, width = len(sites), 2 * flank + 1
    observed, expected = np.asarray(observed), np.asarray(expected)
    if observed.shape != (n, width) or expected.shape != observed.shape or len(valid) != n:
        raise ValueError("profile shape mismatch")
    if not np.isfinite(observed).all() or not np.isfinite(expected).all():
        raise ValueError("nonfinite profile")
    sites = sites.reset_index(drop=True).copy()
    sites["row"] = np.arange(n)
    sites["center"] = (sites.TFBS_start + sites.TFBS_end) // 2
    rows = []
    for chrom, task in sites.loc[np.asarray(valid, bool)].groupby("TFBS_chr"):
        unique = task.sort_values(["center", "tf", "TFBS_start", "TFBS_end", "TFBS_strand"], kind="stable").drop_duplicates("center")
        pairs = list(unique.itertuples(index=False))
        for left, right in zip(pairs, pairs[1:]):
            shift = int(right.center - left.center)
            if shift >= width:
                continue
            obs_left = observed[left.row][::-1] if left.TFBS_strand == "-" else observed[left.row]
            obs_right = observed[right.row][::-1] if right.TFBS_strand == "-" else observed[right.row]
            exp_left = expected[left.row][::-1] if left.TFBS_strand == "-" else expected[left.row]
            exp_right = expected[right.row][::-1] if right.TFBS_strand == "-" else expected[right.row]
            obs_delta = np.abs(obs_left[shift:] - obs_right[:width-shift])
            exp_delta = np.abs(exp_left[shift:] - exp_right[:width-shift])
            disagreement = ~np.isclose(exp_left[shift:], exp_right[:width-shift], rtol=1e-8, atol=1e-10)
            rows.append(dict(chromosome=chrom, left_center=int(left.center), right_center=int(right.center),
                left_tf=left.tf, right_tf=right.tf, overlap_positions=width-shift,
                observed_disagreements=int((obs_delta != 0).sum()), expected_disagreements=int(disagreement.sum()),
                maximum_expected_difference=float(exp_delta.max())))
    return pd.DataFrame(rows, columns=["chromosome", "left_center", "right_center", "left_tf", "right_tf",
        "overlap_positions", "observed_disagreements", "expected_disagreements", "maximum_expected_difference"])


def freeze(root, out):
    checked_path(dict(path=str(PROFILES / "configuration.freeze.json"), sha256=PROFILE_SHA), ROOT)
    if out.exists():
        raise ValueError("use a fresh output directory")
    inputs, stages = [source.safety.record(PROFILES / "configuration.freeze.json")], {}
    for cell in ("HepG2", "K562"):
        path = PROFILES / f"{cell}.validation.json"
        stages[cell] = source.safety.record(path)
        stage = json.loads(path.read_text())
        inputs.append(stages[cell])
        inputs.extend(source.safety.record(checked_path(item, root)) for item in stage["outputs"].values())
    original = json.loads((PROFILES / "configuration.freeze.json").read_text())
    code_paths = {entry["path"] for entry in original["code"]}
    code_paths.update([str(Path(__file__)), source.__file__, source.safety.__file__])
    config = dict(schema="fp-tools-site-conditioned-track-audit-v1", stages=stages, inputs=inputs,
        code=[source.safety.record(Path(path)) for path in sorted(code_paths)],
        runtime=source.builder.extension.runtime(), chromosomes=TRAIN,
        selection="valid_training_sites_first_sorted_tf_at_each_center_then_adjacent_distinct_centers",
        expected_tolerance=dict(rtol=1e-8, atol=1e-10), promotion_eligible=False,
        model_refitted=False, output_operator_changed=False, new_holdout_labels_read=False)
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    source.safety.write_immutable(out / "configuration.freeze.json", config)
    print("Frozen training-only output-contract audit", flush=True)


def execute(root, out, verify_only=False):
    config = source.builder.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    if config["chromosomes"] != TRAIN or config["runtime"] != source.builder.extension.runtime():
        raise ValueError("scope or runtime changed")
    for entry in config["inputs"] + config["code"]:
        checked_path(entry, root)
    result_path = out / "contract_manifest.json"
    if result_path.exists():
        source.builder.base.verify_stage(result_path, config["configuration_id"], root)
        print("Verified output-contract audit", flush=True)
        return
    if verify_only or list(out.glob("*.tsv")):
        raise ValueError("missing result or partial outputs preserved")
    rows, summaries = [], []
    for cell, stage in config["stages"].items():
        sites, arrays = source.read_profiles(stage, root)
        selection = sites.TFBS_chr.isin(TRAIN).to_numpy()
        pairs = adjacent_overlaps(sites.loc[selection], arrays["observed"][selection], arrays["expected"][selection], arrays["valid"][selection])
        pairs["cell"] = cell
        rows.append(pairs)
        if pairs.observed_disagreements.sum():
            raise ValueError("observed cuts disagree across shared genomic positions")
        summaries.append(dict(cell=cell, training_sites=int(selection.sum()), adjacent_overlapping_pairs=len(pairs),
            pairs_with_expected_conflicts=int(pairs.expected_disagreements.gt(0).sum()),
            expected_conflicts=int(pairs.expected_disagreements.sum()), observed_conflicts=0,
            single_track_representation_disproved=bool(pairs.expected_disagreements.gt(0).any()),
            full_output_contract_passed=False))
    outputs = {}
    for name, table in (("overlaps", pd.concat(rows, ignore_index=True)), ("summary", pd.DataFrame(summaries))):
        path = out / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        outputs[name] = source.safety.record(path)
    result = dict(configuration_id=config["configuration_id"], outputs=outputs, promotion_eligible=False,
        evidence="training_only_output_representation_audit", output_operator_changed=False,
        limitation="Adjacent-window subset can disprove consistency but cannot certify every output product.")
    result["stage_id"] = canonical_id(result)
    source.safety.write_immutable(result_path, result)
    print("Completed training-only output-contract audit", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--freeze-only", action="store_true")
    modes.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("use research results")
    if args.freeze_only:
        freeze(args.data_root.resolve(), out)
    else:
        with source.builder.extension.exclusive_run(out):
            execute(args.data_root.resolve(), out, args.verify_only)
