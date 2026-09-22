#!/usr/bin/env python3
"""Describe the locked 96-run development result; never certify promotion.

Reads only checksummed summary tables. Nested scientific-input verification is
a separate prerequisite, performed by extend_direct_bias_stability.py.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path


COHORT = "motif_accessibility_observed_total"
CANDIDATE = "frozen_parametric_bias_conventional_geometry"
BASELINE = "DWM_sample_matched_geometry"
TASKS = {
    "HepG2": "ARID3A CTCF FOXA1 HMGA1 MEF2A MEF2D MYC TBP TCF7L2 ZNF362 ZNF384 ZNF558".split(),
    "K562": "ARID3A CTCF FOXA1 MEF2A MEF2D MYC TBP TCF7L2 ZNF384".split(),
}
KEYS = ("cell", "tf", "sample", "depth", "seed")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def key(row):
    return tuple(str(row[k]) for k in KEYS)


def number(row, name):
    value = float(row[name])
    if not math.isfinite(value):
        raise ValueError(f"nonfinite {name} at {key(row)}")
    return value


def truth(value):
    if str(value) not in ("True", "False"):
        raise ValueError("invalid protection boolean")
    return str(value) == "True"


def indexed(rows):
    result = {}
    for row in rows:
        k = key(row)
        if k in result:
            raise ValueError(f"duplicate run: {k}")
        result[k] = row
    return result


def expected_keys():
    return {(cell, tf, f"{cell}_rep{rep}", depth, str(seed))
            for cell, tfs in TASKS.items() for tf in tfs for rep in (1, 2, 3)
            for depth in ("10m", "25m", "50m", "full")
            for seed in ([2026] if depth == "full" else range(2026, 2031))}


def summarize(tables):
    selected = {name: [r for r in rows if r["cohort"] == COHORT]
                for name, rows in tables.items()}
    coverage = indexed(selected["coverage"])
    if set(coverage) != expected_keys():
        raise ValueError("coverage must retain all 21 tasks and 48 runs per task")
    candidate = indexed(r for r in selected["metrics"] if r["method"] == CANDIDATE)
    baseline = indexed(r for r in selected["metrics"] if r["method"] == BASELINE)
    shape = indexed(r for r in selected["shapes"] if r["method"] == CANDIDATE)
    intervals = indexed(r for r in selected["bootstrap"]
                        if r["method"] == CANDIDATE and r["baseline"] == BASELINE)
    scorable = {k for k, r in coverage.items()
                if min(int(r["n_positive"]), int(r["n_negative"])) >= 2}
    if any(set(table) != scorable for table in (candidate, baseline, shape)):
        raise ValueError("candidate, comparator and shape rows must cover every scorable run")
    eligible = {k for k in scorable if coverage[k]["status"] == "eligible"}
    required_intervals = {k for k in eligible if k[3] == "full" or (k[3], k[4]) == ("10m", "2026")}
    if set(intervals) != required_intervals:
        raise ValueError("principal interval rows differ from powered prespecified runs")
    ledger, exceptions = [], []
    for cell, tfs in TASKS.items():
        for tf in tfs:
            keys = sorted(k for k in coverage if k[:2] == (cell, tf))
            scored = [k for k in keys if k in scorable]
            ci = [k for k in keys if k in intervals]
            gains, pr_gains, separations = [], [], []
            for k in scored:
                row, ref = candidate[k], baseline[k]
                if any(int(row[n]) != int(ref[n]) or int(row[n]) != int(coverage[k][n])
                       for n in ("n_positive", "n_negative")):
                    raise ValueError("candidate/comparator cohort counts differ")
                gain = number(row, "auroc") - number(ref, "auroc")
                pr_gain = number(row, "auprc") / number(ref, "auprc") - 1
                if not (math.isclose(gain, number(row, "auroc_gain_over_sample_dwm"), abs_tol=1e-10)
                        and math.isclose(pr_gain, number(row, "relative_auprc_gain_over_sample_dwm"), abs_tol=1e-10)):
                    raise ValueError("saved gains do not match paired comparator")
                gains.append(gain)
                pr_gains.append(pr_gain)
                ref_sep = number(ref, "functional_separation")
                if ref_sep > 0:
                    separations.append(number(row, "functional_separation") / ref_sep - 1)
                if not truth(shape[k]["canonical_protection"]):
                    exceptions.append({**dict(zip(KEYS, k)), "reason": "canonical_protection_failed",
                                       **{n: shape[k][n] for n in ("center_contrast", "left_shoulder", "right_shoulder", "significant_center_points")}})
            lower_roc = [number(intervals[k], "auroc_gain_lower_95") for k in ci]
            lower_pr = [number(intervals[k], "relative_auprc_gain_lower_95") for k in ci]
            def bounds(values):
                return [min(values), max(values)] if values else [None, None]
            smd = [abs(number(coverage[k], "observed_total_smd")) for k in keys
                   if coverage[k]["observed_total_smd"] != ""]
            ledger.append({"cell": cell, "tf": tf, "expected_runs": 48,
                           "scored_runs": len(scored), "powered_runs": sum(k in eligible for k in keys),
                           "minimum_sites_per_class": min(min(int(coverage[k][n]) for n in ("n_positive", "n_negative")) for k in keys),
                           "canonical_protection_runs": sum(truth(shape[k]["canonical_protection"]) for k in scored),
                           "auroc_gain_range": bounds(gains), "relative_auprc_gain_range": bounds(pr_gains),
                           "relative_separation_gain_range": bounds(separations),
                           "separation_denominator_positive_runs": len(separations),
                           "principal_intervals_measured": len(ci),
                           "principal_intervals_both_positive": sum(a > 0 and b > 0 for a, b in zip(lower_roc, lower_pr)),
                           "principal_auroc_lower_bound_range": bounds(lower_roc),
                           "principal_relative_auprc_lower_bound_range": bounds(lower_pr),
                           "maximum_absolute_current_count_smd": max(smd) if smd else None,
                           "safety_status": "not_assessed_by_this_ledger",
                           "resource_status": "end_to_end_unmeasured",
                           "uncertainty_status": "functional_coverage_uncalibrated",
                           "promotion_eligible": False})
    return ledger, exceptions


def build(manifest_path, expected_sha, outdir):
    if digest(manifest_path) != expected_sha:
        raise ValueError("manifest checksum mismatch")
    manifest = json.loads(Path(manifest_path).read_text())
    if (manifest.get("schema") != "fp-tools-direct-bias-stability-96-v1"
            or manifest.get("evidence") != "post_hoc_development_complete_replicate_depth_stability"
            or any(manifest.get(k) is not False for k in
                   ("candidate_retuned", "new_holdout_labels_read", "final_untouched_control_read", "promotion_eligible"))):
        raise ValueError("unexpected evidence scope")
    tables, inputs = {}, {}
    for name in ("metrics", "shapes", "coverage", "bootstrap"):
        record = manifest["outputs"][name]
        if digest(record["path"]) != record["sha256"]:
            raise ValueError(f"{name} checksum mismatch")
        with gzip.open(record["path"], "rt") as stream:
            tables[name] = list(csv.DictReader(stream, delimiter="\t"))
        inputs[name] = record
    ledger, exceptions = summarize(tables)
    document = {"schema": "fp-tools-direct-bias-gate-ledger-v1",
                "evidence": manifest["evidence"], "source_manifest_sha256": expected_sha,
                "source_result_id": manifest["stage_id"], "inputs": inputs,
                "script_sha256": digest(__file__), "task_count": len(ledger),
                "promotion_eligible": False,
                "scope": "Descriptive summary only; no new scoring, thresholds, selection or nested-input verification. Per-task ranges are not the locked mean-gain gate. Existing exact-operator safety is separate evidence.",
                "tasks": ledger, "shape_exceptions": exceptions}
    content = json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n"
    outdir = Path(outdir)
    target = outdir / "gate_ledger.json"
    if target.exists():
        if target.read_text() != content:
            raise ValueError("refusing to overwrite different ledger")
    else:
        outdir.mkdir(parents=True, exist_ok=True)
        with target.open("x") as stream:
            stream.write(content)
    print(f"Verified ledger: {len(ledger)} tasks; {len(exceptions)} shape exceptions; no promotion")
    return document


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    build(args.manifest, args.manifest_sha256, args.outdir)
