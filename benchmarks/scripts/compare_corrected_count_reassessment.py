#!/usr/bin/env python3
"""Compare historical and corrected-count summaries without changing selection."""

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reassess_corrected_strand_counts as reassessment
from audit_chrombpnet_comparison import canonical_id, check_identity, checked_path
from audit_direct_bias_safety import record, write_immutable

KEYS = ["cell", "tf", "sample", "depth", "seed", "cohort"]
CANDIDATE = "frozen_parametric_bias_conventional_geometry"
COHORT = "motif_accessibility_observed_total"


def pair_frames(old, new, keys):
    if old.duplicated(keys).any() or new.duplicated(keys).any():
        raise ValueError("duplicate comparison identity")
    paired = old.merge(new, on=keys, how="outer", suffixes=("_old", "_new"),
                       indicator=True, validate="one_to_one")
    if not paired["_merge"].eq("both").all():
        raise ValueError("before/after task or run identity differs")
    return paired.drop(columns="_merge")


def boolean(series):
    if not series.astype(str).isin(["True", "False"]).all():
        raise ValueError("shape boolean is missing or invalid")
    return series.astype(str).eq("True")


def comparisons(old, new):
    pairs = {name: pair_frames(old[name], new[name], KEYS +
             (["method", "baseline"] if name == "bootstrap" else ["method"] if name in ("metrics", "shapes") else []))
             for name in ("metrics", "shapes", "coverage", "bootstrap")}
    for name in ("metrics", "coverage"):
        frame = pairs[name]
        for field in ("n_positive", "n_negative", "status"):
            if not frame[field + "_old"].equals(frame[field + "_new"]):
                raise ValueError("site class counts or eligibility changed")
    metric = pairs["metrics"]
    for field in ("auroc", "auprc", "functional_separation", "auroc_gain_over_sample_dwm", "relative_auprc_gain_over_sample_dwm"):
        metric[field + "_change"] = metric[field + "_new"] - metric[field + "_old"]
    shape = pairs["shapes"]
    shape["protection_lost"] = boolean(shape.canonical_protection_old) & ~boolean(shape.canonical_protection_new)
    shape["protection_gained"] = ~boolean(shape.canonical_protection_old) & boolean(shape.canonical_protection_new)
    summary = []
    # Coverage contains tasks too small to score; use it as the denominator.
    for (cell, tf), coverage in pairs["coverage"].loc[pairs["coverage"].cohort.eq(COHORT)].groupby(["cell", "tf"], sort=True):
        task = metric.loc[metric.cell.eq(cell) & metric.tf.eq(tf) & metric.cohort.eq(COHORT) & metric.method.eq(CANDIDATE)]
        shapes = shape.loc[shape.cell.eq(cell) & shape.tf.eq(tf) & shape.cohort.eq(COHORT) & shape.method.eq(CANDIDATE)]
        row = {"cell": cell, "tf": tf, "cohort": COHORT, "expected_runs": len(coverage),
               "scored_runs": len(task), "powered_runs": int(coverage.status_new.eq("eligible").sum()),
               "protection_lost_runs": int(shapes.protection_lost.sum()),
               "protection_gained_runs": int(shapes.protection_gained.sum()),
               "new_canonical_runs": int(boolean(shapes.canonical_protection_new).sum()),
               "promotion_eligible": False}
        for version in ("old", "new"):
            full = task.loc[task.depth.eq("full")]
            row[f"full_depth_mean_auroc_gain_vs_sample_dwm_{version}"] = full[f"auroc_gain_over_sample_dwm_{version}"].mean()
            row[f"full_depth_mean_relative_auprc_gain_vs_sample_dwm_{version}"] = full[f"relative_auprc_gain_over_sample_dwm_{version}"].mean()
        row["largest_absolute_auroc_change"] = task.auroc_change.abs().max()
        row["largest_absolute_auprc_change"] = task.auprc_change.abs().max()
        summary.append(row)
    return pairs, pd.DataFrame(summary)


def build(root, corrected, output):
    config_path = corrected / "configuration.freeze.json"
    config = reassessment.extension.read_document(config_path, "configuration_id")
    reassessment.validate_configuration(config)
    new_path = corrected / "reassessment_manifest.json"
    new = reassessment.base.verify_stage(new_path, config["configuration_id"], root)
    if (new.get("schema") != reassessment.SCHEMA or len(new.get("stages", [])) != 96
            or any(new.get(k) is not False for k in ("promotion_eligible", "candidate_retuned", "new_holdout_labels_read", "final_untouched_control_read"))):
        raise ValueError("corrected result is not the complete non-promoted reassessment")
    # Bind every stage to its run and validate its output hashes. Full input
    # verification remains the reassessment runner's separate --verify-only.
    runs = set()
    for item in new["stages"]:
        extraction = json.loads(checked_path(item["extraction"], root).read_text())
        check_identity(extraction, "stage_id")
        stage = reassessment.base.verify_stage(checked_path(item["evaluation"], root), config["configuration_id"], root,
                                               expected_run=extraction["run"])
        run = tuple(stage["run"][k] for k in reassessment.base.RUN_KEYS)
        if run in runs:
            raise ValueError("duplicate corrected run")
        runs.add(run)
    if runs != {tuple(row[k] for k in reassessment.base.RUN_KEYS) for row in config["matrix"]}:
        raise ValueError("corrected run matrix differs")
    old = reassessment.source_documents()["stability_manifest.json"]
    before, after = {}, {}
    for name in ("metrics", "shapes", "coverage", "bootstrap"):
        before[name] = pd.read_csv(checked_path(old["outputs"][name], root), sep="\t")
        after[name] = pd.read_csv(checked_path(new["outputs"][name], root), sep="\t")
    pairs, summary = comparisons(before, after)
    if len(summary) != 21:
        raise ValueError("all 21 development tasks must remain visible")
    output.mkdir(parents=True, exist_ok=False)
    outputs = {}
    for name, table in {**pairs, "per_tf_summary": summary}.items():
        path = output / f"{name}.tsv.gz"
        table.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
        outputs[name] = record(path)
    profiles = pd.read_csv(checked_path(new["outputs"]["profiles"], root), sep="\t")
    pdf = output / "corrected_aggregate_profiles.pdf"
    reassessment.extension.render_stability(after["metrics"], profiles, pdf)
    outputs["aggregate_pdf"] = record(pdf)
    result = {"schema": "fp-tools-corrected-count-comparison-v1",
              "source_configuration": record(config_path), "corrected_manifest": record(new_path),
              "historical_manifest_sha256": reassessment.SOURCE_HASHES["stability_manifest.json"],
              "reporter": record(Path(__file__)), "outputs": outputs,
              "scope": "post_hoc_fixed_candidate_implementation_reassessment",
              "limitations": ["No new method, fitting, selection or untouched evaluation.",
                              "Changes are numerical differences, not confidence intervals for a reader effect.",
                              "Safety and complete-workflow resource gates remain separate and unestablished.",
                              "Source-input verification and visual plot inspection must be reported separately."],
              "promotion_eligible": False}
    result["stage_id"] = canonical_id(result)
    write_immutable(output / "comparison.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--corrected", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if not args.outdir.resolve().is_relative_to(reassessment.ROOT / "benchmarks/results"):
        parser.error("report output must remain in research benchmark results")
    build(args.data_root.resolve(), args.corrected.resolve(), args.outdir.resolve())
