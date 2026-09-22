#!/usr/bin/env python3
"""Enumerate tuning chromosome-block resamples for every frozen comparison.

Only three blocks are available: intervals are descriptive sensitivity, not
calibrated uncertainty or an untouched principal promotion test.
"""

import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reassess_fixed_shape_scores as source
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = source.ROOT
SOURCE = ROOT / "benchmarks/results/fixed_shape_tuning_20260922_v1"
SOURCE_SHA = "24096520aaab92947a0351a9208b32575392b5512d0a5aaab6ddf5ac71b5b65d"
BLOCKS = ("chr16", "chr17", "chr18")
SITE_KEYS = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand", "chip_label"]
KEYS = ["cell", "tf", "correction", "detector"]


def paired_scores(candidate, baseline):
    if candidate.duplicated(SITE_KEYS).any() or baseline.duplicated(SITE_KEYS).any():
        raise ValueError("duplicate sites")
    paired = candidate[SITE_KEYS + ["score"]].merge(baseline[SITE_KEYS + ["score"]],
        on=SITE_KEYS, how="outer", suffixes=("_candidate", "_baseline"), indicator=True, validate="one_to_one")
    if not paired._merge.eq("both").all():
        raise ValueError("candidate and baseline site or label support differs")
    if not paired.TFBS_chr.isin(BLOCKS).all():
        raise ValueError("only tuning chromosomes are allowed")
    if not paired.chip_label.isin([0, 1]).all() or not np.isfinite(paired[["score_candidate", "score_baseline"]]).all().all():
        raise ValueError("invalid scores or binary labels")
    return paired


def delta(frame, weights):
    weights = np.asarray(weights, dtype=float)
    labels = frame.chip_label.to_numpy()
    if weights.shape != labels.shape or not np.isfinite(weights).all() or (weights < 0).any():
        raise ValueError("invalid resample weights")
    if any(weights[labels == value].sum() <= 0 for value in (0, 1)):
        return dict(auroc_gain=np.nan, auprc_gain=np.nan, relative_auprc_gain=np.nan, valid=False)
    # Remove zero-weight sites explicitly to avoid undefined sklearn boundary
    # cases when a chromosome contains only one class.
    keep = weights > 0
    y, w = labels[keep], weights[keep]
    candidate, baseline = frame.score_candidate.to_numpy()[keep], frame.score_baseline.to_numpy()[keep]
    auc = roc_auc_score(y, candidate, sample_weight=w) - roc_auc_score(y, baseline, sample_weight=w)
    ap0 = average_precision_score(y, baseline, sample_weight=w)
    ap1 = average_precision_score(y, candidate, sample_weight=w)
    return dict(auroc_gain=float(auc), auprc_gain=float(ap1-ap0), relative_auprc_gain=float(ap1/ap0-1), valid=True)


def resamples(frame):
    rows = []
    chrom = frame.TFBS_chr.to_numpy()
    for i, selected in enumerate(itertools.product(BLOCKS, repeat=len(BLOCKS))):
        weights = sum((chrom == block).astype(int) for block in selected)
        rows.append(dict(resample=i, blocks=",".join(selected), **delta(frame, weights)))
    return pd.DataFrame(rows)


def freeze(root, out):
    checked_path(dict(path=str(SOURCE / "reassessment_manifest.json"), sha256=SOURCE_SHA), ROOT)
    source.execute(root, SOURCE, verify_only=True)
    if out.exists():
        raise ValueError("use a fresh output directory")
    manifest = json.loads((SOURCE / "reassessment_manifest.json").read_text())
    config = dict(schema="fp-tools-shape-chromosome-sensitivity-v1", source=record(SOURCE / "reassessment_manifest.json"),
        inputs=list(manifest["outputs"].values()), code=[record(Path(__file__))],
        runtime=source.fixed.parent.base.source.builder.extension.runtime(), blocks=list(BLOCKS),
        resampling="all_27_ordered_three_chromosome_draws_with_replacement", quantiles=[.025, .975],
        comparisons=["same_correction_conventional", "DWM_conventional"], retain_all_189_metric_rows=True,
        candidate_selected=False, model_refitted=False, new_holdout_labels_read=False,
        promotion_eligible=False, calibrated_uncertainty=False)
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / "configuration.freeze.json", config)
    print("Frozen diagnostic sensitivity configuration", flush=True)


def execute(root, out, verify_only=False):
    builder = source.fixed.parent.base.source.builder
    config = builder.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    if (config["source"]["sha256"] != SOURCE_SHA or config["blocks"] != list(BLOCKS)
            or config["runtime"] != builder.extension.runtime()
            or any(config[k] is not False for k in ("candidate_selected", "model_refitted", "new_holdout_labels_read", "promotion_eligible", "calibrated_uncertainty"))):
        raise ValueError("sensitivity scope changed")
    for item in config["inputs"] + config["code"] + [config["source"]]:
        checked_path(item, root)
    result_path = out / "sensitivity_manifest.json"
    if result_path.exists():
        result = builder.base.verify_stage(result_path, config["configuration_id"], root)
        if result.get("promotion_eligible") is not False or result.get("calibrated_uncertainty") is not False:
            raise ValueError("result overstates evidence")
        print("Verified completed chromosome sensitivity diagnostics", flush=True)
        return
    if verify_only or list(out.glob("*.tsv")):
        raise ValueError("missing result or partial outputs preserved")
    manifest = json.loads(checked_path(config["source"], root).read_text())
    metrics = pd.read_csv(checked_path(manifest["outputs"]["metrics"], root), sep="\t")
    scores = pd.read_csv(checked_path(manifest["outputs"]["scores"], root), sep="\t")
    if len(metrics) != 189 or metrics.duplicated(KEYS).any() or len(metrics[["cell", "tf"]].drop_duplicates()) != 21:
        raise ValueError("task denominator changed")
    groups = {key: group for key, group in scores.groupby(KEYS)}
    summaries, draws, chromosome_rows = [], [], []
    for row in metrics.itertuples(index=False):
        key = tuple(getattr(row, name) for name in KEYS)
        candidate = groups.get(key)
        for comparison in config["comparisons"]:
            baseline_key = (row.cell, row.tf, row.correction if comparison == "same_correction_conventional" else "DWM_expected", "conventional")
            baseline = groups.get(baseline_key)
            identity = dict(zip(KEYS, key), comparison=comparison)
            summary = dict(identity, source_status=row.status, status="unavailable", valid_resamples=0)
            if candidate is not None and baseline is not None:
                paired = paired_scores(candidate, baseline)
                distribution = resamples(paired)
                summary.update(delta(paired, np.ones(len(paired))), sites=len(paired),
                    observed_blocks=paired.TFBS_chr.nunique(), valid_resamples=int(distribution.valid.sum()),
                    status="self_comparison" if key == baseline_key else "diagnostic")
                for metric in ("auroc_gain", "auprc_gain", "relative_auprc_gain"):
                    # Undefined draws make the interval unavailable, rather
                    # than selecting only favorable resamples.
                    interval = np.quantile(distribution[metric], config["quantiles"]) if distribution.valid.all() else [np.nan, np.nan]
                    summary[metric + "_lower"], summary[metric + "_upper"] = interval
                for name, value in identity.items():
                    distribution[name] = value
                draws.append(distribution)
                for block in BLOCKS:
                    chromosome_rows.append(dict(identity, chromosome=block,
                        **delta(paired, paired.TFBS_chr.eq(block).to_numpy(int))))
            summaries.append(summary)
    outputs = {}
    for name, table in (("summary", pd.DataFrame(summaries)), ("resamples", pd.concat(draws, ignore_index=True)),
                        ("chromosomes", pd.DataFrame(chromosome_rows))):
        path = out / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        outputs[name] = record(path)
    result = dict(configuration_id=config["configuration_id"], outputs=outputs,
        promotion_eligible=False, calibrated_uncertainty=False, evidence="post_hoc_three_tuning_chromosome_sensitivity",
        limitations="Three chromosomes cannot establish calibrated interval coverage or the principal untouched promotion gate.")
    result["stage_id"] = canonical_id(result)
    write_immutable(result_path, result)
    print("Completed all 378 comparisons without model or task selection", flush=True)


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
        with source.fixed.parent.base.source.builder.extension.exclusive_run(out):
            execute(args.data_root.resolve(), out, args.verify_only)
