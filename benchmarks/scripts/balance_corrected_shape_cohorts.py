#!/usr/bin/env python3
"""Select a covariate-balanced prefix of frozen distance-ordered matched pairs."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rematch_corrected_shape_cohorts as parent
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = parent.ROOT
SOURCE = ROOT / "benchmarks/results/corrected_shape_cohorts_20260922_v2"
SOURCE_HASH = "39300d43e6cee4da314da08999729ffd6c48202dc83667957d42c6411cec023d"
RULE = "largest_balanced_prefix_sorted_by_distance_then_pair_id"


def balanced_prefix(positive, negative, distance, pair_ids, maximum_smd):
    positive, negative = np.asarray(positive, float), np.asarray(negative, float)
    distance, pair_ids = np.asarray(distance, float), np.asarray(pair_ids)
    if (positive.ndim != 2 or positive.shape != negative.shape or positive.shape[1] != 2
            or distance.shape != (len(positive),) or pair_ids.shape != distance.shape
            or len(np.unique(pair_ids)) != len(pair_ids)
            or not all(np.isfinite(x).all() for x in (positive, negative, distance))
            or np.any(distance < 0) or not 0 < maximum_smd <= .1):
        raise ValueError("invalid pair covariates or balance threshold")
    order = np.lexsort((pair_ids, distance))
    # This is the largest balanced prefix, not a claim of globally optimal
    # cardinality among all possible balanced subsets.
    for n in range(len(order), 1, -1):
        ix = order[:n]
        smd = parent.base.smd(np.concatenate([positive[ix], negative[ix]]), np.repeat([1, 0], n))
        if np.isfinite(smd).all() and np.max(np.abs(smd)) <= maximum_smd:
            return ix, smd
    return np.array([], dtype=int), np.full(2, np.nan)


def balance(sites, counts, pairs, diagnostics, protocol):
    features = np.column_stack([sites.motif_score.to_numpy(float), np.log1p(counts)])
    output, rows = [], []
    for prior in diagnostics.to_dict("records"):
        group = pairs.loc[pairs.cell.eq(prior["cell"]) & pairs.tf.eq(prior["tf"]) & pairs.split.eq(prior["split"])].copy()
        row = dict(prior)
        row["previous_matched_pairs"] = prior["matched_pairs"]
        if group.empty:
            row.update(matched_pairs=0, status="missing_class" if prior["status"] == "missing_class" else "underpowered")
            rows.append(row)
            continue
        if group.site_index.duplicated().any() or not group.chip_label.isin([0,1]).all():
            raise ValueError("invalid frozen pair membership")
        positive = group[group.chip_label.eq(1)].sort_values("pair_id")
        negative = group[group.chip_label.eq(0)].sort_values("pair_id")
        if not np.array_equal(positive.pair_id, negative.pair_id) or len(positive)*2 != len(group):
            raise ValueError("every frozen pair must retain both classes")
        for part in (positive, negative):
            indexes = part.site_index.to_numpy(int)
            if np.any(indexes < 0) or np.any(indexes >= len(sites)):
                raise ValueError("frozen pair index outside site support")
            subset = sites.iloc[indexes]
            expected_chr = protocol["train_chromosomes"] if prior["split"] == "train" else protocol["tune_chromosomes"]
            if (not subset.chip_label.to_numpy().tolist() == part.chip_label.tolist()
                    or not subset.cell.eq(prior["cell"]).all() or not subset.tf.eq(prior["tf"]).all()
                    or not subset.TFBS_chr.isin(expected_chr).all()):
                raise ValueError("frozen pair identity or split changed")
        if not np.array_equal(positive.match_distance.to_numpy(), negative.match_distance.to_numpy()):
            raise ValueError("paired distances differ")
        ix, smd = balanced_prefix(features[positive.site_index], features[negative.site_index],
                                  positive.match_distance, positive.pair_id,
                                  protocol["matching"]["maximum_absolute_standardized_mean_difference"])
        status = "eligible" if len(ix) >= protocol["minimum_sites_per_class"] else "underpowered"
        row.update(matched_pairs=len(ix), status=status, after_motif_smd=smd[0], after_total_smd=smd[1])
        rows.append(row)
        selected = pd.concat([positive.iloc[ix], negative.iloc[ix]], ignore_index=True)
        selected["status"] = status
        output.append(selected)
    return pd.concat(output, ignore_index=True) if output else pairs.iloc[:0].copy(), pd.DataFrame(rows)


def freeze(root, out):
    checked_path({"path": str(SOURCE / "cohorts_manifest.json"), "sha256": SOURCE_HASH}, ROOT)
    parent.execute(root, SOURCE, verify_only=True)
    if out.exists():
        raise ValueError("use a fresh output directory")
    old = json.loads((SOURCE / "configuration.freeze.json").read_text())
    previous = json.loads((SOURCE / "cohorts_manifest.json").read_text())
    protocol = old["protocol"]
    protocol["schema"] = "fp-tools-corrected-shape-development-balanced-v3"
    protocol["matching"]["balance_subset_rule"] = RULE
    config = dict(protocol=protocol, source_configuration=record(SOURCE / "configuration.freeze.json"),
                  source_manifest=record(SOURCE / "cohorts_manifest.json"),
                  inputs=old["inputs"] + list(previous["outputs"].values()),
                  code=old["code"]+[record(Path(__file__))], runtime=parent.base.source.builder.extension.runtime(),
                  candidates_fitted=False, tuning_scores_read=False, promotion_eligible=False)
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / "configuration.freeze.json", config)
    print("Frozen covariate-only balance subset rule", flush=True)


def execute(root, out, verify_only=False):
    config = parent.base.source.builder.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    if config["runtime"] != parent.base.source.builder.extension.runtime() or any(config[k] is not False for k in ("candidates_fitted", "tuning_scores_read", "promotion_eligible")):
        raise ValueError("balance scope or runtime changed")
    for item in config["inputs"]+config["code"]+[config["source_configuration"], config["source_manifest"]]:
        checked_path(item, root)
    if config["source_manifest"]["sha256"] != SOURCE_HASH:
        raise ValueError("source cohort identity changed")
    expected = json.loads(checked_path(config["source_configuration"], root).read_text())["protocol"]
    expected["schema"] = "fp-tools-corrected-shape-development-balanced-v3"
    expected["matching"]["balance_subset_rule"] = RULE
    if config["protocol"] != expected:
        raise ValueError("only the balance subset rule may change")
    path = out / "cohorts_manifest.json"
    if path.exists():
        parent.base.source.builder.base.verify_stage(path, config["configuration_id"], root)
        print("Verified balanced corrected cohorts", flush=True)
        return
    if verify_only or list(out.glob("*.tsv")):
        raise ValueError("missing result or partial cohort outputs preserved")
    previous = json.loads(checked_path(config["source_manifest"], root).read_text())
    outputs = {}
    for entry in config["protocol"]["profile_inputs"]:
        cell = entry["cell"]
        sites, arrays = parent.base.source.read_profiles(entry, ROOT)
        pairs, diagnostics = [pd.read_csv(checked_path(previous["outputs"][cell + "_" + name], root), sep="\t") for name in ("pairs", "diagnostics")]
        pairs, rows = balance(sites, arrays["observed"].sum(axis=1), pairs, diagnostics, config["protocol"])
        for name, table in (("pairs", pairs), ("diagnostics", rows)):
            p = out / f"{cell}.{name}.tsv"
            table.to_csv(p, sep="\t", index=False)
            outputs[f"{cell}_{name}"] = record(p)
        outputs[cell + "_dwm"] = previous["outputs"][cell + "_dwm"]
        print("Prepared balanced cohorts:", cell, flush=True)
    result = dict(configuration_id=config["configuration_id"], outputs=outputs, candidates_fitted=False, promotion_eligible=False)
    result["stage_id"] = canonical_id(result)
    write_immutable(path, result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--freeze-only", action="store_true")
    mode.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("use research results")
    if args.freeze_only:
        freeze(args.data_root.resolve(), out)
    else:
        with parent.base.source.builder.extension.exclusive_run(out):
            execute(args.data_root.resolve(), out, args.verify_only)
