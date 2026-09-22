#!/usr/bin/env python3
"""Preserve failed v1 cohorts and rematch using a prespecified in-caliper objective."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_corrected_shape_cohorts as base
import caliper_matching
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = base.ROOT
SOURCE = ROOT / "benchmarks/results/corrected_shape_cohorts_20260922_v1"
SOURCE_HASH = "1593a42b7d9ca3ec90b6508a2dbb88627a0732481ff319da14350fe155e7ceef"


def rematch(sites, counts, support, protocol):
    sites = sites.reset_index(drop=True)
    counts, support = np.asarray(counts, dtype=float), np.asarray(support, dtype=bool)
    if counts.shape != (len(sites),) or support.shape != counts.shape or not np.isfinite(counts).all() or np.any(counts < 0):
        raise ValueError("counts and support must match finite nonnegative site totals")
    if (not set(sites.TFBS_chr).issubset(set(protocol["train_chromosomes"] + protocol["tune_chromosomes"]))
            or not sites.chip_label.isin([0, 1]).all() or not set(sites.cell).issubset({"HepG2", "K562"})):
        raise ValueError("non-development chromosomes")
    features = np.column_stack([sites.motif_score.to_numpy(float), np.log1p(counts)])
    if not np.isfinite(features[support]).all():
        raise ValueError("supported matching features must be finite")
    split = np.where(sites.TFBS_chr.isin(protocol["train_chromosomes"]), "train", "tune")
    pairs, rows = [], []
    spec = protocol["matching"]
    for (cell, tf), task in sites.groupby(["cell", "tf"], sort=True):
        for part in ("train", "tune"):
            all_ix = task.index[split[task.index] == part].to_numpy()
            ix = all_ix[support[all_ix]]
            y = sites.iloc[ix].chip_label.to_numpy(int)
            pos, neg = ix[y == 1], ix[y == 0]
            row = dict(cell=cell, tf=tf, split=part, original_sites=len(all_ix), supported_sites=len(ix),
                       supported_positive=len(pos), supported_negative=len(neg))
            if not len(pos) or not len(neg):
                rows.append(dict(row, matched_pairs=0, status="missing_class"))
                continue
            if len(pos) > len(neg):
                pos = np.sort(np.random.default_rng(spec["positive_subsampling_seed"]).choice(pos, len(neg), replace=False))
            scale = features[np.concatenate([pos, neg])].std(axis=0)
            scale[scale == 0] = 1
            distance = cdist(features[pos]/scale, features[neg]/scale, metric="sqeuclidean")
            pr, nr = caliper_matching.caliper_assignment(distance, spec["maximum_standardized_pair_distance"])
            pi, ni = pos[pr], neg[nr]
            selected = np.concatenate([pi, ni])
            after = base.smd(features[selected], sites.iloc[selected].chip_label.to_numpy())
            balanced = bool(len(pi) and np.isfinite(after).all() and np.max(np.abs(after)) <= spec["maximum_absolute_standardized_mean_difference"])
            status = "eligible" if len(pi) >= protocol["minimum_sites_per_class"] and balanced else "underpowered" if len(pi) < protocol["minimum_sites_per_class"] else "unbalanced"
            rows.append(dict(row, matched_pairs=len(pi), status=status, after_motif_smd=after[0], after_total_smd=after[1]))
            for number, (p, n, d) in enumerate(zip(pi, ni, np.sqrt(distance[pr, nr]))):
                for index in (p, n):
                    pairs.append(dict(cell=cell, tf=tf, split=part, pair_id=number, site_index=int(index),
                                      chip_label=int(sites.iloc[index].chip_label), match_distance=float(d), status=status))
    return pd.DataFrame(pairs, columns=["cell", "tf", "split", "pair_id", "site_index", "chip_label", "match_distance", "status"]), pd.DataFrame(rows)


def freeze(root, out):
    checked_path({"path": str(SOURCE / "cohorts_manifest.json"), "sha256": SOURCE_HASH}, ROOT)
    base.execute(root, SOURCE, verify_only=True)
    if out.exists():
        raise ValueError("use a fresh output directory")
    old = json.loads((SOURCE / "configuration.freeze.json").read_text())
    previous = json.loads((SOURCE / "cohorts_manifest.json").read_text())
    protocol = old["protocol"]
    protocol["schema"] = "fp-tools-corrected-shape-development-caliper-v2"
    protocol["matching"]["method"] = "maximum_in_caliper_cardinality_then_minimum_squared_distance"
    config = dict(protocol=protocol, source_configuration=record(SOURCE / "configuration.freeze.json"),
                  source_manifest=record(SOURCE / "cohorts_manifest.json"),
                  inputs=old["inputs"] + list(previous["outputs"].values()),
                  code=old["code"] + [record(Path(__file__)), record(Path(caliper_matching.__file__))],
                  runtime=base.source.builder.extension.runtime(),
                  candidates_fitted=False, tuning_scores_read=False, promotion_eligible=False,
                  reason="v1 forced full assignment discarded feasible pairs before fitting any candidate")
    config["configuration_id"] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / "configuration.freeze.json", config)
    print("Frozen constrained matching with unchanged distance, balance and power gates", flush=True)


def execute(root, out, verify_only=False):
    config = base.source.builder.extension.read_document(out / "configuration.freeze.json", "configuration_id")
    if config["runtime"] != base.source.builder.extension.runtime() or any(config[k] is not False for k in ("candidates_fitted", "tuning_scores_read", "promotion_eligible")):
        raise ValueError("matching scope or runtime changed")
    for item in config["inputs"] + config["code"] + [config["source_configuration"], config["source_manifest"]]:
        checked_path(item, root)
    if config["source_manifest"]["sha256"] != SOURCE_HASH:
        raise ValueError("source cohort identity changed")
    expected = json.loads(checked_path(config["source_configuration"], root).read_text())["protocol"]
    expected["schema"] = "fp-tools-corrected-shape-development-caliper-v2"
    expected["matching"]["method"] = "maximum_in_caliper_cardinality_then_minimum_squared_distance"
    if config["protocol"] != expected:
        raise ValueError("only the matching objective may change from v1")
    path = out / "cohorts_manifest.json"
    if path.exists():
        base.source.builder.base.verify_stage(path, config["configuration_id"], root)
        print("Verified constrained corrected cohorts", flush=True)
        return
    if verify_only or list(out.glob("*.tsv")):
        raise ValueError("missing result or partial cohort outputs preserved")
    previous = json.loads(checked_path(config["source_manifest"], root).read_text())
    outputs = {}
    for entry in config["protocol"]["profile_inputs"]:
        cell = entry["cell"]
        sites, arrays = base.source.read_profiles(entry, ROOT)
        with np.load(checked_path(previous["outputs"][cell + "_dwm"], root), allow_pickle=False) as data:
            if not np.array_equal(data["site_hash"], arrays["site_hash"]):
                raise ValueError("common-support site identity changed")
            support = data["support"]
        pairs, rows = rematch(sites, arrays["observed"].sum(axis=1), support, config["protocol"])
        for name, table in (("pairs", pairs), ("diagnostics", rows)):
            p = out / f"{cell}.{name}.tsv"
            table.to_csv(p, sep="\t", index=False)
            outputs[f"{cell}_{name}"] = record(p)
        outputs[cell + "_dwm"] = previous["outputs"][cell + "_dwm"]
        print("Prepared constrained cohorts:", cell, flush=True)
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
        with base.source.builder.extension.exclusive_run(out):
            execute(args.data_root.resolve(), out, args.verify_only)
