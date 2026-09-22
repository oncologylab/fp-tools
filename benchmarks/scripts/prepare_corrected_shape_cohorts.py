#!/usr/bin/env python3
"""Freeze common-support matched training/tuning cohorts before detector fitting."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import reassess_corrected_safety as source
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

PROTOCOL = ROOT / "benchmarks/manifests/corrected_shape_development_v1.json"


def smd(values, labels):
    first, second = values[labels == 1], values[labels == 0]
    if not len(first) or not len(second):
        return np.full(values.shape[1], np.nan)
    delta = first.mean(axis=0) - second.mean(axis=0)
    scale = np.sqrt((first.var(axis=0) + second.var(axis=0)) / 2)
    return np.divide(delta, scale, out=np.where(delta == 0, 0.0, np.copysign(np.inf, delta)), where=scale > 0)


def match_cohorts(sites, observed, expected, dwm, valid, dwm_valid, protocol):
    sites = sites.reset_index(drop=True).copy()
    permitted = set(protocol["train_chromosomes"] + protocol["tune_chromosomes"])
    if (sites.empty or not set(sites.TFBS_chr).issubset(permitted)
            or not sites.chip_label.isin([0, 1]).all()
            or not set(sites.cell).issubset({"HepG2", "K562"})):
        raise ValueError("cohort preparation accepts development chromosomes and cells only")
    arrays = [np.asarray(a, dtype=float) for a in (observed, expected, dwm)]
    if any(a.shape != (len(sites), 201) for a in arrays):
        raise ValueError("profile dimensions do not match sites")
    if np.asarray(valid).shape != (len(sites),) or np.asarray(dwm_valid).shape != (len(sites),):
        raise ValueError("validity dimensions do not match sites")
    observed, expected, dwm = arrays
    finite = np.logical_and.reduce([np.isfinite(a).all(axis=1) & (a >= 0).all(axis=1) for a in arrays])
    total = observed.sum(axis=1)
    support = finite & valid & dwm_valid & (dwm.sum(axis=1) > 0) & ((total == 0) | (expected.sum(axis=1) > 0))
    features = np.column_stack([sites.motif_score.to_numpy(dtype=float), np.log1p(np.maximum(total, 0))])
    support &= np.isfinite(features).all(axis=1)
    sites["split"] = np.where(sites.TFBS_chr.isin(protocol["train_chromosomes"]), "train", "tune")
    pairs, rows = [], []
    spec = protocol["matching"]
    for (cell, tf), task in sites.groupby(["cell", "tf"], sort=True):
        for split in ("train", "tune"):
            all_ix = task.index[task.split.eq(split)].to_numpy()
            ix = all_ix[support[all_ix]]
            labels = sites.iloc[ix].chip_label.to_numpy(dtype=int)
            pos, neg = ix[labels == 1], ix[labels == 0]
            row = dict(cell=cell, tf=tf, split=split, original_sites=len(all_ix), supported_sites=len(ix),
                       original_positive=int(sites.iloc[all_ix].chip_label.eq(1).sum()),
                       supported_positive=len(pos), supported_negative=len(neg))
            if not len(pos) or not len(neg):
                rows.append(dict(row, matched_pairs=0, status="missing_class"))
                continue
            if len(pos) > len(neg):
                pos = np.sort(np.random.default_rng(spec["positive_subsampling_seed"]).choice(pos, len(neg), replace=False))
            fit_ix = np.concatenate([pos, neg])
            scale = features[fit_ix].std(axis=0)
            scale[scale == 0] = 1
            distance = cdist(features[pos] / scale, features[neg] / scale, metric="sqeuclidean")
            pr, nr = linear_sum_assignment(distance)
            lengths = np.sqrt(distance[pr, nr])
            keep = lengths <= spec["maximum_standardized_pair_distance"]
            pi, ni, lengths = pos[pr[keep]], neg[nr[keep]], lengths[keep]
            selected = np.concatenate([pi, ni])
            after = smd(features[selected], sites.iloc[selected].chip_label.to_numpy())
            balanced = bool(len(pi) and np.isfinite(after).all() and np.max(np.abs(after)) <= spec["maximum_absolute_standardized_mean_difference"])
            status = "eligible" if len(pi) >= protocol["minimum_sites_per_class"] and balanced else "underpowered" if len(pi) < protocol["minimum_sites_per_class"] else "unbalanced"
            before = smd(features[ix], labels)
            rows.append(dict(row, matched_pairs=len(pi), status=status,
                             before_motif_smd=before[0], before_total_smd=before[1],
                             after_motif_smd=after[0], after_total_smd=after[1]))
            for number, (p, n, distance) in enumerate(zip(pi, ni, lengths)):
                for index in (p, n):
                    pairs.append(dict(cell=cell, tf=tf, split=split, pair_id=number, site_index=int(index),
                                      chip_label=int(sites.iloc[index].chip_label), match_distance=float(distance), status=status))
    columns = ["cell", "tf", "split", "pair_id", "site_index", "chip_label", "match_distance", "status"]
    return pd.DataFrame(pairs, columns=columns), pd.DataFrame(rows), support


def freeze(root, output):
    protocol = json.loads(PROTOCOL.read_text())
    if output.exists():
        raise ValueError("freeze requires a new directory")
    old = json.loads(checked_path(protocol["dwm_source_configuration"], ROOT).read_text())
    inputs = [record(PROTOCOL), record(checked_path(protocol["dwm_source_configuration"], ROOT))]
    for entry in protocol["profile_inputs"]:
        path = checked_path(entry, ROOT)
        doc = json.loads(path.read_text())
        inputs.append(record(path))
        inputs.extend(record(checked_path(x, root)) for x in doc["outputs"].values())
        inputs.append(record(checked_path(old["cells"][entry["cell"]]["dwm_track"], root)))
    paths = [Path(__file__), Path(source.__file__), Path(source.safety.__file__),
             ROOT / "benchmarks/scripts/evaluate_chrombpnet_reference.py",
             ROOT / "benchmarks/scripts/audit_chrombpnet_comparison.py"]
    config = dict(protocol=protocol, source_configuration=old, inputs=inputs, code=[record(p) for p in paths],
                  runtime=source.builder.extension.runtime(), candidates_fitted=False, tuning_scores_read=False,
                  new_holdout_labels_read=False, promotion_eligible=False)
    config["configuration_id"] = canonical_id(config)
    output.mkdir(parents=True)
    write_immutable(output / "configuration.freeze.json", config)


def execute(root, output, verify_only=False):
    config = source.builder.extension.read_document(output / "configuration.freeze.json", "configuration_id")
    if config["runtime"] != source.builder.extension.runtime() or config["protocol"] != json.loads(PROTOCOL.read_text()):
        raise ValueError("cohort configuration or runtime changed")
    for key in ("candidates_fitted", "tuning_scores_read", "new_holdout_labels_read", "promotion_eligible"):
        if config.get(key) is not False:
            raise ValueError("invalid cohort evidence scope")
    for item in config["inputs"] + config["code"]:
        checked_path(item, root)
    manifest = output / "cohorts_manifest.json"
    if manifest.exists():
        source.builder.base.verify_stage(manifest, config["configuration_id"], root)
        print("Verified corrected shape cohorts", flush=True)
        return
    if verify_only:
        raise ValueError("missing cohort manifest")
    if list(output.glob("*.tsv")) or list(output.glob("*.npz")):
        raise ValueError("partial cohort outputs preserved")
    outputs = {}
    for entry in config["protocol"]["profile_inputs"]:
        cell = entry["cell"]
        sites, arrays = source.read_profiles(entry, ROOT)
        dwm, valid = source.safety.extract_prediction_profiles(sites,
            checked_path(config["source_configuration"]["cells"][cell]["dwm_track"], root), 100, require_dense=False)
        pairs, diagnostics, support = match_cohorts(sites, arrays["observed"], arrays["expected"],
            np.nan_to_num(dwm), arrays["valid"], valid, config["protocol"])
        for name, table in (("pairs", pairs), ("diagnostics", diagnostics)):
            path = output / f"{cell}.{name}.tsv"
            table.to_csv(path, sep="\t", index=False)
            outputs[f"{cell}_{name}"] = record(path)
        path = output / f"{cell}.dwm.npz"
        np.savez_compressed(path, expected=np.nan_to_num(dwm), support=support, site_hash=arrays["site_hash"])
        outputs[f"{cell}_dwm"] = record(path)
        print("Prepared common cohorts:", cell, flush=True)
    result = dict(configuration_id=config["configuration_id"], outputs=outputs, candidates_fitted=False, promotion_eligible=False)
    result["stage_id"] = canonical_id(result)
    write_immutable(manifest, result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--freeze-only", action="store_true")
    mode.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    output = args.outdir.resolve()
    if not output.is_relative_to(ROOT / "benchmarks/results"):
        parser.error("output must stay in research results")
    if args.freeze_only:
        freeze(args.data_root.resolve(), output)
    else:
        with source.builder.extension.exclusive_run(output):
            execute(args.data_root.resolve(), output, args.verify_only)
