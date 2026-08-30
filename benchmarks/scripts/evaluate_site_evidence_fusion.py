#!/usr/bin/env python3
"""Evaluate a frozen label-free footprint/PWM site-evidence fusion candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from fp_tools.tools.tfbs_model import fuse_ranked_evidence


REQUIRED_COLUMNS = {
    "cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end",
    "chip_label", "footprint_score", "pwm_score",
}


def file_record(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {
        "path": str(path.resolve()),
        "bytes": int(path.stat().st_size),
        "sha256": digest.hexdigest(),
    }


def load_study(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_site_scores(paths: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path, sep="\t") for path in paths]
    frame = pd.concat(frames, ignore_index=True)
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError("site-score table is missing columns: " + ", ".join(missing))
    key = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end"]
    duplicates = frame.duplicated(key, keep=False)
    if duplicates.any():
        examples = frame.loc[duplicates, key].head(3).astype(str).agg(":".join, axis=1)
        raise ValueError("duplicate site rows across inputs: " + ", ".join(examples))
    frame["chip_label"] = pd.to_numeric(frame["chip_label"], errors="raise").astype(int)
    if not set(frame["chip_label"].unique()).issubset({0, 1}):
        raise ValueError("chip_label must contain only 0 and 1")
    return frame


def attach_design(frame: pd.DataFrame, study: dict[str, object]) -> pd.DataFrame:
    output = frame.copy()
    chromosome_split = study["chromosome_split"]
    split_by_chromosome = {
        chromosome: split
        for split, chromosomes in chromosome_split.items()
        for chromosome in chromosomes
    }
    output["chromosome_split"] = output["TFBS_chr"].map(split_by_chromosome).fillna("excluded")
    task_metadata = {
        (str(task["cell"]), str(task["tf"])): task
        for task in study["tasks"]
    }
    output["cell_split"] = [
        str(task_metadata.get((str(cell), str(tf)), {}).get("split", "supplemental"))
        for cell, tf in zip(output["cell"], output["tf"])
    ]
    output["role"] = [
        str(task_metadata.get((str(cell), str(tf)), {}).get("role", "supplemental"))
        for cell, tf in zip(output["cell"], output["tf"])
    ]
    output["motif_family"] = [
        str(task_metadata.get((str(cell), str(tf)), {}).get("motif_family", tf))
        for cell, tf in zip(output["cell"], output["tf"])
    ]
    return output


def add_candidate_scores(frame: pd.DataFrame) -> pd.DataFrame:
    return fuse_ranked_evidence(
        frame,
        ["footprint_score", "pwm_score"],
        group_columns=["cell", "tf"],
        output_column="evidence_fusion_score",
    )


def select_evaluation_rows(
    frame: pd.DataFrame,
    *,
    unlock_development_test: bool,
    unlock_holdout: bool,
) -> pd.DataFrame:
    allowed = (
        (frame["cell_split"].isin(["development", "supplemental"]))
        & (frame["chromosome_split"] == "validation")
    )
    if unlock_development_test:
        allowed |= (
            frame["cell_split"].isin(["development", "supplemental"])
            & (frame["chromosome_split"] == "test")
        )
    if unlock_holdout:
        allowed |= (
            (frame["cell_split"] == "locked_holdout")
            & (frame["chromosome_split"] == "test")
        )
    return frame.loc[allowed].copy()


def task_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cell_split, chromosome_split, cell, tf, role, family), group in frame.groupby(
        ["cell_split", "chromosome_split", "cell", "tf", "role", "motif_family"],
        sort=True,
        dropna=False,
    ):
        labels = group["chip_label"].to_numpy(dtype=int)
        if len(np.unique(labels)) != 2:
            continue
        for method, column in (
            ("fp-tools footprint", "footprint_score"),
            ("PWM", "pwm_score"),
            ("fp-tools evidence fusion", "evidence_fusion_score"),
        ):
            scores = group[column].to_numpy(dtype=float)
            rows.append(
                {
                    "cell_split": cell_split,
                    "chromosome_split": chromosome_split,
                    "cell": cell,
                    "tf": tf,
                    "role": role,
                    "motif_family": family,
                    "method": method,
                    "n_sites": int(len(group)),
                    "positive_sites": int(labels.sum()),
                    "unique_scores": int(pd.Series(scores).nunique(dropna=True)),
                    "auroc": float(roc_auc_score(labels, scores)),
                    "auprc": float(average_precision_score(labels, scores)),
                }
            )
    return pd.DataFrame(rows)


def paired_task_deltas(metrics: pd.DataFrame) -> pd.DataFrame:
    identifiers = [
        "cell_split", "chromosome_split", "cell", "tf", "role", "motif_family",
        "n_sites", "positive_sites",
    ]
    wide = metrics.pivot(
        index=identifiers,
        columns="method",
        values=["auroc", "auprc", "unique_scores"],
    )
    wide.columns = [f"{metric}__{method}" for metric, method in wide.columns]
    wide = wide.reset_index()
    for metric in ("auroc", "auprc"):
        wide[f"delta_{metric}"] = (
            wide[f"{metric}__fp-tools evidence fusion"]
            - wide[f"{metric}__fp-tools footprint"]
        )
        wide[f"pwm_minus_footprint_{metric}"] = (
            wide[f"{metric}__PWM"] - wide[f"{metric}__fp-tools footprint"]
        )
    return wide


def paired_chromosome_bootstrap(
    frame: pd.DataFrame,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    """Bootstrap paired candidate-minus-baseline metrics by chromosome."""

    if n_bootstrap <= 0:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    rows = []
    group_columns = [
        "cell_split", "chromosome_split", "cell", "tf", "role", "motif_family",
    ]
    for key, group in frame.groupby(group_columns, sort=True, dropna=False):
        blocks = [block.index.to_numpy() for _, block in group.groupby("TFBS_chr", sort=True)]
        if not blocks:
            continue
        distributions = {"auroc": [], "auprc": []}
        for _ in range(n_bootstrap):
            selected_blocks = rng.integers(0, len(blocks), size=len(blocks))
            indices = np.concatenate([blocks[index] for index in selected_blocks])
            sample = frame.loc[indices]
            labels = sample["chip_label"].to_numpy(dtype=int)
            if len(np.unique(labels)) != 2:
                continue
            baseline = sample["footprint_score"].to_numpy(dtype=float)
            candidate = sample["evidence_fusion_score"].to_numpy(dtype=float)
            distributions["auroc"].append(
                roc_auc_score(labels, candidate) - roc_auc_score(labels, baseline)
            )
            distributions["auprc"].append(
                average_precision_score(labels, candidate)
                - average_precision_score(labels, baseline)
            )
        base = dict(zip(group_columns, key))
        for metric, values in distributions.items():
            array = np.asarray(values, dtype=float)
            rows.append(
                {
                    **base,
                    "metric": metric,
                    "ci_low": float(np.quantile(array, 0.025)) if len(array) else np.nan,
                    "ci_high": float(np.quantile(array, 0.975)) if len(array) else np.nan,
                    "probability_delta_gt_zero": float(np.mean(array > 0)) if len(array) else np.nan,
                    "successful_bootstraps": int(len(array)),
                    "requested_bootstraps": int(n_bootstrap),
                    "resampling_unit": "chromosome",
                    "n_blocks": int(len(blocks)),
                }
            )
    return pd.DataFrame(rows)


def attach_bootstrap(deltas: pd.DataFrame, bootstrap: pd.DataFrame) -> pd.DataFrame:
    if bootstrap.empty:
        return deltas
    identifiers = [
        "cell_split", "chromosome_split", "cell", "tf", "role", "motif_family",
    ]
    wide = bootstrap.pivot(
        index=identifiers,
        columns="metric",
        values=["ci_low", "ci_high", "probability_delta_gt_zero"],
    )
    wide.columns = [f"{value}_{metric}" for value, metric in wide.columns]
    return deltas.merge(wide.reset_index(), on=identifiers, how="left", validate="one_to_one")


def summary_payload(deltas: pd.DataFrame, study: dict[str, object]) -> dict[str, object]:
    gates = study["promotion_gates"]
    summaries = []
    for (cell_split, chromosome_split), group in deltas.groupby(
        ["cell_split", "chromosome_split"], sort=True
    ):
        baseline_auprc = group["auprc__fp-tools footprint"].mean()
        candidate_auprc = group["auprc__fp-tools evidence fusion"].mean()
        positive_controls = group[group["role"] == "positive_control"]
        maximum_control_loss = float(
            np.maximum(0.0, -positive_controls["delta_auroc"]).max()
        ) if len(positive_controls) else None
        probability_column = "probability_delta_gt_zero_auroc"
        has_bootstrap = probability_column in group
        confident = (
            group[probability_column] >= gates["minimum_detectability_probability"]
            if has_bootstrap
            else pd.Series(False, index=group.index)
        )
        summaries.append(
            {
                "cell_split": str(cell_split),
                "chromosome_split": str(chromosome_split),
                "tasks": int(len(group)),
                "cells": int(group["cell"].nunique()),
                "families_improved_auroc": int(group.loc[group["delta_auroc"] > 0, "motif_family"].nunique()),
                "families_improved_with_bootstrap_support": (
                    int(group.loc[(group["delta_auroc"] > 0) & confident, "motif_family"].nunique())
                    if has_bootstrap
                    else None
                ),
                "mean_baseline_auroc": float(group["auroc__fp-tools footprint"].mean()),
                "mean_candidate_auroc": float(group["auroc__fp-tools evidence fusion"].mean()),
                "mean_delta_auroc": float(group["delta_auroc"].mean()),
                "mean_baseline_auprc": float(baseline_auprc),
                "mean_candidate_auprc": float(candidate_auprc),
                "relative_mean_auprc_gain": float(candidate_auprc / baseline_auprc - 1.0),
                "maximum_positive_control_auroc_loss": maximum_control_loss,
                "passes_mean_auroc_gate": bool(group["delta_auroc"].mean() >= gates["minimum_auroc_gain"]),
                "passes_relative_auprc_gate": bool(candidate_auprc / baseline_auprc - 1.0 >= gates["minimum_relative_auprc_gain"]),
                "passes_positive_control_gate": bool(
                    maximum_control_loss is None
                    or maximum_control_loss <= gates["maximum_strong_control_auroc_loss"]
                ),
            }
        )
    return {
        "candidate": "label-free within-task percentile soft-OR of footprint and PWM evidence",
        "default_changed": False,
        "summaries": summaries,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-scores", nargs="+", type=Path, required=True)
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_detectability_v1.spec.json"),
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--unlock-development-test", action="store_true")
    parser.add_argument("--unlock-holdout", action="store_true")
    parser.add_argument("--write-site-predictions", action="store_true")
    parser.add_argument("--bootstrap", type=int, default=0)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    study = load_study(args.study)
    frame = add_candidate_scores(attach_design(load_site_scores(args.site_scores), study))
    selected = select_evaluation_rows(
        frame,
        unlock_development_test=args.unlock_development_test,
        unlock_holdout=args.unlock_holdout,
    )
    if selected.empty:
        raise ValueError("no site rows matched the requested evaluation scope")
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics = task_metrics(selected)
    deltas = paired_task_deltas(metrics)
    bootstrap = paired_chromosome_bootstrap(selected, args.bootstrap, args.seed)
    deltas = attach_bootstrap(deltas, bootstrap)
    metrics.to_csv(args.outdir / "site_evidence_fusion_metrics.tsv", sep="\t", index=False)
    deltas.to_csv(args.outdir / "site_evidence_fusion_deltas.tsv", sep="\t", index=False)
    if not bootstrap.empty:
        bootstrap.to_csv(
            args.outdir / "site_evidence_fusion_bootstrap.tsv",
            sep="\t",
            index=False,
        )
    payload = summary_payload(deltas, study)
    payload["study"] = file_record(args.study)
    payload["site_score_inputs"] = [file_record(path) for path in args.site_scores]
    payload["evaluation_scope"] = {
        "unlock_development_test": bool(args.unlock_development_test),
        "unlock_holdout": bool(args.unlock_holdout),
        "bootstrap": int(args.bootstrap),
        "seed": int(args.seed),
    }
    (args.outdir / "site_evidence_fusion_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if args.write_site_predictions:
        selected.to_csv(
            args.outdir / "site_evidence_fusion_predictions.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
