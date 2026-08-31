#!/usr/bin/env python3
"""Hold TF geometry fixed and transfer it across raw/PWM/DWM signal arms."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from search_tf_footprint_models import binary_metrics, candidate_from_row, score_candidate


def evaluate_corrections(
    profiles_by_correction: dict[str, np.ndarray],
    labels: np.ndarray,
    candidate,
) -> pd.DataFrame:
    rows = []
    for correction, profiles in sorted(profiles_by_correction.items()):
        transferred = replace(candidate, correction=correction)
        rows.append(
            {
                "correction": correction,
                "candidate": transferred.identifier,
                **binary_metrics(labels, score_candidate(profiles, transferred)),
            }
        )
    result = pd.DataFrame(rows)
    result["auroc_rank"] = result["auroc"].rank(method="min", ascending=False).astype(int)
    result["auprc_rank"] = result["auprc"].rank(method="min", ascending=False).astype(int)
    return result


def evaluate(
    sites: pd.DataFrame,
    winners: pd.DataFrame,
    cache_dir: Path,
    flank: int,
    split: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    parts = []
    summaries = []
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        caches = {
            path.name.split(".")[1]: np.load(path)["profiles"]
            for path in cache_dir.glob(f"{cell}.*.flank{flank}.npz")
        }
        for winner in winners[winners["cell"] == cell].itertuples(index=False):
            positions = np.flatnonzero(
                (cell_sites["tf"].to_numpy() == str(winner.tf))
                & (cell_sites["chromosome_split"].to_numpy() == split)
            )
            labels = cell_sites.iloc[positions]["chip_label"].to_numpy(dtype=int)
            candidate = candidate_from_row(winner)
            metrics = evaluate_corrections(
                {correction: profiles[positions] for correction, profiles in caches.items()},
                labels,
                candidate,
            )
            metrics.insert(0, "cell", str(cell))
            metrics.insert(1, "tf", str(winner.tf))
            metrics["chromosome_split"] = split
            metrics["geometry_selected_with"] = candidate.correction
            metrics["is_selected_correction"] = metrics["correction"] == candidate.correction
            parts.append(metrics)
            chosen = metrics[metrics["is_selected_correction"]].iloc[0]
            summaries.append(
                {
                    "cell": str(cell), "tf": str(winner.tf), "chromosome_split": split,
                    "selected_correction": candidate.correction,
                    "best_test_auroc_correction_posthoc": str(metrics.loc[metrics["auroc"].idxmax(), "correction"]),
                    "best_test_auprc_correction_posthoc": str(metrics.loc[metrics["auprc"].idxmax(), "correction"]),
                    "selected_auroc": float(chosen.auroc),
                    "selected_auprc": float(chosen.auprc),
                    "correction_auroc_range": float(metrics["auroc"].max() - metrics["auroc"].min()),
                    "correction_auprc_range": float(metrics["auprc"].max() - metrics["auprc"].min()),
                }
            )
    return pd.concat(parts, ignore_index=True), pd.DataFrame(summaries)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--split", default="test")
    args = parser.parse_args(argv)
    metrics, summary = evaluate(
        pd.read_csv(args.sites, sep="\t"),
        pd.read_csv(args.winners, sep="\t"),
        args.cache_dir,
        args.flank,
        args.split,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.output, sep="\t", index=False)
    summary.to_csv(args.summary, sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
