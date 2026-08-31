#!/usr/bin/env python3
"""Plot bound/unbound aggregate profiles for frozen TF-specific candidates."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from search_tf_footprint_models import candidate_from_row, score_candidate


def normalize_profiles_for_display(profiles: np.ndarray, outer_width: int = 20) -> np.ndarray:
    outer = np.concatenate([profiles[:, :outer_width], profiles[:, -outer_width:]], axis=1)
    center = np.mean(outer, axis=1, keepdims=True)
    residual = profiles - center
    rms = np.sqrt(np.mean(np.square(outer - center), axis=1))
    positive = rms[np.isfinite(rms) & (rms > 0)]
    floor = float(np.quantile(positive, 0.25)) if len(positive) else 1.0
    return residual / np.maximum(rms, max(floor, 1e-6))[:, None]


def mean_ci(profiles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(profiles, axis=0)
    error = 1.96 * np.std(profiles, axis=0, ddof=1) / math.sqrt(max(1, len(profiles)))
    return mean, error


def plot_profiles(
    sites: pd.DataFrame,
    winners: pd.DataFrame,
    cache_dir: Path,
    comparison: pd.DataFrame,
    outdir: Path,
    flank: int,
) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    stats = []
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        cell_winners = winners[winners["cell"] == cell].sort_values("tf")
        columns = 3
        rows = math.ceil(len(cell_winners) / columns)
        figure, axes = plt.subplots(rows, columns, figsize=(5.2 * columns, 3.8 * rows), squeeze=False)
        x = np.arange(-flank, flank + 1)
        for axis, winner in zip(axes.flat, cell_winners.itertuples(index=False)):
            candidate = candidate_from_row(winner)
            profiles = np.load(cache_dir / f"{cell}.{candidate.correction}.flank{flank}.npz")["profiles"]
            positions = np.flatnonzero(cell_sites["tf"].to_numpy() == str(winner.tf))
            labels = cell_sites.iloc[positions]["chip_label"].to_numpy(dtype=int)
            selected = normalize_profiles_for_display(profiles[positions])
            colors = {0: "#4C78A8", 1: "#E45756"}
            names = {0: "ChIP-negative", 1: "ChIP-positive"}
            for label in (0, 1):
                values = selected[labels == label]
                mean, error = mean_ci(values)
                axis.plot(x, mean, color=colors[label], lw=2, label=f"{names[label]} (n={len(values)})")
                axis.fill_between(x, mean - error, mean + error, color=colors[label], alpha=0.16)
            half = candidate.center_width // 2
            axis.axvspan(-half, half, color="#777777", alpha=0.09)
            axis.axvline(0, color="#333333", lw=0.7, alpha=0.5)
            metric = comparison[(comparison["cell"] == cell) & (comparison["tf"] == winner.tf)]
            delta = float(metric.iloc[0].delta_auroc) if not metric.empty else np.nan
            auroc = float(metric.iloc[0].candidate_auroc) if not metric.empty else np.nan
            axis.set_title(f"{winner.tf} · {candidate.correction}\nAUROC {auroc:.3f}; Δ {delta:+.3f}")
            axis.set_xlabel("Distance from motif center (bp)")
            axis.set_ylabel("Outer-flank standardized cut signal")
            axis.legend(frameon=False, fontsize=8)
            axis.spines[["top", "right"]].set_visible(False)

            candidate_scores = score_candidate(profiles[positions], candidate)
            stats.append(
                {
                    "cell": cell,
                    "tf": str(winner.tf),
                    "correction": candidate.correction,
                    "candidate": candidate.identifier,
                    "positive_sites": int((labels == 1).sum()),
                    "negative_sites": int((labels == 0).sum()),
                    "positive_mean_score": float(np.mean(candidate_scores[labels == 1])),
                    "negative_mean_score": float(np.mean(candidate_scores[labels == 0])),
                    "mean_score_difference": float(
                        np.mean(candidate_scores[labels == 1]) - np.mean(candidate_scores[labels == 0])
                    ),
                }
            )
        for axis in axes.flat[len(cell_winners):]:
            axis.set_visible(False)
        figure.suptitle(f"{cell}: frozen TF-specific profiles on untouched test chromosomes", fontsize=15)
        figure.tight_layout()
        figure.savefig(outdir / f"{cell}_frozen_tf_profiles.png", dpi=180, bbox_inches="tight")
        figure.savefig(outdir / f"{cell}_frozen_tf_profiles.pdf", bbox_inches="tight")
        plt.close(figure)
    return pd.DataFrame(stats)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    args = parser.parse_args(argv)
    stats = plot_profiles(
        pd.read_csv(args.sites, sep="\t"),
        pd.read_csv(args.winners, sep="\t"),
        args.cache_dir,
        pd.read_csv(args.comparison, sep="\t"),
        args.outdir,
        args.flank,
    )
    stats.to_csv(args.outdir / "frozen_tf_profile_stats.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
