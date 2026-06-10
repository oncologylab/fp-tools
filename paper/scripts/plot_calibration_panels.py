#!/usr/bin/env python3
"""Create multi-panel calibration figures from fp-tools reliability tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_calibration(bins: pd.DataFrame, summary: pd.DataFrame, out_prefix: str | Path, title: str = "fp-tools calibration") -> list[Path]:
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    groups = [group for group in summary["group"].tolist() if group != "global"][:8]
    if not groups and "global" in set(summary["group"]):
        groups = ["global"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    fig.suptitle(title)

    ax = axes[0]
    ax.plot([0, 1], [0, 1], color="0.4", linestyle="--", linewidth=1)
    for group in groups:
        group_bins = bins[(bins["group"] == group) & (bins["n"] > 0)].sort_values("bin")
        ax.plot(group_bins["mean_score"], group_bins["observed_rate"], marker="o", linewidth=1.5, label=group)
    ax.set_xlabel("Mean predicted score")
    ax.set_ylabel("Observed positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="best")

    ax = axes[1]
    plot_summary = summary[summary["group"].isin(groups)].copy()
    ax.bar(range(len(plot_summary)), plot_summary["ece"])
    ax.set_title("Expected calibration error")
    ax.set_xticks(range(len(plot_summary)))
    ax.set_xticklabels(plot_summary["group"], rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("ECE")
    ax.grid(axis="y", alpha=0.25)

    outputs = []
    for suffix in ("pdf", "svg", "png"):
        path = out_prefix.with_suffix(f".{suffix}")
        fig.savefig(path, dpi=300 if suffix == "png" else None)
        outputs.append(path)
    plt.close(fig)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins", required=True, help="Reliability-bin TSV from compute_calibration.py.")
    parser.add_argument("--summary", required=True, help="Calibration summary TSV from compute_calibration.py.")
    parser.add_argument("--out-prefix", default="paper/figures/figure_calibration")
    parser.add_argument("--title", default="fp-tools calibration")
    args = parser.parse_args()

    bins = pd.read_csv(args.bins, sep="	")
    summary = pd.read_csv(args.summary, sep="	")
    outputs = plot_calibration(bins, summary, args.out_prefix, title=args.title)
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
