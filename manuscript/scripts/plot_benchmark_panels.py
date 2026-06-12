#!/usr/bin/env python3
"""Create a compact multi-panel benchmark summary figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_metrics(metrics: pd.DataFrame, out_prefix: str | Path, title: str = "fp-tools benchmark summary") -> list[Path]:
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    df = metrics.copy()
    df = df[df["group"] != "global"].head(20) if "global" in set(df["group"]) else df.head(20)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    fig.suptitle(title)
    panels = [
        ("auroc", "AUROC"),
        ("auprc", "AUPRC"),
        ("recall_at_5pct_fdr", "Recall at 5% FDR"),
        ("brier", "Brier score"),
    ]
    for ax, (column, label) in zip(axes.flat, panels):
        plot_df = df[["group", column]].dropna()
        ax.bar(range(len(plot_df)), plot_df[column])
        ax.set_title(label)
        ax.set_xticks(range(len(plot_df)))
        ax.set_xticklabels(plot_df["group"], rotation=75, ha="right", fontsize=7)
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
    parser.add_argument("--metrics", required=True, help="Metrics TSV from compute_binary_metrics.py.")
    parser.add_argument("--out-prefix", default="manuscript/figures/figure_benchmark_summary")
    parser.add_argument("--title", default="fp-tools benchmark summary")
    args = parser.parse_args()

    metrics = pd.read_csv(args.metrics, sep="	")
    outputs = plot_metrics(metrics, args.out_prefix, title=args.title)
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
