#!/usr/bin/env python3
"""Plot a paired method comparison (e.g. accessibility vs motif) per TF.

Reads a binary-metrics TSV whose ``group`` column encodes ``<tf>/<method>`` (the
default produced by run_benchmark_pipeline with ``--group-cols tf method``) and
draws grouped bars of AUROC and AUPRC, one bar pair per TF, plus the per-method
mean. Saves PDF/SVG/PNG.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _split_groups(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = metrics[metrics["group"].str.contains("/", na=False)].copy()
    rows[["tf", "method"]] = rows["group"].str.split("/", n=1, expand=True)
    return rows


def plot_method_comparison(
    metrics: pd.DataFrame,
    out_prefix: str | Path,
    baseline_method: str = "accessibility",
    title: str = "Motif-aware vs accessibility baseline",
) -> list[Path]:
    """Draw grouped AUROC/AUPRC bars per TF for a baseline vs improved method."""

    rows = _split_groups(metrics)
    methods = list(dict.fromkeys(rows["method"]))
    # Put the baseline first, then the other method(s).
    methods.sort(key=lambda m: (m != baseline_method, m))
    tfs = list(dict.fromkeys(rows["tf"]))

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), constrained_layout=True)
    fig.suptitle(title)
    colors = {methods[0]: "0.6"}
    palette = ["tab:blue", "tab:green", "tab:purple"]
    for i, m in enumerate(methods[1:]):
        colors[m] = palette[i % len(palette)]

    width = 0.8 / max(len(methods), 1)
    for ax, metric, label in ((axes[0], "auroc", "AUROC"), (axes[1], "auprc", "AUPRC")):
        x = np.arange(len(tfs))
        for j, method in enumerate(methods):
            vals = [
                float(rows[(rows["tf"] == tf) & (rows["method"] == method)][metric].iloc[0])
                if not rows[(rows["tf"] == tf) & (rows["method"] == method)].empty
                else np.nan
                for tf in tfs
            ]
            ax.bar(x + j * width, vals, width=width, label=method, color=colors.get(method))
        if metric == "auroc":
            ax.axhline(0.5, color="red", linewidth=0.8, linestyle="--", label="chance")
        ax.set_xticks(x + width * (len(methods) - 1) / 2)
        ax.set_xticklabels(tfs, rotation=0)
        ax.set_ylabel(label)
        ax.set_ylim(0, 1)
        ax.set_title(f"{label} by transcription factor")
        ax.grid(axis="y", alpha=0.25)
        if metric == "auroc":
            ax.legend(fontsize=8, ncol=2)

    outputs = []
    for suffix in (".pdf", ".svg", ".png"):
        path = out_prefix.with_suffix(suffix)
        fig.savefig(path, dpi=300 if suffix == ".png" else None)
        outputs.append(path)
    plt.close(fig)
    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True, help="binary_metrics.tsv with <tf>/<method> groups.")
    parser.add_argument("--out-prefix", required=True, help="Output figure path prefix.")
    parser.add_argument("--baseline-method", default="accessibility")
    parser.add_argument("--title", default="Motif-aware vs accessibility baseline")
    args = parser.parse_args(argv)

    metrics = pd.read_csv(args.metrics, sep="\t")
    outputs = plot_method_comparison(metrics, args.out_prefix, baseline_method=args.baseline_method, title=args.title)
    for path in outputs:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
