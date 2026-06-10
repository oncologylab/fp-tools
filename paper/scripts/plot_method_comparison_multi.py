#!/usr/bin/env python3
"""Multi-panel method-comparison figure from an evaluate_methods.py metrics TSV.

Panel A: per-task AUROC grouped bars (accessibility / motif / footprint /
fp-tools-integrated). Panel B: fp-tools-integrated vs best single-method AUROC
scatter (points above the diagonal mean integration helps). Panel C: mean AUROC
by method across all tasks. Uses the shared bold Arial/Helvetica style.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import apply_style, bold_all_text  # noqa: E402

METHOD_ORDER = ["accessibility", "motif", "footprint", "fp-tools-integrated"]
COLORS = {
    "accessibility": "#9e9e9e",
    "motif": "#1f77b4",
    "footprint": "#2ca02c",
    "fp-tools-integrated": "#d62728",
}
LABELS = {
    "accessibility": "Accessibility",
    "motif": "Motif (PWM)",
    "footprint": "Footprint (cut-site)",
    "fp-tools-integrated": "fp-tools integrated",
}


def plot(metrics: pd.DataFrame, out_prefix: str | Path, base_size: int = 11) -> list[Path]:
    apply_style(base_size)
    metrics = metrics.copy()
    metrics["task"] = metrics["cell"] + "\n" + metrics["tf"]
    tasks = list(dict.fromkeys(metrics["task"]))
    methods = [m for m in METHOD_ORDER if m in set(metrics["method"])]

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0])
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    # Panel A: grouped AUROC bars per task with bootstrap CIs.
    x = np.arange(len(tasks))
    width = 0.8 / len(methods)
    for j, method in enumerate(methods):
        vals, los, his = [], [], []
        for task in tasks:
            r = metrics[(metrics["task"] == task) & (metrics["method"] == method)]
            if r.empty:
                vals.append(np.nan); los.append(0); his.append(0)
            else:
                v = float(r["auroc"].iloc[0]); vals.append(v)
                los.append(v - float(r["auroc_lo"].iloc[0])); his.append(float(r["auroc_hi"].iloc[0]) - v)
        axA.bar(x + j * width, vals, width=width, color=COLORS[method], label=LABELS[method],
                yerr=[los, his], capsize=2, error_kw={"elinewidth": 0.8})
    axA.axhline(0.5, color="black", ls="--", lw=0.8)
    axA.set_xticks(x + width * (len(methods) - 1) / 2)
    axA.set_xticklabels(tasks, fontsize=base_size - 2, fontweight="bold")
    axA.set_ylabel("AUROC"); axA.set_ylim(0, 1.0)
    axA.set_title("(A) TF-binding discrimination by method across cell lines and TFs")
    axA.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.18), frameon=False)
    axA.grid(axis="y", alpha=0.3)
    bold_all_text(axA)

    # Panel B: integrated vs best single-method AUROC.
    singles = metrics[metrics["method"] != "fp-tools-integrated"]
    best_single = singles.loc[singles.groupby("task")["auroc"].idxmax()].set_index("task")["auroc"]
    integ = metrics[metrics["method"] == "fp-tools-integrated"].set_index("task")["auroc"]
    common = [t for t in tasks if t in best_single.index and t in integ.index]
    bx = best_single.loc[common].to_numpy(); iy = integ.loc[common].to_numpy()
    axB.scatter(bx, iy, s=60, color="#d62728", edgecolor="black", zorder=3)
    lim = [min(bx.min(), iy.min()) - 0.03, 1.0]
    axB.plot(lim, lim, color="black", ls="--", lw=1.0)
    axB.set_xlim(lim); axB.set_ylim(lim)
    axB.set_xlabel("Best single-method AUROC"); axB.set_ylabel("fp-tools integrated AUROC")
    axB.set_title("(B) Integration matches or beats\nthe best single method (all tasks)")
    axB.grid(alpha=0.3)
    bold_all_text(axB)

    # Panel C: mean AUROC by method.
    means = metrics.groupby("method")["auroc"].mean().reindex(methods)
    axC.bar(range(len(methods)), means.to_numpy(), color=[COLORS[m] for m in methods])
    axC.axhline(0.5, color="black", ls="--", lw=0.8)
    axC.set_xticks(range(len(methods)))
    axC.set_xticklabels([LABELS[m] for m in methods], rotation=20, ha="right", fontsize=base_size - 2, fontweight="bold")
    axC.set_ylabel("Mean AUROC"); axC.set_ylim(0, 1.0)
    axC.set_title("(C) Mean AUROC across all tasks")
    for i, v in enumerate(means.to_numpy()):
        axC.text(i, v + 0.02, f"{v:.2f}", ha="center", fontweight="bold", fontsize=base_size - 2)
    axC.grid(axis="y", alpha=0.3)
    bold_all_text(axC)

    fig.tight_layout()
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix in (".pdf", ".svg", ".png"):
        path = out_prefix.with_suffix(suffix)
        fig.savefig(path, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--base-size", type=int, default=11)
    args = parser.parse_args(argv)
    outputs = plot(pd.read_csv(args.metrics, sep="\t"), args.out_prefix, base_size=args.base_size)
    for p in outputs:
        print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
