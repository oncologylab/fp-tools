#!/usr/bin/env python3
"""Calibration figure for the integrated fp-tools model on genome-wide tables.

For each cell-TF feature table it computes chromosome-held-out out-of-fold
integrated probabilities, then renders three panels:

* (A) a pooled reliability diagram (predicted vs observed binding frequency) with
  the pooled expected calibration error (ECE);
* (B) per-task ECE bars;
* (C) the pooled predicted-probability histogram.

Probabilities are the same chromosome-held-out integrated scores used in the
benchmark, so the calibration shown is the held-out calibration, not a train fit.
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
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import GroupKFold, StratifiedKFold

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "benchmarks" / "scripts"))
from figure_style import apply_style, bold_all_text  # noqa: E402
from evaluate_methods import integrated_oof_scores  # noqa: E402


def crossfit_isotonic(probs: np.ndarray, labels: np.ndarray, groups: np.ndarray | None) -> np.ndarray:
    """Leak-free isotonic recalibration of out-of-fold probabilities.

    Isotonic regression is fit on training folds and applied to the held-out fold,
    grouped by chromosome when available, so calibrated probabilities are never fit
    on their own data.
    """
    out = np.full(len(probs), np.nan)
    if groups is not None and len(np.unique(groups)) >= 2:
        splitter = GroupKFold(n_splits=min(5, len(np.unique(groups))))
        split = splitter.split(probs, labels, groups)
    else:
        splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=2026)
        split = splitter.split(probs, labels)
    for tr, te in split:
        iso = IsotonicRegression(out_of_bounds="clip")
        iso.fit(probs[tr], labels[tr])
        out[te] = iso.predict(probs[te])
    return out


def expected_calibration_error(labels: np.ndarray, probs: np.ndarray, n_bins: int = 10) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.clip(np.digitize(probs, edges[1:-1]), 0, n_bins - 1)
    conf = np.full(n_bins, np.nan)
    acc = np.full(n_bins, np.nan)
    weight = np.zeros(n_bins)
    n = len(probs)
    ece = 0.0
    for b in range(n_bins):
        m = idx == b
        if not m.any():
            continue
        conf[b] = probs[m].mean()
        acc[b] = labels[m].mean()
        weight[b] = m.sum() / n
        ece += weight[b] * abs(acc[b] - conf[b])
    return float(ece), conf, acc, weight


def collect(tables: list[tuple[str, str, str]], seed: int = 2026):
    per_task = {}
    pooled_labels, pooled_raw, pooled_cal = [], [], []
    for cell, tf, path in tables:
        frame = pd.read_csv(path, sep="\t")
        has_fp = "footprint" in frame.columns and frame["footprint"].abs().sum() > 0
        feats = ["accessibility", "motif", "gc"] + (["footprint"] if has_fp else [])
        groups = frame["chrom"].to_numpy() if "chrom" in frame.columns else None
        probs = integrated_oof_scores(frame, feats, seed=seed, groups=groups)
        labels = frame["label"].to_numpy(dtype=int)
        mask = ~np.isnan(probs)
        if mask.sum() < 20:
            continue
        cal = crossfit_isotonic(probs[mask], labels[mask], groups[mask] if groups is not None else None)
        cmask = ~np.isnan(cal)
        ece_raw, *_ = expected_calibration_error(labels[mask], probs[mask])
        ece_cal, *_ = expected_calibration_error(labels[mask][cmask], cal[cmask])
        per_task[f"{cell}\n{tf}"] = (ece_raw, ece_cal)
        pooled_labels.append(labels[mask][cmask])
        pooled_raw.append(probs[mask][cmask])
        pooled_cal.append(cal[cmask])
    return per_task, np.concatenate(pooled_labels), np.concatenate(pooled_raw), np.concatenate(pooled_cal)


def plot(tables: list[tuple[str, str, str]], out_prefix: str | Path, base_size: int = 11, seed: int = 2026) -> list[Path]:
    apply_style(base_size)
    per_task, labels, raw, cal = collect(tables, seed=seed)
    ece_raw, conf_r, acc_r, _ = expected_calibration_error(labels, raw, n_bins=10)
    ece_cal, conf_c, acc_c, _ = expected_calibration_error(labels, cal, n_bins=10)

    fig = plt.figure(figsize=(14, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.4, 1.0])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])

    # (A) reliability diagram: raw (over-confident) vs isotonic-recalibrated
    axA.plot([0, 1], [0, 1], color="black", ls="--", lw=1.0, label="Perfect")
    vr = ~np.isnan(conf_r); vc = ~np.isnan(conf_c)
    axA.plot(conf_r[vr], acc_r[vr], "o-", color="#d62728", lw=1.6, markersize=5,
             label=f"Raw (ECE {ece_raw:.3f})")
    axA.plot(conf_c[vc], acc_c[vc], "s-", color="#1f77b4", lw=1.6, markersize=5,
             label=f"Isotonic (ECE {ece_cal:.3f})")
    axA.set_xlabel("Mean predicted probability")
    axA.set_ylabel("Observed binding fraction")
    axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    axA.set_title("(A) Reliability (chrom-held-out)")
    axA.legend(frameon=False, loc="upper left", fontsize=base_size - 3)
    axA.grid(alpha=0.3)
    bold_all_text(axA)

    # (B) per-task ECE: raw vs calibrated
    items = sorted(per_task.items(), key=lambda kv: kv[1][0])
    names = [k for k, _ in items]
    raw_v = [v[0] for _, v in items]
    cal_v = [v[1] for _, v in items]
    x = np.arange(len(names))
    axB.bar(x - 0.2, raw_v, width=0.4, color="#d62728", label="Raw")
    axB.bar(x + 0.2, cal_v, width=0.4, color="#1f77b4", label="Isotonic")
    axB.set_xticks(x)
    axB.set_xticklabels(names, rotation=90, fontsize=base_size - 4, fontweight="bold")
    axB.set_ylabel("ECE")
    axB.set_title("(B) Per-task calibration error (raw vs isotonic)")
    axB.legend(frameon=False, fontsize=base_size - 2)
    axB.grid(axis="y", alpha=0.3)
    bold_all_text(axB)

    # (C) calibrated probability histogram
    axC.hist(cal, bins=20, color="#9e9e9e", edgecolor="black", linewidth=0.5)
    axC.set_xlabel("Calibrated probability")
    axC.set_ylabel("Sites")
    axC.set_title("(C) Calibrated probability\ndistribution (pooled)")
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
    parser.add_argument("--tables", nargs="+", required=True, help="cell:tf:path triples.")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--base-size", type=int, default=11)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    tables = [tuple(spec.split(":", 2)) for spec in args.tables]
    outputs = plot(tables, args.out_prefix, base_size=args.base_size, seed=args.seed)
    for p in outputs:
        print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
