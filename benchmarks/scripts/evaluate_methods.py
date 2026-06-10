#!/usr/bin/env python3
"""Evaluate competing TF-binding scoring strategies on shared feature tables.

For each per-peak feature table (one cell/TF task) this scores four strategies on
identical inputs and reports AUROC/AUPRC with bootstrap confidence intervals:

* ``accessibility`` -- rank by raw accessibility magnitude (peak-caller baseline);
* ``motif``         -- best PWM log-odds (FIMO-style sequence scan);
* ``footprint``     -- Tn5 cut-site footprint-occupancy (TOBIAS-style, if present);
* ``fp-tools-integrated`` -- a cross-validated logistic model over accessibility,
  motif, GC, and footprint features (the fp-tools supervised integration).

The integrated score uses out-of-fold predictions from stratified k-fold CV so the
reported metrics are not optimistically biased.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from compute_binary_metrics import bootstrap_confidence_intervals  # noqa: E402


def integrated_oof_scores(frame: pd.DataFrame, features: list[str], seed: int = 2026) -> np.ndarray:
    X = frame[features].to_numpy(dtype=float)
    y = frame["label"].to_numpy(dtype=int)
    if y.sum() < 5 or (len(y) - y.sum()) < 5:
        return np.full(len(y), np.nan)
    model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, class_weight="balanced"))
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    return cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]


def evaluate_table(path: str | Path, cell: str, tf: str, n_bootstrap: int = 300, seed: int = 2026) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    has_fp = "footprint" in frame.columns and frame["footprint"].abs().sum() > 0
    methods = {
        "accessibility": frame["accessibility"],
        "motif": frame["motif"],
    }
    if has_fp:
        methods["footprint"] = frame["footprint"]
    feats = ["accessibility", "motif", "gc"] + (["footprint"] if has_fp else [])
    methods["fp-tools-integrated"] = pd.Series(integrated_oof_scores(frame, feats, seed=seed), index=frame.index)

    rows = []
    for method, score in methods.items():
        sub = pd.DataFrame({"label": frame["label"], "score": score}).dropna()
        if sub.empty:
            continue
        b = bootstrap_confidence_intervals(sub, "label", "score", [], n_bootstrap=n_bootstrap, seed=seed)
        g = b[b["group"] == "global"]
        au = g[g["metric"] == "auroc"].iloc[0]
        ap = g[g["metric"] == "auprc"].iloc[0]
        rows.append({
            "cell": cell, "tf": tf, "method": method,
            "n": int(len(sub)), "positives": int(sub["label"].sum()),
            "auroc": au["estimate"], "auroc_lo": au["ci_low"], "auroc_hi": au["ci_high"],
            "auprc": ap["estimate"], "auprc_lo": ap["ci_low"], "auprc_hi": ap["ci_high"],
        })
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", nargs="+", required=True,
                        help="Feature tables as cell:tf:path triples.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--bootstrap", type=int, default=300)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    frames = []
    for spec in args.tables:
        cell, tf, path = spec.split(":", 2)
        frames.append(evaluate_table(path, cell, tf, n_bootstrap=args.bootstrap, seed=args.seed))
    out = pd.concat(frames, ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\nwrote {len(out)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
