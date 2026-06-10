#!/usr/bin/env python3
"""Compute calibration bins and expected calibration error for fp-tools predictions."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss


BIN_COLUMNS = [
    "group",
    "bin",
    "bin_low",
    "bin_high",
    "n",
    "positives",
    "mean_score",
    "observed_rate",
    "abs_gap",
]
SUMMARY_COLUMNS = ["group", "n", "positives", "ece", "mce", "brier"]


def _probability_frame(df: pd.DataFrame, label_col: str, score_col: str) -> pd.DataFrame:
    out = df.copy()
    out[label_col] = pd.to_numeric(out[label_col], errors="coerce")
    out[score_col] = pd.to_numeric(out[score_col], errors="coerce")
    out = out.dropna(subset=[label_col, score_col])
    return out[(out[score_col] >= 0.0) & (out[score_col] <= 1.0)].copy()


def calibration_for_group(df: pd.DataFrame, label_col: str, score_col: str, group: str, bins: int = 10) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    """Compute reliability bins and ECE for one probability-like prediction group."""

    if df.empty:
        empty = pd.DataFrame(columns=BIN_COLUMNS)
        return empty, {"group": group, "n": 0, "positives": 0, "ece": np.nan, "mce": np.nan, "brier": np.nan}

    labels = df[label_col].astype(float)
    scores = df[score_col].astype(float)
    edges = np.linspace(0.0, 1.0, bins + 1)
    bin_ids = np.digitize(scores.to_numpy(), edges[1:-1], right=True)
    rows = []
    total = len(df)
    ece = 0.0
    mce = 0.0
    for bin_id in range(bins):
        mask = bin_ids == bin_id
        n = int(mask.sum())
        positives = int(labels.to_numpy()[mask].sum()) if n else 0
        mean_score = float(scores.to_numpy()[mask].mean()) if n else np.nan
        observed_rate = float(labels.to_numpy()[mask].mean()) if n else np.nan
        abs_gap = abs(mean_score - observed_rate) if n else np.nan
        if n:
            ece += (n / total) * abs_gap
            mce = max(mce, abs_gap)
        rows.append(
            {
                "group": group,
                "bin": bin_id + 1,
                "bin_low": float(edges[bin_id]),
                "bin_high": float(edges[bin_id + 1]),
                "n": n,
                "positives": positives,
                "mean_score": mean_score,
                "observed_rate": observed_rate,
                "abs_gap": abs_gap,
            }
        )
    summary = {
        "group": group,
        "n": int(total),
        "positives": int(labels.sum()),
        "ece": float(ece),
        "mce": float(mce),
        "brier": float(brier_score_loss(labels, scores)),
    }
    return pd.DataFrame(rows, columns=BIN_COLUMNS), summary


def compute_calibration(df: pd.DataFrame, label_col: str, score_col: str, group_cols: list[str], bins: int = 10) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute calibration bins and summaries globally and by optional groups."""

    prob_df = _probability_frame(df, label_col, score_col)
    bin_frames = []
    summaries = []
    bins_df, summary = calibration_for_group(prob_df, label_col, score_col, "global", bins=bins)
    bin_frames.append(bins_df)
    summaries.append(summary)

    if group_cols and not prob_df.empty:
        for key, group_df in prob_df.groupby(group_cols, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            group_name = "/".join(str(item) for item in key)
            bins_df, summary = calibration_for_group(group_df, label_col, score_col, group_name, bins=bins)
            bin_frames.append(bins_df)
            summaries.append(summary)

    return pd.concat(bin_frames, ignore_index=True), pd.DataFrame(summaries, columns=SUMMARY_COLUMNS)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, help="TSV with binary labels and probability-like scores.")
    parser.add_argument("--out-bins", required=True, help="Output reliability-bin TSV.")
    parser.add_argument("--out-summary", required=True, help="Output ECE/Brier summary TSV.")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--score-col", default="score")
    parser.add_argument("--group-cols", nargs="*", default=["tf", "cell", "method"])
    parser.add_argument("--bins", type=int, default=10)
    args = parser.parse_args()

    df = pd.read_csv(args.predictions, sep="	")
    group_cols = [column for column in args.group_cols if column in df.columns]
    bins_df, summary = compute_calibration(df, args.label_col, args.score_col, group_cols, bins=args.bins)
    out_bins = Path(args.out_bins)
    out_summary = Path(args.out_summary)
    out_bins.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    bins_df.to_csv(out_bins, sep="	", index=False)
    summary.to_csv(out_summary, sep="	", index=False)
    print(f"wrote {len(bins_df)} calibration-bin rows to {out_bins}")
    print(f"wrote {len(summary)} calibration-summary rows to {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
