#!/usr/bin/env python3
"""Compute binary-classification metrics for fp-tools benchmarks."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, precision_recall_curve, roc_auc_score


SCORE_METRICS = ["auroc", "auprc", "recall_at_1pct_fdr", "recall_at_5pct_fdr", "recall_at_10pct_fdr", "brier"]


METRIC_COLUMNS = [
    "group",
    "n",
    "positives",
    "auroc",
    "auprc",
    "recall_at_1pct_fdr",
    "recall_at_5pct_fdr",
    "recall_at_10pct_fdr",
    "brier",
]


def recall_at_fdr(y_true: np.ndarray, y_prob: np.ndarray, fdr_cutoff: float) -> float:
    precision, recall, _ = precision_recall_curve(y_true, y_prob)
    fdr = 1.0 - precision
    ok = recall[fdr <= fdr_cutoff]
    return float(np.max(ok)) if len(ok) else 0.0


def summarize(y_true: np.ndarray, y_score: np.ndarray, group: str = "global") -> dict[str, float | int | str]:
    positives = int(np.sum(y_true))
    has_both_classes = len(set(y_true.tolist())) == 2
    score_min = float(np.nanmin(y_score)) if len(y_score) else np.nan
    score_max = float(np.nanmax(y_score)) if len(y_score) else np.nan
    is_probability_like = score_min >= 0.0 and score_max <= 1.0
    return {
        "group": group,
        "n": int(len(y_true)),
        "positives": positives,
        "auroc": float(roc_auc_score(y_true, y_score)) if has_both_classes else np.nan,
        "auprc": float(average_precision_score(y_true, y_score)) if positives else np.nan,
        "recall_at_1pct_fdr": recall_at_fdr(y_true, y_score, 0.01) if positives else np.nan,
        "recall_at_5pct_fdr": recall_at_fdr(y_true, y_score, 0.05) if positives else np.nan,
        "recall_at_10pct_fdr": recall_at_fdr(y_true, y_score, 0.10) if positives else np.nan,
        "brier": float(brier_score_loss(y_true, y_score)) if is_probability_like else np.nan,
    }


def compute_metrics(df: pd.DataFrame, label_col: str, score_col: str, group_cols: list[str]) -> pd.DataFrame:
    rows = [summarize(df[label_col].to_numpy(), df[score_col].to_numpy(), group="global")]
    if group_cols:
        grouped = df.groupby(group_cols, dropna=False)
        for key, group_df in grouped:
            if not isinstance(key, tuple):
                key = (key,)
            group_name = "/".join(str(item) for item in key)
            rows.append(summarize(group_df[label_col].to_numpy(), group_df[score_col].to_numpy(), group=group_name))
    return pd.DataFrame(rows, columns=METRIC_COLUMNS)


def bootstrap_confidence_intervals(
    df: pd.DataFrame,
    label_col: str,
    score_col: str,
    group_cols: list[str],
    n_bootstrap: int = 1000,
    seed: int = 2026,
    ci: float = 0.95,
) -> pd.DataFrame:
    """Compute bootstrap confidence intervals for binary benchmark metrics."""

    rng = np.random.default_rng(seed)
    alpha = (1.0 - ci) / 2.0
    rows = []

    def add_group(group_name: str, group_df: pd.DataFrame) -> None:
        y_true = group_df[label_col].to_numpy()
        y_score = group_df[score_col].to_numpy()
        estimate = summarize(y_true, y_score, group=group_name)
        boot_values: dict[str, list[float]] = {metric: [] for metric in SCORE_METRICS}
        n = len(group_df)
        if n == 0:
            return
        for _ in range(int(n_bootstrap)):
            indices = rng.integers(0, n, size=n)
            boot = summarize(y_true[indices], y_score[indices], group=group_name)
            for metric in SCORE_METRICS:
                value = boot[metric]
                if not pd.isna(value):
                    boot_values[metric].append(float(value))
        for metric in SCORE_METRICS:
            values = np.asarray(boot_values[metric], dtype=float)
            rows.append(
                {
                    "group": group_name,
                    "metric": metric,
                    "estimate": estimate[metric],
                    "ci_low": float(np.quantile(values, alpha)) if len(values) else np.nan,
                    "ci_high": float(np.quantile(values, 1.0 - alpha)) if len(values) else np.nan,
                    "n_bootstrap": int(n_bootstrap),
                    "successful_bootstraps": int(len(values)),
                }
            )

    add_group("global", df)
    if group_cols:
        for key, group_df in df.groupby(group_cols, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            add_group("/".join(str(item) for item in key), group_df)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, help="TSV with binary labels and scores.")
    parser.add_argument("--out", required=True, help="Output metrics TSV.")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--score-col", default="score")
    parser.add_argument("--group-cols", nargs="*", default=["tf", "cell", "method"])
    parser.add_argument("--bootstrap", type=int, default=0, help="If >0, compute bootstrap confidence intervals with this many resamples.")
    parser.add_argument("--out-bootstrap", help="Optional output TSV for long-form bootstrap confidence intervals.")
    parser.add_argument("--seed", type=int, default=2026, help="Random seed for bootstrap confidence intervals.")
    args = parser.parse_args()

    df = pd.read_csv(args.predictions, sep="	")
    group_cols = [col for col in args.group_cols if col in df.columns]
    metrics = compute_metrics(df, args.label_col, args.score_col, group_cols)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out, sep="	", index=False)
    print(f"wrote {len(metrics)} metric rows to {out}")
    if args.bootstrap > 0:
        if not args.out_bootstrap:
            raise SystemExit("--out-bootstrap is required when --bootstrap > 0")
        bootstrap = bootstrap_confidence_intervals(
            df,
            args.label_col,
            args.score_col,
            group_cols,
            n_bootstrap=args.bootstrap,
            seed=args.seed,
        )
        out_bootstrap = Path(args.out_bootstrap)
        out_bootstrap.parent.mkdir(parents=True, exist_ok=True)
        bootstrap.to_csv(out_bootstrap, sep="	", index=False)
        print(f"wrote {len(bootstrap)} bootstrap CI rows to {out_bootstrap}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
