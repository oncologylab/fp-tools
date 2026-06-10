#!/usr/bin/env python3
"""Build motif-removal benchmark tables from scored candidate sites.

The output is a long-form prediction table compatible with
``compute_binary_metrics.py``.  A removed motif or motif family can be scored
with zeroed baseline motif scores and one or more recovery scores from
motif-relaxed, motif-free, supervised, or reranked fp-tools outputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


DEFAULT_METADATA_COLUMNS = [
    "tf",
    "cell",
    "motif_id",
    "motif_family",
    "chrom",
    "start",
    "end",
    "name",
]
DEFAULT_RECOVERY_SCORE_COLUMNS = ["rank_score", "binding_probability", "candidate_score", "score"]


def _existing(columns: Iterable[str], df: pd.DataFrame) -> list[str]:
    return [column for column in columns if column in df.columns]


def _as_numeric_score(values: pd.Series, column: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.isna().all():
        raise ValueError(f"score column {column!r} has no numeric values")
    return numeric.fillna(float(numeric.min()))


def build_motif_removal_table(
    predictions: pd.DataFrame,
    remove_col: str,
    remove_values: Iterable[str],
    label_col: str = "label",
    baseline_score_col: str | None = None,
    recovery_score_cols: Iterable[str] | None = None,
    metadata_cols: Iterable[str] = DEFAULT_METADATA_COLUMNS,
    include_controls: bool = False,
    zero_baseline: bool = True,
) -> pd.DataFrame:
    """Return long-form scores for a motif-removal benchmark.

    Parameters
    ----------
    predictions:
        Candidate-level table with labels, motif/family annotations, and one or
        more score columns.
    remove_col:
        Column containing the motif ID or motif family to remove.
    remove_values:
        Motif IDs or families to treat as removed from the known-motif catalog.
    label_col:
        Binary label column used by downstream metric scripts.
    baseline_score_col:
        Optional strict motif score column.  When ``zero_baseline`` is true,
        removed rows receive score 0 to model a strict known-motif workflow after
        the selected motif/family has been removed.
    recovery_score_cols:
        Candidate, model, or reranked score columns to evaluate as recovery
        methods.  Missing columns are ignored.
    metadata_cols:
        Metadata columns copied to the long-form output when present.
    include_controls:
        Include non-removed rows as controls.  By default, only rows carrying a
        removed motif/family are retained for a focused recovery benchmark.
    zero_baseline:
        Set the baseline method score to zero for removed rows.
    """
    if remove_col not in predictions.columns:
        raise ValueError(f"missing remove column: {remove_col}")
    if label_col not in predictions.columns:
        raise ValueError(f"missing label column: {label_col}")

    remove_values = [str(value) for value in remove_values]
    if not remove_values:
        raise ValueError("at least one removed motif or family is required")

    df = predictions.copy()
    df["removed"] = df[remove_col].astype(str).isin(remove_values)
    if not include_controls:
        df = df[df["removed"]].copy()
    if df.empty:
        raise ValueError("no rows matched the requested motif-removal benchmark")

    recovery_score_cols = list(recovery_score_cols or DEFAULT_RECOVERY_SCORE_COLUMNS)
    score_columns = _existing(recovery_score_cols, df)
    if baseline_score_col:
        if baseline_score_col not in df.columns:
            raise ValueError(f"missing baseline score column: {baseline_score_col}")
        score_columns = [column for column in score_columns if column != baseline_score_col]

    if not score_columns and not baseline_score_col:
        raise ValueError("no recovery score columns were found")

    base_columns = _existing(metadata_cols, df)
    rows: list[pd.DataFrame] = []
    removal_target = ",".join(remove_values)

    if baseline_score_col:
        baseline = df[base_columns + [label_col, "removed"]].copy()
        baseline["score"] = _as_numeric_score(df[baseline_score_col], baseline_score_col)
        if zero_baseline:
            baseline.loc[baseline["removed"], "score"] = 0.0
        baseline["method"] = "motif_removed_baseline" if zero_baseline else baseline_score_col
        baseline["source_score_column"] = baseline_score_col
        baseline["removal_column"] = remove_col
        baseline["removal_target"] = removal_target
        rows.append(baseline)

    for column in score_columns:
        method = column
        method_df = df[base_columns + [label_col, "removed"]].copy()
        method_df["score"] = _as_numeric_score(df[column], column)
        method_df["method"] = method
        method_df["source_score_column"] = column
        method_df["removal_column"] = remove_col
        method_df["removal_target"] = removal_target
        rows.append(method_df)

    out = pd.concat(rows, ignore_index=True)
    out = out.rename(columns={label_col: "label"}) if label_col != "label" else out
    ordered = [
        "label",
        "score",
        "method",
        "removed",
        "removal_column",
        "removal_target",
        "source_score_column",
    ]
    ordered += [column for column in base_columns if column not in ordered]
    return out[ordered]


def summarize_removal_table(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, group in table.groupby(["removal_target", "method"], dropna=False):
        target, method = key
        positives = int(pd.to_numeric(group["label"], errors="coerce").fillna(0).sum())
        rows.append(
            {
                "removal_target": target,
                "method": method,
                "n": int(len(group)),
                "positives": positives,
                "removed_rows": int(group["removed"].sum()),
                "score_mean": float(group["score"].mean()),
                "score_median": float(group["score"].median()),
                "positive_score_mean": float(group.loc[group["label"].astype(int) == 1, "score"].mean()) if positives else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, help="TSV with labels, motif annotations, and score columns.")
    parser.add_argument("--out-long", required=True, help="Long-form output TSV for compute_binary_metrics.py.")
    parser.add_argument("--out-summary", help="Optional compact method summary TSV.")
    parser.add_argument("--remove-col", default="motif_family", help="Motif ID or family column to remove.")
    parser.add_argument("--remove-values", nargs="+", required=True, help="Motif IDs or families removed from the catalog.")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--baseline-score-col", help="Strict motif score column to zero for removed rows.")
    parser.add_argument("--recovery-score-cols", nargs="*", default=DEFAULT_RECOVERY_SCORE_COLUMNS)
    parser.add_argument("--metadata-cols", nargs="*", default=DEFAULT_METADATA_COLUMNS)
    parser.add_argument("--include-controls", action="store_true", help="Keep non-removed rows as controls.")
    parser.add_argument("--keep-baseline-scores", action="store_true", help="Do not zero baseline scores for removed rows.")
    args = parser.parse_args()

    predictions = pd.read_csv(args.predictions, sep="	")
    table = build_motif_removal_table(
        predictions,
        remove_col=args.remove_col,
        remove_values=args.remove_values,
        label_col=args.label_col,
        baseline_score_col=args.baseline_score_col,
        recovery_score_cols=args.recovery_score_cols,
        metadata_cols=args.metadata_cols,
        include_controls=args.include_controls,
        zero_baseline=not args.keep_baseline_scores,
    )

    out_long = Path(args.out_long)
    out_long.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_long, sep="	", index=False)
    print(f"wrote {len(table)} motif-removal prediction rows to {out_long}")

    if args.out_summary:
        summary = summarize_removal_table(table)
        out_summary = Path(args.out_summary)
        out_summary.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_summary, sep="	", index=False)
        print(f"wrote {len(summary)} motif-removal summary rows to {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
