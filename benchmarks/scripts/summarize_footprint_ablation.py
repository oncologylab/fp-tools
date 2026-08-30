#!/usr/bin/env python3
"""Summarize depth, correction, and scorer ablations without pooling sites."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


TASK_COLUMNS = ["cell", "tf", "motif_id"]
REQUIRED_COLUMNS = TASK_COLUMNS + ["correction", "method", "depth", "seed", "auroc", "auprc"]


def depth_value(value: object) -> float:
    if str(value).strip().lower() == "full":
        return np.inf
    return float(value)


def validate_metrics(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"ablation metrics are missing columns: {', '.join(missing)}")
    duplicate = frame.duplicated(TASK_COLUMNS + ["correction", "method", "depth", "seed"])
    if duplicate.any():
        raise ValueError("ablation metrics contain duplicate task/arm/depth/seed rows")


def aggregate_randomizations(frame: pd.DataFrame) -> pd.DataFrame:
    validate_metrics(frame)
    metadata = [
        column
        for column in ("motif_family", "role", "split", "positive_sites", "coverage_pass", "protein_supported", "motif_ambiguous", "bias_residual")
        if column in frame
    ]
    aggregations: dict[str, tuple[str, str]] = {
        "auroc": ("auroc", "median"),
        "auroc_seed_sd": ("auroc", "std"),
        "auprc": ("auprc", "median"),
        "auprc_seed_sd": ("auprc", "std"),
        "randomizations": ("seed", "nunique"),
    }
    aggregations.update({column: (column, "first") for column in metadata})
    grouped = (
        frame.groupby(TASK_COLUMNS + ["correction", "method", "depth"], as_index=False, dropna=False)
        .agg(**aggregations)
    )
    grouped["depth_value"] = grouped["depth"].map(depth_value)
    grouped["method_arm"] = grouped["correction"].astype(str) + " / " + grouped["method"].astype(str)
    return grouped.sort_values(TASK_COLUMNS + ["method_arm", "depth_value"]).reset_index(drop=True)


def depth_diagnostics(aggregated: pd.DataFrame, plateau_delta: float = 0.01) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = TASK_COLUMNS + ["correction", "method", "method_arm"]
    for key, group in aggregated.groupby(group_columns, sort=True, dropna=False):
        ordered = group.sort_values("depth_value")
        top = ordered.iloc[-1]
        prior = ordered.iloc[-2] if len(ordered) >= 2 else None
        third = ordered.iloc[-3] if len(ordered) >= 3 else None
        last_gain = float(top.auroc - prior.auroc) if prior is not None else np.nan
        prior_gain = float(prior.auroc - third.auroc) if third is not None else np.nan
        plateau = bool(
            third is not None
            and np.isfinite(last_gain)
            and np.isfinite(prior_gain)
            and abs(last_gain) <= plateau_delta
            and abs(prior_gain) <= 2 * plateau_delta
        )
        rows.append(
            {
                **dict(zip(group_columns, key)),
                "depth_levels": len(ordered),
                "highest_depth": top.depth,
                "highest_depth_auroc": top.auroc,
                "highest_depth_auprc": top.auprc,
                "previous_depth": prior.depth if prior is not None else "",
                "last_depth_gain_auroc": last_gain,
                "previous_depth_gain_auroc": prior_gain,
                "depth_plateau": plateau,
                "plateau_delta": plateau_delta,
            }
        )
    return pd.DataFrame(rows)


def correction_diagnostics(
    aggregated: pd.DataFrame,
    current_correction: str = "fp_tools_dwm",
    raw_correction: str = "raw",
) -> pd.DataFrame:
    highest = (
        aggregated.sort_values("depth_value")
        .groupby(TASK_COLUMNS + ["correction", "method"], as_index=False, dropna=False)
        .tail(1)
    )
    rows: list[dict[str, object]] = []
    for key, group in highest.groupby(TASK_COLUMNS + ["method"], sort=True, dropna=False):
        indexed = group.set_index("correction")
        if current_correction not in indexed.index or raw_correction not in indexed.index:
            continue
        current = indexed.loc[current_correction]
        raw = indexed.loc[raw_correction]
        best = group.sort_values(["auroc", "auprc", "correction"], ascending=[False, False, True]).iloc[0]
        rows.append(
            {
                **dict(zip(TASK_COLUMNS + ["method"], key)),
                "current_correction": current_correction,
                "raw_correction": raw_correction,
                "current_auroc": current.auroc,
                "raw_auroc": raw.auroc,
                "correction_gain_auroc": float(current.auroc - raw.auroc),
                "current_auprc": current.auprc,
                "raw_auprc": raw.auprc,
                "correction_gain_auprc": float(current.auprc - raw.auprc),
                "best_correction": best.correction,
                "best_correction_auroc": best.auroc,
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--plateau-delta", type=float, default=0.01)
    parser.add_argument("--current-correction", default="fp_tools_dwm")
    parser.add_argument("--raw-correction", default="raw")
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.metrics, sep="\t")
    aggregated = aggregate_randomizations(frame)
    depth = depth_diagnostics(aggregated, plateau_delta=args.plateau_delta)
    correction = correction_diagnostics(
        aggregated,
        current_correction=args.current_correction,
        raw_correction=args.raw_correction,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    aggregated.to_csv(args.outdir / "ablation_task_metrics.tsv", sep="\t", index=False)
    depth.to_csv(args.outdir / "depth_diagnostics.tsv", sep="\t", index=False)
    correction.to_csv(args.outdir / "correction_diagnostics.tsv", sep="\t", index=False)
    print(f"wrote {len(aggregated):,} task-depth summaries to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
