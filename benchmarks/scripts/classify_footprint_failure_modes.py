#!/usr/bin/env python3
"""Classify matched-label footprint benchmark tasks into diagnostic modes.

This helper consumes long-form, task-level metrics produced by correction and
scoring ablations.  Its labels are operational benchmark diagnoses, not claims
about transcription-factor biology.  In particular, ``atac_information_limited``
requires orthogonal occupancy labels, adequate coverage, protein support, and a
depth plateau; it must not be inferred from a low footprint score alone.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


TASK_COLUMNS = ("cell", "tf", "motif_id")
REQUIRED_COLUMNS = TASK_COLUMNS + ("method", "auroc")


@dataclass(frozen=True)
class DiagnosticThresholds:
    minimum_positive_sites: int = 500
    correction_delta: float = 0.03
    scorer_delta: float = 0.05
    information_limit_auroc: float = 0.65
    weak_auroc: float = 0.70
    detectable_auroc: float = 0.75


def as_bool(value: object, default: bool = False) -> bool:
    if pd.isna(value):
        return default
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float)):
        return bool(value)
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "yes", "y", "1"}:
        return True
    if normalized in {"false", "f", "no", "n", "0", ""}:
        return False
    raise ValueError(f"Cannot interpret boolean value: {value}")


def _consistent_value(group: pd.DataFrame, column: str, default: object) -> object:
    if column not in group:
        return default
    values = group[column].dropna().unique()
    if not len(values):
        return default
    if len(values) > 1:
        task = "/".join(str(group.iloc[0][item]) for item in TASK_COLUMNS)
        raise ValueError(f"Inconsistent {column} values for task {task}")
    return values[0]


def validate_metrics(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"Metrics table is missing required columns: {', '.join(missing)}")
    duplicated = frame.duplicated([*TASK_COLUMNS, "method"])
    if duplicated.any():
        row = frame.loc[duplicated].iloc[0]
        task = "/".join(str(row[column]) for column in TASK_COLUMNS)
        raise ValueError(f"Duplicate method {row['method']} for task {task}")


def normalize_metrics_schema(frame: pd.DataFrame) -> pd.DataFrame:
    """Accept the canonical schema and aliases used by existing benchmark tables."""

    normalized = frame.copy()
    aliases = {
        "motif": "motif_id",
        "chip_positive_sites": "positive_sites",
    }
    for source, target in aliases.items():
        if target not in normalized and source in normalized:
            normalized[target] = normalized[source]
    if "auprc" not in normalized:
        normalized["auprc"] = np.nan
    return normalized


def classify_task(
    group: pd.DataFrame,
    current_method: str,
    raw_method: str,
    thresholds: DiagnosticThresholds,
) -> dict[str, object]:
    first = group.iloc[0]
    result: dict[str, object] = {column: first[column] for column in TASK_COLUMNS}
    indexed = group.set_index("method")
    if current_method not in indexed.index:
        task = "/".join(str(first[column]) for column in TASK_COLUMNS)
        raise ValueError(f"Current method {current_method!r} is absent for task {task}")

    current = indexed.loc[current_method]
    valid = group.dropna(subset=["auroc"])
    task = "/".join(str(first[column]) for column in TASK_COLUMNS)
    if pd.isna(current["auroc"]):
        raise ValueError(f"Current method {current_method!r} has no AUROC for task {task}")
    if valid.empty:
        raise ValueError(f"No finite AUROC values are available for task {task}")
    best = valid.sort_values(["auroc", "auprc", "method"], ascending=[False, False, True]).iloc[0]
    raw = indexed.loc[raw_method] if raw_method in indexed.index else None
    positive_sites = int(_consistent_value(group, "positive_sites", 0))
    coverage_pass = as_bool(_consistent_value(group, "coverage_pass", True), default=True)
    depth_plateau = as_bool(_consistent_value(group, "depth_plateau", False))
    protein_supported = as_bool(_consistent_value(group, "protein_supported", False))
    motif_ambiguous = as_bool(_consistent_value(group, "motif_ambiguous", False))
    bias_residual = as_bool(_consistent_value(group, "bias_residual", False))

    current_auroc = float(current["auroc"])
    current_auprc = float(current["auprc"])
    best_auroc = float(best["auroc"])
    raw_auroc = float(raw["auroc"]) if raw is not None and pd.notna(raw["auroc"]) else np.nan
    correction_gain = current_auroc - raw_auroc if np.isfinite(raw_auroc) else np.nan
    scorer_gain = best_auroc - current_auroc

    if positive_sites < thresholds.minimum_positive_sites:
        status = "insufficient_orthogonal_labels"
    elif not coverage_pass:
        status = "not_callable_low_coverage"
    elif motif_ambiguous:
        status = "not_callable_motif_ambiguous"
    elif bias_residual:
        status = "undercorrection_bias_residual"
    elif np.isfinite(correction_gain) and correction_gain <= -thresholds.correction_delta:
        status = "correction_sensitive"
    elif current_auroc < thresholds.weak_auroc and scorer_gain >= thresholds.scorer_delta:
        status = "scorer_limited"
    elif (
        best_auroc < thresholds.information_limit_auroc
        and depth_plateau
        and protein_supported
    ):
        status = "atac_information_limited"
    elif current_auroc >= thresholds.detectable_auroc:
        status = "detectable"
    elif current_auroc < thresholds.weak_auroc:
        status = "weak_site_discrimination_unresolved"
    else:
        status = "intermediate"

    result.update(
        {
            "diagnostic_status": status,
            "positive_sites": positive_sites,
            "current_method": current_method,
            "current_auroc": current_auroc,
            "current_auprc": current_auprc,
            "raw_method": raw_method if raw is not None else "",
            "raw_auroc": raw_auroc,
            "correction_gain_auroc": correction_gain,
            "best_method": str(best["method"]),
            "best_auroc": best_auroc,
            "best_auprc": float(best["auprc"]),
            "best_method_gain_auroc": scorer_gain,
            "coverage_pass": coverage_pass,
            "depth_plateau": depth_plateau,
            "protein_supported": protein_supported,
            "motif_ambiguous": motif_ambiguous,
            "bias_residual": bias_residual,
            "interpretation": (
                "Operational matched-label diagnosis; do not infer TF absence or a biological mechanism from this status alone."
            ),
        }
    )
    return result


def classify_failure_modes(
    frame: pd.DataFrame,
    current_method: str,
    raw_method: str,
    thresholds: DiagnosticThresholds = DiagnosticThresholds(),
) -> pd.DataFrame:
    frame = normalize_metrics_schema(frame)
    validate_metrics(frame)
    rows = [
        classify_task(group, current_method, raw_method, thresholds)
        for _, group in frame.groupby(list(TASK_COLUMNS), sort=True, dropna=False)
    ]
    return pd.DataFrame(rows).sort_values(
        ["diagnostic_status", "current_auroc", *TASK_COLUMNS]
    ).reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--current-method", default="fp-tools footprint")
    parser.add_argument("--raw-method", default="raw footprint")
    parser.add_argument("--minimum-positive-sites", type=int, default=500)
    parser.add_argument("--correction-delta", type=float, default=0.03)
    parser.add_argument("--scorer-delta", type=float, default=0.05)
    parser.add_argument("--information-limit-auroc", type=float, default=0.65)
    parser.add_argument("--weak-auroc", type=float, default=0.70)
    parser.add_argument("--detectable-auroc", type=float, default=0.75)
    args = parser.parse_args(argv)

    thresholds = DiagnosticThresholds(
        minimum_positive_sites=args.minimum_positive_sites,
        correction_delta=args.correction_delta,
        scorer_delta=args.scorer_delta,
        information_limit_auroc=args.information_limit_auroc,
        weak_auroc=args.weak_auroc,
        detectable_auroc=args.detectable_auroc,
    )
    frame = pd.read_csv(args.metrics, sep="\t")
    output = classify_failure_modes(frame, args.current_method, args.raw_method, thresholds)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out, sep="\t", index=False)
    display_columns = [*TASK_COLUMNS, "diagnostic_status", "current_auroc", "best_method", "best_auroc"]
    print(output[display_columns].to_string(index=False))
    print(f"\nwrote {len(output)} diagnostic tasks to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
