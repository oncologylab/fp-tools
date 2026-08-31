#!/usr/bin/env python3
"""Assemble one auditable result row for every tested cell/TF task."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


KEY = ["cell", "tf"]


def _prefixed(frame: pd.DataFrame, prefix: str, keep: tuple[str, ...] = ()) -> pd.DataFrame:
    keep_columns = set(KEY).union(keep)
    return frame.rename(
        columns={column: f"{prefix}{column}" for column in frame.columns if column not in keep_columns}
    )


def evidence_status(row: pd.Series, minimum_positives: int = 500, balance_smd: float = 0.25) -> str:
    if int(row.get("match_positive_sites", 0)) < minimum_positives:
        return "underpowered"
    if abs(float(row.get("match_after_smd_accessibility", np.nan))) > balance_smd:
        return "accessibility_confounded"
    quick_robust = (
        float(row.get("bootstrap_auroc_ci_low", -np.inf)) > 0
        and float(row.get("bootstrap_auprc_ci_low", -np.inf)) > 0
        and float(row.get("quick_delta_auroc", -np.inf)) > 0
        and float(row.get("quick_delta_auprc", -np.inf)) > 0
    )
    if quick_robust:
        return "geometry_rescued"
    classifier_gain = float(row.get("classifier_delta_auroc", -np.inf))
    classifier_auroc = float(row.get("classifier_test_auroc", -np.inf))
    if classifier_gain >= 0.03 and classifier_auroc >= 0.65:
        return "learned_shape_candidate"
    best = np.nanmax(
        [
            float(row.get("quick_candidate_auroc", np.nan)),
            float(row.get("depth_auroc_full_replicate_mean", np.nan)),
            classifier_auroc,
        ]
    )
    return "not_detected" if best < 0.60 else "weak_or_context_dependent"


def likely_drivers(row: pd.Series) -> str:
    drivers = []
    if abs(float(row.get("match_after_smd_motif_score", 0.0))) >= 0.25:
        drivers.append("motif_score_imbalanced")
    if abs(float(row.get("correction_correction_auroc_range", 0.0))) >= 0.05:
        drivers.append("correction_sensitive")
    if float(row.get("depth_delta_auroc", 0.0)) >= 0.03:
        drivers.append("depth_limited")
    if float(row.get("depth_delta_auroc", 0.0)) <= -0.03:
        drivers.append("replicate_or_depth_instability")
    if float(row.get("classifier_delta_auroc", 0.0)) >= 0.03:
        drivers.append("learnable_shape_missed_by_kernel")
    if float(row.get("quick_delta_auroc", 0.0)) >= 0.03:
        drivers.append("fixed_geometry_mismatch")
    if not drivers:
        drivers.append("weak_or_absent_profile_information")
    return ";".join(drivers)


def assemble(
    matching: pd.DataFrame,
    quick: pd.DataFrame,
    bootstrap: pd.DataFrame,
    correction: pd.DataFrame,
    depth: pd.DataFrame,
    full_grid: pd.DataFrame,
    classifiers: pd.DataFrame,
) -> pd.DataFrame:
    boot = bootstrap.pivot_table(
        index=KEY,
        columns="metric",
        values=["ci_low", "ci_high", "probability_delta_gt_zero"],
        aggfunc="first",
    )
    boot.columns = [f"bootstrap_{metric}_{stat}" for stat, metric in boot.columns]
    boot = boot.reset_index()
    output = _prefixed(matching, "match_").merge(
        _prefixed(quick, "quick_"), on=KEY, how="outer", validate="one_to_one"
    )
    for frame in (
        boot,
        _prefixed(correction, "correction_"),
        _prefixed(depth, "depth_"),
        _prefixed(full_grid, "full_grid_"),
        _prefixed(classifiers, "classifier_"),
    ):
        output = output.merge(frame, on=KEY, how="left", validate="one_to_one")
    output["classifier_delta_auroc"] = (
        output["classifier_test_auroc"] - output["quick_baseline_auroc"]
    )
    output["classifier_delta_auprc"] = (
        output["classifier_test_auprc"] - output["quick_baseline_auprc"]
    )
    output["evidence_status"] = output.apply(evidence_status, axis=1)
    output["likely_drivers"] = output.apply(likely_drivers, axis=1)
    return output.sort_values(KEY, kind="mergesort").reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("matching", "quick", "bootstrap", "correction", "depth", "full-grid", "classifiers"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = assemble(
        *(pd.read_csv(getattr(args, name.replace("-", "_")), sep="\t") for name in (
            "matching", "quick", "bootstrap", "correction", "depth", "full_grid", "classifiers"
        ))
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    columns = [
        "cell", "tf", "match_positive_sites", "evidence_status", "likely_drivers",
        "quick_delta_auroc", "depth_delta_auroc", "classifier_delta_auroc",
    ]
    print(result[columns].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
