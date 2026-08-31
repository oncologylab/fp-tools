#!/usr/bin/env python3
"""Summarize frozen TF tests into detectable, weak, and abstention states."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def classify_row(
    row: pd.Series,
    minimum_positives: int,
    maximum_accessibility_smd: float,
    maximum_motif_smd: float,
    information_limit_auroc: float,
    weak_auroc: float,
    detectable_auroc: float,
) -> str:
    if int(row.positive_sites) < minimum_positives:
        return "underpowered"
    if abs(float(row.after_smd_accessibility)) > maximum_accessibility_smd:
        return "accessibility_confounded"
    if abs(float(row.after_smd_motif_score)) > maximum_motif_smd:
        return "motif_score_confounded"
    if float(row.candidate_auroc) < information_limit_auroc:
        return "not_detected"
    if float(row.candidate_auroc) < weak_auroc:
        return "weak"
    if float(row.candidate_auroc) < detectable_auroc:
        return "detectable"
    return "strong"


def summarize(
    comparison: pd.DataFrame,
    matching: pd.DataFrame,
    study: dict[str, object],
    maximum_accessibility_smd: float = 0.25,
    maximum_motif_smd: float = 0.5,
) -> pd.DataFrame:
    matching = matching.rename(
        columns={
            "positive_sites": "matched_positive_sites",
            "negative_sites": "matched_negative_sites",
        }
    )
    joined = comparison.merge(
        matching,
        on=["cell", "tf", "chromosome_split"],
        how="left",
        validate="one_to_one",
    )
    thresholds = study["diagnostic_thresholds"]
    joined["detectability_status"] = joined.apply(
        classify_row,
        axis=1,
        minimum_positives=int(thresholds["minimum_positive_sites"]),
        maximum_accessibility_smd=maximum_accessibility_smd,
        maximum_motif_smd=maximum_motif_smd,
        information_limit_auroc=float(thresholds["information_limit_auroc"]),
        weak_auroc=float(thresholds["weak_auroc"]),
        detectable_auroc=float(thresholds["detectable_auroc"]),
    )
    gates = study["promotion_gates"]
    relative_auprc = joined["delta_auprc"] / joined["baseline_auprc"].clip(lower=1e-9)
    joined["passes_point_gain_gate"] = (
        (joined["delta_auroc"] >= float(gates["minimum_auroc_gain"]))
        & (relative_auprc >= float(gates["minimum_relative_auprc_gain"]))
    )
    joined["recommended_action"] = joined["detectability_status"].map(
        {
            "underpowered": "collect more orthogonal positive sites or aggregate replicates",
            "accessibility_confounded": "improve common-support matching before scorer evaluation",
            "motif_score_confounded": "improve motif-score common support or family mapping",
            "not_detected": "abstain; test depth/replicates and orthogonal occupancy",
            "weak": "retain as experimental and test depth/replicate stability",
            "detectable": "confirm in new cell lines before routing",
            "strong": "confirm non-regression and naked-DNA specificity",
        }
    )
    return joined.sort_values(["cell", "tf"]).reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--matching", type=Path, required=True)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-accessibility-smd", type=float, default=0.25)
    parser.add_argument("--max-motif-smd", type=float, default=0.5)
    args = parser.parse_args(argv)
    result = summarize(
        pd.read_csv(args.comparison, sep="\t"),
        pd.read_csv(args.matching, sep="\t"),
        json.loads(args.study.read_text(encoding="utf-8")),
        args.max_accessibility_smd,
        args.max_motif_smd,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(result[["cell", "tf", "positive_sites", "candidate_auroc", "delta_auroc", "detectability_status", "passes_point_gain_gate"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
