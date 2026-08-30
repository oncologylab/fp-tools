#!/usr/bin/env python3
"""Prioritize nutrient-stress footprint changes through replication tiers.

Input is a canonical long TSV. ``cohort`` is one of ``local``,
``external_pdac``, or ``external_mechanistic``. ``contrast_class`` is
``stress``, ``recovery``, or ``occupancy``. Positive delta values must always
mean greater activity/occupancy in the stressed state; recovery rows are
therefore expected to have the opposite sign.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
    "cohort",
    "cell",
    "contrast_class",
    "motif_id",
    "tf",
    "delta_footprint",
    "fdr",
]


def as_numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    frame = frame.copy()
    for column in columns:
        if column not in frame:
            frame[column] = np.nan
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def validate_evidence(frame: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in frame]
    if missing:
        raise ValueError(f"nutrient evidence is missing columns: {', '.join(missing)}")
    allowed_cohorts = {"local", "external_pdac", "external_mechanistic"}
    unknown = set(frame["cohort"].dropna()) - allowed_cohorts
    if unknown:
        raise ValueError(f"unknown nutrient cohorts: {', '.join(sorted(unknown))}")
    allowed_contrasts = {"stress", "recovery", "occupancy"}
    unknown = set(frame["contrast_class"].dropna()) - allowed_contrasts
    if unknown:
        raise ValueError(f"unknown contrast classes: {', '.join(sorted(unknown))}")


def sign_or_zero(value: float) -> int:
    return int(np.sign(value)) if np.isfinite(value) else 0


def evaluate_replication(frame: pd.DataFrame, rules: dict) -> pd.DataFrame:
    validate_evidence(frame)
    frame = as_numeric(
        frame,
        ["delta_footprint", "fdr", "rna_log2fc", "rna_fdr", "occupancy_log2fc", "occupancy_fdr"],
    )
    local = frame[(frame["cohort"] == "local") & (frame["contrast_class"] == "stress")].copy()
    if local.empty:
        raise ValueError("nutrient evidence contains no local stress rows")
    local["absolute_delta_footprint"] = local["delta_footprint"].abs()
    local["absolute_rank_fraction"] = local.groupby("cell")["absolute_delta_footprint"].rank(
        method="min", ascending=False, pct=True
    )
    fdr_threshold = float(rules["external_fdr"])
    rows: list[dict[str, object]] = []
    for (motif_id, tf), local_group in local.groupby(["motif_id", "tf"], sort=True):
        finite = local_group.dropna(subset=["delta_footprint"])
        direction = sign_or_zero(float(finite["delta_footprint"].median())) if len(finite) else 0
        directional = finite[
            (np.sign(finite["delta_footprint"]) == direction)
            & (finite["fdr"] <= fdr_threshold)
        ]
        top = directional[
            directional["absolute_rank_fraction"] <= float(rules["top_fraction_absolute_change"])
        ]
        rna = directional[
            (np.sign(directional["rna_log2fc"]) == direction)
            & (directional["rna_fdr"] <= fdr_threshold)
        ]
        local_pass = bool(
            direction
            and directional["cell"].nunique() >= int(rules["required_directional_cell_lines"])
            and top["cell"].nunique() >= int(rules["minimum_top_cell_lines"])
            and rna["cell"].nunique() >= int(rules["minimum_rna_concordant_cell_lines"])
        )

        external = frame[
            (frame["cohort"] == "external_pdac")
            & (frame["motif_id"] == motif_id)
        ]
        stress = external[external["contrast_class"] == "stress"]
        recovery = external[external["contrast_class"] == "recovery"]
        stress_agree = stress[
            (np.sign(stress["delta_footprint"]) == direction)
            & (stress["fdr"] <= fdr_threshold)
        ]
        recovery_agree = recovery[
            (np.sign(recovery["delta_footprint"]) == -direction)
            & (recovery["fdr"] <= fdr_threshold)
        ]
        recovery_fraction = len(recovery_agree) / len(recovery) if len(recovery) else 0.0
        external_pass = bool(
            local_pass
            and len(stress_agree)
            and recovery_fraction >= float(rules["minimum_external_reversal_fraction"])
        )

        mechanistic = frame[
            (frame["cohort"] == "external_mechanistic")
            & (frame["motif_id"] == motif_id)
            & (frame["contrast_class"] == "occupancy")
        ]
        occupancy_agree = mechanistic[
            (np.sign(mechanistic["occupancy_log2fc"]) == direction)
            & (mechanistic["occupancy_fdr"] <= fdr_threshold)
        ]
        mechanism_pass = bool(external_pass and len(occupancy_agree))
        tier = (
            "mechanism_supported"
            if mechanism_pass
            else "external_replicated"
            if external_pass
            else "local_reproducible"
            if local_pass
            else "not_prioritized"
        )
        rows.append(
            {
                "motif_id": motif_id,
                "tf": tf,
                "direction": direction,
                "directional_local_cells": int(directional["cell"].nunique()),
                "top_local_cells": int(top["cell"].nunique()),
                "rna_concordant_local_cells": int(rna["cell"].nunique()),
                "external_stress_rows": int(len(stress)),
                "external_stress_agree": int(len(stress_agree)),
                "external_recovery_rows": int(len(recovery)),
                "external_recovery_agree": int(len(recovery_agree)),
                "external_recovery_fraction": float(recovery_fraction),
                "mechanistic_occupancy_rows": int(len(mechanistic)),
                "mechanistic_occupancy_agree": int(len(occupancy_agree)),
                "local_pass": local_pass,
                "external_pass": external_pass,
                "mechanism_pass": mechanism_pass,
                "evidence_tier": tier,
            }
        )
    order = {
        "mechanism_supported": 0,
        "external_replicated": 1,
        "local_reproducible": 2,
        "not_prioritized": 3,
    }
    result = pd.DataFrame(rows)
    result["_order"] = result["evidence_tier"].map(order)
    return result.sort_values(["_order", "rna_concordant_local_cells", "tf"], ascending=[True, False, True]).drop(columns="_order").reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--study", type=Path, default=Path("benchmarks/manifests/footprint_detectability_v1.spec.json"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    study = json.loads(args.study.read_text(encoding="utf-8"))
    result = evaluate_replication(
        pd.read_csv(args.evidence, sep="\t"),
        study["nutrient_application"]["candidate_rules"],
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    print(result.to_string(index=False))
    print(f"\nwrote {len(result)} nutrient candidates to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
