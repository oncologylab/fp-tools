#!/usr/bin/env python3
"""Calibrate footprint calls against enzyme and cellular-ATAC nulls.

The primary null is naked DNA.  The secondary null contains motif-misaligned
cellular ATAC profiles and therefore captures accessibility and nucleosomal
structure absent from naked DNA.  A call must exceed the conservative
threshold learned from each source.  Occupancy labels are forbidden.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from calibrate_naked_dna_posteriors import (  # noqa: E402
    GROUP_COLUMNS,
    _validate_scores,
    conservative_threshold,
    file_sha256,
)


def _informative(frame: pd.DataFrame) -> np.ndarray:
    return (
        frame["valid"].astype(bool).to_numpy()
        & frame["informative"].astype(bool).to_numpy()
        & np.isfinite(frame["binding_probability"].to_numpy(dtype=float))
    )


def calibrate_dual_null(
    primary: pd.DataFrame,
    secondary: pd.DataFrame,
    validation: pd.DataFrame,
    *,
    primary_alpha: float,
    secondary_alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit per-group thresholds and apply them to an independent null panel."""

    grouped = [
        frame.groupby(GROUP_COLUMNS, sort=True, dropna=False)
        for frame in (primary, secondary, validation)
    ]
    key_sets = [set(value.groups) for value in grouped]
    if key_sets[0] != key_sets[1] or key_sets[0] != key_sets[2]:
        raise ValueError("primary, secondary, and validation groups differ")
    rows: list[dict] = []
    calls: list[pd.DataFrame] = []
    for key in sorted(key_sets[0]):
        primary_group, secondary_group, validation_group = [
            value.get_group(key).copy() for value in grouped
        ]
        primary_mask = _informative(primary_group)
        secondary_mask = _informative(secondary_group)
        validation_mask = _informative(validation_group)
        primary_threshold = conservative_threshold(
            primary_group.loc[primary_mask, "binding_probability"].to_numpy(float),
            primary_alpha,
        )
        secondary_threshold = conservative_threshold(
            secondary_group.loc[secondary_mask, "binding_probability"].to_numpy(float),
            secondary_alpha,
        )
        threshold = max(primary_threshold, secondary_threshold)
        primary_calls = primary_mask & (
            primary_group["binding_probability"].to_numpy(float) >= threshold
        )
        secondary_calls = secondary_mask & (
            secondary_group["binding_probability"].to_numpy(float) >= threshold
        )
        validation_calls = validation_mask & (
            validation_group["binding_probability"].to_numpy(float) >= threshold
        )
        validation_valid = (
            validation_group["valid"].astype(bool).to_numpy()
            & np.isfinite(validation_group["binding_probability"].to_numpy(float))
        )
        validation_group["dual_null_threshold"] = threshold
        validation_group["dual_null_call"] = validation_calls
        calls.append(validation_group)
        row = dict(zip(GROUP_COLUMNS, key))
        row.update(
            {
                "primary_alpha": float(primary_alpha),
                "secondary_alpha": float(secondary_alpha),
                "primary_threshold": primary_threshold,
                "secondary_threshold": secondary_threshold,
                "dual_null_threshold": threshold,
                "threshold_source": (
                    "primary_naked_dna"
                    if primary_threshold >= secondary_threshold
                    else "secondary_shifted_atac"
                ),
                "primary_informative_sites": int(primary_mask.sum()),
                "primary_calls": int(primary_calls.sum()),
                "primary_informative_rate": (
                    float(primary_calls.sum() / primary_mask.sum())
                    if primary_mask.any()
                    else np.nan
                ),
                "secondary_informative_sites": int(secondary_mask.sum()),
                "secondary_calls": int(secondary_calls.sum()),
                "secondary_informative_rate": (
                    float(secondary_calls.sum() / secondary_mask.sum())
                    if secondary_mask.any()
                    else np.nan
                ),
                "validation_valid_sites": int(validation_valid.sum()),
                "validation_informative_sites": int(validation_mask.sum()),
                "validation_calls": int(validation_calls.sum()),
                "validation_all_site_rate": (
                    float(validation_calls.sum() / validation_valid.sum())
                    if validation_valid.any()
                    else np.nan
                ),
                "validation_informative_rate": (
                    float(validation_calls.sum() / validation_mask.sum())
                    if validation_mask.any()
                    else np.nan
                ),
            }
        )
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values(GROUP_COLUMNS, kind="mergesort")
    site_calls = pd.concat(calls, ignore_index=True).sort_values(
        GROUP_COLUMNS + ["TFBS_chr", "TFBS_start"], kind="mergesort"
    )
    return summary.reset_index(drop=True), site_calls.reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--naked-calibration-scores", type=Path, required=True)
    parser.add_argument("--shifted-atac-scores", type=Path, required=True)
    parser.add_argument("--naked-validation-scores", type=Path, required=True)
    parser.add_argument("--naked-alpha", type=float, default=0.025)
    parser.add_argument("--shifted-atac-alpha", type=float, default=0.025)
    parser.add_argument("--validation-fpr-limit", type=float, default=0.05)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    paths = (
        args.naked_calibration_scores,
        args.shifted_atac_scores,
        args.naked_validation_scores,
    )
    frames = [pd.read_csv(path, sep="\t") for path in paths]
    for frame, path in zip(frames, paths):
        _validate_scores(frame, path)
    summary, calls = calibrate_dual_null(
        *frames,
        primary_alpha=args.naked_alpha,
        secondary_alpha=args.shifted_atac_alpha,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "dual_null_calibration.tsv"
    calls_path = args.outdir / "dual_null_naked_validation_calls.tsv.gz"
    summary.to_csv(summary_path, sep="\t", index=False)
    calls.to_csv(calls_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-dual-null-calibration-v1",
        "labels_used": False,
        "naked_alpha": float(args.naked_alpha),
        "shifted_atac_alpha": float(args.shifted_atac_alpha),
        "validation_fpr_limit": float(args.validation_fpr_limit),
        "groups": int(len(summary)),
        "all_validation_groups_pass": bool(
            summary["validation_all_site_rate"]
            .le(args.validation_fpr_limit)
            .fillna(False)
            .all()
        ),
        "all_validation_informative_groups_pass": bool(
            summary["validation_informative_rate"]
            .le(args.validation_fpr_limit)
            .fillna(False)
            .all()
        ),
        "inputs": [
            {"path": str(path), "sha256": file_sha256(path)} for path in paths
        ],
        "outputs": {
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "calls": {"path": str(calls_path), "sha256": file_sha256(calls_path)},
        },
    }
    (args.outdir / "dual_null_calibration_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
