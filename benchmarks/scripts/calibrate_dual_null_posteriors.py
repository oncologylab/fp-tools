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


def empirical_upper_tail_pvalues(
    null_values: np.ndarray, observed_values: np.ndarray
) -> np.ndarray:
    """Return finite-sample-corrected upper-tail empirical p-values."""

    null = np.asarray(null_values, dtype=float)
    null = np.sort(null[np.isfinite(null)])
    observed = np.asarray(observed_values, dtype=float)
    output = np.full(observed.shape, np.nan, dtype=float)
    finite = np.isfinite(observed)
    if not len(null):
        return output
    counts = len(null) - np.searchsorted(null, observed[finite], side="left")
    output[finite] = (counts + 1.0) / (len(null) + 1.0)
    return output


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Adjust one family of p-values while retaining NaN positions."""

    values = np.asarray(pvalues, dtype=float)
    output = np.full(values.shape, np.nan, dtype=float)
    finite_indexes = np.flatnonzero(np.isfinite(values))
    if not len(finite_indexes):
        return output
    finite = values[finite_indexes]
    order = np.argsort(finite, kind="mergesort")
    ranked = finite[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    output[finite_indexes] = restored
    return output


def calibrate_dual_null(
    primary: pd.DataFrame,
    secondary: pd.DataFrame,
    validation: pd.DataFrame,
    *,
    primary_alpha: float,
    secondary_alpha: float,
    fdr: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit per-group thresholds and apply them to an independent null panel."""

    if not 0 < fdr < 1:
        raise ValueError("fdr must be between zero and one")
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
        validation_probabilities = validation_group["binding_probability"].to_numpy(
            float
        )
        primary_null = primary_group.loc[
            primary_mask, "binding_probability"
        ].to_numpy(float)
        secondary_null = secondary_group.loc[
            secondary_mask, "binding_probability"
        ].to_numpy(float)
        primary_pvalue = np.full(len(validation_group), np.nan, dtype=float)
        secondary_pvalue = np.full(len(validation_group), np.nan, dtype=float)
        primary_pvalue[validation_mask] = empirical_upper_tail_pvalues(
            primary_null, validation_probabilities[validation_mask]
        )
        secondary_pvalue[validation_mask] = empirical_upper_tail_pvalues(
            secondary_null, validation_probabilities[validation_mask]
        )
        dual_pvalue = np.fmax(primary_pvalue, secondary_pvalue)
        dual_qvalue = benjamini_hochberg(dual_pvalue)
        fdr_calls = validation_mask & (dual_qvalue <= fdr)
        validation_valid = (
            validation_group["valid"].astype(bool).to_numpy()
            & np.isfinite(validation_group["binding_probability"].to_numpy(float))
        )
        validation_group["dual_null_threshold"] = threshold
        validation_group["dual_null_call"] = validation_calls
        validation_group["naked_dna_pvalue"] = primary_pvalue
        validation_group["shifted_atac_pvalue"] = secondary_pvalue
        validation_group["dual_null_pvalue"] = dual_pvalue
        validation_group["dual_null_qvalue"] = dual_qvalue
        validation_group["dual_null_fdr_call"] = fdr_calls
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
                "fdr": float(fdr),
                "validation_fdr_calls": int(fdr_calls.sum()),
                "validation_fdr_all_site_rate": (
                    float(fdr_calls.sum() / validation_valid.sum())
                    if validation_valid.any()
                    else np.nan
                ),
                "validation_fdr_informative_rate": (
                    float(fdr_calls.sum() / validation_mask.sum())
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
    parser.add_argument("--fdr", type=float, default=0.05)
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
        fdr=args.fdr,
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
        "fdr": float(args.fdr),
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
        "all_validation_fdr_groups_pass": bool(
            summary["validation_fdr_all_site_rate"]
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
