#!/usr/bin/env python3
"""Calibrate footprint calls against an independent naked-DNA null panel.

Posterior probabilities from unsupervised two-state models are useful ranks,
but 0.5 is not a valid universal calling threshold.  This helper chooses a
conservative per-TF threshold on one label-free naked-DNA panel and applies it
unchanged to a second panel.  It never reads occupancy labels.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


GROUP_COLUMNS = [
    "cell",
    "tf",
    "motif_family",
    "method",
    "candidate_id",
    "bias_configuration",
]


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def conservative_threshold(values: np.ndarray, alpha: float) -> float:
    """Return the most permissive observed threshold with empirical FPR <= alpha."""

    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    allowed = int(np.floor(alpha * len(finite)))
    if not len(finite) or allowed == 0:
        return float("inf")
    ordered = np.sort(finite)[::-1]
    threshold = float(ordered[allowed - 1])
    if int(np.sum(finite >= threshold)) > allowed:
        threshold = float(np.nextafter(threshold, np.inf))
    return threshold


def _validate_scores(frame: pd.DataFrame, source: Path) -> None:
    required = set(GROUP_COLUMNS) | {
        "site_hash",
        "binding_probability",
        "valid",
        "informative",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{source} is missing columns: " + ", ".join(sorted(missing)))
    leaked = [
        column
        for column in frame.columns
        if "label" in column.lower() or "chip" in column.lower()
    ]
    if leaked:
        raise ValueError(f"control score table contains label columns: {', '.join(leaked)}")


def calibrate(
    calibration: pd.DataFrame,
    validation: pd.DataFrame,
    *,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    thresholds: list[dict] = []
    validation_parts: list[pd.DataFrame] = []
    calibration_groups = calibration.groupby(GROUP_COLUMNS, sort=True, dropna=False)
    validation_groups = validation.groupby(GROUP_COLUMNS, sort=True, dropna=False)
    calibration_keys = set(calibration_groups.groups)
    validation_keys = set(validation_groups.groups)
    if calibration_keys != validation_keys:
        missing_calibration = validation_keys.difference(calibration_keys)
        missing_validation = calibration_keys.difference(validation_keys)
        raise ValueError(
            "calibration/validation groups differ; "
            f"missing calibration={len(missing_calibration)}, "
            f"missing validation={len(missing_validation)}"
        )
    for key in sorted(calibration_keys):
        calibration_group = calibration_groups.get_group(key).copy()
        validation_group = validation_groups.get_group(key).copy()
        calibration_informative = (
            calibration_group["valid"].astype(bool)
            & calibration_group["informative"].astype(bool)
            & np.isfinite(calibration_group["binding_probability"].to_numpy(float))
        )
        threshold = conservative_threshold(
            calibration_group.loc[calibration_informative, "binding_probability"].to_numpy(
                float
            ),
            alpha,
        )
        validation_valid = (
            validation_group["valid"].astype(bool)
            & np.isfinite(validation_group["binding_probability"].to_numpy(float))
        )
        validation_informative = (
            validation_valid & validation_group["informative"].astype(bool)
        )
        calibration_calls = calibration_informative & (
            calibration_group["binding_probability"].to_numpy(float) >= threshold
        )
        validation_calls = validation_informative & (
            validation_group["binding_probability"].to_numpy(float) >= threshold
        )
        validation_group["null_calibrated_threshold"] = threshold
        validation_group["null_calibrated_call"] = validation_calls
        validation_parts.append(validation_group)
        row = dict(zip(GROUP_COLUMNS, key))
        row.update(
            {
                "alpha": float(alpha),
                "threshold": threshold,
                "calibration_valid_sites": int(calibration_group["valid"].astype(bool).sum()),
                "calibration_informative_sites": int(calibration_informative.sum()),
                "calibration_calls": int(calibration_calls.sum()),
                "calibration_informative_fpr": (
                    float(calibration_calls.sum() / calibration_informative.sum())
                    if calibration_informative.any()
                    else np.nan
                ),
                "validation_valid_sites": int(validation_valid.sum()),
                "validation_informative_sites": int(validation_informative.sum()),
                "validation_calls": int(validation_calls.sum()),
                "validation_all_site_fpr": (
                    float(validation_calls.sum() / validation_valid.sum())
                    if validation_valid.any()
                    else np.nan
                ),
                "validation_informative_fpr": (
                    float(validation_calls.sum() / validation_informative.sum())
                    if validation_informative.any()
                    else np.nan
                ),
            }
        )
        thresholds.append(row)
    summary = pd.DataFrame(thresholds).sort_values(GROUP_COLUMNS, kind="mergesort")
    site_calls = pd.concat(validation_parts, ignore_index=True).sort_values(
        GROUP_COLUMNS + ["TFBS_chr", "TFBS_start"], kind="mergesort"
    )
    return summary.reset_index(drop=True), site_calls.reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration-scores", type=Path, required=True)
    parser.add_argument("--validation-scores", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument(
        "--validation-fpr-limit",
        type=float,
        default=0.05,
        help="Independent-panel false-positive ceiling (default: 0.05).",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if not 0 < args.validation_fpr_limit < 1:
        raise SystemExit("--validation-fpr-limit must be between zero and one")
    calibration_frame = pd.read_csv(args.calibration_scores, sep="\t")
    validation_frame = pd.read_csv(args.validation_scores, sep="\t")
    _validate_scores(calibration_frame, args.calibration_scores)
    _validate_scores(validation_frame, args.validation_scores)
    summary, site_calls = calibrate(
        calibration_frame, validation_frame, alpha=args.alpha
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "naked_dna_null_calibration.tsv"
    calls_path = args.outdir / "naked_dna_null_calibrated_calls.tsv.gz"
    summary.to_csv(summary_path, sep="\t", index=False)
    site_calls.to_csv(calls_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-naked-dna-null-calibration-v1",
        "labels_used": False,
        "alpha": float(args.alpha),
        "validation_fpr_limit": float(args.validation_fpr_limit),
        "calibration_scores": str(args.calibration_scores),
        "calibration_scores_sha256": file_sha256(args.calibration_scores),
        "validation_scores": str(args.validation_scores),
        "validation_scores_sha256": file_sha256(args.validation_scores),
        "groups": int(len(summary)),
        "all_validation_groups_pass": bool(
            (
                summary["validation_all_site_fpr"]
                <= args.validation_fpr_limit
            ).fillna(False).all()
        ),
        "all_validation_informative_groups_pass": bool(
            (
                summary["validation_informative_fpr"]
                <= args.validation_fpr_limit
            ).fillna(False).all()
        ),
        "outputs": {
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "site_calls": {"path": str(calls_path), "sha256": file_sha256(calls_path)},
        },
    }
    (args.outdir / "naked_dna_null_calibration_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
