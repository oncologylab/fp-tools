#!/usr/bin/env python3
"""Screen frozen factorization residuals on independent naked-DNA profiles.

Residual cutoffs are calibrated once from matched ChIP-negative sites on the
development validation chromosomes.  The fixed cutoffs are then applied to a
label-free naked-DNA library.  Naked-DNA data never influence a cutoff, the
factorization fit, or the continuous residual scores.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import sqrt
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_parametric_factorization import (  # noqa: E402
    RESIDUALS,
    load_profiles,
    parse_name_path,
    residual_score,
    sha256_file,
)
from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenParametricFactorization,
)


SCHEMA = "fp-tools-parametric-factorization-residual-safety-v1"


def validate_label_free_sites(sites: pd.DataFrame, path: Path) -> None:
    """Fail closed if an independent-control site table contains labels."""

    forbidden = sorted(
        column
        for column in sites.columns
        if "label" in column.lower() or "chip" in column.lower()
    )
    if forbidden:
        raise ValueError(
            f"naked-DNA site table {path} contains forbidden columns: "
            + ", ".join(forbidden)
        )


def conformal_upper_threshold(scores: np.ndarray, false_positive_rate: float) -> float:
    """Return a deterministic upper-tail cutoff with finite-sample control.

    Calls use ``score > cutoff``.  Selecting the ceil((n + 1) * (1-alpha))-th
    order statistic makes ties conservative and avoids interpolation-dependent
    thresholds.
    """

    if not 0.0 < false_positive_rate < 1.0:
        raise ValueError("false_positive_rate must be between zero and one")
    values = np.sort(np.asarray(scores, dtype=np.float64))
    values = values[np.isfinite(values)]
    if len(values) == 0:
        raise ValueError("cannot calibrate a threshold without finite scores")
    rank = int(np.ceil((len(values) + 1) * (1.0 - false_positive_rate)))
    if rank > len(values):
        return float("inf")
    return float(values[max(rank - 1, 0)])


def wilson_interval(
    successes: int,
    total: int,
    z: float = 1.959963984540054,
) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    spread = (
        z
        * sqrt(
            proportion * (1.0 - proportion) / total
            + z * z / (4.0 * total * total)
        )
        / denominator
    )
    return max(0.0, center - spread), min(1.0, center + spread)


def _profile_scores(
    model: FrozenParametricFactorization,
    arrays: dict[str, np.ndarray],
    sites: pd.DataFrame,
    residuals: Iterable[str],
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    counts = np.asarray(
        arrays["plus_observed"] + arrays["minus_observed"], dtype=np.float64
    )
    log_bias = np.asarray(arrays["combined_log_bias"], dtype=np.float64)
    if counts.shape[1] != len(model.positions):
        raise ValueError("profile width does not match the frozen factorization")
    if "cell" not in sites or "tf" not in sites:
        raise ValueError("profile sites must contain cell and tf columns")
    result = model.predict(
        counts,
        log_bias,
        sites["cell"].astype(str),
        sites["tf"].astype(str),
    )
    scores = {
        residual: residual_score(
            counts,
            result.expected_unbound,
            model.positions,
            residual,
            model.total_dispersion_,
        )[0]
        for residual in residuals
    }
    valid = (
        np.asarray(arrays["valid"], dtype=bool)
        & np.isfinite(log_bias).all(axis=1)
        & np.isfinite(result.expected_unbound).all(axis=1)
    )
    return scores, valid, counts.sum(axis=1)


def calibrate_thresholds(
    model: FrozenParametricFactorization,
    validation_artifacts: dict[str, Path],
    *,
    residuals: Iterable[str] = RESIDUALS,
    target_false_positive_rate: float = 0.05,
) -> tuple[pd.DataFrame, list[Path]]:
    rows: list[dict[str, object]] = []
    inputs: list[Path] = []
    for cell, prefix in sorted(validation_artifacts.items()):
        arrays, sites, _document = load_profiles(prefix, require_log_bias=True)
        inputs.extend(Path(str(prefix) + suffix) for suffix in (".npz", ".json", ".sites.tsv.gz"))
        if set(sites["cell"].astype(str)) != {cell}:
            raise ValueError(f"validation artifact is not exclusive to {cell}: {prefix}")
        required = {"chip_label", "chromosome_split", "tf"}
        missing = required.difference(sites.columns)
        if missing:
            raise ValueError(
                "validation site table is missing columns: "
                + ", ".join(sorted(missing))
            )
        scores, valid, totals = _profile_scores(model, arrays, sites, residuals)
        validation = sites["chromosome_split"].astype(str).eq("validation").to_numpy()
        negatives = sites["chip_label"].to_numpy(dtype=int) == 0
        for tf in sorted(sites.loc[validation, "tf"].astype(str).unique()):
            task = sites["tf"].astype(str).eq(tf).to_numpy()
            support = valid & validation & negatives & task & (totals > 0)
            for residual, values in scores.items():
                finite = support & np.isfinite(values)
                rows.append(
                    {
                        "cell": cell,
                        "tf": tf,
                        "residual": residual,
                        "validation_negative_support": int(np.sum(finite)),
                        "target_false_positive_rate": float(
                            target_false_positive_rate
                        ),
                        "score_threshold": conformal_upper_threshold(
                            values[finite], target_false_positive_rate
                        )
                        if np.any(finite)
                        else np.nan,
                    }
                )
    thresholds = pd.DataFrame(rows)
    if thresholds.duplicated(["cell", "tf", "residual"]).any():
        raise ValueError("validation artifacts contain duplicate cell/TF tasks")
    return thresholds, inputs


def apply_naked_dna_safety(
    model: FrozenParametricFactorization,
    naked_artifacts: dict[str, Path],
    thresholds: pd.DataFrame,
    *,
    residuals: Iterable[str] = RESIDUALS,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    detail_rows: list[dict[str, object]] = []
    score_frames: list[pd.DataFrame] = []
    inputs: list[Path] = []
    lookup = thresholds.set_index(["cell", "tf", "residual"])
    for cell, prefix in sorted(naked_artifacts.items()):
        arrays, sites, _document = load_profiles(prefix, require_log_bias=True)
        paths = [Path(str(prefix) + suffix) for suffix in (".npz", ".json", ".sites.tsv.gz")]
        inputs.extend(paths)
        validate_label_free_sites(sites, paths[-1])
        if set(sites["cell"].astype(str)) != {cell}:
            raise ValueError(f"naked-DNA artifact is not exclusive to {cell}: {prefix}")
        scores, valid, totals = _profile_scores(model, arrays, sites, residuals)
        hashes = np.asarray(arrays["site_hash"], dtype=np.uint64)
        for tf in sorted(sites["tf"].astype(str).unique()):
            task = sites["tf"].astype(str).eq(tf).to_numpy()
            for residual, values in scores.items():
                key = (cell, tf, residual)
                if key not in lookup.index:
                    raise ValueError(f"no frozen validation threshold for {key}")
                threshold_row = lookup.loc[key]
                threshold = float(threshold_row["score_threshold"])
                finite = task & valid & np.isfinite(values)
                informative = finite & (totals > 0)
                calls = informative & (values > threshold)
                n_finite = int(np.sum(finite))
                n_informative = int(np.sum(informative))
                n_calls = int(np.sum(calls))
                low, high = wilson_interval(n_calls, n_finite)
                informative_low, informative_high = wilson_interval(
                    n_calls, n_informative
                )
                detail_rows.append(
                    {
                        "cell": cell,
                        "tf": tf,
                        "residual": residual,
                        "validation_negative_support": int(
                            threshold_row["validation_negative_support"]
                        ),
                        "score_threshold": threshold,
                        "naked_sites": int(np.sum(task)),
                        "finite_support": n_finite,
                        "informative_support": n_informative,
                        "false_positive_calls": n_calls,
                        "false_positive_rate": n_calls / n_finite
                        if n_finite
                        else np.nan,
                        "false_positive_rate_lower_95": low,
                        "false_positive_rate_upper_95": high,
                        "informative_false_positive_rate": n_calls / n_informative
                        if n_informative
                        else np.nan,
                        "informative_false_positive_rate_lower_95": informative_low,
                        "informative_false_positive_rate_upper_95": informative_high,
                    }
                )
                selected = np.flatnonzero(task)
                score_frames.append(
                    pd.DataFrame(
                        {
                            "cell": cell,
                            "tf": tf,
                            "residual": residual,
                            "site_hash": hashes[selected],
                            "total_cuts": totals[selected],
                            "score": values[selected],
                            "score_threshold": threshold,
                            "finite": finite[selected],
                            "informative": informative[selected],
                            "false_positive_call": calls[selected],
                        }
                    )
                )
    return pd.DataFrame(detail_rows), pd.concat(score_frames, ignore_index=True), inputs


def summarize_safety(
    detail: pd.DataFrame,
    *,
    maximum_false_positive_rate: float,
    minimum_validation_negatives: int,
) -> pd.DataFrame:
    rows = []
    for residual, group in detail.groupby("residual", sort=True):
        finite_rates = group["false_positive_rate"].dropna()
        informative_rates = group["informative_false_positive_rate"].dropna()
        complete = len(finite_rates) == len(group) and len(informative_rates) == len(group)
        minimum_support = int(group["validation_negative_support"].min())
        maximum_rate = float(finite_rates.max()) if len(finite_rates) else np.nan
        maximum_informative = (
            float(informative_rates.max()) if len(informative_rates) else np.nan
        )
        rows.append(
            {
                "residual": residual,
                "task_count": int(len(group)),
                "minimum_validation_negative_support": minimum_support,
                "pooled_false_positive_calls": int(group["false_positive_calls"].sum()),
                "pooled_finite_support": int(group["finite_support"].sum()),
                "pooled_informative_support": int(group["informative_support"].sum()),
                "maximum_false_positive_rate": maximum_rate,
                "maximum_informative_false_positive_rate": maximum_informative,
                "maximum_false_positive_rate_upper_95": float(
                    group["false_positive_rate_upper_95"].max()
                ),
                "maximum_informative_false_positive_rate_upper_95": float(
                    group["informative_false_positive_rate_upper_95"].max()
                ),
                "complete_support": bool(complete),
                "passes_minimum_validation_support": bool(
                    minimum_support >= minimum_validation_negatives
                ),
                "passes_naked_dna_safety": bool(
                    complete
                    and minimum_support >= minimum_validation_negatives
                    and maximum_rate <= maximum_false_positive_rate
                    and maximum_informative <= maximum_false_positive_rate
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("residual", kind="mergesort")


def _artifact_record(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": sha256_file(path)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factorization-model", type=Path, required=True)
    parser.add_argument(
        "--validation-candidate", type=parse_name_path, action="append", required=True
    )
    parser.add_argument(
        "--naked-candidate", type=parse_name_path, action="append", required=True
    )
    parser.add_argument("--target-false-positive-rate", type=float, default=0.05)
    parser.add_argument("--maximum-false-positive-rate", type=float, default=0.05)
    parser.add_argument("--minimum-validation-negatives", type=int, default=200)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    validation = dict(args.validation_candidate)
    naked = dict(args.naked_candidate)
    if set(validation) != set(naked):
        raise ValueError("validation and naked-DNA cells must match")
    if len(validation) != len(args.validation_candidate):
        raise ValueError("duplicate validation cell")
    if len(naked) != len(args.naked_candidate):
        raise ValueError("duplicate naked-DNA cell")
    model = FrozenParametricFactorization.load(args.factorization_model)
    thresholds, validation_inputs = calibrate_thresholds(
        model,
        validation,
        target_false_positive_rate=args.target_false_positive_rate,
    )
    detail, site_scores, naked_inputs = apply_naked_dna_safety(
        model,
        naked,
        thresholds,
    )
    summary = summarize_safety(
        detail,
        maximum_false_positive_rate=args.maximum_false_positive_rate,
        minimum_validation_negatives=args.minimum_validation_negatives,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    threshold_path = args.outdir / "factorization_residual_thresholds.tsv"
    detail_path = args.outdir / "factorization_naked_dna_safety.tsv"
    summary_path = args.outdir / "factorization_residual_safety_summary.tsv"
    scores_path = args.outdir / "factorization_naked_dna_site_scores.tsv.gz"
    thresholds.to_csv(threshold_path, sep="\t", index=False)
    detail.to_csv(detail_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    site_scores.to_csv(scores_path, sep="\t", index=False)
    metadata_path = args.factorization_model.with_suffix(".json")
    document = {
        "schema": SCHEMA,
        "factorization_model": _artifact_record(args.factorization_model),
        "factorization_model_metadata": _artifact_record(metadata_path),
        "validation_inputs": [_artifact_record(path) for path in validation_inputs],
        "naked_dna_inputs": [_artifact_record(path) for path in naked_inputs],
        "target_false_positive_rate": float(args.target_false_positive_rate),
        "maximum_false_positive_rate": float(args.maximum_false_positive_rate),
        "minimum_validation_negatives": int(args.minimum_validation_negatives),
        "thresholds": _artifact_record(threshold_path),
        "detail": _artifact_record(detail_path),
        "summary": _artifact_record(summary_path),
        "site_scores": _artifact_record(scores_path),
        "labels_used_for_threshold_calibration": True,
        "naked_dna_labels_used": False,
        "passing_residuals": summary.loc[
            summary["passes_naked_dna_safety"], "residual"
        ].astype(str).tolist(),
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["safety_id"] = sha256(canonical.encode()).hexdigest()
    (args.outdir / "factorization_residual_safety.json").write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
