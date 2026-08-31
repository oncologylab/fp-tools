#!/usr/bin/env python3
"""Apply frozen per-TF dual-null thresholds to independent control scores."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


JOIN_COLUMNS = (
    "cell",
    "tf",
    "method",
    "candidate_id",
    "bias_configuration",
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _boolean(values: pd.Series, column: str) -> np.ndarray:
    if pd.api.types.is_bool_dtype(values):
        return values.to_numpy(dtype=bool)
    normalized = values.astype(str).str.strip().str.lower()
    invalid = ~normalized.isin({"true", "false", "1", "0"})
    if invalid.any():
        raise ValueError(f"{column} contains non-boolean values")
    return normalized.isin({"true", "1"}).to_numpy(dtype=bool)


def wilson_interval(
    successes: int,
    total: int,
    z: float = 1.959963984540054,
) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    radius = z * np.sqrt(
        proportion * (1.0 - proportion) / total
        + z * z / (4.0 * total * total)
    ) / denominator
    lower = 0.0 if successes == 0 else max(0.0, center - radius)
    return float(lower), float(min(1.0, center + radius))


def apply_thresholds(
    scores: pd.DataFrame,
    calibration: pd.DataFrame,
    *,
    method: str = "frozen_policy_candidate",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return site calls and per-TF rates without refitting any threshold."""

    score_required = set(JOIN_COLUMNS).union(
        {
            "motif_family",
            "replicate",
            "binding_probability",
            "total_signal",
            "valid",
        }
    )
    calibration_required = set(JOIN_COLUMNS).union({"dual_null_threshold"})
    missing_scores = score_required.difference(scores.columns)
    missing_calibration = calibration_required.difference(calibration.columns)
    if missing_scores:
        raise ValueError(
            "score table lacks columns: " + ", ".join(sorted(missing_scores))
        )
    if missing_calibration:
        raise ValueError(
            "calibration table lacks columns: "
            + ", ".join(sorted(missing_calibration))
        )

    selected_scores = scores[scores["method"].astype(str).eq(method)].copy()
    thresholds = calibration[
        calibration["method"].astype(str).eq(method)
    ][[*JOIN_COLUMNS, "dual_null_threshold"]].copy()
    if selected_scores.empty:
        raise ValueError(f"score table contains no {method} rows")
    if thresholds.empty:
        raise ValueError(f"calibration table contains no {method} rows")
    duplicated = thresholds.duplicated(list(JOIN_COLUMNS), keep=False)
    if duplicated.any():
        raise ValueError("calibration contains duplicate frozen thresholds")
    thresholds["dual_null_threshold"] = pd.to_numeric(
        thresholds["dual_null_threshold"], errors="raise"
    )
    if (
        ~np.isfinite(thresholds["dual_null_threshold"])
        | ~thresholds["dual_null_threshold"].between(0.0, 1.0, inclusive="neither")
    ).any():
        raise ValueError("dual-null thresholds must be finite and between zero and one")

    calls = selected_scores.merge(
        thresholds,
        on=list(JOIN_COLUMNS),
        how="left",
        validate="many_to_one",
    )
    if calls["dual_null_threshold"].isna().any():
        missing = calls.loc[
            calls["dual_null_threshold"].isna(), list(JOIN_COLUMNS)
        ].drop_duplicates()
        raise ValueError(
            "no frozen threshold for score routes: "
            + "; ".join(",".join(map(str, row)) for row in missing.itertuples(index=False, name=None))
        )

    probabilities = pd.to_numeric(
        calls["binding_probability"], errors="coerce"
    ).to_numpy(dtype=float)
    total_signal = pd.to_numeric(calls["total_signal"], errors="coerce").to_numpy(
        dtype=float
    )
    valid = _boolean(calls["valid"], "valid") & np.isfinite(probabilities)
    informative = valid & np.isfinite(total_signal) & (total_signal > 0)
    dual_null_call = informative & (
        probabilities >= calls["dual_null_threshold"].to_numpy(dtype=float)
    )
    calls["valid"] = valid
    calls["informative"] = informative
    calls["dual_null_call"] = dual_null_call

    group_columns = [
        "cell",
        "tf",
        "motif_family",
        "replicate",
        "method",
        "candidate_id",
        "bias_configuration",
        "dual_null_threshold",
    ]
    rows = []
    for keys, group in calls.groupby(group_columns, sort=True, dropna=False):
        n_valid = int(group["valid"].sum())
        n_informative = int(group["informative"].sum())
        n_calls = int(group["dual_null_call"].sum())
        all_lower, all_upper = wilson_interval(n_calls, n_valid)
        informative_lower, informative_upper = wilson_interval(
            n_calls, n_informative
        )
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "valid": n_valid,
                "informative": n_informative,
                "calls": n_calls,
                "all_site_rate": n_calls / n_valid if n_valid else np.nan,
                "all_site_rate_lower_95": all_lower,
                "all_site_rate_upper_95": all_upper,
                "informative_rate": (
                    n_calls / n_informative if n_informative else np.nan
                ),
                "informative_rate_lower_95": informative_lower,
                "informative_rate_upper_95": informative_upper,
                "mean_signal": float(
                    pd.to_numeric(group["total_signal"], errors="coerce").mean()
                ),
            }
        )
    return calls, pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--method", default="frozen_policy_candidate")
    parser.add_argument("--maximum-rate", type=float, default=0.05)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if not 0 <= args.maximum_rate <= 1:
        raise SystemExit("--maximum-rate must be between zero and one")

    scores = pd.read_csv(args.scores, sep="\t")
    calibration = pd.read_csv(args.calibration, sep="\t")
    calls, rates = apply_thresholds(scores, calibration, method=args.method)
    args.outdir.mkdir(parents=True, exist_ok=True)
    calls_path = args.outdir / "independent_dual_null_site_calls.tsv.gz"
    rates_path = args.outdir / "independent_dual_null_rates.tsv"
    calls.to_csv(calls_path, sep="\t", index=False)
    rates.to_csv(rates_path, sep="\t", index=False)

    maximum_all = float(rates["all_site_rate"].max())
    maximum_informative = float(rates["informative_rate"].max())
    maximum_all_upper = float(rates["all_site_rate_upper_95"].max())
    maximum_informative_upper = float(
        rates["informative_rate_upper_95"].max()
    )
    point_estimate_passes = bool(
        maximum_all <= args.maximum_rate
        and maximum_informative <= args.maximum_rate
    )
    gate = {
        "schema": "fp-tools-independent-dual-null-gate-v1",
        "method": args.method,
        "threshold_policy": "frozen_per_tf_dual_null",
        "maximum_allowed_rate": float(args.maximum_rate),
        "maximum_all_site_rate": maximum_all,
        "maximum_informative_rate": maximum_informative,
        "maximum_all_site_upper_95": maximum_all_upper,
        "maximum_informative_upper_95": maximum_informative_upper,
        "confidence_bound_passes": bool(
            maximum_all_upper <= args.maximum_rate
            and maximum_informative_upper <= args.maximum_rate
        ),
        "point_estimate_passes": point_estimate_passes,
        "passes": point_estimate_passes,
        "inputs": {
            "scores": {"path": str(args.scores), "sha256": file_sha256(args.scores)},
            "calibration": {
                "path": str(args.calibration),
                "sha256": file_sha256(args.calibration),
            },
        },
        "outputs": {
            calls_path.name: {"path": str(calls_path), "sha256": file_sha256(calls_path)},
            rates_path.name: {"path": str(rates_path), "sha256": file_sha256(rates_path)},
        },
    }
    manifest_path = args.outdir / "independent_dual_null_gate.json"
    manifest_path.write_text(
        json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(rates.to_string(index=False))
    print(json.dumps(gate, indent=2, sort_keys=True))
    return 0 if gate["passes"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
