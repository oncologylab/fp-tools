from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks" / "scripts"))

import evaluate_factorization_residual_safety as safety  # noqa: E402


def test_conformal_cutoff_is_deterministic_and_conservative_with_ties() -> None:
    scores = np.arange(200, dtype=float)
    threshold = safety.conformal_upper_threshold(scores, 0.05)
    assert threshold == 190.0
    assert np.mean(scores > threshold) == 0.045

    tied = np.concatenate([np.zeros(190), np.ones(10)])
    tied_threshold = safety.conformal_upper_threshold(tied, 0.05)
    assert tied_threshold == 1.0
    assert not np.any(tied > tied_threshold)


def test_naked_site_validation_fails_closed_on_chip_or_label_columns(
    tmp_path: Path,
) -> None:
    clean = pd.DataFrame({"cell": ["K562"], "tf": ["CTCF"]})
    safety.validate_label_free_sites(clean, tmp_path / "clean.tsv")
    for forbidden in ("chip_label", "bindingLabel"):
        sites = clean.assign(**{forbidden: [0]})
        with pytest.raises(ValueError, match="forbidden columns"):
            safety.validate_label_free_sites(sites, tmp_path / "unsafe.tsv")


def test_safety_summary_gates_overall_and_informative_rates() -> None:
    detail = pd.DataFrame(
        [
            {
                "residual": "deviance",
                "validation_negative_support": 200,
                "false_positive_calls": 2,
                "finite_support": 200,
                "informative_support": 100,
                "false_positive_rate": 0.01,
                "informative_false_positive_rate": 0.02,
                "false_positive_rate_upper_95": 0.035,
                "informative_false_positive_rate_upper_95": 0.07,
            },
            {
                "residual": "deviance",
                "validation_negative_support": 220,
                "false_positive_calls": 1,
                "finite_support": 200,
                "informative_support": 100,
                "false_positive_rate": 0.005,
                "informative_false_positive_rate": 0.01,
                "false_positive_rate_upper_95": 0.028,
                "informative_false_positive_rate_upper_95": 0.054,
            },
            {
                "residual": "pearson",
                "validation_negative_support": 200,
                "false_positive_calls": 6,
                "finite_support": 200,
                "informative_support": 100,
                "false_positive_rate": 0.03,
                "informative_false_positive_rate": 0.06,
                "false_positive_rate_upper_95": 0.064,
                "informative_false_positive_rate_upper_95": 0.125,
            },
        ]
    )
    summary = safety.summarize_safety(
        detail,
        maximum_false_positive_rate=0.05,
        minimum_validation_negatives=200,
    ).set_index("residual")
    assert bool(summary.loc["deviance", "passes_naked_dna_safety"])
    assert not bool(summary.loc["pearson", "passes_naked_dna_safety"])
    assert summary.loc["pearson", "maximum_false_positive_rate"] == 0.03
    assert (
        summary.loc["pearson", "maximum_informative_false_positive_rate"]
        == 0.06
    )


def test_wilson_interval_handles_empty_and_bounds_rate() -> None:
    low, high = safety.wilson_interval(3, 200)
    assert 0.0 <= low < 3 / 200 < high <= 1.0
    empty = safety.wilson_interval(0, 0)
    assert all(np.isnan(value) for value in empty)


def test_profile_scoring_excludes_nonfinite_bias_before_prediction() -> None:
    class FiniteOnlyModel:
        positions = np.arange(-100, 101, dtype=float)
        total_dispersion_ = 0.0

        def predict(self, counts, log_bias, samples, tfs):
            assert np.isfinite(log_bias).all()
            assert len(counts) == 1
            return SimpleNamespace(expected_unbound=np.ones_like(counts))

    arrays = {
        "plus_observed": np.ones((2, 201)),
        "minus_observed": np.ones((2, 201)),
        "combined_log_bias": np.asarray(
            [[0.0] * 201, [0.0, np.nan] + [0.0] * 199]
        ),
        "valid": np.ones(2, dtype=bool),
    }
    sites = pd.DataFrame({"cell": ["K562"] * 2, "tf": ["CTCF"] * 2})
    scores, valid, totals = safety._profile_scores(
        FiniteOnlyModel(), arrays, sites, ["deviance"]
    )
    assert valid.tolist() == [True, False]
    assert np.isfinite(scores["deviance"][0])
    assert np.isnan(scores["deviance"][1])
    assert totals.tolist() == [402.0, 402.0]
