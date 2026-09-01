from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import evaluate_frozen_functional_consensus as consensus  # noqa: E402


def test_empirical_negative_cdf_is_monotone_and_frozen_to_negatives() -> None:
    negative = np.asarray([0.0, 1.0, 2.0])
    values = np.asarray([-1.0, 0.0, 0.5, 2.0, 3.0, np.nan])
    observed = consensus.empirical_negative_cdf(values, negative)
    np.testing.assert_allclose(observed[:-1], [0.125, 0.375, 0.375, 0.875, 0.875])
    assert np.isnan(observed[-1])
    with pytest.raises(ValueError, match="sorted"):
        consensus.empirical_negative_cdf(values, negative[::-1])


def test_consensus_operators_have_predeclared_geometry() -> None:
    count = np.asarray([0.8])
    fda = np.asarray([0.2])
    operators = {record["operator"]: record for record in consensus.OPERATORS}
    assert consensus.consensus_score(count, fda, operators["rank_min"])[0] == 0.2
    assert np.isclose(
        consensus.consensus_score(count, fda, operators["rank_geometric_mean"])[0],
        0.4,
    )
    assert np.isclose(
        consensus.consensus_score(
            count,
            fda,
            operators["rank_mean_count_0p75"],
        )[0],
        0.65,
    )


def _validation_frame(score: list[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell": ["CellA", "CellA"],
            "tf": ["TF1", "TF1"],
            "motif_family": ["F1", "F1"],
            "bias_configuration": ["LOG21", "LOG21"],
            "candidate_id": ["candidate", "candidate"],
            "artifact_index": [7, 9],
            "site_hash": np.asarray([17, 91], dtype=np.uint64),
            "TFBS_chr": ["chr16", "chr17"],
            "label": [0, 1],
            "candidate_score": score,
            "dwm_score": [0.1, 0.5],
        }
    )


def test_validation_alignment_requires_identical_sites() -> None:
    result = consensus.align_validation_scores(
        _validation_frame([0.2, 0.8]),
        _validation_frame([0.3, 0.9]),
    )
    assert result["count_score"].tolist() == [0.2, 0.8]
    assert result["fda_score"].tolist() == [0.3, 0.9]
    changed = _validation_frame([0.3, 0.9])
    changed.loc[1, "site_hash"] = 123
    with pytest.raises(ValueError, match="site_hash"):
        consensus.align_validation_scores(_validation_frame([0.2, 0.8]), changed)


def test_primary_selection_uses_one_se_and_simplicity_order() -> None:
    summary = pd.DataFrame(
        {
            "operator": ["count_only", "rank_min", "fda_only"],
            "simplicity_rank": [0, 1, 6],
            "mean_selection_score_gain_over_dwm": [0.05, 0.08, 0.10],
            "selection_score_gain_standard_error": [0.01, 0.01, 0.03],
            "passes_ctcf_nonregression": [True, True, True],
        }
    )
    selected, cutoff, best = consensus.select_primary_operator(summary)
    assert best == "fda_only"
    assert np.isclose(cutoff, 0.07)
    assert selected == "rank_min"


def test_naked_dna_safety_retains_zero_cut_sites_in_denominator() -> None:
    rate, calls = consensus.safety_record(
        candidate_score=np.asarray([0.99, 0.99, 0.1, 0.1]),
        candidate_valid=np.ones(4, dtype=bool),
        candidate_informative=np.asarray([False, True, False, False]),
        candidate_threshold=0.9,
        dwm_score=np.asarray([0.1, 0.1, 0.1, 0.1]),
        dwm_valid=np.ones(4, dtype=bool),
        dwm_informative=np.asarray([False, True, False, False]),
        dwm_threshold=0.9,
    )
    assert rate["valid_sites"] == 4
    assert rate["informative_sites"] == 1
    assert rate["calls"] == 1
    assert rate["false_positive_rate"] == 0.25
    assert calls.tolist() == [False, True, False, False]


def test_raw_guardrail_alignment_uses_common_support() -> None:
    scored = pd.DataFrame(
        {
            "cell": ["CellA", "CellA"],
            "tf": ["TF1", "TF1"],
            "TFBS_chr": ["chr19", "chr20"],
            "TFBS_start": [10, 20],
            "TFBS_end": [15, 25],
            "TFBS_strand": ["+", "-"],
            "label": [0, 1],
            "dwm_score": [0.1, 0.5],
        }
    )
    raw = scored.rename(columns={"label": "chip_label"}).copy()
    raw["raw_score"] = [0.2, 0.7]
    aligned, coverage = consensus.align_raw_test_scores(scored, raw)
    assert aligned["raw_score"].tolist() == [0.2, 0.7]
    assert coverage.loc[0, "common_fraction"] == 1.0

    incomplete = raw.iloc[:1]
    with pytest.raises(ValueError, match="common support"):
        consensus.align_raw_test_scores(scored, incomplete)


def test_task_metric_record_reports_raw_guardrail_deltas() -> None:
    frame = pd.DataFrame(
        {
            "cell": ["CellA"] * 4,
            "tf": ["TF1"] * 4,
            "motif_family": ["F1"] * 4,
            "label": [0, 0, 1, 1],
            "dwm_score": [0.4, 0.3, 0.2, 0.1],
            "raw_score": [0.2, 0.1, 0.3, 0.4],
        }
    )
    record = consensus.task_metric_record(
        frame,
        np.asarray([0.1, 0.2, 0.8, 0.9]),
        operator="candidate",
        minimum_sites_per_class=2,
    )
    assert record["raw_auroc"] == 1.0
    assert record["auroc_gain_over_raw"] == 0.0
    assert record["relative_auprc_gain_over_raw"] == 0.0
