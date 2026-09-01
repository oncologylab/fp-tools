from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fp_tools.tools.frozen_bias_evaluation import (
    TobiasDwmReferenceModel,
    conditional_control_scores,
    motif_residual_effect,
    negative_control_safety,
    paired_block_bootstrap_gain,
    retain_control_candidates,
    wilson_interval,
)


class FixedModel:
    def __init__(self, probabilities: np.ndarray):
        self.values = np.asarray(probabilities, dtype=float)

    def probabilities(self, contexts: np.ndarray) -> np.ndarray:
        assert contexts.shape[:2] == self.values.shape
        return self.values


def test_conditional_scores_preserve_per_window_likelihood() -> None:
    contexts = np.zeros((3, 4, 3), dtype=np.uint8)
    counts = np.asarray([[8, 1, 1, 0], [0, 7, 2, 1], [1, 1, 1, 7]], dtype=float)
    probabilities = (counts + 1) / (counts + 1).sum(axis=1, keepdims=True)
    scores = conditional_control_scores(FixedModel(probabilities), contexts, counts)
    direct = np.sum(counts * np.log(probabilities), axis=1)
    assert np.allclose(scores.log_likelihood, direct)
    assert scores.conditional_nll == pytest.approx(-direct.sum() / counts.sum())
    assert scores.nll_gain > 0
    assert scores.deviance_per_cut >= 0


def test_paired_block_bootstrap_detects_consistent_gain() -> None:
    reference = np.full(60, -30.0)
    candidate = reference + 2.0
    totals = np.full(60, 10.0)
    blocks = np.repeat(["chr1", "chr2", "chr3"], 20)
    result = paired_block_bootstrap_gain(
        candidate, reference, totals, blocks, bootstraps=100, seed=8
    )
    assert result["paired_log_likelihood_gain_per_cut"] == pytest.approx(0.2)
    assert result["paired_gain_lower_95"] > 0


def test_retention_uses_one_se_and_requires_significant_larger_gain() -> None:
    frame = pd.DataFrame(
        [
            ("SELMA10", 10, 1.01, 0.02, 0.1, 0.01, np.nan),
            ("LOG21", 21, 1.00, 0.01, 0.2, 0.02, 0.01),
            ("LOG81", 81, 0.99, 0.02, 0.2, 0.08, 0.02),
        ],
        columns=[
            "candidate_id",
            "context_length",
            "mean_conditional_nll",
            "standard_error_conditional_nll",
            "minimum_library_nll_gain",
            "model_size_mb",
            "gain_over_smallest_lower_95",
        ],
    )
    selected = retain_control_candidates(frame)
    retained = selected[selected["retained"]]
    assert set(retained["candidate_id"]) == {"SELMA10", "LOG81"}


def test_motif_residual_flag_requires_large_confident_structure() -> None:
    positions = np.arange(-100, 101)
    expected = np.full((300, len(positions)), 5.0)
    observed = expected.copy()
    observed[:, np.abs(positions) <= 15] = 2.0
    result = motif_residual_effect(
        observed, expected, positions, bootstraps=100, seed=7, threshold=0.25
    )
    assert result["motif_residual_flag"]
    assert result["observed_minus_predicted_center_flank_effect"] > 0.25


def test_wilson_and_negative_control_safety_use_finite_support() -> None:
    lower, upper = wilson_interval(0, 1000)
    assert lower == 0
    assert 0 < upper < 0.005
    safe = negative_control_safety(
        np.zeros(1000, dtype=bool),
        np.zeros(1000, dtype=bool),
    )
    assert safe["passed_negative_control_safety"]
    unsafe = negative_control_safety(
        np.r_[np.ones(60, dtype=bool), np.zeros(940, dtype=bool)],
        np.zeros(1000, dtype=bool),
    )
    assert not unsafe["passed_negative_control_safety"]


def test_tobias_dwm_reference_learns_bias_and_roundtrips_safely(tmp_path) -> None:
    rng = np.random.default_rng(19)
    contexts = rng.integers(0, 4, size=(120, 20, 11), dtype=np.uint8)
    counts = rng.poisson(1.0, size=(120, 20)).astype(float)
    favored = contexts[:, :, 5] == 0
    counts[favored] += rng.poisson(5.0, size=np.sum(favored))
    model = TobiasDwmReferenceModel().fit(contexts, counts)
    training_scores = conditional_control_scores(model, contexts, counts)
    assert training_scores.nll_gain > 0
    assert 0 < model.score_scale < 2
    probabilities = model.probabilities(contexts[:4])
    assert probabilities.shape == (4, 20)
    assert np.allclose(probabilities.sum(axis=1), 1.0)
    favored_probability = probabilities[contexts[:4, :, 5] == 0].mean()
    other_probability = probabilities[contexts[:4, :, 5] != 0].mean()
    assert favored_probability > other_probability
    npz_path, _json_path = model.save(
        tmp_path / "dwm", metadata={"read_shift": [4, -4]}
    )
    restored = TobiasDwmReferenceModel.load(npz_path)
    assert restored.score_scale == pytest.approx(model.score_scale)
    assert np.allclose(restored.probabilities(contexts[:4]), probabilities)
    with np.load(npz_path, allow_pickle=False) as arrays:
        assert all(arrays[name].dtype != object for name in arrays.files)
