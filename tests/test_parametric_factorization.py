from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.special import logsumexp
from scipy.stats import multinomial
from sklearn.metrics import roc_auc_score

from fp_tools.tools.parametric_bias import BiasFeatureSpec
from fp_tools.tools.parametric_factorization import (
    FrozenBiasStrengthCalibrator,
    FrozenParametricFactorization,
    conditional_profile_log_likelihood,
    expected_profile_counts,
    fit_flank_accessibility_background,
    footprint_taper,
    natural_cubic_spline_basis,
)


def _multinomial_profiles(
    logits: np.ndarray,
    totals: np.ndarray,
    *,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    probabilities = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
    return np.vstack(
        [rng.multinomial(int(total), probability) for total, probability in zip(totals, probabilities)]
    ).astype(float)


def _synthetic_factorization_data(
    *,
    sites_per_tf: int = 60,
    footprint: bool = True,
    seed: int = 17,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    positions = np.arange(-100, 101, dtype=float)
    count = 2 * sites_per_tf
    samples = np.repeat("sample", count)
    tfs = np.repeat(["TF_A", "TF_B"], sites_per_tf)
    families = np.repeat(["family_1", "family_2"], sites_per_tf)
    log_bias = rng.normal(0.0, 0.35, size=(count, len(positions)))
    slopes = rng.normal(0.0, 0.25, size=count)
    broad = slopes[:, None] * positions[None, :] / 100.0
    labels = np.zeros(count, dtype=int)
    labels[np.arange(count) % 2 == 0] = 1
    effect = np.zeros_like(log_bias)
    if footprint:
        protection = -1.35 * np.exp(-0.5 * np.square(positions / 7.0))
        shoulders = 0.35 * (
            np.exp(-0.5 * np.square((positions - 20.0) / 5.0))
            + np.exp(-0.5 * np.square((positions + 20.0) / 5.0))
        )
        effect = labels[:, None] * (protection + shoulders)[None, :]
    logits = 1.2 * log_bias + broad + effect
    totals = np.repeat(500, count)
    counts = _multinomial_profiles(logits, totals, seed=seed + 1)
    return positions, counts, log_bias, samples, tfs, families, labels


def test_frozen_grid_exposes_all_receptive_fields() -> None:
    specs = [
        BiasFeatureSpec.selma10(),
        BiasFeatureSpec.loglinear21(),
        BiasFeatureSpec.loglinear41(),
        BiasFeatureSpec.loglinear81(),
    ]
    assert [spec.context_length for spec in specs] == [10, 21, 41, 81]
    assert specs[0].pair_distances == (1,)
    assert all(spec.pair_distances == (1, 5, 9) for spec in specs[1:])


def test_vectorized_conditional_likelihood_matches_scipy() -> None:
    counts = np.asarray([[2.0, 1.0, 4.0], [0.0, 3.0, 2.0]])
    logits = np.asarray([[0.2, -0.1, 0.8], [1.0, 0.0, -0.2]])
    # The public helper requires five positions; pad two zero-count categories.
    counts = np.pad(counts, ((0, 0), (0, 2)))
    logits = np.pad(logits, ((0, 0), (0, 2)), constant_values=-0.4)
    actual = conditional_profile_log_likelihood(counts, logits, include_constant=True)
    probabilities = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
    expected = np.asarray(
        [multinomial.logpmf(row, n=int(row.sum()), p=probability) for row, probability in zip(counts, probabilities)]
    )
    assert np.allclose(actual, expected)


def test_conditional_total_invariance_and_expected_total_preservation() -> None:
    counts = np.asarray([[1.0, 3.0, 2.0, 7.0, 1.0], [4.0, 0.0, 2.0, 1.0, 5.0]])
    logits = np.asarray([[0.2, 0.4, -0.3, 1.0, 0.1], [-0.4, 0.5, 0.2, 0.7, 0.0]])
    first = conditional_profile_log_likelihood(counts, logits).sum() / counts.sum()
    second = conditional_profile_log_likelihood(7.0 * counts, logits).sum() / (7.0 * counts.sum())
    assert first == pytest.approx(second)
    expected = expected_profile_counts(counts, logits)
    assert np.allclose(expected.sum(axis=1), counts.sum(axis=1))


def test_bias_strength_recovery_and_checksum_roundtrip(tmp_path: Path) -> None:
    rng = np.random.default_rng(5)
    log_bias = rng.normal(0.0, 0.8, size=(360, 41))
    samples = np.repeat(["A", "B"], 180)
    strengths = np.where(samples == "A", 0.45, 1.55)
    counts = _multinomial_profiles(
        strengths[:, None] * log_bias,
        np.repeat(80, len(samples)),
        seed=6,
    )
    calibrator = FrozenBiasStrengthCalibrator().fit(counts, log_bias, samples)
    assert calibrator.strength("A") == pytest.approx(0.45, abs=0.06)
    assert calibrator.strength("B") == pytest.approx(1.55, abs=0.06)
    assert all(estimate.nll_gain > 0 for estimate in calibrator.estimates.values())
    npz_path, _json_path = calibrator.save(tmp_path / "calibration", {"split": "chr1-15"})
    loaded = FrozenBiasStrengthCalibrator.load(npz_path)
    assert loaded.strength("B") == pytest.approx(calibrator.strength("B"))
    assert loaded.metadata["split"] == "chr1-15"
    with npz_path.open("ab") as handle:
        handle.write(b"broken")
    with pytest.raises(ValueError, match="checksum"):
        FrozenBiasStrengthCalibrator.load(npz_path)


def test_flank_only_natural_spline_recovers_broad_background() -> None:
    positions = np.arange(-100, 101, dtype=float)
    log_bias = np.zeros((3, len(positions)))
    broad = np.vstack(
        [
            0.2 + 0.6 * positions / 100.0,
            -0.4 - 0.3 * positions / 100.0,
            0.1 + 0.1 * positions / 100.0,
        ]
    )
    probabilities = np.exp(broad - logsumexp(broad, axis=1, keepdims=True))
    counts = 100000.0 * probabilities
    fitted = fit_flank_accessibility_background(
        counts,
        log_bias,
        positions,
        1.0,
        outer_start=50,
        outer_end=100,
    )
    centered_truth = broad - broad.mean(axis=1, keepdims=True)
    centered_fit = fitted.log_background - fitted.log_background.mean(axis=1, keepdims=True)
    assert np.corrcoef(centered_truth.ravel(), centered_fit.ravel())[0, 1] > 0.999
    assert np.max(np.abs(centered_truth - centered_fit)) < 0.02
    assert not fitted.outer_mask[np.argmin(np.abs(positions))]
    assert np.allclose(fitted.expected.sum(axis=1), counts.sum(axis=1))


def test_natural_spline_and_footprint_taper_contracts() -> None:
    positions = np.arange(-100, 101, dtype=float)
    basis = natural_cubic_spline_basis(positions, df=5)
    assert basis.matrix.shape == (201, 5)
    taper = footprint_taper(positions, limit=50)
    assert np.all(taper[np.abs(positions) >= 50] == 0)
    assert np.all(taper[np.abs(positions) <= 40] == 1)


def test_factorization_recovers_protection_without_erasing_total() -> None:
    positions, counts, log_bias, samples, tfs, families, labels = _synthetic_factorization_data()
    model = FrozenParametricFactorization(
        positions,
        family_shrinkage=5.0,
        tf_shrinkage=5.0,
        use_total_component=True,
    ).fit(
        counts,
        log_bias,
        samples,
        tfs,
        families,
        {"sample": 1.2},
        max_iter=12,
    )
    result = model.predict(counts, log_bias, samples, tfs)
    assert roc_auc_score(labels, result.posterior_bound) > 0.85
    center = np.abs(positions) <= 5
    outer = np.abs(positions) >= 60
    bound_curve = np.mean(result.footprint_log_effect[labels == 1], axis=0)
    assert np.mean(bound_curve[center]) < np.mean(bound_curve[outer]) - 0.25
    assert np.all(result.footprint_log_effect[:, np.abs(positions) >= 50] == 0)
    assert np.allclose(result.expected_unbound.sum(axis=1), counts.sum(axis=1))
    assert np.allclose(result.expected_bound.sum(axis=1), counts.sum(axis=1))
    for mode in ("difference", "pearson", "deviance", "log-ratio"):
        assert result.residual(counts, mode, dispersion=model.total_dispersion_).shape == counts.shape
    assert result.residual(counts, "negative-binomial").shape == (len(counts),)


def test_bias_only_profiles_do_not_create_a_footprint() -> None:
    positions, _counts, log_bias, samples, tfs, families, _labels = _synthetic_factorization_data(
        footprint=False,
        sites_per_tf=30,
    )
    broad = np.zeros_like(log_bias)
    logits = 1.2 * log_bias + broad
    probabilities = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
    # Exact expectation removes sampling-driven apparent protection.
    counts = 5000.0 * probabilities
    model = FrozenParametricFactorization(
        positions,
        family_shrinkage=5.0,
        tf_shrinkage=5.0,
        use_total_component=False,
    ).fit(
        counts,
        log_bias,
        samples,
        tfs,
        families,
        {"sample": 1.2},
        max_iter=6,
    )
    result = model.predict(counts, log_bias, samples, tfs)
    center = np.abs(positions) <= 15
    assert np.max(np.abs(result.footprint_log_effect[:, center])) < 0.08


def test_factorization_is_deterministic_and_checksum_safe(tmp_path: Path) -> None:
    positions, counts, log_bias, samples, tfs, families, _labels = _synthetic_factorization_data(
        sites_per_tf=24,
    )
    kwargs = dict(
        family_shrinkage=5.0,
        tf_shrinkage=5.0,
        use_total_component=False,
        seed=29,
    )
    first = FrozenParametricFactorization(positions, **kwargs).fit(
        counts,
        log_bias,
        samples,
        tfs,
        families,
        {"sample": 1.2},
        max_iter=7,
    )
    second = FrozenParametricFactorization(positions, **kwargs).fit(
        counts,
        log_bias,
        samples,
        tfs,
        families,
        {"sample": 1.2},
        max_iter=7,
    )
    first_result = first.predict(counts, log_bias, samples, tfs)
    second_result = second.predict(counts, log_bias, samples, tfs)
    assert np.array_equal(first_result.posterior_bound, second_result.posterior_bound)
    npz_path, _json_path = first.save(tmp_path / "factorization", {"configuration_hash": "abc"})
    loaded = FrozenParametricFactorization.load(npz_path)
    loaded_result = loaded.predict(counts, log_bias, samples, tfs)
    assert np.allclose(first_result.posterior_bound, loaded_result.posterior_bound)
    assert loaded.metadata["configuration_hash"] == "abc"
    with npz_path.open("ab") as handle:
        handle.write(b"broken")
    with pytest.raises(ValueError, match="checksum"):
        FrozenParametricFactorization.load(npz_path)
