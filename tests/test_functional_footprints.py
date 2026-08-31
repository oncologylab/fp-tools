from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

from fp_tools.tools.functional_footprints import (
    BiasAwareFunctionalMixture,
    ExactAdditiveGPSmoother,
    FdaMixtureModel,
    FunctionalPCA,
    HybridFdaGpModel,
    PenalizedSplineSmoother,
    SparseAdditiveGPSmoother,
    deviance_profiles,
    functional_differential_test,
    orient_profiles,
    profile_descriptors,
)


def _positions(width: int = 101) -> np.ndarray:
    return np.arange(width, dtype=float) - width // 2


def _footprint_shape(x: np.ndarray) -> np.ndarray:
    center = -0.8 * np.exp(-0.5 * np.square(x / 6.0))
    shoulders = 0.25 * (
        np.exp(-0.5 * np.square((x - 20.0) / 5.0))
        + np.exp(-0.5 * np.square((x + 20.0) / 5.0))
    )
    return center + shoulders


def test_profile_orientation_and_descriptors() -> None:
    x = _positions()
    shape = _footprint_shape(x)
    profiles = np.vstack([shape, shape[::-1]])
    oriented = orient_profiles(profiles, ["+", "-"])
    assert np.allclose(oriented[0], oriented[1])
    descriptors = profile_descriptors(shape, x)
    assert descriptors.depletion > 0.5
    assert 5 <= descriptors.width <= 30
    assert 10 <= descriptors.shoulder_distance <= 30
    assert abs(descriptors.asymmetry) < 0.05


def test_weighted_functional_pca_roundtrip(tmp_path: Path) -> None:
    rng = np.random.default_rng(4)
    x = _positions()
    shape = _footprint_shape(x)
    slope = x / np.max(np.abs(x))
    coefficients = rng.normal(size=(300, 2))
    profiles = (
        coefficients[:, [0]] * shape[None, :]
        + coefficients[:, [1]] * slope[None, :]
        + rng.normal(scale=0.03, size=(300, len(x)))
    )
    weights = rng.uniform(0.5, 3.0, size=len(profiles))
    model = FunctionalPCA(variance_threshold=0.95, max_components=20, seed=9)
    scores = model.fit_transform(profiles, sample_weight=weights)
    reconstructed = model.inverse_transform(scores)
    assert scores.shape[1] <= 20
    assert np.mean(np.square(reconstructed - profiles)) < 0.003
    model.save(tmp_path / "functional")
    loaded = FunctionalPCA.load(tmp_path / "functional.npz")
    assert np.allclose(loaded.transform(profiles[:10]), model.transform(profiles[:10]))


def test_sparse_gp_tracks_exact_reference_and_spline() -> None:
    rng = np.random.default_rng(8)
    x = _positions()
    truth = _footprint_shape(x) + 0.1 * np.sin(x / 25.0)
    noisy = truth + rng.normal(scale=0.12, size=len(x))
    weights = np.full(len(x), 1.0 / 0.12**2)
    sparse = SparseAdditiveGPSmoother(x, inducing_points=25).fit(noisy, weights)
    exact = ExactAdditiveGPSmoother(x, noise=0.12**2).fit(noisy, weights)
    spline = PenalizedSplineSmoother(x, n_basis=25, penalty=10.0).fit(noisy, weights)
    assert np.corrcoef(sparse.mean, exact.mean)[0, 1] > 0.95
    assert np.mean(np.square(sparse.mean - truth)) < np.mean(np.square(noisy - truth))
    assert np.mean(np.square(spline.mean - truth)) < np.mean(np.square(noisy - truth))
    assert np.all(np.isfinite(sparse.standard_error))


def _synthetic_counts(seed: int = 12) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = _positions()
    sites = 500
    labels = np.repeat([0, 1], sites // 2)
    rng.shuffle(labels)
    sequence_bias = np.exp(0.25 * np.sin(x / 5.0) + 0.15 * np.cos(x / 11.0))
    sequence_bias /= sequence_bias.sum()
    footprint = _footprint_shape(x)
    expected = np.empty((sites, len(x)))
    observed = np.empty_like(expected)
    for index, label in enumerate(labels):
        total = int(rng.integers(130, 260))
        expected[index] = total * sequence_bias
        probability = sequence_bias * np.exp(footprint * label)
        probability /= probability.sum()
        observed[index] = rng.multinomial(total, probability)
    motif_score = rng.normal(loc=labels * 0.35, scale=1.0, size=sites)
    return observed, expected, labels, motif_score


@pytest.mark.parametrize("smoother", ["spline", "gp"])
def test_bias_aware_functional_mixture_detects_bound_profiles(
    smoother: str,
    tmp_path: Path,
) -> None:
    observed, expected, labels, motif_score = _synthetic_counts()
    x = _positions()
    model = BiasAwareFunctionalMixture(
        x,
        smoother=smoother,
        dispersion=0.02,
        max_iter=60,
        tolerance=1e-6,
    )
    result = model.fit(observed, expected, motif_score=motif_score, accessibility=observed.sum(axis=1))
    assert roc_auc_score(labels, result.posterior) > 0.80
    assert result.descriptors.depletion > 0.25
    assert result.iterations <= 60
    model.save(tmp_path / f"{smoother}_model")
    loaded = BiasAwareFunctionalMixture.load(tmp_path / f"{smoother}_model.npz")
    assert np.allclose(
        loaded.predict(observed[:20], expected[:20], motif_score=motif_score[:20], accessibility=observed[:20].sum(axis=1)),
        model.predict(observed[:20], expected[:20], motif_score=motif_score[:20], accessibility=observed[:20].sum(axis=1)),
    )


def test_fda_and_hybrid_models_detect_shape() -> None:
    observed, expected, labels, _motif_score = _synthetic_counts()
    x = _positions()
    residuals = deviance_profiles(observed, expected, dispersion=0.02)
    fda = FdaMixtureModel(max_components=12, seed=6).fit(residuals, positions=x)
    hybrid = HybridFdaGpModel(x, max_components=12, seed=6).fit(residuals)
    assert roc_auc_score(labels, fda.predict_proba(residuals)) > 0.70
    assert roc_auc_score(labels, hybrid.predict_proba(residuals)) > 0.70


def test_replicate_level_functional_differential_test() -> None:
    rng = np.random.default_rng(23)
    x = _positions()
    condition_effect = 0.55 * _footprint_shape(x)
    profiles = []
    conditions = []
    replicates = []
    for condition in ("control", "stress"):
        for replicate in range(4):
            replicate_offset = rng.normal(scale=0.03, size=len(x))
            for _site in range(40):
                value = replicate_offset + rng.normal(scale=0.12, size=len(x))
                if condition == "stress":
                    value = value + condition_effect
                profiles.append(value)
                conditions.append(condition)
                replicates.append(f"{condition}_{replicate}")
    result = functional_differential_test(
        np.asarray(profiles),
        conditions,
        replicates,
        ("stress", "control"),
        positions=x,
        n_bootstrap=300,
        seed=5,
    )
    assert result.unit == "replicate"
    assert result.global_pvalue < 0.1
    assert result.descriptor_change.depletion > 0.2
    assert result.difference.shape == x.shape
    assert np.all(result.simultaneous_lower <= result.simultaneous_upper)

