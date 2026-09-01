from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

from fp_tools.tools.functional_footprints import (
    BiasAwareFunctionalMixture,
    ConditionalMultinomialMixture,
    CovariateAnchoredFdaModel,
    CovariateResidualizedFdaModel,
    construct_strand_functional_profiles,
    ExactAdditiveGPSmoother,
    FdaMixtureModel,
    FunctionalTemplateDetector,
    MultichannelFunctionalTemplateDetector,
    FunctionalPCA,
    HybridFdaGpModel,
    PenalizedSplineSmoother,
    SparseAdditiveGPSmoother,
    deviance_profiles,
    functional_differential_test,
    orient_profiles,
    profile_descriptors,
    site_accessibility_background,
    standardized_functional_separation,
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


def test_standardized_functional_separation_detects_curve_difference() -> None:
    rng = np.random.default_rng(41)
    x = _positions()
    labels = np.repeat([0, 1], 100)
    profiles = rng.normal(scale=0.3, size=(len(labels), len(x)))
    profiles[labels == 1] += _footprint_shape(x)
    separated = standardized_functional_separation(profiles, labels, x)
    absent = standardized_functional_separation(
        rng.normal(scale=0.3, size=profiles.shape),
        labels,
        x,
    )
    assert separated > 0.5
    assert separated > absent * 3


@pytest.mark.parametrize("smoother", ["spline", "gp"])
def test_shape_only_template_detector_transfers_footprint(smoother: str) -> None:
    rng = np.random.default_rng(92)
    x = _positions()
    shape = _footprint_shape(x)

    def sample(n: int) -> tuple[np.ndarray, np.ndarray]:
        labels = rng.integers(0, 2, size=n)
        broad = rng.normal(scale=0.2, size=(n, 1)) * (x / 50.0)[None, :]
        profiles = broad + rng.normal(scale=0.45, size=(n, len(x)))
        profiles += labels[:, None] * shape
        return profiles, labels

    train, train_labels = sample(500)
    transfer, transfer_labels = sample(300)
    model = FunctionalTemplateDetector(
        x,
        smoother=smoother,
        window_limit=45,
    ).fit(train, train_labels)
    scores = model.decision_function(transfer)
    assert roc_auc_score(transfer_labels, scores) > 0.85
    assert model.positive_sites_ > 100
    assert model.negative_sites_ > 100
    assert np.all(model.footprint_template_[np.abs(x) >= 45] == 0)


def test_template_detector_supports_hierarchical_prior() -> None:
    rng = np.random.default_rng(17)
    x = _positions()
    labels = np.repeat([0, 1], 50)
    profiles = rng.normal(scale=0.8, size=(len(labels), len(x)))
    profiles[labels == 1] += _footprint_shape(x)
    prior = 2.0 * _footprint_shape(x)
    model = FunctionalTemplateDetector(x).fit(
        profiles,
        labels,
        prior_template=prior,
        prior_strength=500,
    )
    assert np.corrcoef(model.raw_template_, prior)[0, 1] > 0.95


def test_multichannel_template_recovers_strand_specific_shape() -> None:
    rng = np.random.default_rng(123)
    x = _positions()
    shape = _footprint_shape(x)

    def sample(n: int) -> tuple[np.ndarray, np.ndarray]:
        labels = rng.integers(0, 2, size=n)
        values = rng.normal(scale=0.55, size=(n, 3, len(x)))
        values[:, 0, :] += labels[:, None] * 0.35 * shape
        values[:, 2, :] += labels[:, None] * 1.2 * shape
        return values, labels

    train, train_labels = sample(500)
    validation, validation_labels = sample(250)
    model = MultichannelFunctionalTemplateDetector(
        x,
        smoother="gp",
        window_limit=45,
    ).fit(train, train_labels)
    assert roc_auc_score(validation_labels, model.decision_function(validation)) > 0.9
    assert len(model.channel_models_) == 3
    assert np.argmax(np.abs(model.channel_discriminant_)) == 2


def test_strand_functional_profiles_reverse_and_swap_channels() -> None:
    plus = np.tile(np.arange(21, dtype=float), (2, 1))
    minus = np.tile(np.arange(100, 121, dtype=float), (2, 1))
    plus_expected = plus + 1.0
    minus_expected = minus + 1.0
    profiles = construct_strand_functional_profiles(
        plus,
        minus,
        plus_expected,
        minus_expected,
        ["+", "-"],
        dispersion=0.05,
    )
    assert np.array_equal(profiles.plus_observed[0], plus[0])
    assert np.array_equal(profiles.minus_observed[0], minus[0])
    assert np.array_equal(profiles.plus_observed[1], minus[1, ::-1])
    assert np.array_equal(profiles.minus_observed[1], plus[1, ::-1])
    assert profiles.combined_residual.shape == plus.shape
    assert profiles.shared_strand_residual.shape == plus.shape
    assert profiles.antisymmetric_strand_residual.shape == plus.shape


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


@pytest.mark.parametrize("method", ["linear", "quadratic", "gp-long"])
def test_site_accessibility_background_preserves_injected_footprint(method: str) -> None:
    rng = np.random.default_rng(31)
    x = np.arange(-100, 101, dtype=float)
    sites = 120
    labels = np.repeat([0, 1], sites // 2)
    slopes = rng.normal(scale=0.65, size=sites)
    curvature = rng.normal(loc=-0.35, scale=0.12, size=sites)
    broad = np.exp(
        slopes[:, None] * x[None, :] / 100.0
        + curvature[:, None] * np.square(x[None, :] / 100.0)
    )
    broad *= 800.0 / broad.sum(axis=1, keepdims=True)
    footprint = _footprint_shape(x)
    observed_probability = broad * np.exp(labels[:, None] * footprint[None, :])
    observed = observed_probability * 800.0 / observed_probability.sum(axis=1, keepdims=True)
    sequence_only = np.full_like(observed, 800.0 / len(x))

    adjusted = site_accessibility_background(
        observed,
        sequence_only,
        x,
        method=method,
        exclusion=50.0,
        ridge=3.0,
    )
    assert np.allclose(adjusted.sum(axis=1), observed.sum(axis=1))
    assert np.mean(np.square(adjusted - broad)) < np.mean(np.square(sequence_only - broad))
    center = np.abs(x) <= 6
    residual = np.log((observed + 0.5) / (adjusted + 0.5))
    assert np.mean(residual[labels == 1][:, center]) < -0.15
    assert abs(np.mean(residual[labels == 0][:, center])) < 0.15


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


def test_bias_aware_mixture_can_anchor_binding_prior_direction(tmp_path: Path) -> None:
    observed, expected, labels, motif_score = _synthetic_counts(seed=19)
    model = BiasAwareFunctionalMixture(
        _positions(),
        smoother="gp",
        dispersion=0.02,
        prior_constraint="motif-accessibility",
        accessibility_background="linear",
        background_exclusion=30.0,
        profile_inner_limit=30.0,
        profile_outer_limit=40.0,
        likelihood_limit=30.0,
        max_iter=60,
    )
    result = model.fit(
        observed,
        expected,
        motif_score=motif_score,
        accessibility=observed.sum(axis=1),
    )
    assert np.all(result.prior_coefficients[1:] >= 0)
    assert roc_auc_score(labels, result.posterior) > 0.75
    assert np.allclose(result.footprint_profile[np.abs(_positions()) >= 40], 0.0)
    model.save(tmp_path / "anchored")
    loaded = BiasAwareFunctionalMixture.load(tmp_path / "anchored.npz")
    assert loaded.prior_constraint == "motif-accessibility"
    assert loaded.accessibility_background == "linear"
    assert loaded.profile_outer_limit == 40.0
    assert loaded.likelihood_limit == 30.0
    shape_log_odds, prior_log_odds = loaded.predict_log_odds_components(
        observed[:25],
        expected[:25],
        motif_score=motif_score[:25],
        accessibility=observed[:25].sum(axis=1),
    )
    reconstructed = 1.0 / (1.0 + np.exp(-(shape_log_odds + prior_log_odds)))
    assert np.allclose(
        reconstructed,
        loaded.predict(
            observed[:25],
            expected[:25],
            motif_score=motif_score[:25],
            accessibility=observed[:25].sum(axis=1),
        ),
    )


@pytest.mark.parametrize("smoother", ["spline", "gp"])
def test_conditional_multinomial_mixture_detects_profile_not_total(
    smoother: str,
    tmp_path: Path,
) -> None:
    observed, expected, labels, _motif_score = _synthetic_counts(seed=44)
    model = ConditionalMultinomialMixture(
        _positions(),
        smoother=smoother,
        max_iter=60,
        tolerance=1e-6,
        profile_outer_limit=50.0,
        likelihood_limit=50.0,
    )
    result = model.fit(observed, expected)
    probability = model.predict(observed, expected)
    assert result.converged
    assert roc_auc_score(labels, probability) > 0.90
    assert result.descriptors.depletion > 0.25
    assert np.allclose(
        model._conditional_probabilities(expected),
        model._conditional_probabilities(7.0 * expected),
    )
    model.save(tmp_path / f"conditional_{smoother}")
    loaded = ConditionalMultinomialMixture.load(
        tmp_path / f"conditional_{smoother}.npz"
    )
    assert np.allclose(
        loaded.predict(observed[:30], expected[:30]),
        probability[:30],
    )


def test_conditional_multinomial_bias_only_profiles_remain_shallow() -> None:
    rng = np.random.default_rng(31)
    positions = _positions()
    probability = np.exp(
        0.25 * np.sin(positions / 5.0) + 0.15 * np.cos(positions / 11.0)
    )
    probability /= probability.sum()
    totals = rng.integers(130, 260, size=400)
    expected = totals[:, None] * probability[None, :]
    observed = np.vstack([rng.multinomial(int(total), probability) for total in totals])
    result = ConditionalMultinomialMixture(
        positions,
        smoother="spline",
        max_iter=60,
        tolerance=1e-6,
    ).fit(observed, expected)
    assert result.converged
    assert result.descriptors.depletion < 0.15


def test_fda_and_hybrid_models_detect_shape(tmp_path: Path) -> None:
    observed, expected, labels, _motif_score = _synthetic_counts()
    x = _positions()
    residuals = deviance_profiles(observed, expected, dispersion=0.02)
    fda = FdaMixtureModel(max_components=12, seed=6).fit(residuals, positions=x)
    hybrid = HybridFdaGpModel(x, max_components=12, seed=6).fit(residuals)
    fda_probability = fda.predict_proba(residuals)
    hybrid_probability = hybrid.predict_proba(residuals)
    assert roc_auc_score(labels, fda_probability) > 0.70
    assert roc_auc_score(labels, hybrid_probability) > 0.70

    fda.save(tmp_path / "fda", metadata={"labels_used": False})
    hybrid.save(tmp_path / "hybrid", metadata={"labels_used": False})
    loaded_fda = FdaMixtureModel.load(tmp_path / "fda.npz")
    loaded_hybrid = HybridFdaGpModel.load(tmp_path / "hybrid.npz")
    assert np.allclose(loaded_fda.predict_proba(residuals), fda_probability)
    assert np.allclose(loaded_fda.component_profiles(), fda.component_profiles())
    assert np.allclose(loaded_hybrid.predict_proba(residuals), hybrid_probability)


def test_fda_model_rejects_checksum_mismatch(tmp_path: Path) -> None:
    observed, expected, _labels, _motif_score = _synthetic_counts(seed=42)
    residuals = deviance_profiles(observed, expected, dispersion=0.02)
    model_path, _metadata_path = FdaMixtureModel(seed=3).fit(
        residuals,
        positions=_positions(),
    ).save(tmp_path / "fda")
    with model_path.open("ab") as handle:
        handle.write(b"modified")
    with pytest.raises(ValueError, match="checksum"):
        FdaMixtureModel.load(model_path)


def test_covariate_anchored_fda_separates_shape_from_prior(tmp_path: Path) -> None:
    observed, expected, labels, motif_score = _synthetic_counts(seed=27)
    residuals = deviance_profiles(observed, expected, dispersion=0.02)
    accessibility = observed.sum(axis=1)
    model = CovariateAnchoredFdaModel(
        max_components=12,
        anchor_strength=0.7,
        seed=11,
    ).fit(
        residuals,
        motif_score=motif_score,
        accessibility=accessibility,
        positions=_positions(),
        sample_weight=np.sqrt(accessibility),
    )
    shape, prior = model.predict_log_odds_components(
        residuals,
        motif_score=motif_score,
        accessibility=accessibility,
    )
    assert model.converged_
    assert roc_auc_score(labels, shape) > 0.70
    assert shape.shape == prior.shape == labels.shape
    assert model.profile_difference().shape == _positions().shape
    model.save(tmp_path / "anchored")
    loaded = CovariateAnchoredFdaModel.load(tmp_path / "anchored.npz")
    loaded_shape, loaded_prior = loaded.predict_log_odds_components(
        residuals,
        motif_score=motif_score,
        accessibility=accessibility,
    )
    assert np.allclose(loaded_shape, shape)
    assert np.allclose(loaded_prior, prior)


def test_covariate_residualized_fda_removes_accessibility_pc_trend(
    tmp_path: Path,
) -> None:
    observed, expected, labels, motif_score = _synthetic_counts(seed=31)
    residuals = deviance_profiles(observed, expected, dispersion=0.02)
    accessibility = observed.sum(axis=1)
    model = CovariateResidualizedFdaModel(
        max_components=12,
        covariate_ridge=1.0,
        seed=9,
    ).fit(
        residuals,
        motif_score=motif_score,
        accessibility=accessibility,
        positions=_positions(),
        sample_weight=np.sqrt(accessibility),
    )
    scores = model.transform_residual_scores(
        residuals,
        motif_score=motif_score,
        accessibility=accessibility,
    )
    standardized_accessibility = (
        np.log1p(accessibility) - np.mean(np.log1p(accessibility))
    ) / np.std(np.log1p(accessibility))
    correlations = [
        abs(float(np.corrcoef(scores[:, index], standardized_accessibility)[0, 1]))
        for index in range(scores.shape[1])
    ]
    probabilities = model.predict_proba(
        residuals,
        motif_score=motif_score,
        accessibility=accessibility,
    )
    assert max(correlations) < 0.1
    assert roc_auc_score(labels, probabilities) > 0.65
    assert model.profile_difference().shape == _positions().shape
    model.save(tmp_path / "residualized")
    loaded = CovariateResidualizedFdaModel.load(tmp_path / "residualized.npz")
    assert np.allclose(
        loaded.predict_proba(
            residuals,
            motif_score=motif_score,
            accessibility=accessibility,
        ),
        probabilities,
    )


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
