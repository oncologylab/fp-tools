#!/usr/bin/env python
"""Functional-data and Gaussian-process models for motif-centered footprints.

The functions in this module are deliberately independent of the CLI.  They
support the locked footprint-improvement benchmark and can later be promoted
to command options only if the preregistered scientific and performance gates
pass.  All deployable implementations use NumPy/SciPy/scikit-learn already
required by fp-tools; no deep-learning framework is introduced.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.interpolate import BSpline
from scipy.linalg import solve_triangular
from scipy.special import expit, logsumexp
from sklearn.mixture import GaussianMixture
from sklearn.utils.extmath import randomized_svd

from fp_tools.tools.parametric_bias import (
    calibrated_residuals,
    negative_binomial_log_likelihood,
)


FUNCTIONAL_SCHEMA = "fp-tools-functional-footprint-v1"


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_profiles(profiles: np.ndarray, *, nonnegative: bool = False) -> np.ndarray:
    values = np.asarray(profiles, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] < 1 or values.shape[1] < 5:
        raise ValueError("profiles must have shape (sites, positions) with at least five positions")
    if nonnegative and np.any(values[np.isfinite(values)] < 0):
        raise ValueError("count profiles must be non-negative")
    return values


def orient_profiles(profiles: np.ndarray, strands: Iterable[str | int | bool]) -> np.ndarray:
    """Orient motif-centered profiles so every motif points in one direction."""

    values = _validate_profiles(profiles)
    raw_strands = list(strands)
    if len(raw_strands) != len(values):
        raise ValueError("strands must contain one value per profile")
    reverse = np.asarray(
        [value in ("-", -1, "reverse", True) for value in raw_strands],
        dtype=bool,
    )
    output = values.copy()
    output[reverse] = output[reverse, ::-1]
    return output


@dataclass(frozen=True)
class StrandFunctionalProfiles:
    """Orientation-aligned strand channels for one motif-site collection."""

    plus_observed: np.ndarray
    minus_observed: np.ndarray
    plus_expected: np.ndarray
    minus_expected: np.ndarray
    combined_residual: np.ndarray
    shared_strand_residual: np.ndarray
    antisymmetric_strand_residual: np.ndarray


def construct_strand_functional_profiles(
    plus_observed: np.ndarray,
    minus_observed: np.ndarray,
    plus_expected: np.ndarray,
    minus_expected: np.ndarray,
    strands: Iterable[str | int | bool],
    *,
    dispersion: float = 0.0,
) -> StrandFunctionalProfiles:
    """Orient cut channels and construct shared/antisymmetric residuals.

    On reverse-strand motifs genomic coordinates are reversed and the genomic
    plus/minus tracks are swapped. This preserves biological strand in a
    common motif orientation rather than merely reversing a combined track.
    """

    arrays = [
        _validate_profiles(values, nonnegative=True)
        for values in (plus_observed, minus_observed, plus_expected, minus_expected)
    ]
    if any(values.shape != arrays[0].shape for values in arrays[1:]):
        raise ValueError("all strand profiles must have equal shape")
    raw_strands = list(strands)
    if len(raw_strands) != len(arrays[0]):
        raise ValueError("strands must contain one value per profile")
    reverse = np.asarray(
        [value in ("-", -1, "reverse", True) for value in raw_strands],
        dtype=bool,
    )
    plus_observed_out, minus_observed_out, plus_expected_out, minus_expected_out = (
        values.copy() for values in arrays
    )
    plus_observed_out[reverse] = arrays[1][reverse, ::-1]
    minus_observed_out[reverse] = arrays[0][reverse, ::-1]
    plus_expected_out[reverse] = arrays[3][reverse, ::-1]
    minus_expected_out[reverse] = arrays[2][reverse, ::-1]

    plus_residual = calibrated_residuals(
        plus_observed_out,
        plus_expected_out,
        "deviance",
        dispersion=dispersion,
    )
    minus_residual = calibrated_residuals(
        minus_observed_out,
        minus_expected_out,
        "deviance",
        dispersion=dispersion,
    )
    combined_residual = calibrated_residuals(
        plus_observed_out + minus_observed_out,
        plus_expected_out + minus_expected_out,
        "deviance",
        dispersion=dispersion,
    )
    normalization = np.sqrt(2.0)
    return StrandFunctionalProfiles(
        plus_observed=plus_observed_out,
        minus_observed=minus_observed_out,
        plus_expected=plus_expected_out,
        minus_expected=minus_expected_out,
        combined_residual=combined_residual,
        shared_strand_residual=(plus_residual + minus_residual) / normalization,
        antisymmetric_strand_residual=(plus_residual - minus_residual) / normalization,
    )


def normalize_functional_profiles(
    profiles: np.ndarray,
    positions: np.ndarray | None = None,
    *,
    outer_flank_start: float = 60.0,
    clip: float = 12.0,
) -> np.ndarray:
    """Remove broad per-site trends and scale by outer-flank RMS.

    Functional mixtures should distinguish localized protection shapes rather
    than rediscover total accessibility.  A line is fitted only to the outer
    flanks, extrapolated across the window, and the detrended curve is divided
    by its outer-flank RMS.  This transformation is label-free and symmetric.
    """

    values = _validate_profiles(profiles)
    x = (
        np.asarray(positions, dtype=np.float64)
        if positions is not None
        else np.arange(values.shape[1], dtype=float) - values.shape[1] // 2
    )
    if x.shape != (values.shape[1],):
        raise ValueError("positions must match profile width")
    effective_outer_start = min(float(outer_flank_start), 0.6 * float(np.max(np.abs(x))))
    outer = np.abs(x) >= effective_outer_start
    if outer.sum() < 4:
        raise ValueError("outer_flank_start leaves fewer than four baseline positions")
    finite = np.isfinite(values)
    column_fill = np.nanmedian(np.where(finite, values, np.nan), axis=0)
    column_fill = np.nan_to_num(column_fill, nan=0.0)
    filled = np.where(finite, values, column_fill)
    design_outer = np.column_stack([np.ones(outer.sum()), x[outer]])
    projection = np.linalg.pinv(design_outer)
    coefficients = filled[:, outer] @ projection.T
    baseline = coefficients[:, [0]] + coefficients[:, [1]] * x[None, :]
    detrended = filled - baseline
    scale = np.sqrt(np.mean(np.square(detrended[:, outer]), axis=1))
    positive = scale[scale > 1e-8]
    fallback = float(np.median(positive)) if len(positive) else 1.0
    scale = np.where(scale > 1e-8, scale, fallback)
    normalized = detrended / scale[:, None]
    return np.clip(normalized, -abs(float(clip)), abs(float(clip)))


@dataclass(frozen=True)
class ProfileDescriptors:
    center: float
    shoulders: float
    depletion: float
    width: float
    shoulder_distance: float
    asymmetry: float
    periodicity: float


def profile_descriptors(profile: np.ndarray, positions: np.ndarray | None = None) -> ProfileDescriptors:
    """Summarize a functional footprint without treating it as occupancy proof."""

    values = np.asarray(profile, dtype=np.float64)
    if values.ndim != 1 or len(values) < 21:
        raise ValueError("profile must be a one-dimensional array with at least 21 positions")
    x = (
        np.asarray(positions, dtype=np.float64)
        if positions is not None
        else np.arange(len(values), dtype=np.float64) - len(values) // 2
    )
    if x.shape != values.shape:
        raise ValueError("positions must match profile length")
    finite = np.isfinite(values)
    if not finite.any():
        return ProfileDescriptors(*(float("nan") for _ in range(7)))
    filled = np.where(finite, values, np.nanmedian(values[finite]))
    center_mask = np.abs(x) <= 5
    left_shoulder = (x >= -40) & (x <= -15)
    right_shoulder = (x >= 15) & (x <= 40)
    shoulder_mask = left_shoulder | right_shoulder
    center = float(np.mean(filled[center_mask]))
    shoulders = float(np.mean(filled[shoulder_mask]))
    depletion = shoulders - center

    threshold = center + 0.5 * depletion
    protected = filled <= threshold if depletion >= 0 else filled >= threshold
    midpoint = int(np.argmin(np.abs(x)))
    left = midpoint
    right = midpoint
    while left > 0 and protected[left - 1]:
        left -= 1
    while right + 1 < len(values) and protected[right + 1]:
        right += 1
    width = float(abs(x[right] - x[left]) + 1.0)

    left_index = np.flatnonzero(left_shoulder)
    right_index = np.flatnonzero(right_shoulder)
    if depletion >= 0:
        left_peak = left_index[np.argmax(filled[left_index])]
        right_peak = right_index[np.argmax(filled[right_index])]
    else:
        left_peak = left_index[np.argmin(filled[left_index])]
        right_peak = right_index[np.argmin(filled[right_index])]
    shoulder_distance = float((abs(x[left_peak]) + abs(x[right_peak])) / 2.0)
    denominator = abs(float(np.mean(filled[shoulder_mask]))) + np.finfo(float).eps
    asymmetry = float(
        (np.mean(filled[right_shoulder]) - np.mean(filled[left_shoulder])) / denominator
    )

    centered = filled - np.mean(filled)
    spectrum = np.abs(np.fft.rfft(centered))
    frequencies = np.fft.rfftfreq(len(centered), d=1.0)
    periodic_band = (frequencies >= 1.0 / 13.0) & (frequencies <= 1.0 / 8.0)
    periodicity = float(np.sqrt(np.sum(np.square(spectrum[periodic_band]))) / (np.linalg.norm(spectrum) + np.finfo(float).eps))
    return ProfileDescriptors(center, shoulders, depletion, width, shoulder_distance, asymmetry, periodicity)


def standardized_functional_separation(
    profiles: np.ndarray,
    labels: Iterable[int | bool],
    positions: np.ndarray | None = None,
    *,
    limit: float = 50.0,
) -> float:
    """RMS standardized bound/unbound curve difference within a fixed window."""

    values = _validate_profiles(profiles)
    group = np.asarray(list(labels), dtype=bool)
    if group.shape != (len(values),):
        raise ValueError("labels must contain one value per profile")
    if np.sum(group) < 2 or np.sum(~group) < 2:
        return float("nan")
    x = (
        np.asarray(positions, dtype=float)
        if positions is not None
        else np.arange(values.shape[1], dtype=float) - values.shape[1] // 2
    )
    if x.shape != (values.shape[1],) or limit <= 0:
        raise ValueError("positions must match profiles and limit must be positive")
    selected = np.abs(x) <= float(limit)
    positive = values[group][:, selected]
    negative = values[~group][:, selected]
    difference = np.nanmean(positive, axis=0) - np.nanmean(negative, axis=0)
    positive_variance = np.nanvar(positive, axis=0, ddof=1)
    negative_variance = np.nanvar(negative, axis=0, ddof=1)
    pooled = (
        (len(positive) - 1) * positive_variance
        + (len(negative) - 1) * negative_variance
    ) / max(len(positive) + len(negative) - 2, 1)
    finite_positive = pooled[np.isfinite(pooled) & (pooled > 0)]
    variance_floor = (
        0.01 * float(np.median(finite_positive)) if len(finite_positive) else 1.0
    )
    standardized = difference / np.sqrt(np.maximum(pooled, variance_floor))
    return float(np.sqrt(np.nanmean(np.square(standardized))))


class FunctionalPCA:
    """Coverage-weighted functional PCA with deterministic component signs."""

    def __init__(self, variance_threshold: float = 0.95, max_components: int = 20, seed: int = 2026):
        if not 0 < variance_threshold <= 1:
            raise ValueError("variance_threshold must be in (0, 1]")
        if max_components < 1:
            raise ValueError("max_components must be positive")
        self.variance_threshold = float(variance_threshold)
        self.max_components = int(max_components)
        self.seed = int(seed)
        self.mean_: np.ndarray | None = None
        self.components_: np.ndarray | None = None
        self.explained_variance_ratio_: np.ndarray | None = None
        self.impute_: np.ndarray | None = None

    def fit(self, profiles: np.ndarray, sample_weight: np.ndarray | None = None) -> "FunctionalPCA":
        values = _validate_profiles(profiles)
        weights = np.ones(len(values), dtype=np.float64) if sample_weight is None else np.asarray(sample_weight, dtype=np.float64)
        if weights.shape != (len(values),) or np.any(weights < 0) or not np.any(weights > 0):
            raise ValueError("sample_weight must be non-negative with one positive value per profile")
        finite = np.isfinite(values)
        weighted_counts = np.sum(weights[:, None] * finite, axis=0)
        weighted_sums = np.sum(weights[:, None] * np.where(finite, values, 0.0), axis=0)
        impute = weighted_sums / np.maximum(weighted_counts, np.finfo(float).eps)
        filled = np.where(finite, values, impute)
        mean = np.average(filled, axis=0, weights=weights)
        centered = (filled - mean) * np.sqrt(weights[:, None] / np.mean(weights))
        rank = min(self.max_components, centered.shape[0] - 1, centered.shape[1])
        if rank < 1:
            raise ValueError("at least two profiles are required for functional PCA")
        _u, singular_values, components = randomized_svd(
            centered,
            n_components=rank,
            random_state=self.seed,
        )
        total_variance = float(np.sum(np.square(centered)))
        ratios = np.square(singular_values) / max(total_variance, np.finfo(float).eps)
        cumulative = np.cumsum(ratios)
        selected = int(np.searchsorted(cumulative, self.variance_threshold) + 1)
        selected = min(selected, rank)
        components = components[:selected]
        for index in range(len(components)):
            pivot = int(np.argmax(np.abs(components[index])))
            if components[index, pivot] < 0:
                components[index] *= -1.0
        self.mean_ = mean
        self.components_ = components
        self.explained_variance_ratio_ = ratios[:selected]
        self.impute_ = impute
        return self

    def _check_fitted(self) -> None:
        if self.mean_ is None or self.components_ is None or self.impute_ is None:
            raise ValueError("functional PCA has not been fitted")

    def transform(self, profiles: np.ndarray) -> np.ndarray:
        self._check_fitted()
        values = _validate_profiles(profiles)
        if values.shape[1] != len(self.mean_):
            raise ValueError("profile width does not match fitted functional PCA")
        filled = np.where(np.isfinite(values), values, self.impute_)
        return (filled - self.mean_) @ self.components_.T

    def inverse_transform(self, scores: np.ndarray) -> np.ndarray:
        self._check_fitted()
        values = np.asarray(scores, dtype=np.float64)
        if values.shape[-1] != len(self.components_):
            raise ValueError("score width does not match fitted functional PCA")
        return values @ self.components_ + self.mean_

    def fit_transform(self, profiles: np.ndarray, sample_weight: np.ndarray | None = None) -> np.ndarray:
        return self.fit(profiles, sample_weight=sample_weight).transform(profiles)

    def save(self, path: str | Path, metadata: dict[str, Any] | None = None) -> tuple[Path, Path]:
        self._check_fitted()
        npz_path = Path(path).with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            mean=self.mean_,
            components=self.components_,
            explained_variance_ratio=self.explained_variance_ratio_,
            impute=self.impute_,
        )
        document = {
            "schema": FUNCTIONAL_SCHEMA,
            "model_type": "functional_pca",
            "npz_sha256": _sha256_file(npz_path),
            "variance_threshold": self.variance_threshold,
            "max_components": self.max_components,
            "seed": self.seed,
            "metadata": dict(metadata or {}),
        }
        json_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "FunctionalPCA":
        npz_path = Path(path).with_suffix(".npz")
        document = json.loads(npz_path.with_suffix(".json").read_text(encoding="utf-8"))
        if document.get("schema") != FUNCTIONAL_SCHEMA or document.get("model_type") != "functional_pca":
            raise ValueError("unsupported functional PCA model")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("functional PCA checksum does not match its metadata")
        model = cls(
            variance_threshold=float(document["variance_threshold"]),
            max_components=int(document["max_components"]),
            seed=int(document["seed"]),
        )
        with np.load(npz_path, allow_pickle=False) as arrays:
            model.mean_ = np.asarray(arrays["mean"], dtype=np.float64)
            model.components_ = np.asarray(arrays["components"], dtype=np.float64)
            model.explained_variance_ratio_ = np.asarray(arrays["explained_variance_ratio"], dtype=np.float64)
            model.impute_ = np.asarray(arrays["impute"], dtype=np.float64)
        return model


@dataclass(frozen=True)
class SmoothResult:
    mean: np.ndarray
    standard_error: np.ndarray
    effective_parameters: int


def _bspline_basis(positions: np.ndarray, n_basis: int = 25, degree: int = 3) -> np.ndarray:
    x = np.asarray(positions, dtype=np.float64)
    n_basis = min(int(n_basis), len(x))
    if n_basis <= degree:
        raise ValueError("n_basis must exceed spline degree")
    internal_count = n_basis - degree - 1
    internal = (
        np.linspace(x.min(), x.max(), internal_count + 2)[1:-1]
        if internal_count > 0
        else np.array([], dtype=float)
    )
    knots = np.concatenate(
        [np.repeat(x.min(), degree + 1), internal, np.repeat(x.max(), degree + 1)]
    )
    return np.asarray(BSpline(knots, np.eye(n_basis), degree, extrapolate=False)(x), dtype=np.float64)


class PenalizedSplineSmoother:
    """Fast 25-basis penalized spline baseline."""

    def __init__(self, positions: np.ndarray, n_basis: int = 25, penalty: float = 10.0):
        self.positions = np.asarray(positions, dtype=np.float64)
        self.basis = _bspline_basis(self.positions, n_basis=n_basis)
        self.penalty = float(penalty)
        differences = np.diff(np.eye(self.basis.shape[1]), n=2, axis=0)
        self.precision = differences.T @ differences

    def fit(self, values: np.ndarray, weights: np.ndarray | None = None) -> SmoothResult:
        y = np.asarray(values, dtype=np.float64)
        if y.shape != self.positions.shape:
            raise ValueError("values must match smoother positions")
        w = np.ones_like(y) if weights is None else np.asarray(weights, dtype=np.float64)
        if w.shape != y.shape or np.any(w < 0):
            raise ValueError("weights must be non-negative and match values")
        finite = np.isfinite(y) & np.isfinite(w) & (w > 0)
        if finite.sum() < 4:
            return SmoothResult(np.nan_to_num(y), np.full_like(y, np.nan), 0)
        basis = self.basis[finite]
        weighted = basis * w[finite, None]
        system = basis.T @ weighted + self.penalty * self.precision + 1e-8 * np.eye(basis.shape[1])
        right = basis.T @ (w[finite] * y[finite])
        coefficients = np.linalg.solve(system, right)
        mean = self.basis @ coefficients
        residual = y[finite] - basis @ coefficients
        variance = float(np.average(np.square(residual), weights=w[finite]))
        inverse = np.linalg.inv(system)
        standard_error = np.sqrt(np.maximum(0.0, variance * np.einsum("ij,jk,ik->i", self.basis, inverse, self.basis)))
        return SmoothResult(mean, standard_error, basis.shape[1])


def matern32_kernel(
    first: np.ndarray,
    second: np.ndarray,
    length_scale: float,
    variance: float = 1.0,
) -> np.ndarray:
    """Matérn-3/2 covariance kernel."""

    if length_scale <= 0 or variance < 0:
        raise ValueError("length_scale must be positive and variance non-negative")
    distance = np.abs(np.subtract.outer(np.asarray(first, dtype=float), np.asarray(second, dtype=float)))
    scaled = np.sqrt(3.0) * distance / float(length_scale)
    return float(variance) * (1.0 + scaled) * np.exp(-scaled)


def _cosine_taper(positions: np.ndarray, limit: float) -> np.ndarray:
    distance = np.abs(np.asarray(positions, dtype=float))
    taper = np.zeros_like(distance)
    inside = distance < limit
    taper[inside] = 0.5 * (1.0 + np.cos(np.pi * distance[inside] / limit))
    return taper


def _flat_top_taper(
    positions: np.ndarray,
    inner_limit: float,
    outer_limit: float,
) -> np.ndarray:
    distance = np.abs(np.asarray(positions, dtype=float))
    if inner_limit < 0 or outer_limit <= inner_limit:
        raise ValueError("taper limits must satisfy 0 <= inner < outer")
    taper = np.ones_like(distance)
    taper[distance >= outer_limit] = 0.0
    transition = (distance > inner_limit) & (distance < outer_limit)
    phase = (distance[transition] - inner_limit) / (outer_limit - inner_limit)
    taper[transition] = 0.5 * (1.0 + np.cos(np.pi * phase))
    return taper


class SparseAdditiveGPSmoother:
    """Nyström GP-equivalent smoother with a fixed small inducing grid."""

    def __init__(
        self,
        positions: np.ndarray,
        *,
        inducing_points: int = 25,
        long_length_scale: float = 50.0,
        short_length_scale: float = 10.0,
        short_taper: float = 50.0,
        ridge: float = 1.0,
    ):
        self.positions = np.asarray(positions, dtype=np.float64)
        self.inducing = np.linspace(self.positions.min(), self.positions.max(), min(inducing_points, len(self.positions)))
        self.long_length_scale = float(long_length_scale)
        self.short_length_scale = float(short_length_scale)
        self.short_taper = float(short_taper)
        self.ridge = float(ridge)
        taper_x = _cosine_taper(self.positions, self.short_taper)
        taper_u = _cosine_taper(self.inducing, self.short_taper)
        cross = matern32_kernel(self.positions, self.inducing, self.long_length_scale)
        cross += np.outer(taper_x, taper_u) * matern32_kernel(
            self.positions, self.inducing, self.short_length_scale
        )
        inducing_covariance = matern32_kernel(self.inducing, self.inducing, self.long_length_scale)
        inducing_covariance += np.outer(taper_u, taper_u) * matern32_kernel(
            self.inducing, self.inducing, self.short_length_scale
        )
        cholesky = np.linalg.cholesky(inducing_covariance + 1e-7 * np.eye(len(self.inducing)))
        self.features = solve_triangular(cholesky, cross.T, lower=True).T

    def fit(self, values: np.ndarray, weights: np.ndarray | None = None) -> SmoothResult:
        y = np.asarray(values, dtype=np.float64)
        if y.shape != self.positions.shape:
            raise ValueError("values must match smoother positions")
        w = np.ones_like(y) if weights is None else np.asarray(weights, dtype=np.float64)
        if w.shape != y.shape or np.any(w < 0):
            raise ValueError("weights must be non-negative and match values")
        finite = np.isfinite(y) & np.isfinite(w) & (w > 0)
        if finite.sum() < 4:
            return SmoothResult(np.nan_to_num(y), np.full_like(y, np.nan), 0)
        features = self.features[finite]
        scale = w[finite] / max(float(np.mean(w[finite])), np.finfo(float).eps)
        system = features.T @ (features * scale[:, None]) + self.ridge * np.eye(features.shape[1])
        right = features.T @ (scale * y[finite])
        coefficients = np.linalg.solve(system, right)
        mean = self.features @ coefficients
        residual = y[finite] - features @ coefficients
        variance = float(np.average(np.square(residual), weights=scale))
        inverse = np.linalg.inv(system)
        standard_error = np.sqrt(
            np.maximum(0.0, variance * np.einsum("ij,jk,ik->i", self.features, inverse, self.features))
        )
        return SmoothResult(mean, standard_error, features.shape[1])


class ExactAdditiveGPSmoother:
    """Exact 201-position GP reference used to validate the sparse approximation."""

    def __init__(
        self,
        positions: np.ndarray,
        *,
        long_length_scale: float = 50.0,
        short_length_scale: float = 10.0,
        short_taper: float = 50.0,
        noise: float = 0.2,
    ):
        self.positions = np.asarray(positions, dtype=np.float64)
        taper = _cosine_taper(self.positions, short_taper)
        self.kernel = matern32_kernel(self.positions, self.positions, long_length_scale)
        self.kernel += np.outer(taper, taper) * matern32_kernel(
            self.positions, self.positions, short_length_scale
        )
        self.noise = float(noise)

    def fit(self, values: np.ndarray, weights: np.ndarray | None = None) -> SmoothResult:
        y = np.asarray(values, dtype=np.float64)
        w = np.ones_like(y) if weights is None else np.asarray(weights, dtype=np.float64)
        if y.shape != self.positions.shape or w.shape != y.shape or np.any(w < 0):
            raise ValueError("values and non-negative weights must match positions")
        finite = np.isfinite(y) & np.isfinite(w) & (w > 0)
        indexes = np.flatnonzero(finite)
        if len(indexes) < 4:
            return SmoothResult(np.nan_to_num(y), np.full_like(y, np.nan), 0)
        covariance = self.kernel[np.ix_(indexes, indexes)]
        observation_noise = self.noise / np.maximum(w[indexes] / np.mean(w[indexes]), 1e-6)
        system = covariance + np.diag(observation_noise) + 1e-8 * np.eye(len(indexes))
        cross = self.kernel[:, indexes]
        solved = np.linalg.solve(system, y[indexes])
        mean = cross @ solved
        posterior = self.kernel - cross @ np.linalg.solve(system, cross.T)
        standard_error = np.sqrt(np.maximum(0.0, np.diag(posterior)))
        return SmoothResult(mean, standard_error, len(indexes))


class FunctionalTemplateDetector:
    """Shape-only diagonal-LDA detector with spline or sparse-GP templates.

    This model is primarily an information-ceiling and transfer diagnostic. It
    deliberately excludes motif score and accessibility. Broad per-site trends
    are removed before a positive-minus-negative functional template is learned,
    and a pooled diagonal covariance converts the template to a fast matched
    filter. An optional prior template supports TF-to-family-to-global shrinkage.
    """

    def __init__(
        self,
        positions: np.ndarray,
        *,
        smoother: str = "spline",
        window_limit: float = 50.0,
        spline_penalty: float = 10.0,
        long_length_scale: float = 50.0,
        short_length_scale: float = 10.0,
        gp_ridge: float = 1.0,
        variance_shrinkage: float = 0.5,
        variance_floor: float = 0.05,
    ):
        self.positions = np.asarray(positions, dtype=np.float64)
        if self.positions.ndim != 1 or len(self.positions) < 5:
            raise ValueError("positions must be one-dimensional with at least five values")
        if window_limit <= 0 or window_limit > float(np.max(np.abs(self.positions))):
            raise ValueError("window_limit must be within the profile coordinates")
        if not 0 <= variance_shrinkage <= 1 or variance_floor <= 0:
            raise ValueError("variance shrinkage and floor are invalid")
        if smoother == "spline":
            self.smoother = PenalizedSplineSmoother(
                self.positions,
                n_basis=min(25, len(self.positions)),
                penalty=spline_penalty,
            )
        elif smoother == "gp":
            self.smoother = SparseAdditiveGPSmoother(
                self.positions,
                inducing_points=min(25, len(self.positions)),
                long_length_scale=long_length_scale,
                short_length_scale=short_length_scale,
                short_taper=window_limit,
                ridge=gp_ridge,
            )
        else:
            raise ValueError("smoother must be spline or gp")
        self.smoother_name = smoother
        self.window_limit = float(window_limit)
        self.variance_shrinkage = float(variance_shrinkage)
        self.variance_floor = float(variance_floor)

    @staticmethod
    def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
        return np.average(values, axis=0, weights=weights)

    @staticmethod
    def _effective_sample_size(weights: np.ndarray) -> float:
        total = float(np.sum(weights))
        squared = float(np.sum(np.square(weights)))
        return total * total / max(squared, np.finfo(float).eps)

    def fit(
        self,
        profiles: np.ndarray,
        labels: Iterable[int | bool],
        *,
        sample_weight: np.ndarray | None = None,
        prior_template: np.ndarray | None = None,
        prior_strength: float = 0.0,
    ) -> "FunctionalTemplateDetector":
        values = normalize_functional_profiles(profiles, self.positions)
        group = np.asarray(list(labels), dtype=bool)
        if group.shape != (len(values),) or np.sum(group) < 2 or np.sum(~group) < 2:
            raise ValueError("labels must define at least two sites in each class")
        weights = (
            np.ones(len(values), dtype=np.float64)
            if sample_weight is None
            else np.asarray(sample_weight, dtype=np.float64)
        )
        if weights.shape != (len(values),) or np.any(~np.isfinite(weights)) or np.any(weights <= 0):
            raise ValueError("sample weights must be finite and positive")
        positive_weights = weights[group]
        negative_weights = weights[~group]
        positive_mean = self._weighted_mean(values[group], positive_weights)
        negative_mean = self._weighted_mean(values[~group], negative_weights)
        positive_variance = self._weighted_mean(
            np.square(values[group] - positive_mean), positive_weights
        )
        negative_variance = self._weighted_mean(
            np.square(values[~group] - negative_mean), negative_weights
        )
        n_positive = self._effective_sample_size(positive_weights)
        n_negative = self._effective_sample_size(negative_weights)
        pooled_variance = (
            max(n_positive - 1.0, 1.0) * positive_variance
            + max(n_negative - 1.0, 1.0) * negative_variance
        ) / max(n_positive + n_negative - 2.0, 2.0)
        finite_variance = pooled_variance[np.isfinite(pooled_variance) & (pooled_variance > 0)]
        variance_center = float(np.median(finite_variance)) if len(finite_variance) else 1.0
        pooled_variance = (
            (1.0 - self.variance_shrinkage) * np.nan_to_num(pooled_variance, nan=variance_center)
            + self.variance_shrinkage * variance_center
        )
        pooled_variance = np.maximum(pooled_variance, self.variance_floor * variance_center)

        raw_template = positive_mean - negative_mean
        effective_sites = 2.0 / (1.0 / n_positive + 1.0 / n_negative)
        if prior_template is not None:
            prior = np.asarray(prior_template, dtype=np.float64)
            if prior.shape != self.positions.shape or prior_strength < 0:
                raise ValueError("prior template must match positions and strength must be non-negative")
            fraction = effective_sites / max(effective_sites + float(prior_strength), np.finfo(float).eps)
            raw_template = fraction * raw_template + (1.0 - fraction) * prior
        standard_error = np.sqrt(
            pooled_variance * (1.0 / max(n_positive, 1.0) + 1.0 / max(n_negative, 1.0))
        )
        smooth_weights = 1.0 / np.maximum(np.square(standard_error), 1e-6)
        smoothed = self.smoother.fit(raw_template, smooth_weights)
        inner = max(0.0, self.window_limit - min(10.0, self.window_limit / 2.0))
        taper = _flat_top_taper(self.positions, inner, self.window_limit)
        template = smoothed.mean * taper
        discriminant = template / pooled_variance
        norm = float(np.sqrt(max(np.dot(template, discriminant), 0.0)))
        if not np.isfinite(norm) or norm <= 1e-10:
            raise ValueError("training profiles do not define a non-zero functional template")
        discriminant /= norm
        midpoint = 0.5 * float(np.dot(positive_mean + negative_mean, discriminant))
        training_scores = values @ discriminant - midpoint
        score_scale = float(np.std(training_scores))
        if not np.isfinite(score_scale) or score_scale <= 1e-8:
            score_scale = 1.0

        self.positive_mean_ = positive_mean
        self.negative_mean_ = negative_mean
        self.raw_template_ = raw_template
        self.footprint_template_ = template
        self.template_standard_error_ = smoothed.standard_error
        self.pooled_variance_ = pooled_variance
        self.discriminant_ = discriminant
        self.midpoint_ = midpoint
        self.score_scale_ = score_scale
        self.positive_sites_ = int(np.sum(group))
        self.negative_sites_ = int(np.sum(~group))
        self.effective_sites_ = float(effective_sites)
        self.prior_strength_ = float(prior_strength)
        return self

    def decision_function(self, profiles: np.ndarray) -> np.ndarray:
        if not hasattr(self, "discriminant_"):
            raise RuntimeError("fit must be called before prediction")
        values = normalize_functional_profiles(profiles, self.positions)
        return (values @ self.discriminant_ - self.midpoint_) / self.score_scale_

    def predict_proba(self, profiles: np.ndarray) -> np.ndarray:
        return expit(self.decision_function(profiles))


class MultichannelFunctionalTemplateDetector:
    """Combine independently smoothed strand/channel templates shape-only.

    Each functional channel is detrended and scored by a
    :class:`FunctionalTemplateDetector`. A small shrinkage-LDA layer then
    combines the channel scores. This preserves channel-specific smoothness and
    avoids smoothing across artificial boundaries in a concatenated vector.
    """

    def __init__(
        self,
        positions: np.ndarray,
        *,
        smoother: str = "spline",
        window_limit: float = 50.0,
        covariance_shrinkage: float = 0.5,
        **detector_options: Any,
    ):
        if not 0 <= covariance_shrinkage <= 1:
            raise ValueError("covariance_shrinkage must be between zero and one")
        self.positions = np.asarray(positions, dtype=np.float64)
        self.smoother = smoother
        self.window_limit = float(window_limit)
        self.covariance_shrinkage = float(covariance_shrinkage)
        self.detector_options = dict(detector_options)

    @staticmethod
    def _validate_channels(profiles: np.ndarray, positions: np.ndarray) -> np.ndarray:
        values = np.asarray(profiles, dtype=np.float64)
        if values.ndim != 3 or values.shape[2] != len(positions) or values.shape[1] < 1:
            raise ValueError("profiles must have shape (sites, channels, positions)")
        return values

    def fit(
        self,
        profiles: np.ndarray,
        labels: Iterable[int | bool],
        *,
        sample_weight: np.ndarray | None = None,
    ) -> "MultichannelFunctionalTemplateDetector":
        values = self._validate_channels(profiles, self.positions)
        group = np.asarray(list(labels), dtype=bool)
        if group.shape != (len(values),) or np.sum(group) < 2 or np.sum(~group) < 2:
            raise ValueError("labels must define at least two sites in each class")
        self.channel_models_ = []
        channel_scores = []
        for channel in range(values.shape[1]):
            model = FunctionalTemplateDetector(
                self.positions,
                smoother=self.smoother,
                window_limit=self.window_limit,
                **self.detector_options,
            ).fit(values[:, channel, :], group, sample_weight=sample_weight)
            self.channel_models_.append(model)
            channel_scores.append(model.decision_function(values[:, channel, :]))
        score_matrix = np.column_stack(channel_scores)
        positive_mean = np.mean(score_matrix[group], axis=0)
        negative_mean = np.mean(score_matrix[~group], axis=0)
        positive_covariance = np.atleast_2d(np.cov(score_matrix[group], rowvar=False))
        negative_covariance = np.atleast_2d(np.cov(score_matrix[~group], rowvar=False))
        pooled = 0.5 * (positive_covariance + negative_covariance)
        diagonal = np.diag(np.diag(pooled))
        covariance = (
            (1.0 - self.covariance_shrinkage) * pooled
            + self.covariance_shrinkage * diagonal
            + 1e-4 * np.eye(values.shape[1])
        )
        discriminant = np.linalg.solve(covariance, positive_mean - negative_mean)
        norm = float(np.sqrt(max(np.dot(discriminant, covariance @ discriminant), 0.0)))
        if not np.isfinite(norm) or norm <= 1e-10:
            raise ValueError("channels do not define a non-zero discriminant")
        discriminant /= norm
        midpoint = 0.5 * float(np.dot(positive_mean + negative_mean, discriminant))
        combined = score_matrix @ discriminant - midpoint
        scale = float(np.std(combined))
        self.channel_discriminant_ = discriminant
        self.channel_midpoint_ = midpoint
        self.score_scale_ = scale if np.isfinite(scale) and scale > 1e-8 else 1.0
        self.channel_positive_mean_ = positive_mean
        self.channel_negative_mean_ = negative_mean
        self.channel_covariance_ = covariance
        return self

    def decision_function(self, profiles: np.ndarray) -> np.ndarray:
        if not hasattr(self, "channel_models_"):
            raise RuntimeError("fit must be called before prediction")
        values = self._validate_channels(profiles, self.positions)
        if values.shape[1] != len(self.channel_models_):
            raise ValueError("prediction channel count differs from fitted model")
        scores = np.column_stack(
            [
                model.decision_function(values[:, channel, :])
                for channel, model in enumerate(self.channel_models_)
            ]
        )
        return (scores @ self.channel_discriminant_ - self.channel_midpoint_) / self.score_scale_

    def predict_proba(self, profiles: np.ndarray) -> np.ndarray:
        return expit(self.decision_function(profiles))


def site_accessibility_background(
    observed_profiles: np.ndarray,
    expected_profiles: np.ndarray,
    positions: np.ndarray,
    *,
    method: str = "none",
    exclusion: float = 50.0,
    ridge: float = 10.0,
    length_scale: float = 80.0,
    inducing_points: int = 9,
    pseudocount: float = 0.5,
) -> np.ndarray:
    """Add a broad, site-specific accessibility trend to sequence bias.

    The trend is learned only outside ``+/- exclusion`` and extrapolated
    through the motif center. This accounts for position within a broad ATAC
    peak without fitting the central protection/shoulder signal that the
    detector is intended to discover. Adjusted profiles are normalized back
    to their observed totals, retaining the conditional-profile formulation.
    """

    observed = _validate_profiles(observed_profiles, nonnegative=True)
    expected = _validate_profiles(expected_profiles, nonnegative=True)
    x = np.asarray(positions, dtype=np.float64)
    if observed.shape != expected.shape or x.shape != (observed.shape[1],):
        raise ValueError("observed/expected profiles and positions must agree")
    name = str(method).lower().replace("_", "-")
    if name not in {"none", "linear", "quadratic", "gp-long"}:
        raise ValueError("method must be none, linear, quadratic, or gp-long")

    totals = observed.sum(axis=1)
    expected_total = expected.sum(axis=1)
    uniform = np.full_like(expected, 1.0 / expected.shape[1])
    probabilities = np.divide(
        expected,
        expected_total[:, None],
        out=uniform,
        where=expected_total[:, None] > 0,
    )
    baseline = probabilities * totals[:, None]
    if name == "none":
        return baseline
    if exclusion <= 0 or ridge < 0 or length_scale <= 0 or inducing_points < 3:
        raise ValueError("background parameters must be positive")
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive")
    outer = np.abs(x) >= float(exclusion)
    if outer.sum() < 8:
        raise ValueError("background exclusion leaves fewer than eight outer positions")

    log_ratio = np.log((observed + pseudocount) / (baseline + pseudocount))
    scaled_x = x / max(float(np.max(np.abs(x))), 1.0)
    if name == "linear":
        design = np.column_stack([np.ones(len(x)), scaled_x])
        penalty = np.diag([0.0, float(ridge)])
    elif name == "quadratic":
        design = np.column_stack([np.ones(len(x)), scaled_x, np.square(scaled_x)])
        penalty = np.diag([0.0, float(ridge), float(ridge)])
    else:
        inducing = np.linspace(x.min(), x.max(), min(int(inducing_points), len(x)))
        inducing_covariance = matern32_kernel(inducing, inducing, length_scale)
        cholesky = np.linalg.cholesky(
            inducing_covariance + 1e-7 * np.eye(len(inducing))
        )
        cross = matern32_kernel(x, inducing, length_scale)
        features = solve_triangular(cholesky, cross.T, lower=True).T
        design = np.column_stack([np.ones(len(x)), features])
        penalty = np.diag([0.0] + [float(ridge)] * features.shape[1])

    outer_design = design[outer]
    system = outer_design.T @ outer_design + penalty + 1e-8 * np.eye(design.shape[1])
    projection = outer_design @ np.linalg.inv(system)
    coefficients = log_ratio[:, outer] @ projection
    trend = coefficients @ design.T
    # Low-depth profiles receive stronger shrinkage toward sequence bias.
    trend *= (totals / (totals + 50.0))[:, None]
    adjusted = baseline * np.exp(np.clip(trend, -2.0, 2.0))
    adjusted_total = adjusted.sum(axis=1)
    return np.divide(
        adjusted * totals[:, None],
        adjusted_total[:, None],
        out=baseline.copy(),
        where=adjusted_total[:, None] > 0,
    )


def _standardize(
    values: np.ndarray | None,
    length: int,
    *,
    location: float | None = None,
    scale: float | None = None,
) -> tuple[np.ndarray, float, float]:
    if values is None:
        return np.zeros(length, dtype=float), 0.0, 1.0
    array = np.asarray(values, dtype=float)
    if array.shape != (length,):
        raise ValueError("site covariates must contain one value per profile")
    finite = np.isfinite(array)
    fill = np.median(array[finite]) if finite.any() else 0.0
    array = np.where(finite, array, fill)
    fitted_location = float(np.mean(array)) if location is None else float(location)
    fitted_scale = float(np.std(array)) if scale is None else float(scale)
    fitted_scale = fitted_scale if fitted_scale > 0 else 1.0
    return (array - fitted_location) / fitted_scale, fitted_location, fitted_scale


def _fit_fractional_logistic(
    design: np.ndarray,
    response: np.ndarray,
    initial: np.ndarray,
    penalty: float = 1.0,
    iterations: int = 20,
    nonnegative_indexes: tuple[int, ...] = (),
) -> np.ndarray:
    coefficients = initial.copy()
    regularizer = np.eye(design.shape[1]) * penalty
    regularizer[0, 0] = 1e-8
    for _ in range(iterations):
        probabilities = expit(design @ coefficients)
        weights = np.maximum(probabilities * (1.0 - probabilities), 1e-5)
        adjusted = design @ coefficients + (response - probabilities) / weights
        system = design.T @ (design * weights[:, None]) + regularizer
        update = np.linalg.solve(system, design.T @ (weights * adjusted))
        if nonnegative_indexes:
            update[np.asarray(nonnegative_indexes, dtype=int)] = np.maximum(
                update[np.asarray(nonnegative_indexes, dtype=int)],
                0.0,
            )
        if np.max(np.abs(update - coefficients)) < 1e-7:
            coefficients = update
            break
        coefficients = update
    return coefficients


@dataclass
class FunctionalMixtureResult:
    posterior: np.ndarray
    footprint_profile: np.ndarray
    standard_error: np.ndarray
    prior_coefficients: np.ndarray
    converged: bool
    iterations: int
    log_likelihood: float
    descriptors: ProfileDescriptors


class BiasAwareFunctionalMixture:
    """Label-free bound/unbound mixture using expected Tn5 cuts as background."""

    def __init__(
        self,
        positions: np.ndarray,
        *,
        smoother: str = "spline",
        dispersion: float = 0.1,
        max_iter: int = 100,
        tolerance: float = 1e-5,
        shrinkage: float = 50.0,
        long_length_scale: float = 50.0,
        short_length_scale: float = 10.0,
        spline_penalty: float = 10.0,
        inducing_points: int = 25,
        gp_ridge: float = 1.0,
        accessibility_background: str = "none",
        background_exclusion: float = 50.0,
        background_ridge: float = 10.0,
        background_length_scale: float = 80.0,
        prior_constraint: str = "none",
        profile_inner_limit: float = 40.0,
        profile_outer_limit: float | None = None,
        likelihood_limit: float | None = None,
    ):
        self.positions = np.asarray(positions, dtype=float)
        self.smoother_name = str(smoother)
        self.dispersion = float(dispersion)
        self.max_iter = int(max_iter)
        self.tolerance = float(tolerance)
        self.shrinkage = float(shrinkage)
        self.long_length_scale = float(long_length_scale)
        self.short_length_scale = float(short_length_scale)
        self.spline_penalty = float(spline_penalty)
        self.inducing_points = int(inducing_points)
        self.gp_ridge = float(gp_ridge)
        self.accessibility_background = str(accessibility_background)
        self.background_exclusion = float(background_exclusion)
        self.background_ridge = float(background_ridge)
        self.background_length_scale = float(background_length_scale)
        self.prior_constraint = str(prior_constraint)
        if self.prior_constraint not in {"none", "motif", "motif-accessibility"}:
            raise ValueError(
                "prior_constraint must be none, motif, or motif-accessibility"
            )
        self.profile_inner_limit = float(profile_inner_limit)
        self.profile_outer_limit = (
            None if profile_outer_limit is None else float(profile_outer_limit)
        )
        self.likelihood_limit = None if likelihood_limit is None else float(likelihood_limit)
        if self.profile_outer_limit is not None:
            self.profile_taper = _flat_top_taper(
                self.positions,
                self.profile_inner_limit,
                self.profile_outer_limit,
            )
        else:
            self.profile_taper = np.ones_like(self.positions)
        if self.likelihood_limit is not None and self.likelihood_limit <= 0:
            raise ValueError("likelihood_limit must be positive")
        self.likelihood_mask = (
            np.ones(len(self.positions), dtype=bool)
            if self.likelihood_limit is None
            else np.abs(self.positions) <= self.likelihood_limit
        )
        if self.likelihood_mask.sum() < 5:
            raise ValueError("likelihood_limit leaves fewer than five positions")
        if smoother == "spline":
            self.smoother = PenalizedSplineSmoother(
                self.positions,
                penalty=self.spline_penalty,
            )
        elif smoother == "gp":
            self.smoother = SparseAdditiveGPSmoother(
                self.positions,
                inducing_points=self.inducing_points,
                long_length_scale=long_length_scale,
                short_length_scale=short_length_scale,
                ridge=self.gp_ridge,
            )
        elif smoother == "exact-gp":
            self.smoother = ExactAdditiveGPSmoother(
                self.positions,
                long_length_scale=long_length_scale,
                short_length_scale=short_length_scale,
            )
        else:
            raise ValueError("smoother must be spline, gp, or exact-gp")
        self.result_: FunctionalMixtureResult | None = None
        self.motif_location_: float = 0.0
        self.motif_scale_: float = 1.0
        self.accessibility_location_: float = 0.0
        self.accessibility_scale_: float = 1.0

    def _background(self, observed: np.ndarray, expected: np.ndarray) -> np.ndarray:
        return site_accessibility_background(
            observed,
            np.maximum(expected, 0.0),
            self.positions,
            method=self.accessibility_background,
            exclusion=self.background_exclusion,
            ridge=self.background_ridge,
            length_scale=self.background_length_scale,
        )

    def _bound_mean(self, background: np.ndarray, profile: np.ndarray) -> np.ndarray:
        weighted = background * np.exp(np.clip(profile, -8.0, 8.0))[None, :]
        total = background.sum(axis=1)
        weighted_total = weighted.sum(axis=1)
        return np.divide(
            weighted * total[:, None],
            weighted_total[:, None],
            out=background.copy(),
            where=weighted_total[:, None] > 0,
        )

    def fit(
        self,
        observed_profiles: np.ndarray,
        expected_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
        prior_profile: np.ndarray | None = None,
    ) -> FunctionalMixtureResult:
        observed = _validate_profiles(observed_profiles, nonnegative=True)
        expected = _validate_profiles(expected_profiles, nonnegative=True)
        if observed.shape != expected.shape or observed.shape[1] != len(self.positions):
            raise ValueError("observed/expected profiles and positions must agree")
        observed = np.nan_to_num(observed, nan=0.0, posinf=0.0, neginf=0.0)
        expected = np.nan_to_num(expected, nan=0.0, posinf=0.0, neginf=0.0)
        background = self._background(observed, expected)
        motif, self.motif_location_, self.motif_scale_ = _standardize(motif_score, len(observed))
        access, self.accessibility_location_, self.accessibility_scale_ = _standardize(
            accessibility, len(observed)
        )
        design = np.column_stack([np.ones(len(observed)), motif, access])
        prior_coefficients = np.zeros(3, dtype=float)

        if prior_profile is None:
            center = np.exp(-0.5 * np.square(self.positions / 7.0))
            left = np.exp(-0.5 * np.square((self.positions + 20.0) / 7.0))
            right = np.exp(-0.5 * np.square((self.positions - 20.0) / 7.0))
            profile = -0.35 * center + 0.08 * (left + right)
            prior = np.zeros_like(profile)
            prior_weight = 0.0
        else:
            prior = np.asarray(prior_profile, dtype=float)
            if prior.shape != self.positions.shape:
                raise ValueError("prior_profile must match positions")
            profile = prior.copy()
            prior_weight = self.shrinkage
        profile = profile * self.profile_taper

        previous = -np.inf
        converged = False
        posterior = np.full(len(observed), 0.5, dtype=float)
        smooth = SmoothResult(profile, np.full_like(profile, np.nan), 0)
        for iteration in range(1, self.max_iter + 1):
            bound = self._bound_mean(background, profile)
            unbound_ll = negative_binomial_log_likelihood(
                observed[:, self.likelihood_mask],
                background[:, self.likelihood_mask] + 1e-8,
                self.dispersion,
            ).sum(axis=1)
            bound_ll = negative_binomial_log_likelihood(
                observed[:, self.likelihood_mask],
                bound[:, self.likelihood_mask] + 1e-8,
                self.dispersion,
            ).sum(axis=1)
            log_prior = design @ prior_coefficients
            posterior = expit(np.clip(bound_ll - unbound_ll + log_prior, -40.0, 40.0))
            posterior = np.clip(posterior, 1e-5, 1.0 - 1e-5)
            nonnegative_indexes = {
                "none": (),
                "motif": (1,),
                "motif-accessibility": (1, 2),
            }[self.prior_constraint]
            prior_coefficients = _fit_fractional_logistic(
                design,
                posterior,
                prior_coefficients,
                penalty=1.0,
                nonnegative_indexes=nonnegative_indexes,
            )

            numerator = np.sum(posterior[:, None] * observed, axis=0)
            denominator = np.sum(posterior[:, None] * background, axis=0)
            target = np.log((numerator + 0.5) / (denominator + 0.5))
            weights = denominator + 1.0
            if prior_weight > 0:
                target = (weights * target + prior_weight * prior) / (weights + prior_weight)
                weights = weights + prior_weight
            smooth = self.smoother.fit(target, weights)
            profile = smooth.mean * self.profile_taper
            smooth = SmoothResult(
                profile,
                smooth.standard_error * self.profile_taper,
                smooth.effective_parameters,
            )

            log_likelihood = float(
                np.sum(
                    logsumexp(
                        np.column_stack(
                            [
                                unbound_ll - np.logaddexp(0.0, log_prior),
                                bound_ll - np.logaddexp(0.0, -log_prior),
                            ]
                        ),
                        axis=1,
                    )
                )
            )
            if np.isfinite(previous) and abs(log_likelihood - previous) <= self.tolerance * (
                1.0 + abs(previous)
            ):
                converged = True
                break
            previous = log_likelihood

        result = FunctionalMixtureResult(
            posterior=posterior,
            footprint_profile=profile,
            standard_error=smooth.standard_error,
            prior_coefficients=prior_coefficients,
            converged=converged,
            iterations=iteration,
            log_likelihood=log_likelihood,
            descriptors=profile_descriptors(profile, self.positions),
        )
        self.result_ = result
        return result

    def predict(
        self,
        observed_profiles: np.ndarray,
        expected_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
    ) -> np.ndarray:
        likelihood_ratio, log_prior = self.predict_log_odds_components(
            observed_profiles,
            expected_profiles,
            motif_score=motif_score,
            accessibility=accessibility,
        )
        return expit(np.clip(likelihood_ratio + log_prior, -40.0, 40.0))

    def predict_log_odds_components(
        self,
        observed_profiles: np.ndarray,
        expected_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return profile-likelihood and prior contributions separately."""

        if self.result_ is None:
            raise ValueError("functional mixture has not been fitted")
        observed = _validate_profiles(observed_profiles, nonnegative=True)
        expected = _validate_profiles(expected_profiles, nonnegative=True)
        if observed.shape != expected.shape or observed.shape[1] != len(self.positions):
            raise ValueError("observed/expected profiles and positions must agree")
        background = self._background(np.nan_to_num(observed), np.nan_to_num(expected))
        bound = self._bound_mean(background, self.result_.footprint_profile)
        unbound_ll = negative_binomial_log_likelihood(
            observed[:, self.likelihood_mask],
            background[:, self.likelihood_mask] + 1e-8,
            self.dispersion,
        ).sum(axis=1)
        bound_ll = negative_binomial_log_likelihood(
            observed[:, self.likelihood_mask],
            bound[:, self.likelihood_mask] + 1e-8,
            self.dispersion,
        ).sum(axis=1)
        motif, _location, _scale = _standardize(
            motif_score,
            len(observed),
            location=self.motif_location_,
            scale=self.motif_scale_,
        )
        access, _location, _scale = _standardize(
            accessibility,
            len(observed),
            location=self.accessibility_location_,
            scale=self.accessibility_scale_,
        )
        design = np.column_stack([np.ones(len(observed)), motif, access])
        return bound_ll - unbound_ll, design @ self.result_.prior_coefficients

    def save(self, path: str | Path, metadata: dict[str, Any] | None = None) -> tuple[Path, Path]:
        if self.result_ is None:
            raise ValueError("functional mixture has not been fitted")
        npz_path = Path(path).with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            positions=self.positions,
            footprint_profile=self.result_.footprint_profile,
            standard_error=self.result_.standard_error,
            prior_coefficients=self.result_.prior_coefficients,
        )
        document = {
            "schema": FUNCTIONAL_SCHEMA,
            "model_type": "bias_aware_mixture",
            "npz_sha256": _sha256_file(npz_path),
            "smoother": self.smoother_name,
            "dispersion": self.dispersion,
            "max_iter": self.max_iter,
            "tolerance": self.tolerance,
            "shrinkage": self.shrinkage,
            "long_length_scale": self.long_length_scale,
            "short_length_scale": self.short_length_scale,
            "spline_penalty": self.spline_penalty,
            "inducing_points": self.inducing_points,
            "gp_ridge": self.gp_ridge,
            "accessibility_background": self.accessibility_background,
            "background_exclusion": self.background_exclusion,
            "background_ridge": self.background_ridge,
            "background_length_scale": self.background_length_scale,
            "prior_constraint": self.prior_constraint,
            "profile_inner_limit": self.profile_inner_limit,
            "profile_outer_limit": self.profile_outer_limit,
            "likelihood_limit": self.likelihood_limit,
            "motif_location": self.motif_location_,
            "motif_scale": self.motif_scale_,
            "accessibility_location": self.accessibility_location_,
            "accessibility_scale": self.accessibility_scale_,
            "converged": self.result_.converged,
            "iterations": self.result_.iterations,
            "descriptors": asdict(self.result_.descriptors),
            "metadata": dict(metadata or {}),
        }
        json_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "BiasAwareFunctionalMixture":
        npz_path = Path(path).with_suffix(".npz")
        document = json.loads(npz_path.with_suffix(".json").read_text(encoding="utf-8"))
        if document.get("schema") != FUNCTIONAL_SCHEMA or document.get("model_type") != "bias_aware_mixture":
            raise ValueError("unsupported functional mixture model")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("functional mixture checksum does not match its metadata")
        with np.load(npz_path, allow_pickle=False) as arrays:
            positions = np.asarray(arrays["positions"], dtype=np.float64)
            model = cls(
                positions,
                smoother=str(document["smoother"]),
                dispersion=float(document["dispersion"]),
                max_iter=int(document["max_iter"]),
                tolerance=float(document["tolerance"]),
                shrinkage=float(document["shrinkage"]),
                long_length_scale=float(document.get("long_length_scale", 50.0)),
                short_length_scale=float(document.get("short_length_scale", 10.0)),
                spline_penalty=float(document.get("spline_penalty", 10.0)),
                inducing_points=int(document.get("inducing_points", 25)),
                gp_ridge=float(document.get("gp_ridge", 1.0)),
                accessibility_background=str(document.get("accessibility_background", "none")),
                background_exclusion=float(document.get("background_exclusion", 50.0)),
                background_ridge=float(document.get("background_ridge", 10.0)),
                background_length_scale=float(document.get("background_length_scale", 80.0)),
                prior_constraint=str(document.get("prior_constraint", "none")),
                profile_inner_limit=float(document.get("profile_inner_limit", 40.0)),
                profile_outer_limit=(
                    None
                    if document.get("profile_outer_limit") is None
                    else float(document["profile_outer_limit"])
                ),
                likelihood_limit=(
                    None
                    if document.get("likelihood_limit") is None
                    else float(document["likelihood_limit"])
                ),
            )
            profile = np.asarray(arrays["footprint_profile"], dtype=np.float64)
            standard_error = np.asarray(arrays["standard_error"], dtype=np.float64)
            prior_coefficients = np.asarray(arrays["prior_coefficients"], dtype=np.float64)
        descriptors = profile_descriptors(profile, positions)
        model.result_ = FunctionalMixtureResult(
            posterior=np.array([], dtype=float),
            footprint_profile=profile,
            standard_error=standard_error,
            prior_coefficients=prior_coefficients,
            converged=bool(document.get("converged", False)),
            iterations=int(document.get("iterations", 0)),
            log_likelihood=float("nan"),
            descriptors=descriptors,
        )
        model.motif_location_ = float(document.get("motif_location", 0.0))
        model.motif_scale_ = float(document.get("motif_scale", 1.0))
        model.accessibility_location_ = float(document.get("accessibility_location", 0.0))
        model.accessibility_scale_ = float(document.get("accessibility_scale", 1.0))
        return model


class ConditionalMultinomialMixture(BiasAwareFunctionalMixture):
    """Label-free footprint mixture using a conditional profile likelihood.

    The expected Tn5 signal defines the unbound positional probabilities.  A
    smooth multiplicative footprint modifies those probabilities for the
    bound state.  Conditioning on the cut total prevents the profile head from
    learning library-size or local-accessibility differences; those can be
    evaluated independently as a separate count component.
    """

    def __init__(self, positions: np.ndarray, **kwargs: Any):
        super().__init__(positions, **kwargs)
        self.evidence_temperature_: float = 1.0

    def _conditional_probabilities(
        self,
        background: np.ndarray,
        profile: np.ndarray | None = None,
    ) -> np.ndarray:
        values = np.asarray(background, dtype=float)[:, self.likelihood_mask]
        log_probability = np.log(np.maximum(values, 1e-8))
        if profile is not None:
            footprint = np.asarray(profile, dtype=float)
            if footprint.shape != self.positions.shape:
                raise ValueError("footprint profile must match positions")
            log_probability = log_probability + footprint[self.likelihood_mask][
                None, :
            ]
        log_probability -= logsumexp(log_probability, axis=1, keepdims=True)
        return np.exp(log_probability)

    @staticmethod
    def _profile_log_likelihood(
        observed: np.ndarray,
        probability: np.ndarray,
    ) -> np.ndarray:
        counts = np.asarray(observed, dtype=float)
        probabilities = np.asarray(probability, dtype=float)
        if counts.shape != probabilities.shape:
            raise ValueError("observed counts and profile probabilities must agree")
        return np.sum(counts * np.log(np.maximum(probabilities, 1e-300)), axis=1)

    def fit(
        self,
        observed_profiles: np.ndarray,
        expected_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
        prior_profile: np.ndarray | None = None,
    ) -> FunctionalMixtureResult:
        observed = _validate_profiles(observed_profiles, nonnegative=True)
        expected = _validate_profiles(expected_profiles, nonnegative=True)
        if observed.shape != expected.shape or observed.shape[1] != len(self.positions):
            raise ValueError("observed/expected profiles and positions must agree")
        observed = np.nan_to_num(observed, nan=0.0, posinf=0.0, neginf=0.0)
        expected = np.nan_to_num(expected, nan=0.0, posinf=0.0, neginf=0.0)
        background = self._background(observed, expected)
        observed_window = observed[:, self.likelihood_mask]
        totals = observed_window.sum(axis=1)
        unbound_probability = self._conditional_probabilities(background)
        motif, self.motif_location_, self.motif_scale_ = _standardize(
            motif_score, len(observed)
        )
        access, self.accessibility_location_, self.accessibility_scale_ = _standardize(
            accessibility, len(observed)
        )
        design = np.column_stack([np.ones(len(observed)), motif, access])
        prior_coefficients = np.zeros(3, dtype=float)

        if prior_profile is None:
            center = np.exp(-0.5 * np.square(self.positions / 7.0))
            left = np.exp(-0.5 * np.square((self.positions + 20.0) / 7.0))
            right = np.exp(-0.5 * np.square((self.positions - 20.0) / 7.0))
            profile = -0.35 * center + 0.08 * (left + right)
            prior = np.zeros_like(profile)
            prior_weight = 0.0
        else:
            prior = np.asarray(prior_profile, dtype=float)
            if prior.shape != self.positions.shape:
                raise ValueError("prior_profile must match positions")
            profile = prior.copy()
            prior_weight = self.shrinkage
        profile = profile * self.profile_taper

        previous = -np.inf
        converged = False
        posterior = np.full(len(observed), 0.5, dtype=float)
        smooth = SmoothResult(profile, np.full_like(profile, np.nan), 0)
        for iteration in range(1, self.max_iter + 1):
            bound_probability = self._conditional_probabilities(background, profile)
            unbound_ll = self._profile_log_likelihood(
                observed_window, unbound_probability
            )
            bound_ll = self._profile_log_likelihood(
                observed_window, bound_probability
            )
            log_prior = design @ prior_coefficients
            posterior = expit(np.clip(bound_ll - unbound_ll + log_prior, -40.0, 40.0))
            posterior = np.clip(posterior, 1e-5, 1.0 - 1e-5)
            nonnegative_indexes = {
                "none": (),
                "motif": (1,),
                "motif-accessibility": (1, 2),
            }[self.prior_constraint]
            prior_coefficients = _fit_fractional_logistic(
                design,
                posterior,
                prior_coefficients,
                penalty=1.0,
                nonnegative_indexes=nonnegative_indexes,
            )

            weighted_observed = np.sum(
                posterior[:, None] * observed_window,
                axis=0,
            )
            weighted_unbound = np.sum(
                posterior[:, None] * totals[:, None] * unbound_probability,
                axis=0,
            )
            target_window = np.log(
                (weighted_observed + 0.5) / (weighted_unbound + 0.5)
            )
            target = np.zeros_like(self.positions)
            weights = np.full_like(self.positions, 1e-6)
            target[self.likelihood_mask] = target_window
            weights[self.likelihood_mask] = weighted_unbound + 1.0
            if prior_weight > 0:
                target = (weights * target + prior_weight * prior) / (
                    weights + prior_weight
                )
                weights = weights + prior_weight
            smooth = self.smoother.fit(target, weights)
            profile = smooth.mean * self.profile_taper
            smooth = SmoothResult(
                profile,
                smooth.standard_error * self.profile_taper,
                smooth.effective_parameters,
            )

            log_likelihood = float(
                np.sum(
                    logsumexp(
                        np.column_stack(
                            [
                                unbound_ll - np.logaddexp(0.0, log_prior),
                                bound_ll - np.logaddexp(0.0, -log_prior),
                            ]
                        ),
                        axis=1,
                    )
                )
            )
            if np.isfinite(previous) and abs(log_likelihood - previous) <= self.tolerance * (
                1.0 + abs(previous)
            ):
                converged = True
                break
            previous = log_likelihood

        raw_log_ratio = bound_ll - unbound_ll
        finite = np.abs(raw_log_ratio[np.isfinite(raw_log_ratio)])
        robust_range = float(np.quantile(finite, 0.95)) if len(finite) else 1.0
        self.evidence_temperature_ = max(1.0, robust_range / 8.0)
        result = FunctionalMixtureResult(
            posterior=expit(
                np.clip(
                    raw_log_ratio / self.evidence_temperature_
                    + design @ prior_coefficients,
                    -40.0,
                    40.0,
                )
            ),
            footprint_profile=profile,
            standard_error=smooth.standard_error,
            prior_coefficients=prior_coefficients,
            converged=converged,
            iterations=iteration,
            log_likelihood=log_likelihood,
            descriptors=profile_descriptors(profile, self.positions),
        )
        self.result_ = result
        return result

    def predict_log_odds_components(
        self,
        observed_profiles: np.ndarray,
        expected_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        if self.result_ is None:
            raise ValueError("conditional mixture has not been fitted")
        observed = _validate_profiles(observed_profiles, nonnegative=True)
        expected = _validate_profiles(expected_profiles, nonnegative=True)
        if observed.shape != expected.shape or observed.shape[1] != len(self.positions):
            raise ValueError("observed/expected profiles and positions must agree")
        background = self._background(np.nan_to_num(observed), np.nan_to_num(expected))
        unbound = self._conditional_probabilities(background)
        bound = self._conditional_probabilities(
            background, self.result_.footprint_profile
        )
        observed_window = np.nan_to_num(observed)[:, self.likelihood_mask]
        shape = (
            self._profile_log_likelihood(observed_window, bound)
            - self._profile_log_likelihood(observed_window, unbound)
        ) / self.evidence_temperature_
        motif, _location, _scale = _standardize(
            motif_score,
            len(observed),
            location=self.motif_location_,
            scale=self.motif_scale_,
        )
        access, _location, _scale = _standardize(
            accessibility,
            len(observed),
            location=self.accessibility_location_,
            scale=self.accessibility_scale_,
        )
        design = np.column_stack([np.ones(len(observed)), motif, access])
        return shape, design @ self.result_.prior_coefficients

    def save(
        self,
        path: str | Path,
        metadata: dict[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        npz_path, json_path = super().save(path, metadata=metadata)
        document = json.loads(json_path.read_text(encoding="utf-8"))
        document["model_type"] = "conditional_multinomial_mixture"
        document["evidence_temperature"] = self.evidence_temperature_
        json_path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "ConditionalMultinomialMixture":
        npz_path = Path(path).with_suffix(".npz")
        document = json.loads(
            npz_path.with_suffix(".json").read_text(encoding="utf-8")
        )
        if (
            document.get("schema") != FUNCTIONAL_SCHEMA
            or document.get("model_type") != "conditional_multinomial_mixture"
        ):
            raise ValueError("unsupported conditional multinomial mixture")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("conditional mixture checksum does not match its metadata")
        with np.load(npz_path, allow_pickle=False) as arrays:
            positions = np.asarray(arrays["positions"], dtype=np.float64)
            profile = np.asarray(arrays["footprint_profile"], dtype=np.float64)
            standard_error = np.asarray(arrays["standard_error"], dtype=np.float64)
            prior_coefficients = np.asarray(
                arrays["prior_coefficients"], dtype=np.float64
            )
        model = cls(
            positions,
            smoother=str(document["smoother"]),
            dispersion=float(document["dispersion"]),
            max_iter=int(document["max_iter"]),
            tolerance=float(document["tolerance"]),
            shrinkage=float(document["shrinkage"]),
            long_length_scale=float(document.get("long_length_scale", 50.0)),
            short_length_scale=float(document.get("short_length_scale", 10.0)),
            spline_penalty=float(document.get("spline_penalty", 10.0)),
            inducing_points=int(document.get("inducing_points", 25)),
            gp_ridge=float(document.get("gp_ridge", 1.0)),
            accessibility_background=str(
                document.get("accessibility_background", "none")
            ),
            background_exclusion=float(document.get("background_exclusion", 50.0)),
            background_ridge=float(document.get("background_ridge", 10.0)),
            background_length_scale=float(
                document.get("background_length_scale", 80.0)
            ),
            prior_constraint=str(document.get("prior_constraint", "none")),
            profile_inner_limit=float(document.get("profile_inner_limit", 40.0)),
            profile_outer_limit=(
                None
                if document.get("profile_outer_limit") is None
                else float(document["profile_outer_limit"])
            ),
            likelihood_limit=(
                None
                if document.get("likelihood_limit") is None
                else float(document["likelihood_limit"])
            ),
        )
        model.result_ = FunctionalMixtureResult(
            posterior=np.array([], dtype=float),
            footprint_profile=profile,
            standard_error=standard_error,
            prior_coefficients=prior_coefficients,
            converged=bool(document.get("converged", False)),
            iterations=int(document.get("iterations", 0)),
            log_likelihood=float("nan"),
            descriptors=profile_descriptors(profile, positions),
        )
        model.motif_location_ = float(document.get("motif_location", 0.0))
        model.motif_scale_ = float(document.get("motif_scale", 1.0))
        model.accessibility_location_ = float(
            document.get("accessibility_location", 0.0)
        )
        model.accessibility_scale_ = float(
            document.get("accessibility_scale", 1.0)
        )
        model.evidence_temperature_ = float(
            document.get("evidence_temperature", 1.0)
        )
        return model


def _restore_diagonal_gaussian_mixture(
    *,
    weights: np.ndarray,
    means: np.ndarray,
    covariances: np.ndarray,
    seed: int,
    converged: bool,
    iterations: int,
    lower_bound: float,
) -> GaussianMixture:
    """Restore the numeric state needed by a diagonal Gaussian mixture."""

    weights = np.asarray(weights, dtype=np.float64)
    means = np.asarray(means, dtype=np.float64)
    covariances = np.asarray(covariances, dtype=np.float64)
    if (
        weights.ndim != 1
        or means.ndim != 2
        or covariances.shape != means.shape
        or means.shape[0] != len(weights)
        or len(weights) < 2
        or np.any(weights <= 0)
        or np.any(covariances <= 0)
    ):
        raise ValueError("invalid diagonal Gaussian-mixture state")
    mixture = GaussianMixture(
        n_components=len(weights),
        covariance_type="diag",
        random_state=int(seed),
    )
    mixture.weights_ = weights / weights.sum()
    mixture.means_ = means
    mixture.covariances_ = covariances
    mixture.precisions_ = 1.0 / covariances
    mixture.precisions_cholesky_ = np.sqrt(mixture.precisions_)
    mixture.converged_ = bool(converged)
    mixture.n_iter_ = int(iterations)
    mixture.lower_bound_ = float(lower_bound)
    mixture.n_features_in_ = means.shape[1]
    return mixture


class FdaMixtureModel:
    """Two-state Gaussian mixture in coverage-weighted functional-PC space."""

    def __init__(self, *, variance_threshold: float = 0.95, max_components: int = 20, seed: int = 2026):
        self.fpca = FunctionalPCA(variance_threshold, max_components, seed)
        self.seed = int(seed)
        self.mixture: GaussianMixture | None = None
        self.binding_component_: int | None = None
        self.positions_: np.ndarray | None = None
        self.temperature_: float = 1.0

    def _component_log_probabilities(self, scores: np.ndarray) -> np.ndarray:
        if self.mixture is None:
            raise ValueError("FDA mixture has not been fitted")
        covariance = np.maximum(self.mixture.covariances_, 1e-8)
        difference = scores[:, None, :] - self.mixture.means_[None, :, :]
        log_density = -0.5 * np.sum(
            np.square(difference) / covariance[None, :, :]
            + np.log(2.0 * np.pi * covariance)[None, :, :],
            axis=2,
        )
        return log_density + np.log(np.maximum(self.mixture.weights_, 1e-12))[None, :]

    def fit(
        self,
        residual_profiles: np.ndarray,
        *,
        positions: np.ndarray | None = None,
        sample_weight: np.ndarray | None = None,
    ) -> "FdaMixtureModel":
        values = _validate_profiles(residual_profiles)
        x = np.asarray(positions, dtype=float) if positions is not None else np.arange(values.shape[1]) - values.shape[1] // 2
        values = normalize_functional_profiles(values, x)
        scores = self.fpca.fit_transform(values, sample_weight=sample_weight)
        midpoint = values.shape[1] // 2
        center = values[:, midpoint - 5:midpoint + 6].mean(axis=1)
        flanks = np.concatenate(
            [values[:, midpoint - 40:midpoint - 15], values[:, midpoint + 16:midpoint + 41]],
            axis=1,
        ).mean(axis=1)
        shape_score = flanks - center
        lower = scores[shape_score <= np.quantile(shape_score, 0.35)]
        upper = scores[shape_score >= np.quantile(shape_score, 0.65)]
        means_init = None
        if len(lower) and len(upper):
            means_init = np.vstack([np.mean(lower, axis=0), np.mean(upper, axis=0)])
        self.mixture = GaussianMixture(
            n_components=2,
            covariance_type="diag",
            reg_covar=1e-5,
            n_init=5,
            random_state=self.seed,
            means_init=means_init,
        ).fit(scores)
        component_profiles = self.fpca.inverse_transform(self.mixture.means_)
        descriptors = [profile_descriptors(profile, x) for profile in component_profiles]
        self.binding_component_ = int(np.argmax([item.depletion for item in descriptors]))
        self.positions_ = x
        other = 1 - self.binding_component_
        log_probabilities = self._component_log_probabilities(scores)
        log_ratio = log_probabilities[:, self.binding_component_] - log_probabilities[:, other]
        robust_range = float(np.quantile(np.abs(log_ratio[np.isfinite(log_ratio)]), 0.95))
        self.temperature_ = max(1.0, robust_range / 8.0)
        return self

    def predict_proba(self, residual_profiles: np.ndarray) -> np.ndarray:
        if self.mixture is None or self.binding_component_ is None:
            raise ValueError("FDA mixture has not been fitted")
        values = normalize_functional_profiles(residual_profiles, self.positions_)
        scores = self.fpca.transform(values)
        other = 1 - self.binding_component_
        log_probabilities = self._component_log_probabilities(scores)
        log_ratio = log_probabilities[:, self.binding_component_] - log_probabilities[:, other]
        return expit(np.clip(log_ratio / self.temperature_, -30.0, 30.0))

    def component_profiles(self) -> np.ndarray:
        if self.mixture is None:
            raise ValueError("FDA mixture has not been fitted")
        return self.fpca.inverse_transform(self.mixture.means_)

    def save(
        self,
        path: str | Path,
        metadata: dict[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        """Serialize a fitted FDA mixture without executable Python objects."""

        if (
            self.mixture is None
            or self.binding_component_ is None
            or self.positions_ is None
        ):
            raise ValueError("FDA mixture has not been fitted")
        self.fpca._check_fitted()
        npz_path = Path(path).with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            positions=self.positions_,
            fpca_mean=self.fpca.mean_,
            fpca_components=self.fpca.components_,
            fpca_explained_variance_ratio=self.fpca.explained_variance_ratio_,
            fpca_impute=self.fpca.impute_,
            mixture_weights=self.mixture.weights_,
            mixture_means=self.mixture.means_,
            mixture_covariances=self.mixture.covariances_,
        )
        document = {
            "schema": FUNCTIONAL_SCHEMA,
            "model_type": "fda_mixture",
            "npz_sha256": _sha256_file(npz_path),
            "variance_threshold": self.fpca.variance_threshold,
            "max_components": self.fpca.max_components,
            "seed": self.seed,
            "binding_component": self.binding_component_,
            "temperature": self.temperature_,
            "converged": bool(self.mixture.converged_),
            "iterations": int(self.mixture.n_iter_),
            "lower_bound": float(self.mixture.lower_bound_),
            "metadata": dict(metadata or {}),
        }
        json_path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "FdaMixtureModel":
        npz_path = Path(path).with_suffix(".npz")
        document = json.loads(
            npz_path.with_suffix(".json").read_text(encoding="utf-8")
        )
        if (
            document.get("schema") != FUNCTIONAL_SCHEMA
            or document.get("model_type") != "fda_mixture"
        ):
            raise ValueError("unsupported FDA mixture model")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("FDA mixture checksum does not match its metadata")
        model = cls(
            variance_threshold=float(document["variance_threshold"]),
            max_components=int(document["max_components"]),
            seed=int(document["seed"]),
        )
        with np.load(npz_path, allow_pickle=False) as arrays:
            model.positions_ = np.asarray(arrays["positions"], dtype=np.float64)
            model.fpca.mean_ = np.asarray(arrays["fpca_mean"], dtype=np.float64)
            model.fpca.components_ = np.asarray(
                arrays["fpca_components"], dtype=np.float64
            )
            model.fpca.explained_variance_ratio_ = np.asarray(
                arrays["fpca_explained_variance_ratio"], dtype=np.float64
            )
            model.fpca.impute_ = np.asarray(
                arrays["fpca_impute"], dtype=np.float64
            )
            model.mixture = _restore_diagonal_gaussian_mixture(
                weights=np.asarray(arrays["mixture_weights"], dtype=np.float64),
                means=np.asarray(arrays["mixture_means"], dtype=np.float64),
                covariances=np.asarray(
                    arrays["mixture_covariances"], dtype=np.float64
                ),
                seed=model.seed,
                converged=bool(document.get("converged", False)),
                iterations=int(document.get("iterations", 0)),
                lower_bound=float(document.get("lower_bound", float("nan"))),
            )
        model.binding_component_ = int(document["binding_component"])
        if model.binding_component_ not in range(len(model.mixture.weights_)):
            raise ValueError("invalid FDA binding component")
        model.temperature_ = float(document.get("temperature", 1.0))
        return model


class CovariateAnchoredFdaModel:
    """Functional-PC mixture whose component identity is weakly anchored.

    Motif strength and accessibility enter only as a label-free prior during
    EM. The fitted Gaussian shape likelihood can then be evaluated separately
    from that prior, which prevents an accessibility classifier from being
    mistaken for improved footprint detection.
    """

    def __init__(
        self,
        *,
        variance_threshold: float = 0.95,
        max_components: int = 20,
        anchor_strength: float = 1.0,
        covariance_shrinkage: float = 10.0,
        max_iter: int = 100,
        tolerance: float = 1e-5,
        seed: int = 2026,
    ):
        if anchor_strength < 0 or covariance_shrinkage < 0:
            raise ValueError("anchor and covariance shrinkage must be non-negative")
        self.fpca = FunctionalPCA(variance_threshold, max_components, seed)
        self.anchor_strength = float(anchor_strength)
        self.covariance_shrinkage = float(covariance_shrinkage)
        self.max_iter = int(max_iter)
        self.tolerance = float(tolerance)
        self.positions_: np.ndarray | None = None
        self.bound_mean_: np.ndarray | None = None
        self.unbound_mean_: np.ndarray | None = None
        self.bound_variance_: np.ndarray | None = None
        self.unbound_variance_: np.ndarray | None = None
        self.motif_location_: float = 0.0
        self.motif_scale_: float = 1.0
        self.accessibility_location_: float = 0.0
        self.accessibility_scale_: float = 1.0
        self.temperature_: float = 1.0
        self.converged_: bool = False
        self.iterations_: int = 0

    @staticmethod
    def _log_density(
        scores: np.ndarray,
        mean: np.ndarray,
        variance: np.ndarray,
    ) -> np.ndarray:
        return -0.5 * np.sum(
            np.square(scores - mean) / variance + np.log(2.0 * np.pi * variance),
            axis=1,
        )

    def _anchor_log_odds(
        self,
        motif_score: np.ndarray | None,
        accessibility: np.ndarray | None,
        length: int,
        *,
        fit: bool,
    ) -> np.ndarray:
        motif, motif_location, motif_scale = _standardize(
            motif_score,
            length,
            location=None if fit else self.motif_location_,
            scale=None if fit else self.motif_scale_,
        )
        access, access_location, access_scale = _standardize(
            accessibility,
            length,
            location=None if fit else self.accessibility_location_,
            scale=None if fit else self.accessibility_scale_,
        )
        if fit:
            self.motif_location_, self.motif_scale_ = motif_location, motif_scale
            self.accessibility_location_, self.accessibility_scale_ = (
                access_location,
                access_scale,
            )
        return self.anchor_strength * (motif + access) / np.sqrt(2.0)

    def fit(
        self,
        residual_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
        positions: np.ndarray | None = None,
        sample_weight: np.ndarray | None = None,
    ) -> "CovariateAnchoredFdaModel":
        raw = _validate_profiles(residual_profiles)
        x = (
            np.asarray(positions, dtype=float)
            if positions is not None
            else np.arange(raw.shape[1], dtype=float) - raw.shape[1] // 2
        )
        values = normalize_functional_profiles(raw, x)
        weights = (
            np.ones(len(values), dtype=float)
            if sample_weight is None
            else np.asarray(sample_weight, dtype=float)
        )
        if weights.shape != (len(values),) or np.any(weights < 0) or not np.any(weights > 0):
            raise ValueError("sample_weight must be non-negative with one positive value")
        scores = self.fpca.fit_transform(values, sample_weight=weights)
        anchor = self._anchor_log_odds(
            motif_score,
            accessibility,
            len(values),
            fit=True,
        )
        posterior = np.clip(expit(anchor), 0.05, 0.95)
        global_variance = np.maximum(np.var(scores, axis=0), 1e-4)
        previous = -np.inf
        for iteration in range(1, self.max_iter + 1):
            bound_weight = weights * posterior
            unbound_weight = weights * (1.0 - posterior)
            bound_total = max(float(bound_weight.sum()), np.finfo(float).eps)
            unbound_total = max(float(unbound_weight.sum()), np.finfo(float).eps)
            bound_mean = np.sum(bound_weight[:, None] * scores, axis=0) / bound_total
            unbound_mean = np.sum(unbound_weight[:, None] * scores, axis=0) / unbound_total
            bound_variance = (
                np.sum(bound_weight[:, None] * np.square(scores - bound_mean), axis=0)
                + self.covariance_shrinkage * global_variance
            ) / (bound_total + self.covariance_shrinkage)
            unbound_variance = (
                np.sum(unbound_weight[:, None] * np.square(scores - unbound_mean), axis=0)
                + self.covariance_shrinkage * global_variance
            ) / (unbound_total + self.covariance_shrinkage)
            bound_variance = np.maximum(bound_variance, 1e-4)
            unbound_variance = np.maximum(unbound_variance, 1e-4)
            bound_ll = self._log_density(scores, bound_mean, bound_variance)
            unbound_ll = self._log_density(scores, unbound_mean, unbound_variance)
            posterior = expit(np.clip(bound_ll - unbound_ll + anchor, -40.0, 40.0))
            posterior = np.clip(posterior, 1e-5, 1.0 - 1e-5)
            likelihood = float(
                np.sum(
                    weights
                    * logsumexp(
                        np.column_stack(
                            [
                                unbound_ll - np.logaddexp(0.0, anchor),
                                bound_ll - np.logaddexp(0.0, -anchor),
                            ]
                        ),
                        axis=1,
                    )
                )
            )
            if np.isfinite(previous) and abs(likelihood - previous) <= self.tolerance * (
                1.0 + abs(previous)
            ):
                self.converged_ = True
                break
            previous = likelihood

        self.positions_ = x
        self.bound_mean_ = bound_mean
        self.unbound_mean_ = unbound_mean
        self.bound_variance_ = bound_variance
        self.unbound_variance_ = unbound_variance
        self.iterations_ = iteration
        shape_log_odds = bound_ll - unbound_ll
        finite = np.abs(shape_log_odds[np.isfinite(shape_log_odds)])
        robust_range = float(np.quantile(finite, 0.95)) if len(finite) else 1.0
        self.temperature_ = max(1.0, robust_range / 8.0)
        return self

    def predict_log_odds_components(
        self,
        residual_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        if any(
            value is None
            for value in (
                self.positions_,
                self.bound_mean_,
                self.unbound_mean_,
                self.bound_variance_,
                self.unbound_variance_,
            )
        ):
            raise ValueError("anchored FDA model has not been fitted")
        values = normalize_functional_profiles(residual_profiles, self.positions_)
        scores = self.fpca.transform(values)
        bound_ll = self._log_density(scores, self.bound_mean_, self.bound_variance_)
        unbound_ll = self._log_density(scores, self.unbound_mean_, self.unbound_variance_)
        shape = (bound_ll - unbound_ll) / self.temperature_
        anchor = self._anchor_log_odds(
            motif_score,
            accessibility,
            len(values),
            fit=False,
        )
        return shape, anchor

    def predict_proba(
        self,
        residual_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
    ) -> np.ndarray:
        shape, anchor = self.predict_log_odds_components(
            residual_profiles,
            motif_score=motif_score,
            accessibility=accessibility,
        )
        return expit(np.clip(shape + anchor, -40.0, 40.0))

    def profile_difference(self) -> np.ndarray:
        if self.bound_mean_ is None or self.unbound_mean_ is None:
            raise ValueError("anchored FDA model has not been fitted")
        profiles = self.fpca.inverse_transform(
            np.vstack([self.unbound_mean_, self.bound_mean_])
        )
        return profiles[1] - profiles[0]

    def save(
        self,
        path: str | Path,
        metadata: dict[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        if any(
            value is None
            for value in (
                self.positions_,
                self.bound_mean_,
                self.unbound_mean_,
                self.bound_variance_,
                self.unbound_variance_,
            )
        ):
            raise ValueError("anchored FDA model has not been fitted")
        self.fpca._check_fitted()
        npz_path = Path(path).with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            positions=self.positions_,
            bound_mean=self.bound_mean_,
            unbound_mean=self.unbound_mean_,
            bound_variance=self.bound_variance_,
            unbound_variance=self.unbound_variance_,
            fpca_mean=self.fpca.mean_,
            fpca_components=self.fpca.components_,
            fpca_explained_variance_ratio=self.fpca.explained_variance_ratio_,
            fpca_impute=self.fpca.impute_,
        )
        document = {
            "schema": FUNCTIONAL_SCHEMA,
            "model_type": "covariate_anchored_fda",
            "npz_sha256": _sha256_file(npz_path),
            "variance_threshold": self.fpca.variance_threshold,
            "max_components": self.fpca.max_components,
            "seed": self.fpca.seed,
            "anchor_strength": self.anchor_strength,
            "covariance_shrinkage": self.covariance_shrinkage,
            "max_iter": self.max_iter,
            "tolerance": self.tolerance,
            "motif_location": self.motif_location_,
            "motif_scale": self.motif_scale_,
            "accessibility_location": self.accessibility_location_,
            "accessibility_scale": self.accessibility_scale_,
            "temperature": self.temperature_,
            "converged": self.converged_,
            "iterations": self.iterations_,
            "metadata": dict(metadata or {}),
        }
        json_path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "CovariateAnchoredFdaModel":
        npz_path = Path(path).with_suffix(".npz")
        document = json.loads(
            npz_path.with_suffix(".json").read_text(encoding="utf-8")
        )
        if (
            document.get("schema") != FUNCTIONAL_SCHEMA
            or document.get("model_type") != "covariate_anchored_fda"
        ):
            raise ValueError("unsupported covariate-anchored FDA model")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("anchored FDA checksum does not match its metadata")
        model = cls(
            variance_threshold=float(document["variance_threshold"]),
            max_components=int(document["max_components"]),
            anchor_strength=float(document["anchor_strength"]),
            covariance_shrinkage=float(document["covariance_shrinkage"]),
            max_iter=int(document["max_iter"]),
            tolerance=float(document["tolerance"]),
            seed=int(document["seed"]),
        )
        with np.load(npz_path, allow_pickle=False) as arrays:
            model.positions_ = np.asarray(arrays["positions"], dtype=np.float64)
            model.bound_mean_ = np.asarray(arrays["bound_mean"], dtype=np.float64)
            model.unbound_mean_ = np.asarray(
                arrays["unbound_mean"], dtype=np.float64
            )
            model.bound_variance_ = np.asarray(
                arrays["bound_variance"], dtype=np.float64
            )
            model.unbound_variance_ = np.asarray(
                arrays["unbound_variance"], dtype=np.float64
            )
            model.fpca.mean_ = np.asarray(arrays["fpca_mean"], dtype=np.float64)
            model.fpca.components_ = np.asarray(
                arrays["fpca_components"], dtype=np.float64
            )
            model.fpca.explained_variance_ratio_ = np.asarray(
                arrays["fpca_explained_variance_ratio"], dtype=np.float64
            )
            model.fpca.impute_ = np.asarray(
                arrays["fpca_impute"], dtype=np.float64
            )
        model.motif_location_ = float(document.get("motif_location", 0.0))
        model.motif_scale_ = float(document.get("motif_scale", 1.0))
        model.accessibility_location_ = float(
            document.get("accessibility_location", 0.0)
        )
        model.accessibility_scale_ = float(
            document.get("accessibility_scale", 1.0)
        )
        model.temperature_ = float(document.get("temperature", 1.0))
        model.converged_ = bool(document.get("converged", False))
        model.iterations_ = int(document.get("iterations", 0))
        return model


class CovariateResidualizedFdaModel:
    """Functional-PC mixture after removing label-free covariate trends.

    Accessibility and motif strength can induce broad profile modes that an
    unsupervised mixture mistakes for occupancy. This model regresses those
    covariates from the functional-PC scores before fitting the two-state
    mixture. Covariates are used only to remove their fitted contribution;
    they are not added to the binding log odds.
    """

    def __init__(
        self,
        *,
        variance_threshold: float = 0.95,
        max_components: int = 20,
        covariate_ridge: float = 10.0,
        seed: int = 2026,
    ):
        if covariate_ridge < 0:
            raise ValueError("covariate_ridge must be non-negative")
        self.fpca = FunctionalPCA(variance_threshold, max_components, seed)
        self.covariate_ridge = float(covariate_ridge)
        self.seed = int(seed)
        self.mixture: GaussianMixture | None = None
        self.binding_component_: int | None = None
        self.positions_: np.ndarray | None = None
        self.covariate_coefficients_: np.ndarray | None = None
        self.motif_location_: float = 0.0
        self.motif_scale_: float = 1.0
        self.accessibility_location_: float = 0.0
        self.accessibility_scale_: float = 1.0
        self.temperature_: float = 1.0

    def _design(
        self,
        motif_score: np.ndarray | None,
        accessibility: np.ndarray | None,
        length: int,
        *,
        fit: bool,
    ) -> np.ndarray:
        motif, motif_location, motif_scale = _standardize(
            motif_score,
            length,
            location=None if fit else self.motif_location_,
            scale=None if fit else self.motif_scale_,
        )
        log_accessibility = (
            None
            if accessibility is None
            else np.log1p(np.maximum(np.asarray(accessibility, dtype=float), 0.0))
        )
        access, access_location, access_scale = _standardize(
            log_accessibility,
            length,
            location=None if fit else self.accessibility_location_,
            scale=None if fit else self.accessibility_scale_,
        )
        if fit:
            self.motif_location_, self.motif_scale_ = motif_location, motif_scale
            self.accessibility_location_, self.accessibility_scale_ = (
                access_location,
                access_scale,
            )
        return np.column_stack([np.ones(length), motif, access])

    @staticmethod
    def _component_log_probabilities(
        scores: np.ndarray,
        mixture: GaussianMixture,
    ) -> np.ndarray:
        covariance = np.maximum(mixture.covariances_, 1e-8)
        difference = scores[:, None, :] - mixture.means_[None, :, :]
        log_density = -0.5 * np.sum(
            np.square(difference) / covariance[None, :, :]
            + np.log(2.0 * np.pi * covariance)[None, :, :],
            axis=2,
        )
        return log_density + np.log(np.maximum(mixture.weights_, 1e-12))[None, :]

    def fit(
        self,
        residual_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
        positions: np.ndarray | None = None,
        sample_weight: np.ndarray | None = None,
    ) -> "CovariateResidualizedFdaModel":
        raw = _validate_profiles(residual_profiles)
        x = (
            np.asarray(positions, dtype=float)
            if positions is not None
            else np.arange(raw.shape[1], dtype=float) - raw.shape[1] // 2
        )
        values = normalize_functional_profiles(raw, x)
        weights = (
            np.ones(len(values), dtype=float)
            if sample_weight is None
            else np.asarray(sample_weight, dtype=float)
        )
        if weights.shape != (len(values),) or np.any(weights < 0) or not np.any(weights > 0):
            raise ValueError("sample_weight must be non-negative with one positive value")
        regression_weights = weights / np.mean(weights)
        scores = self.fpca.fit_transform(values, sample_weight=weights)
        design = self._design(motif_score, accessibility, len(values), fit=True)
        penalty = np.eye(design.shape[1], dtype=float) * self.covariate_ridge
        penalty[0, 0] = 1e-8
        system = design.T @ (design * regression_weights[:, None]) + penalty
        right = design.T @ (scores * regression_weights[:, None])
        self.covariate_coefficients_ = np.linalg.solve(system, right)
        residual_scores = scores - design @ self.covariate_coefficients_

        residual_profiles_pc = residual_scores @ self.fpca.components_
        position_span = float(np.max(np.abs(x)))
        center_limit = min(5.0, max(1.0, position_span * 0.2))
        flank_start = min(15.0, max(center_limit + 1.0, position_span * 0.3))
        flank_end = min(40.0, max(flank_start + 1.0, position_span * 0.8))
        center_mask = np.abs(x) <= center_limit
        flank_mask = (np.abs(x) >= flank_start) & (np.abs(x) <= flank_end)
        if not np.any(center_mask) or not np.any(flank_mask):
            raise ValueError("positions do not span usable center and flank intervals")
        center = residual_profiles_pc[:, center_mask].mean(axis=1)
        flanks = residual_profiles_pc[:, flank_mask].mean(axis=1)
        shape_score = flanks - center
        lower = residual_scores[shape_score <= np.quantile(shape_score, 0.35)]
        upper = residual_scores[shape_score >= np.quantile(shape_score, 0.65)]
        means_init = None
        if len(lower) and len(upper):
            means_init = np.vstack([np.mean(lower, axis=0), np.mean(upper, axis=0)])
        self.mixture = GaussianMixture(
            n_components=2,
            covariance_type="diag",
            reg_covar=1e-5,
            n_init=5,
            random_state=self.seed,
            means_init=means_init,
        ).fit(residual_scores)
        component_profiles = self.fpca.inverse_transform(self.mixture.means_)
        descriptors = [profile_descriptors(profile, x) for profile in component_profiles]
        self.binding_component_ = int(np.argmax([item.depletion for item in descriptors]))
        self.positions_ = x
        other = 1 - self.binding_component_
        log_probabilities = self._component_log_probabilities(
            residual_scores, self.mixture
        )
        log_ratio = (
            log_probabilities[:, self.binding_component_] - log_probabilities[:, other]
        )
        finite = np.abs(log_ratio[np.isfinite(log_ratio)])
        robust_range = float(np.quantile(finite, 0.95)) if len(finite) else 1.0
        self.temperature_ = max(1.0, robust_range / 8.0)
        return self

    def transform_residual_scores(
        self,
        residual_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
    ) -> np.ndarray:
        if self.positions_ is None or self.covariate_coefficients_ is None:
            raise ValueError("covariate-residualized FDA model has not been fitted")
        values = normalize_functional_profiles(residual_profiles, self.positions_)
        scores = self.fpca.transform(values)
        design = self._design(motif_score, accessibility, len(values), fit=False)
        return scores - design @ self.covariate_coefficients_

    def predict_proba(
        self,
        residual_profiles: np.ndarray,
        *,
        motif_score: np.ndarray | None = None,
        accessibility: np.ndarray | None = None,
    ) -> np.ndarray:
        if self.mixture is None or self.binding_component_ is None:
            raise ValueError("covariate-residualized FDA model has not been fitted")
        scores = self.transform_residual_scores(
            residual_profiles,
            motif_score=motif_score,
            accessibility=accessibility,
        )
        other = 1 - self.binding_component_
        log_probabilities = self._component_log_probabilities(scores, self.mixture)
        log_ratio = (
            log_probabilities[:, self.binding_component_] - log_probabilities[:, other]
        )
        return expit(np.clip(log_ratio / self.temperature_, -30.0, 30.0))

    def profile_difference(self) -> np.ndarray:
        if self.mixture is None or self.binding_component_ is None:
            raise ValueError("covariate-residualized FDA model has not been fitted")
        profiles = self.fpca.inverse_transform(self.mixture.means_)
        return profiles[self.binding_component_] - profiles[1 - self.binding_component_]

    def save(
        self,
        path: str | Path,
        metadata: dict[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        if (
            self.mixture is None
            or self.binding_component_ is None
            or self.positions_ is None
            or self.covariate_coefficients_ is None
        ):
            raise ValueError("covariate-residualized FDA model has not been fitted")
        self.fpca._check_fitted()
        npz_path = Path(path).with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            positions=self.positions_,
            covariate_coefficients=self.covariate_coefficients_,
            fpca_mean=self.fpca.mean_,
            fpca_components=self.fpca.components_,
            fpca_explained_variance_ratio=self.fpca.explained_variance_ratio_,
            fpca_impute=self.fpca.impute_,
            mixture_weights=self.mixture.weights_,
            mixture_means=self.mixture.means_,
            mixture_covariances=self.mixture.covariances_,
        )
        document = {
            "schema": FUNCTIONAL_SCHEMA,
            "model_type": "covariate_residualized_fda",
            "npz_sha256": _sha256_file(npz_path),
            "variance_threshold": self.fpca.variance_threshold,
            "max_components": self.fpca.max_components,
            "covariate_ridge": self.covariate_ridge,
            "seed": self.seed,
            "binding_component": self.binding_component_,
            "motif_location": self.motif_location_,
            "motif_scale": self.motif_scale_,
            "accessibility_location": self.accessibility_location_,
            "accessibility_scale": self.accessibility_scale_,
            "temperature": self.temperature_,
            "converged": bool(self.mixture.converged_),
            "iterations": int(self.mixture.n_iter_),
            "lower_bound": float(self.mixture.lower_bound_),
            "metadata": dict(metadata or {}),
        }
        json_path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "CovariateResidualizedFdaModel":
        npz_path = Path(path).with_suffix(".npz")
        document = json.loads(
            npz_path.with_suffix(".json").read_text(encoding="utf-8")
        )
        if (
            document.get("schema") != FUNCTIONAL_SCHEMA
            or document.get("model_type") != "covariate_residualized_fda"
        ):
            raise ValueError("unsupported covariate-residualized FDA model")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("residualized FDA checksum does not match its metadata")
        model = cls(
            variance_threshold=float(document["variance_threshold"]),
            max_components=int(document["max_components"]),
            covariate_ridge=float(document["covariate_ridge"]),
            seed=int(document["seed"]),
        )
        with np.load(npz_path, allow_pickle=False) as arrays:
            model.positions_ = np.asarray(arrays["positions"], dtype=np.float64)
            model.covariate_coefficients_ = np.asarray(
                arrays["covariate_coefficients"], dtype=np.float64
            )
            model.fpca.mean_ = np.asarray(arrays["fpca_mean"], dtype=np.float64)
            model.fpca.components_ = np.asarray(
                arrays["fpca_components"], dtype=np.float64
            )
            model.fpca.explained_variance_ratio_ = np.asarray(
                arrays["fpca_explained_variance_ratio"], dtype=np.float64
            )
            model.fpca.impute_ = np.asarray(
                arrays["fpca_impute"], dtype=np.float64
            )
            model.mixture = _restore_diagonal_gaussian_mixture(
                weights=np.asarray(arrays["mixture_weights"], dtype=np.float64),
                means=np.asarray(arrays["mixture_means"], dtype=np.float64),
                covariances=np.asarray(
                    arrays["mixture_covariances"], dtype=np.float64
                ),
                seed=model.seed,
                converged=bool(document.get("converged", False)),
                iterations=int(document.get("iterations", 0)),
                lower_bound=float(document.get("lower_bound", float("nan"))),
            )
        model.binding_component_ = int(document["binding_component"])
        model.motif_location_ = float(document.get("motif_location", 0.0))
        model.motif_scale_ = float(document.get("motif_scale", 1.0))
        model.accessibility_location_ = float(
            document.get("accessibility_location", 0.0)
        )
        model.accessibility_scale_ = float(
            document.get("accessibility_scale", 1.0)
        )
        model.temperature_ = float(document.get("temperature", 1.0))
        return model


class HybridFdaGpModel:
    """FDA initialization followed by a GP-smoothed profile likelihood."""

    def __init__(
        self,
        positions: np.ndarray,
        *,
        variance_threshold: float = 0.95,
        max_components: int = 20,
        seed: int = 2026,
    ):
        self.positions = np.asarray(positions, dtype=float)
        self.fda = FdaMixtureModel(
            variance_threshold=variance_threshold,
            max_components=max_components,
            seed=seed,
        )
        self.smoother = SparseAdditiveGPSmoother(self.positions)
        self.unbound_mean_: np.ndarray | None = None
        self.bound_mean_: np.ndarray | None = None
        self.variance_: np.ndarray | None = None
        self.prior_: float | None = None
        self.temperature_: float = 1.0

    def fit(self, residual_profiles: np.ndarray, sample_weight: np.ndarray | None = None) -> "HybridFdaGpModel":
        raw_values = _validate_profiles(residual_profiles)
        self.fda.fit(raw_values, positions=self.positions, sample_weight=sample_weight)
        responsibility = self.fda.predict_proba(raw_values)
        values = normalize_functional_profiles(raw_values, self.positions)
        site_weight = (
            np.ones(len(values), dtype=float)
            if sample_weight is None
            else np.asarray(sample_weight, dtype=float)
        )
        bound_weight = site_weight * (responsibility + 1e-4)
        unbound_weight = site_weight * (1.0 - responsibility + 1e-4)
        bound = np.average(values, axis=0, weights=bound_weight)
        unbound = np.average(values, axis=0, weights=unbound_weight)
        difference = bound - unbound
        information = np.sum(bound_weight[:, None] * np.square(values - bound), axis=0)
        information += np.sum(unbound_weight[:, None] * np.square(values - unbound), axis=0)
        smoothing_weight = (len(values) + 1.0) / np.maximum(information, 1e-3)
        smooth = self.smoother.fit(difference, smoothing_weight)
        self.unbound_mean_ = unbound
        self.bound_mean_ = unbound + smooth.mean
        residual = (
            responsibility[:, None] * np.square(values - self.bound_mean_)
            + (1.0 - responsibility[:, None]) * np.square(values - self.unbound_mean_)
        )
        self.variance_ = np.maximum(np.mean(residual, axis=0), 1e-4)
        self.prior_ = float(np.clip(np.mean(responsibility), 1e-4, 1.0 - 1e-4))
        training_log_ratio = self._log_likelihood_ratio(values)
        robust_range = float(
            np.quantile(np.abs(training_log_ratio[np.isfinite(training_log_ratio)]), 0.95)
        )
        self.temperature_ = max(1.0, robust_range / 8.0)
        return self

    def _log_likelihood_ratio(self, values: np.ndarray) -> np.ndarray:
        if self.unbound_mean_ is None or self.bound_mean_ is None or self.variance_ is None or self.prior_ is None:
            raise ValueError("hybrid FDA-GP model has not been fitted")
        bound_ll = -0.5 * np.sum(
            np.square(values - self.bound_mean_) / self.variance_ + np.log(self.variance_), axis=1
        )
        unbound_ll = -0.5 * np.sum(
            np.square(values - self.unbound_mean_) / self.variance_ + np.log(self.variance_), axis=1
        )
        return bound_ll - unbound_ll + np.log(self.prior_ / (1.0 - self.prior_))

    def predict_proba(self, residual_profiles: np.ndarray) -> np.ndarray:
        if self.unbound_mean_ is None or self.bound_mean_ is None or self.variance_ is None or self.prior_ is None:
            raise ValueError("hybrid FDA-GP model has not been fitted")
        values = normalize_functional_profiles(residual_profiles, self.positions)
        log_ratio = self._log_likelihood_ratio(values)
        return expit(np.clip(log_ratio / self.temperature_, -30.0, 30.0))

    def save(
        self,
        path: str | Path,
        metadata: dict[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        """Serialize the deployable hybrid likelihood and its FDA initializer."""

        if (
            self.unbound_mean_ is None
            or self.bound_mean_ is None
            or self.variance_ is None
            or self.prior_ is None
            or self.fda.mixture is None
            or self.fda.binding_component_ is None
        ):
            raise ValueError("hybrid FDA-GP model has not been fitted")
        self.fda.fpca._check_fitted()
        npz_path = Path(path).with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            positions=self.positions,
            unbound_mean=self.unbound_mean_,
            bound_mean=self.bound_mean_,
            variance=self.variance_,
            fpca_mean=self.fda.fpca.mean_,
            fpca_components=self.fda.fpca.components_,
            fpca_explained_variance_ratio=(
                self.fda.fpca.explained_variance_ratio_
            ),
            fpca_impute=self.fda.fpca.impute_,
            mixture_weights=self.fda.mixture.weights_,
            mixture_means=self.fda.mixture.means_,
            mixture_covariances=self.fda.mixture.covariances_,
        )
        document = {
            "schema": FUNCTIONAL_SCHEMA,
            "model_type": "hybrid_fda_gp",
            "npz_sha256": _sha256_file(npz_path),
            "variance_threshold": self.fda.fpca.variance_threshold,
            "max_components": self.fda.fpca.max_components,
            "seed": self.fda.seed,
            "binding_component": self.fda.binding_component_,
            "fda_temperature": self.fda.temperature_,
            "fda_converged": bool(self.fda.mixture.converged_),
            "fda_iterations": int(self.fda.mixture.n_iter_),
            "fda_lower_bound": float(self.fda.mixture.lower_bound_),
            "prior": self.prior_,
            "temperature": self.temperature_,
            "metadata": dict(metadata or {}),
        }
        json_path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "HybridFdaGpModel":
        npz_path = Path(path).with_suffix(".npz")
        document = json.loads(
            npz_path.with_suffix(".json").read_text(encoding="utf-8")
        )
        if (
            document.get("schema") != FUNCTIONAL_SCHEMA
            or document.get("model_type") != "hybrid_fda_gp"
        ):
            raise ValueError("unsupported hybrid FDA-GP model")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("hybrid FDA-GP checksum does not match its metadata")
        with np.load(npz_path, allow_pickle=False) as arrays:
            positions = np.asarray(arrays["positions"], dtype=np.float64)
            model = cls(
                positions,
                variance_threshold=float(document["variance_threshold"]),
                max_components=int(document["max_components"]),
                seed=int(document["seed"]),
            )
            model.unbound_mean_ = np.asarray(
                arrays["unbound_mean"], dtype=np.float64
            )
            model.bound_mean_ = np.asarray(arrays["bound_mean"], dtype=np.float64)
            model.variance_ = np.asarray(arrays["variance"], dtype=np.float64)
            model.fda.positions_ = positions.copy()
            model.fda.fpca.mean_ = np.asarray(
                arrays["fpca_mean"], dtype=np.float64
            )
            model.fda.fpca.components_ = np.asarray(
                arrays["fpca_components"], dtype=np.float64
            )
            model.fda.fpca.explained_variance_ratio_ = np.asarray(
                arrays["fpca_explained_variance_ratio"], dtype=np.float64
            )
            model.fda.fpca.impute_ = np.asarray(
                arrays["fpca_impute"], dtype=np.float64
            )
            model.fda.mixture = _restore_diagonal_gaussian_mixture(
                weights=np.asarray(arrays["mixture_weights"], dtype=np.float64),
                means=np.asarray(arrays["mixture_means"], dtype=np.float64),
                covariances=np.asarray(
                    arrays["mixture_covariances"], dtype=np.float64
                ),
                seed=model.fda.seed,
                converged=bool(document.get("fda_converged", False)),
                iterations=int(document.get("fda_iterations", 0)),
                lower_bound=float(
                    document.get("fda_lower_bound", float("nan"))
                ),
            )
        model.fda.binding_component_ = int(document["binding_component"])
        model.fda.temperature_ = float(document.get("fda_temperature", 1.0))
        model.prior_ = float(document["prior"])
        if not 0.0 < model.prior_ < 1.0:
            raise ValueError("invalid hybrid FDA-GP prior")
        model.temperature_ = float(document.get("temperature", 1.0))
        return model


@dataclass(frozen=True)
class DifferentialFunctionalResult:
    positions: np.ndarray
    difference: np.ndarray
    pointwise_lower: np.ndarray
    pointwise_upper: np.ndarray
    simultaneous_lower: np.ndarray
    simultaneous_upper: np.ndarray
    global_pvalue: float
    bootstrap_iterations: int
    unit: str
    descriptor_change: ProfileDescriptors


def functional_differential_test(
    profiles: np.ndarray,
    conditions: Iterable[str],
    replicates: Iterable[str],
    contrast: tuple[str, str],
    *,
    positions: np.ndarray | None = None,
    n_bootstrap: int = 1000,
    seed: int = 2026,
) -> DifferentialFunctionalResult:
    """Replicate-level functional difference with simultaneous uncertainty."""

    values = _validate_profiles(profiles)
    condition_array = np.asarray(list(conditions), dtype=str)
    replicate_array = np.asarray(list(replicates), dtype=str)
    if condition_array.shape != (len(values),) or replicate_array.shape != (len(values),):
        raise ValueError("conditions and replicates must contain one value per profile")
    first, second = contrast
    if first == second:
        raise ValueError("contrast conditions must differ")
    curves: list[np.ndarray] = []
    curve_conditions: list[str] = []
    for (condition, replicate), indexes in _group_indexes(condition_array, replicate_array):
        if condition in contrast:
            curves.append(np.nanmean(values[indexes], axis=0))
            curve_conditions.append(condition)
    replicate_curves = np.asarray(curves, dtype=float)
    curve_conditions_array = np.asarray(curve_conditions, dtype=str)
    first_curves = replicate_curves[curve_conditions_array == first]
    second_curves = replicate_curves[curve_conditions_array == second]
    if len(first_curves) < 2 or len(second_curves) < 2:
        raise ValueError("functional differential testing requires at least two replicates per condition")
    difference = np.mean(first_curves, axis=0) - np.mean(second_curves, axis=0)
    rng = np.random.default_rng(seed)
    bootstrap = np.empty((n_bootstrap, values.shape[1]), dtype=float)
    for index in range(n_bootstrap):
        draw_first = first_curves[rng.integers(0, len(first_curves), size=len(first_curves))]
        draw_second = second_curves[rng.integers(0, len(second_curves), size=len(second_curves))]
        bootstrap[index] = np.mean(draw_first, axis=0) - np.mean(draw_second, axis=0)
    lower, upper = np.quantile(bootstrap, [0.025, 0.975], axis=0)
    standard_error = np.std(bootstrap, axis=0, ddof=1)
    finite_se = np.where(standard_error > 1e-12, standard_error, np.inf)
    maximum = np.max(np.abs((bootstrap - difference) / finite_se), axis=1)
    critical = float(np.quantile(maximum, 0.95))
    simultaneous_lower = difference - critical * np.where(np.isfinite(finite_se), standard_error, 0.0)
    simultaneous_upper = difference + critical * np.where(np.isfinite(finite_se), standard_error, 0.0)

    observed_statistic = float(np.sum(np.square(difference)))
    all_curves = np.vstack([first_curves, second_curves])
    n_first = len(first_curves)
    combinations = list(itertools.combinations(range(len(all_curves)), n_first))
    if len(combinations) > 5000:
        combinations = [tuple(rng.choice(len(all_curves), n_first, replace=False)) for _ in range(5000)]
    permutation_statistics = []
    all_indexes = np.arange(len(all_curves))
    for selected in combinations:
        selected_array = np.asarray(selected, dtype=int)
        other = np.setdiff1d(all_indexes, selected_array, assume_unique=False)
        permuted = np.mean(all_curves[selected_array], axis=0) - np.mean(all_curves[other], axis=0)
        permutation_statistics.append(float(np.sum(np.square(permuted))))
    pvalue = (1.0 + np.sum(np.asarray(permutation_statistics) >= observed_statistic)) / (
        1.0 + len(permutation_statistics)
    )
    x = np.asarray(positions, dtype=float) if positions is not None else np.arange(values.shape[1]) - values.shape[1] // 2
    return DifferentialFunctionalResult(
        positions=x,
        difference=difference,
        pointwise_lower=lower,
        pointwise_upper=upper,
        simultaneous_lower=simultaneous_lower,
        simultaneous_upper=simultaneous_upper,
        global_pvalue=float(pvalue),
        bootstrap_iterations=int(n_bootstrap),
        unit="replicate",
        descriptor_change=profile_descriptors(difference, x),
    )


def _group_indexes(
    conditions: np.ndarray,
    replicates: np.ndarray,
) -> Iterable[tuple[tuple[str, str], np.ndarray]]:
    keys = np.asarray([f"{condition}\0{replicate}" for condition, replicate in zip(conditions, replicates)])
    for key in np.unique(keys):
        condition, replicate = key.split("\0", 1)
        yield (condition, replicate), np.flatnonzero(keys == key)


def deviance_profiles(
    observed_profiles: np.ndarray,
    expected_profiles: np.ndarray,
    dispersion: float = 0.0,
) -> np.ndarray:
    """Convenience wrapper for functional signed-deviance profiles."""

    observed = _validate_profiles(observed_profiles, nonnegative=True)
    expected = _validate_profiles(expected_profiles, nonnegative=True)
    if observed.shape != expected.shape:
        raise ValueError("observed and expected profiles must match")
    return calibrated_residuals(observed, expected, "deviance", dispersion=dispersion)
