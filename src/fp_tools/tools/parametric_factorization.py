#!/usr/bin/env python
"""Frozen parametric Tn5-bias factorization research models.

This module is intentionally not wired into the public CLI.  It implements the
CPU-only research arm that separates a frozen sequence-bias propensity from a
sample calibration, broad accessibility, and a partially pooled TF footprint.
The established PWM/DWM correction remains unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, gammaln, logsumexp

from fp_tools.tools.parametric_bias import (
    calibrated_residuals,
    center_flank_likelihood_score,
    estimate_nb_dispersion,
    negative_binomial_log_likelihood,
)


CALIBRATION_SCHEMA = "fp-tools-frozen-bias-calibration-v1"
FACTORIZATION_SCHEMA = "fp-tools-parametric-factorization-v1"


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _npz_json_paths(path: str | Path) -> tuple[Path, Path]:
    npz_path = Path(path)
    if npz_path.suffix != ".npz":
        npz_path = Path(str(npz_path) + ".npz")
    return npz_path, npz_path.with_suffix(".json")


def _validate_profiles(
    counts: np.ndarray,
    logits: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    observed = np.asarray(counts, dtype=np.float64)
    scores = np.asarray(logits, dtype=np.float64)
    if observed.ndim != 2 or observed.shape != scores.shape:
        raise ValueError("counts and logits must be equal two-dimensional arrays")
    if observed.shape[0] < 1 or observed.shape[1] < 5:
        raise ValueError("at least one profile and five positions are required")
    if np.any(~np.isfinite(observed)) or np.any(observed < 0):
        raise ValueError("counts must be finite and non-negative")
    if np.any(~np.isfinite(scores)):
        raise ValueError("logits must be finite")
    return observed, scores


def conditional_profile_log_likelihood(
    counts: np.ndarray,
    logits: np.ndarray,
    *,
    include_constant: bool = False,
) -> np.ndarray:
    """Vectorized conditional multinomial log likelihood for each profile."""

    observed, scores = _validate_profiles(counts, logits)
    totals = observed.sum(axis=1)
    result = np.sum(observed * scores, axis=1) - totals * logsumexp(scores, axis=1)
    if include_constant:
        result += gammaln(totals + 1.0) - np.sum(gammaln(observed + 1.0), axis=1)
    return result


def expected_profile_counts(counts: np.ndarray, logits: np.ndarray) -> np.ndarray:
    """Scale a log-propensity profile to each profile's observed total."""

    observed, scores = _validate_profiles(counts, logits)
    probabilities = np.exp(scores - logsumexp(scores, axis=1, keepdims=True))
    return observed.sum(axis=1, keepdims=True) * probabilities


@dataclass(frozen=True)
class BiasStrengthEstimate:
    """One frozen-bias calibration fitted without occupancy labels."""

    sample: str
    bias_strength: float
    conditional_nll: float
    null_nll: float
    windows: int
    cuts: float
    lower_bound: float = 0.0
    upper_bound: float = 2.0

    @property
    def nll_gain(self) -> float:
        return self.null_nll - self.conditional_nll


class FrozenBiasStrengthCalibrator:
    """Fit one bounded bias-strength scalar per sample on background cuts."""

    def __init__(self, bounds: tuple[float, float] = (0.0, 2.0)):
        lower, upper = (float(value) for value in bounds)
        if not (0.0 <= lower < upper):
            raise ValueError("bias-strength bounds must satisfy 0 <= lower < upper")
        self.bounds = (lower, upper)
        self.estimates: dict[str, BiasStrengthEstimate] = {}
        self.metadata: dict[str, Any] = {}

    @staticmethod
    def _fit_one(
        sample: str,
        counts: np.ndarray,
        log_bias: np.ndarray,
        bounds: tuple[float, float],
    ) -> BiasStrengthEstimate:
        observed, scores = _validate_profiles(counts, log_bias)
        keep = observed.sum(axis=1) > 0
        if not np.any(keep):
            raise ValueError(f"sample {sample!r} has no cuts in calibration windows")
        observed = observed[keep]
        scores = scores[keep]
        total_cuts = float(observed.sum())

        def objective(value: float) -> float:
            ll = conditional_profile_log_likelihood(observed, float(value) * scores)
            return float(-np.sum(ll) / total_cuts)

        fitted = minimize_scalar(
            objective,
            method="bounded",
            bounds=bounds,
            options={"xatol": 1e-8, "maxiter": 256},
        )
        if not fitted.success or not np.isfinite(fitted.fun):
            raise RuntimeError(f"bias-strength optimization failed for {sample}: {fitted.message}")
        null_nll = objective(0.0)
        return BiasStrengthEstimate(
            sample=str(sample),
            bias_strength=float(fitted.x),
            conditional_nll=float(fitted.fun),
            null_nll=float(null_nll),
            windows=int(len(observed)),
            cuts=total_cuts,
            lower_bound=float(bounds[0]),
            upper_bound=float(bounds[1]),
        )

    def fit(
        self,
        counts: np.ndarray,
        log_bias: np.ndarray,
        sample_ids: Iterable[str] | None = None,
    ) -> "FrozenBiasStrengthCalibrator":
        observed, scores = _validate_profiles(counts, log_bias)
        samples = (
            np.repeat("sample", len(observed))
            if sample_ids is None
            else np.asarray(list(sample_ids), dtype=str)
        )
        if samples.shape != (len(observed),):
            raise ValueError("sample_ids must contain one value per profile")
        self.estimates = {}
        for sample in sorted(set(samples.tolist())):
            selected = samples == sample
            self.estimates[sample] = self._fit_one(
                sample,
                observed[selected],
                scores[selected],
                self.bounds,
            )
        self.metadata = {
            "training_windows": int(len(observed)),
            "training_cuts": float(observed.sum()),
            "samples": sorted(self.estimates),
        }
        return self

    def strength(self, sample: str) -> float:
        if sample not in self.estimates:
            raise KeyError(f"no frozen bias calibration for sample {sample!r}")
        return self.estimates[sample].bias_strength

    def save(
        self,
        path: str | Path,
        metadata: Mapping[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        if not self.estimates:
            raise ValueError("bias-strength calibrator has not been fitted")
        npz_path, json_path = _npz_json_paths(path)
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        names = np.asarray(sorted(self.estimates), dtype="U")
        np.savez_compressed(
            npz_path,
            sample=names,
            bias_strength=np.asarray([self.estimates[name].bias_strength for name in names]),
            conditional_nll=np.asarray([self.estimates[name].conditional_nll for name in names]),
            null_nll=np.asarray([self.estimates[name].null_nll for name in names]),
            windows=np.asarray([self.estimates[name].windows for name in names], dtype=np.int64),
            cuts=np.asarray([self.estimates[name].cuts for name in names]),
        )
        document = {
            "schema": CALIBRATION_SCHEMA,
            "npz_sha256": _sha256_file(npz_path),
            "bounds": list(self.bounds),
            "metadata": {**self.metadata, **dict(metadata or {})},
        }
        json_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "FrozenBiasStrengthCalibrator":
        npz_path, json_path = _npz_json_paths(path)
        document = json.loads(json_path.read_text(encoding="utf-8"))
        if document.get("schema") != CALIBRATION_SCHEMA:
            raise ValueError(f"unsupported calibration schema: {document.get('schema')}")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("calibration checksum does not match its JSON metadata")
        model = cls(tuple(float(value) for value in document["bounds"]))
        with np.load(npz_path, allow_pickle=False) as arrays:
            required = {"sample", "bias_strength", "conditional_nll", "null_nll", "windows", "cuts"}
            if not required.issubset(arrays.files):
                raise ValueError("calibration NPZ is missing required arrays")
            names = np.asarray(arrays["sample"], dtype=str)
            fields = {key: np.asarray(arrays[key]) for key in required - {"sample"}}
            if any(value.shape != names.shape for value in fields.values()):
                raise ValueError("calibration arrays have incompatible shapes")
            for index, name in enumerate(names.tolist()):
                model.estimates[name] = BiasStrengthEstimate(
                    sample=name,
                    bias_strength=float(fields["bias_strength"][index]),
                    conditional_nll=float(fields["conditional_nll"][index]),
                    null_nll=float(fields["null_nll"][index]),
                    windows=int(fields["windows"][index]),
                    cuts=float(fields["cuts"][index]),
                    lower_bound=model.bounds[0],
                    upper_bound=model.bounds[1],
                )
        model.metadata = dict(document.get("metadata", {}))
        return model


@dataclass(frozen=True)
class NaturalCubicBasis:
    positions: np.ndarray
    knots: np.ndarray
    matrix: np.ndarray


def natural_cubic_spline_basis(
    positions: np.ndarray,
    *,
    df: int = 5,
    knots: np.ndarray | None = None,
) -> NaturalCubicBasis:
    """Return a numerically scaled restricted-natural-cubic spline basis."""

    x = np.asarray(positions, dtype=np.float64)
    if x.ndim != 1 or len(x) < 5 or np.any(~np.isfinite(x)):
        raise ValueError("positions must be a finite one-dimensional array")
    if np.any(np.diff(x) <= 0):
        raise ValueError("positions must be strictly increasing")
    if df < 3 or df > len(x):
        raise ValueError("df must be between 3 and the number of positions")
    center = 0.5 * (x[0] + x[-1])
    scale = max(0.5 * (x[-1] - x[0]), np.finfo(float).eps)
    normalized = (x - center) / scale
    knot_values = (
        np.linspace(-1.0, 1.0, int(df), dtype=np.float64)
        if knots is None
        else np.asarray(knots, dtype=np.float64)
    )
    if knot_values.shape != (int(df),) or np.any(np.diff(knot_values) <= 0):
        raise ValueError("knots must contain df strictly increasing values")

    def truncated(values: np.ndarray, knot: float) -> np.ndarray:
        return np.maximum(values - knot, 0.0) ** 3

    last = knot_values[-1]
    penultimate = knot_values[-2]
    columns = [np.ones_like(normalized), normalized]
    reference = (truncated(normalized, penultimate) - truncated(normalized, last)) / (
        last - penultimate
    )
    for knot in knot_values[:-2]:
        term = (truncated(normalized, knot) - truncated(normalized, last)) / (last - knot)
        columns.append(term - reference)
    matrix = np.column_stack(columns)
    return NaturalCubicBasis(x.copy(), knot_values.copy(), matrix)


@dataclass(frozen=True)
class FlankAccessibilityFit:
    coefficients: np.ndarray
    log_background: np.ndarray
    expected: np.ndarray
    outer_mask: np.ndarray
    basis: NaturalCubicBasis


def fit_flank_accessibility_background(
    counts: np.ndarray,
    log_bias: np.ndarray,
    positions: np.ndarray,
    bias_strength: float | np.ndarray,
    *,
    outer_start: float = 50.0,
    outer_end: float = 100.0,
    df: int = 5,
    ridge: float = 1e-3,
    pseudocount: float = 0.5,
) -> FlankAccessibilityFit:
    """Fit broad per-site accessibility only on the two outer flanks.

    The fitted target is the log conditional cut proportion after removing the
    frozen sequence-bias offset.  The natural spline is then extrapolated into
    the protected center without observing center counts during fitting.
    """

    observed, bias = _validate_profiles(counts, log_bias)
    x = np.asarray(positions, dtype=np.float64)
    if x.shape != (observed.shape[1],):
        raise ValueError("positions must match profile width")
    if not (0 <= outer_start < outer_end):
        raise ValueError("outer flank limits must satisfy 0 <= start < end")
    outer = (np.abs(x) >= float(outer_start)) & (np.abs(x) <= float(outer_end))
    if int(outer.sum()) < df + 2:
        raise ValueError("outer flank limits leave too few positions for the spline")
    strengths = np.asarray(bias_strength, dtype=np.float64)
    if strengths.ndim == 0:
        strengths = np.repeat(float(strengths), len(observed))
    if strengths.shape != (len(observed),) or np.any(~np.isfinite(strengths)) or np.any(strengths < 0):
        raise ValueError("bias_strength must be non-negative with one value per profile")
    if ridge < 0 or pseudocount <= 0:
        raise ValueError("ridge must be non-negative and pseudocount positive")

    basis = natural_cubic_spline_basis(x, df=df)
    design = basis.matrix
    design_outer = design[outer]
    penalty = np.eye(df, dtype=np.float64)
    penalty[:2, :2] = 0.0
    totals = observed.sum(axis=1)
    denominator = totals + pseudocount * observed.shape[1]
    log_proportions = np.log((observed + pseudocount) / denominator[:, None])
    target = log_proportions - strengths[:, None] * bias
    coefficients = np.zeros((len(observed), df), dtype=np.float64)
    for row in range(len(observed)):
        weights = observed[row, outer] + pseudocount
        lhs = design_outer.T @ (weights[:, None] * design_outer) + ridge * penalty
        rhs = design_outer.T @ (weights * target[row, outer])
        coefficients[row] = np.linalg.solve(lhs + 1e-10 * np.eye(df), rhs)
    background = coefficients @ design.T
    logits = strengths[:, None] * bias + background
    expected = expected_profile_counts(observed, logits)
    return FlankAccessibilityFit(coefficients, background, expected, outer, basis)


def footprint_taper(positions: np.ndarray, *, limit: float = 50.0) -> np.ndarray:
    """Flat-top cosine taper that is exactly zero outside ``limit``."""

    x = np.abs(np.asarray(positions, dtype=np.float64))
    if limit <= 0:
        raise ValueError("footprint taper limit must be positive")
    inner = 0.8 * float(limit)
    weights = np.ones_like(x)
    transition = (x > inner) & (x < limit)
    weights[transition] = 0.5 * (
        1.0 + np.cos(np.pi * (x[transition] - inner) / (float(limit) - inner))
    )
    weights[x >= limit] = 0.0
    return weights


@dataclass(frozen=True)
class FactorizationResult:
    """Predictions from a frozen parametric footprint factorization."""

    posterior_bound: np.ndarray
    expected_unbound: np.ndarray
    expected_bound: np.ndarray
    log_background: np.ndarray
    footprint_log_effect: np.ndarray
    conditional_log_bayes_factor: np.ndarray
    total_log_bayes_factor: np.ndarray

    def residual(self, observed: np.ndarray, mode: str, *, dispersion: float = 0.0) -> np.ndarray:
        name = mode.lower().replace("_", "-")
        if name in {"nb", "nb-center-flank", "negative-binomial"}:
            return center_flank_likelihood_score(
                observed,
                self.expected_unbound,
                dispersion=dispersion,
            )
        return calibrated_residuals(
            observed,
            self.expected_unbound,
            name,
            dispersion=dispersion,
        )


class FrozenParametricFactorization:
    """Label-free mixture of frozen bias, accessibility, and TF protection.

    The sequence coefficients and per-sample bias strengths are external and
    frozen.  Site backgrounds use only positions 50--100 bp from the motif.
    TF footprint curves are fitted in a tapered spline basis, then partially
    pooled through global, family, and TF effects.
    """

    def __init__(
        self,
        positions: np.ndarray,
        *,
        background_df: int = 5,
        footprint_df: int = 15,
        footprint_limit: float = 50.0,
        outer_start: float = 50.0,
        outer_end: float = 100.0,
        background_ridge: float = 1e-3,
        footprint_ridge: float = 1e-2,
        family_shrinkage: float = 25.0,
        tf_shrinkage: float = 50.0,
        use_total_component: bool = True,
        seed: int = 2026,
    ):
        self.positions = np.asarray(positions, dtype=np.float64)
        self.background_df = int(background_df)
        self.footprint_df = int(footprint_df)
        self.footprint_limit = float(footprint_limit)
        self.outer_start = float(outer_start)
        self.outer_end = float(outer_end)
        self.background_ridge = float(background_ridge)
        self.footprint_ridge = float(footprint_ridge)
        self.family_shrinkage = float(family_shrinkage)
        self.tf_shrinkage = float(tf_shrinkage)
        self.use_total_component = bool(use_total_component)
        self.seed = int(seed)
        if self.positions.ndim != 1 or len(self.positions) < 21:
            raise ValueError("positions must contain at least 21 coordinates")
        if np.any(np.diff(self.positions) <= 0):
            raise ValueError("positions must be strictly increasing")
        if self.background_df != 5:
            raise ValueError("the frozen experiment requires a five-df accessibility spline")
        if self.footprint_df < 5 or self.footprint_df > len(self.positions):
            raise ValueError("footprint_df must be between 5 and profile width")
        if min(self.background_ridge, self.footprint_ridge, self.family_shrinkage, self.tf_shrinkage) < 0:
            raise ValueError("regularization and shrinkage values must be non-negative")

        footprint_basis = natural_cubic_spline_basis(self.positions, df=self.footprint_df)
        self.footprint_knots_ = footprint_basis.knots
        self.footprint_basis_ = footprint_basis.matrix * footprint_taper(
            self.positions,
            limit=self.footprint_limit,
        )[:, None]
        self.bias_strengths_: dict[str, float] = {}
        self.tf_family_: dict[str, str] = {}
        self.global_coefficients_: np.ndarray | None = None
        self.family_coefficients_: dict[str, np.ndarray] = {}
        self.tf_coefficients_: dict[str, np.ndarray] = {}
        self.mixing_probabilities_: dict[str, float] = {}
        self.total_unbound_means_: dict[str, float] = {}
        self.total_bound_means_: dict[str, float] = {}
        self.total_dispersion_: float = 0.0
        self.training_history: list[dict[str, float | int]] = []
        self.metadata: dict[str, Any] = {}

    @staticmethod
    def _string_vector(values: Iterable[str], length: int, name: str) -> np.ndarray:
        array = np.asarray(list(values), dtype=str)
        if array.shape != (length,):
            raise ValueError(f"{name} must contain one value per profile")
        if np.any(array == ""):
            raise ValueError(f"{name} cannot contain empty values")
        return array

    @staticmethod
    def _strength_mapping(
        calibration: FrozenBiasStrengthCalibrator | Mapping[str, float],
    ) -> dict[str, float]:
        if isinstance(calibration, FrozenBiasStrengthCalibrator):
            mapping = {
                name: estimate.bias_strength
                for name, estimate in calibration.estimates.items()
            }
        else:
            mapping = {str(name): float(value) for name, value in calibration.items()}
        if not mapping or any(not np.isfinite(value) or value < 0 or value > 2 for value in mapping.values()):
            raise ValueError("calibration must provide finite bias strengths in [0, 2]")
        return mapping

    def _strength_vector(self, sample_ids: np.ndarray) -> np.ndarray:
        missing = sorted(set(sample_ids.tolist()) - set(self.bias_strengths_))
        if missing:
            raise ValueError(f"missing frozen bias strengths for samples: {', '.join(missing)}")
        return np.asarray([self.bias_strengths_[sample] for sample in sample_ids], dtype=np.float64)

    def _fit_background(
        self,
        counts: np.ndarray,
        log_bias: np.ndarray,
        sample_ids: np.ndarray,
    ) -> FlankAccessibilityFit:
        return fit_flank_accessibility_background(
            counts,
            log_bias,
            self.positions,
            self._strength_vector(sample_ids),
            outer_start=self.outer_start,
            outer_end=self.outer_end,
            df=self.background_df,
            ridge=self.background_ridge,
        )

    def _initial_posteriors(
        self,
        counts: np.ndarray,
        baseline_expected: np.ndarray,
        tf_ids: np.ndarray,
    ) -> np.ndarray:
        scores = center_flank_likelihood_score(
            counts,
            baseline_expected,
            center_width=15,
            flank_width=30,
            gap=5,
        )
        posteriors = np.zeros(len(counts), dtype=np.float64)
        for tf in sorted(set(tf_ids.tolist())):
            selected = tf_ids == tf
            values = scores[selected]
            center = float(np.median(values))
            scale = float(np.median(np.abs(values - center)) * 1.4826)
            if not np.isfinite(scale) or scale < 1e-6:
                scale = float(np.std(values))
            if not np.isfinite(scale) or scale < 1e-6:
                scale = 1.0
            posteriors[selected] = np.clip(expit((values - center) / scale), 0.02, 0.98)
        return posteriors

    def _fit_unpooled_curve(
        self,
        counts: np.ndarray,
        baseline_logits: np.ndarray,
        posterior: np.ndarray,
        start: np.ndarray,
    ) -> np.ndarray:
        totals = counts.sum(axis=1)
        basis = self.footprint_basis_
        effective_cuts = max(float(np.sum(posterior * totals)), 1.0)
        second_difference = np.diff(basis, n=2, axis=0)

        def objective(coefficients: np.ndarray) -> tuple[float, np.ndarray]:
            curve = basis @ coefficients
            logits = baseline_logits + curve[None, :]
            probabilities = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
            ll = conditional_profile_log_likelihood(counts, logits)
            penalty_curve = second_difference @ coefficients
            penalty = 0.5 * self.footprint_ridge * float(np.sum(np.square(penalty_curve)))
            value = -float(np.sum(posterior * ll)) / effective_cuts + penalty
            error = posterior[:, None] * (totals[:, None] * probabilities - counts)
            gradient = np.sum(error, axis=0) @ basis / effective_cuts
            gradient += self.footprint_ridge * (second_difference.T @ penalty_curve)
            return value, np.asarray(gradient, dtype=np.float64)

        fitted = minimize(
            objective,
            np.asarray(start, dtype=np.float64),
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": 160, "ftol": 1e-10, "gtol": 1e-7},
        )
        if not fitted.success and not np.isfinite(fitted.fun):
            raise RuntimeError(f"footprint curve optimization failed: {fitted.message}")
        return np.asarray(fitted.x, dtype=np.float64)

    def _pool_curves(
        self,
        raw: Mapping[str, np.ndarray],
        posterior: np.ndarray,
        totals: np.ndarray,
        tf_ids: np.ndarray,
    ) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray]]:
        tf_weights = {
            tf: max(float(np.sum(posterior[tf_ids == tf] * totals[tf_ids == tf])), 1.0)
            for tf in raw
        }
        ordered = sorted(raw)
        weights = np.asarray([tf_weights[tf] for tf in ordered])
        global_coefficients = np.average(
            np.stack([raw[tf] for tf in ordered]),
            axis=0,
            weights=weights,
        )
        families = sorted(set(self.tf_family_[tf] for tf in ordered))
        family_coefficients: dict[str, np.ndarray] = {}
        for family in families:
            members = [tf for tf in ordered if self.tf_family_[tf] == family]
            member_weights = np.asarray([tf_weights[tf] for tf in members])
            raw_mean = np.average(
                np.stack([raw[tf] for tf in members]),
                axis=0,
                weights=member_weights,
            )
            effective = float(member_weights.sum())
            factor = effective / (effective + self.family_shrinkage)
            family_coefficients[family] = factor * (raw_mean - global_coefficients)
        tf_coefficients: dict[str, np.ndarray] = {}
        for tf in ordered:
            family = self.tf_family_[tf]
            residual = raw[tf] - global_coefficients - family_coefficients[family]
            factor = tf_weights[tf] / (tf_weights[tf] + self.tf_shrinkage)
            tf_coefficients[tf] = factor * residual
        return global_coefficients, family_coefficients, tf_coefficients

    def _effective_coefficients(self, tf: str) -> np.ndarray:
        if self.global_coefficients_ is None or tf not in self.tf_coefficients_:
            raise ValueError(f"factorization has no fitted curve for TF {tf!r}")
        family = self.tf_family_[tf]
        return (
            self.global_coefficients_
            + self.family_coefficients_[family]
            + self.tf_coefficients_[tf]
        )

    def _update_total_component(
        self,
        totals: np.ndarray,
        posterior: np.ndarray,
        tf_ids: np.ndarray,
    ) -> None:
        expected = np.zeros(len(totals), dtype=np.float64)
        for tf in sorted(set(tf_ids.tolist())):
            selected = tf_ids == tf
            q = posterior[selected]
            values = totals[selected]
            mean_bound = float(np.sum(q * values) / max(np.sum(q), 1e-8))
            mean_unbound = float(np.sum((1.0 - q) * values) / max(np.sum(1.0 - q), 1e-8))
            pooled = float(np.mean(values)) if len(values) else 1.0
            self.total_bound_means_[tf] = max(mean_bound, 0.05 * max(pooled, 1.0))
            self.total_unbound_means_[tf] = max(mean_unbound, 0.05 * max(pooled, 1.0))
            expected[selected] = (
                q * self.total_bound_means_[tf]
                + (1.0 - q) * self.total_unbound_means_[tf]
            )
        self.total_dispersion_ = estimate_nb_dispersion(totals, expected)

    def _posterior_update(
        self,
        counts: np.ndarray,
        baseline_logits: np.ndarray,
        tf_ids: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ll_unbound = conditional_profile_log_likelihood(counts, baseline_logits)
        ll_bound = np.empty(len(counts), dtype=np.float64)
        total_delta = np.zeros(len(counts), dtype=np.float64)
        totals = counts.sum(axis=1)
        prior_log_odds = np.zeros(len(counts), dtype=np.float64)
        for tf in sorted(set(tf_ids.tolist())):
            selected = tf_ids == tf
            curve = self.footprint_basis_ @ self._effective_coefficients(tf)
            ll_bound[selected] = conditional_profile_log_likelihood(
                counts[selected],
                baseline_logits[selected] + curve[None, :],
            )
            prior = float(np.clip(self.mixing_probabilities_.get(tf, 0.5), 0.01, 0.99))
            prior_log_odds[selected] = np.log(prior / (1.0 - prior))
            if self.use_total_component:
                mean0 = np.repeat(self.total_unbound_means_[tf], int(np.sum(selected)))
                mean1 = np.repeat(self.total_bound_means_[tf], int(np.sum(selected)))
                total_delta[selected] = (
                    negative_binomial_log_likelihood(
                        totals[selected],
                        mean1,
                        dispersion=self.total_dispersion_,
                    )
                    - negative_binomial_log_likelihood(
                        totals[selected],
                        mean0,
                        dispersion=self.total_dispersion_,
                    )
                )
        conditional_delta = ll_bound - ll_unbound
        posterior = np.clip(expit(prior_log_odds + conditional_delta + total_delta), 1e-5, 1.0 - 1e-5)
        return posterior, conditional_delta, total_delta

    def fit(
        self,
        counts: np.ndarray,
        log_bias: np.ndarray,
        sample_ids: Iterable[str],
        tf_ids: Iterable[str],
        family_ids: Iterable[str],
        calibration: FrozenBiasStrengthCalibrator | Mapping[str, float],
        *,
        max_iter: int = 30,
        tolerance: float = 1e-4,
    ) -> "FrozenParametricFactorization":
        observed, bias = _validate_profiles(counts, log_bias)
        if observed.shape[1] != len(self.positions):
            raise ValueError("profile width must match factorization positions")
        samples = self._string_vector(sample_ids, len(observed), "sample_ids")
        tfs = self._string_vector(tf_ids, len(observed), "tf_ids")
        families = self._string_vector(family_ids, len(observed), "family_ids")
        self.bias_strengths_ = self._strength_mapping(calibration)
        self._strength_vector(samples)
        self.tf_family_ = {}
        for tf in sorted(set(tfs.tolist())):
            tf_families = sorted(set(families[tfs == tf].tolist()))
            if len(tf_families) != 1:
                raise ValueError(f"TF {tf!r} maps to more than one motif family")
            self.tf_family_[tf] = tf_families[0]

        background_fit = self._fit_background(observed, bias, samples)
        strengths = self._strength_vector(samples)
        baseline_logits = strengths[:, None] * bias + background_fit.log_background
        posterior = self._initial_posteriors(observed, background_fit.expected, tfs)
        totals = observed.sum(axis=1)
        self.global_coefficients_ = np.zeros(self.footprint_df, dtype=np.float64)
        self.family_coefficients_ = {
            family: np.zeros(self.footprint_df, dtype=np.float64)
            for family in sorted(set(families.tolist()))
        }
        self.tf_coefficients_ = {
            tf: np.zeros(self.footprint_df, dtype=np.float64)
            for tf in sorted(set(tfs.tolist()))
        }
        self.training_history = []

        for iteration in range(1, int(max_iter) + 1):
            raw: dict[str, np.ndarray] = {}
            for tf in sorted(self.tf_coefficients_):
                selected = tfs == tf
                raw[tf] = self._fit_unpooled_curve(
                    observed[selected],
                    baseline_logits[selected],
                    posterior[selected],
                    self._effective_coefficients(tf),
                )
            (
                self.global_coefficients_,
                self.family_coefficients_,
                self.tf_coefficients_,
            ) = self._pool_curves(raw, posterior, totals, tfs)
            for tf in sorted(self.tf_coefficients_):
                selected = tfs == tf
                self.mixing_probabilities_[tf] = float(
                    np.clip((np.sum(posterior[selected]) + 1.0) / (np.sum(selected) + 2.0), 0.01, 0.99)
                )
            self._update_total_component(totals, posterior, tfs)
            updated, conditional_delta, total_delta = self._posterior_update(
                observed,
                baseline_logits,
                tfs,
            )
            change = float(np.max(np.abs(updated - posterior)))
            objective = float(
                np.mean(
                    np.logaddexp(
                        np.log(np.maximum(1.0 - posterior, 1e-12)),
                        np.log(np.maximum(posterior, 1e-12)) + conditional_delta + total_delta,
                    )
                )
            )
            self.training_history.append(
                {
                    "iteration": int(iteration),
                    "max_posterior_change": change,
                    "mean_mixture_evidence": objective,
                    "mean_posterior": float(np.mean(updated)),
                }
            )
            posterior = updated
            if change < float(tolerance):
                break

        self.metadata = {
            "training_sites": int(len(observed)),
            "training_cuts": float(observed.sum()),
            "samples": sorted(set(samples.tolist())),
            "tfs": sorted(self.tf_family_),
            "families": sorted(set(self.tf_family_.values())),
            "iterations_completed": int(len(self.training_history)),
        }
        return self

    def predict(
        self,
        counts: np.ndarray,
        log_bias: np.ndarray,
        sample_ids: Iterable[str],
        tf_ids: Iterable[str],
    ) -> FactorizationResult:
        if self.global_coefficients_ is None:
            raise ValueError("factorization has not been fitted")
        observed, bias = _validate_profiles(counts, log_bias)
        if observed.shape[1] != len(self.positions):
            raise ValueError("profile width must match factorization positions")
        samples = self._string_vector(sample_ids, len(observed), "sample_ids")
        tfs = self._string_vector(tf_ids, len(observed), "tf_ids")
        missing = sorted(set(tfs.tolist()) - set(self.tf_coefficients_))
        if missing:
            raise ValueError(f"factorization has no frozen TF curve for: {', '.join(missing)}")
        background_fit = self._fit_background(observed, bias, samples)
        strengths = self._strength_vector(samples)
        baseline_logits = strengths[:, None] * bias + background_fit.log_background
        footprint_effect = np.vstack(
            [self.footprint_basis_ @ self._effective_coefficients(tf) for tf in tfs]
        )
        bound_logits = baseline_logits + footprint_effect
        expected_unbound = expected_profile_counts(observed, baseline_logits)
        expected_bound = expected_profile_counts(observed, bound_logits)
        posterior, conditional_delta, total_delta = self._posterior_update(
            observed,
            baseline_logits,
            tfs,
        )
        return FactorizationResult(
            posterior_bound=posterior,
            expected_unbound=expected_unbound,
            expected_bound=expected_bound,
            log_background=background_fit.log_background,
            footprint_log_effect=footprint_effect,
            conditional_log_bayes_factor=conditional_delta,
            total_log_bayes_factor=total_delta,
        )

    def save(
        self,
        path: str | Path,
        metadata: Mapping[str, Any] | None = None,
    ) -> tuple[Path, Path]:
        if self.global_coefficients_ is None:
            raise ValueError("factorization has not been fitted")
        npz_path, json_path = _npz_json_paths(path)
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        samples = np.asarray(sorted(self.bias_strengths_), dtype="U")
        families = np.asarray(sorted(self.family_coefficients_), dtype="U")
        tfs = np.asarray(sorted(self.tf_coefficients_), dtype="U")
        np.savez_compressed(
            npz_path,
            positions=self.positions,
            footprint_knots=self.footprint_knots_,
            sample=samples,
            bias_strength=np.asarray([self.bias_strengths_[name] for name in samples]),
            family=families,
            family_coefficients=np.stack([self.family_coefficients_[name] for name in families]),
            tf=tfs,
            tf_family=np.asarray([self.tf_family_[name] for name in tfs], dtype="U"),
            tf_coefficients=np.stack([self.tf_coefficients_[name] for name in tfs]),
            global_coefficients=self.global_coefficients_,
            mixing_probability=np.asarray([self.mixing_probabilities_[name] for name in tfs]),
            total_unbound_mean=np.asarray([self.total_unbound_means_[name] for name in tfs]),
            total_bound_mean=np.asarray([self.total_bound_means_[name] for name in tfs]),
            total_dispersion=np.asarray([self.total_dispersion_]),
        )
        document = {
            "schema": FACTORIZATION_SCHEMA,
            "npz_sha256": _sha256_file(npz_path),
            "configuration": {
                "background_df": self.background_df,
                "footprint_df": self.footprint_df,
                "footprint_limit": self.footprint_limit,
                "outer_start": self.outer_start,
                "outer_end": self.outer_end,
                "background_ridge": self.background_ridge,
                "footprint_ridge": self.footprint_ridge,
                "family_shrinkage": self.family_shrinkage,
                "tf_shrinkage": self.tf_shrinkage,
                "use_total_component": self.use_total_component,
                "seed": self.seed,
            },
            "training_history": self.training_history,
            "metadata": {**self.metadata, **dict(metadata or {})},
        }
        json_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "FrozenParametricFactorization":
        npz_path, json_path = _npz_json_paths(path)
        document = json.loads(json_path.read_text(encoding="utf-8"))
        if document.get("schema") != FACTORIZATION_SCHEMA:
            raise ValueError(f"unsupported factorization schema: {document.get('schema')}")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("factorization checksum does not match its JSON metadata")
        with np.load(npz_path, allow_pickle=False) as arrays:
            required = {
                "positions", "sample", "bias_strength", "family", "family_coefficients",
                "tf", "tf_family", "tf_coefficients", "global_coefficients",
                "mixing_probability", "total_unbound_mean", "total_bound_mean",
                "total_dispersion",
            }
            if not required.issubset(arrays.files):
                raise ValueError("factorization NPZ is missing required arrays")
            configuration = dict(document["configuration"])
            model = cls(np.asarray(arrays["positions"], dtype=np.float64), **configuration)
            samples = np.asarray(arrays["sample"], dtype=str)
            strengths = np.asarray(arrays["bias_strength"], dtype=np.float64)
            families = np.asarray(arrays["family"], dtype=str)
            family_coefficients = np.asarray(arrays["family_coefficients"], dtype=np.float64)
            tfs = np.asarray(arrays["tf"], dtype=str)
            tf_families = np.asarray(arrays["tf_family"], dtype=str)
            tf_coefficients = np.asarray(arrays["tf_coefficients"], dtype=np.float64)
            expected_family_shape = (len(families), model.footprint_df)
            expected_tf_shape = (len(tfs), model.footprint_df)
            if strengths.shape != samples.shape:
                raise ValueError("factorization sample arrays have incompatible shapes")
            if family_coefficients.shape != expected_family_shape:
                raise ValueError("factorization family coefficients have incompatible shape")
            if tf_families.shape != tfs.shape or tf_coefficients.shape != expected_tf_shape:
                raise ValueError("factorization TF arrays have incompatible shapes")
            model.bias_strengths_ = dict(zip(samples.tolist(), strengths.tolist()))
            model.family_coefficients_ = {
                name: family_coefficients[index]
                for index, name in enumerate(families.tolist())
            }
            model.tf_family_ = dict(zip(tfs.tolist(), tf_families.tolist()))
            model.tf_coefficients_ = {
                name: tf_coefficients[index]
                for index, name in enumerate(tfs.tolist())
            }
            model.global_coefficients_ = np.asarray(arrays["global_coefficients"], dtype=np.float64)
            for key, target in (
                ("mixing_probability", model.mixing_probabilities_),
                ("total_unbound_mean", model.total_unbound_means_),
                ("total_bound_mean", model.total_bound_means_),
            ):
                values = np.asarray(arrays[key], dtype=np.float64)
                if values.shape != tfs.shape:
                    raise ValueError(f"factorization {key} has incompatible shape")
                target.update(zip(tfs.tolist(), values.tolist()))
            dispersion = np.asarray(arrays["total_dispersion"], dtype=np.float64)
            if dispersion.shape != (1,):
                raise ValueError("factorization total dispersion has incompatible shape")
            model.total_dispersion_ = float(dispersion[0])
        model.training_history = list(document.get("training_history", []))
        model.metadata = dict(document.get("metadata", {}))
        return model
