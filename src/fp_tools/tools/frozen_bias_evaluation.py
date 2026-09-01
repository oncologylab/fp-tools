"""Statistics for the frozen parametric-bias research experiment.

The functions in this module deliberately have no command-line integration.
They provide leakage-safe control-data comparisons used by the research
benchmark while leaving the established PWM/DWM workflows unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable, Protocol

import numpy as np
import pandas as pd
from scipy.special import rel_entr
from scipy.stats import norm

from fp_tools.tools.parametric_bias import BiasFeatureSpec


TOBIAS_DWM_SCHEMA = "fp-tools-tobias-dwm-control-reference-v1"


def _file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class ConditionalBiasModel(Protocol):
    """Minimal interface shared by parametric and reference bias scorers."""

    def probabilities(self, contexts: np.ndarray) -> np.ndarray: ...


class TobiasDwmReferenceModel:
    """Safe, vectorized reference implementation of the TOBIAS-style DWM.

    The conventional model estimates an all-pairs dinucleotide distribution
    around observed cuts and contrasts it with the corresponding background
    sequence distribution.  Cut multiplicity is capped at ten, matching the
    established ATACorrect estimator.  This class is benchmark-only and does
    not replace or modify the production Cython implementation.
    """

    def __init__(self, context_length: int = 11) -> None:
        if context_length < 3 or context_length % 2 != 1:
            raise ValueError("DWM context length must be an odd integer >=3")
        self.feature_spec = BiasFeatureSpec(
            "tobias_dwm", int(context_length), tuple(range(1, context_length))
        )
        length = int(context_length)
        self.bias_pwm = np.full((4, length), 0.25, dtype=np.float64)
        self.background_pwm = np.full((4, length), 0.25, dtype=np.float64)
        self.bias_dwm = np.full((4, 4, length, length), 1.0 / 16.0)
        self.background_dwm = np.full((4, 4, length, length), 1.0 / 16.0)
        self.metadata: dict[str, object] = {}

    @staticmethod
    def _weighted_counts(
        contexts: np.ndarray,
        weights: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, float]:
        values = np.asarray(contexts, dtype=np.uint8)
        weight = np.asarray(weights, dtype=np.float64)
        if values.ndim != 2 or weight.shape != (len(values),):
            raise ValueError("flattened contexts and weights have incompatible shapes")
        valid = np.all(values < 4, axis=1) & np.isfinite(weight) & (weight >= 0)
        values = values[valid]
        weight = weight[valid]
        if not len(values) or float(np.sum(weight)) <= 0:
            raise ValueError("DWM fitting requires valid positively weighted contexts")
        length = values.shape[1]
        pwm_counts = np.zeros((4, length), dtype=np.float64)
        pair_counts = np.zeros((4, 4, length, length), dtype=np.float64)
        for position in range(length):
            pwm_counts[:, position] = np.bincount(
                values[:, position], weights=weight, minlength=4
            )[:4]
        for first in range(length):
            left = values[:, first].astype(np.int64)
            for second in range(length):
                code = left * 4 + values[:, second].astype(np.int64)
                pair_counts[:, :, first, second] = np.bincount(
                    code, weights=weight, minlength=16
                ).reshape(4, 4)
        return pwm_counts, pair_counts, float(np.sum(weight))

    def fit(
        self,
        contexts: np.ndarray,
        counts: np.ndarray,
        *,
        maximum_cut_multiplicity: float = 10.0,
    ) -> "TobiasDwmReferenceModel":
        values = np.asarray(contexts, dtype=np.uint8)
        observed = np.asarray(counts, dtype=np.float64)
        if values.ndim != 3 or observed.shape != values.shape[:2]:
            raise ValueError("contexts must be 3D and counts must match")
        if np.any(observed < 0) or maximum_cut_multiplicity <= 0:
            raise ValueError("counts must be non-negative and the cut cap positive")
        flat = values.reshape(-1, values.shape[-1])
        cut_weight = np.minimum(observed.reshape(-1), maximum_cut_multiplicity)
        background_weight = np.ones(len(flat), dtype=np.float64)
        bias_pwm_counts, bias_pair_counts, no_bias = self._weighted_counts(
            flat, cut_weight
        )
        background_pwm_counts, background_pair_counts, no_background = (
            self._weighted_counts(flat, background_weight)
        )
        self.bias_pwm = (bias_pwm_counts + 1.0) / (no_bias + 4.0)
        self.background_pwm = (background_pwm_counts + 1.0) / (no_background + 4.0)
        bias_pseudo = max(16.0, no_bias)
        background_pseudo = max(16.0, no_background)
        for first in range(self.feature_spec.context_length):
            for second in range(self.feature_spec.context_length):
                bias_prior = np.outer(self.bias_pwm[:, first], self.bias_pwm[:, second])
                background_prior = np.outer(
                    self.background_pwm[:, first], self.background_pwm[:, second]
                )
                self.bias_dwm[:, :, first, second] = (
                    bias_pair_counts[:, :, first, second] + bias_pseudo * bias_prior
                ) / (no_bias + bias_pseudo)
                self.background_dwm[:, :, first, second] = (
                    background_pair_counts[:, :, first, second]
                    + background_pseudo * background_prior
                ) / (no_background + background_pseudo)
        self.metadata.update(
            {
                "training_contexts": int(len(flat)),
                "training_capped_cuts": float(np.sum(cut_weight)),
                "maximum_cut_multiplicity": float(maximum_cut_multiplicity),
            }
        )
        return self

    @staticmethod
    def _distribution_score(
        contexts: np.ndarray,
        pwm: np.ndarray,
        dwm: np.ndarray,
    ) -> np.ndarray:
        values = np.asarray(contexts, dtype=np.uint8)
        length = values.shape[1]
        log_pwm = np.log2(np.maximum(pwm, np.finfo(float).tiny))
        log_dwm = np.log2(np.maximum(dwm, np.finfo(float).tiny))
        output = np.zeros(len(values), dtype=np.float64)
        row = np.arange(len(values))
        for second in range(length):
            actual = values[:, second].astype(np.int64)
            candidate_logs = np.zeros((len(values), 4), dtype=np.float64)
            for candidate_base in range(4):
                score = np.full(len(values), log_pwm[candidate_base, second])
                for first in range(length):
                    if first == second:
                        continue
                    left = values[:, first].astype(np.int64)
                    score += (
                        log_dwm[left, candidate_base, first, second]
                        - log_pwm[candidate_base, second]
                    )
                candidate_logs[:, candidate_base] = score
            maximum = np.max(candidate_logs, axis=1)
            normalizer = maximum + np.log2(
                np.sum(np.exp2(candidate_logs - maximum[:, None]), axis=1)
            )
            output += candidate_logs[row, actual] - normalizer
        return output

    def log_scores(self, contexts: np.ndarray, batch_size: int = 65536) -> np.ndarray:
        values = np.asarray(contexts, dtype=np.uint8)
        if values.ndim != 3 or values.shape[2] != self.feature_spec.context_length:
            raise ValueError("contexts have an incompatible DWM shape")
        flat = values.reshape(-1, values.shape[-1])
        valid = np.all(flat < 4, axis=1)
        output = np.full(len(flat), -np.inf, dtype=np.float64)
        valid_indexes = np.flatnonzero(valid)
        for start in range(0, len(valid_indexes), batch_size):
            selected = valid_indexes[start : start + batch_size]
            batch = flat[selected]
            bias = self._distribution_score(batch, self.bias_pwm, self.bias_dwm)
            background = self._distribution_score(
                batch, self.background_pwm, self.background_dwm
            )
            output[selected] = bias - background
        return output.reshape(values.shape[:2])

    def probabilities(self, contexts: np.ndarray) -> np.ndarray:
        logits = self.log_scores(contexts)
        finite = np.isfinite(logits)
        maximum = np.max(np.where(finite, logits, -np.inf), axis=1, keepdims=True)
        weights = np.exp(np.where(finite, logits - maximum, -np.inf))
        totals = weights.sum(axis=1, keepdims=True)
        if np.any(totals <= 0):
            raise ValueError("DWM contexts contain no valid positions")
        return weights / totals

    def save(
        self, path: str | Path, metadata: dict[str, object] | None = None
    ) -> tuple[Path, Path]:
        npz_path = Path(path)
        if npz_path.suffix != ".npz":
            npz_path = Path(str(npz_path) + ".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            bias_pwm=self.bias_pwm,
            background_pwm=self.background_pwm,
            bias_dwm=self.bias_dwm,
            background_dwm=self.background_dwm,
        )
        document = {
            "schema": TOBIAS_DWM_SCHEMA,
            "context_length": self.feature_spec.context_length,
            "npz_sha256": _file_sha256(npz_path),
            "metadata": {**self.metadata, **dict(metadata or {})},
        }
        json_path.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "TobiasDwmReferenceModel":
        npz_path = Path(path)
        if npz_path.suffix != ".npz":
            npz_path = Path(str(npz_path) + ".npz")
        document = json.loads(npz_path.with_suffix(".json").read_text(encoding="utf-8"))
        if document.get("schema") != TOBIAS_DWM_SCHEMA:
            raise ValueError("unsupported TOBIAS DWM reference schema")
        if document.get("npz_sha256") != _file_sha256(npz_path):
            raise ValueError("TOBIAS DWM reference checksum mismatch")
        model = cls(int(document["context_length"]))
        with np.load(npz_path, allow_pickle=False) as arrays:
            for name in (
                "bias_pwm",
                "background_pwm",
                "bias_dwm",
                "background_dwm",
            ):
                if name not in arrays:
                    raise ValueError(f"TOBIAS DWM reference is missing {name}")
                setattr(model, name, np.asarray(arrays[name], dtype=np.float64))
        length = model.feature_spec.context_length
        if model.bias_pwm.shape != (4, length) or model.background_pwm.shape != (
            4,
            length,
        ):
            raise ValueError("TOBIAS DWM PWM arrays have invalid shapes")
        expected_pair_shape = (4, 4, length, length)
        if (
            model.bias_dwm.shape != expected_pair_shape
            or model.background_dwm.shape != expected_pair_shape
        ):
            raise ValueError("TOBIAS DWM pair arrays have invalid shapes")
        model.metadata = dict(document.get("metadata", {}))
        return model


@dataclass(frozen=True)
class ConditionalControlScores:
    """Per-window quantities needed for paired control comparisons."""

    log_likelihood: np.ndarray
    null_log_likelihood: np.ndarray
    deviance: np.ndarray
    totals: np.ndarray
    calibration_error: float
    aggregate_jsd: float

    @property
    def conditional_nll(self) -> float:
        return float(-np.sum(self.log_likelihood) / np.sum(self.totals))

    @property
    def nll_gain(self) -> float:
        return float(
            (np.sum(self.log_likelihood) - np.sum(self.null_log_likelihood))
            / np.sum(self.totals)
        )

    @property
    def deviance_per_cut(self) -> float:
        return float(np.sum(self.deviance) / np.sum(self.totals))


def conditional_control_scores(
    model: ConditionalBiasModel,
    contexts: np.ndarray,
    counts: np.ndarray,
) -> ConditionalControlScores:
    """Score each conditional-total window without pooling observations.

    Keeping one likelihood contribution per window allows models to be compared
    with a paired chromosome/library bootstrap.  Invalid sequence positions are
    removed before probabilities are normalized.
    """

    sequence = np.asarray(contexts)
    observed = np.asarray(counts, dtype=np.float64)
    if sequence.ndim != 3 or observed.shape != sequence.shape[:2]:
        raise ValueError(
            "contexts must be 3D and counts must match their first two axes"
        )
    if np.any(observed < 0) or np.any(~np.isfinite(observed)):
        raise ValueError("counts must be finite and non-negative")
    valid = np.all(sequence < 4, axis=2)
    if np.any((observed > 0) & ~valid):
        raise ValueError("positive cuts cannot occur at invalid sequence positions")
    observed = np.where(valid, observed, 0.0)
    totals = observed.sum(axis=1)
    keep = (totals > 0) & valid.any(axis=1)
    if not np.any(keep):
        raise ValueError("at least one valid window containing cuts is required")
    sequence = sequence[keep]
    observed = observed[keep]
    valid = valid[keep]
    totals = totals[keep]

    probabilities = np.asarray(model.probabilities(sequence), dtype=np.float64)
    if probabilities.shape != observed.shape:
        raise ValueError("model probabilities have an incompatible shape")
    probabilities = np.where(valid, probabilities, 0.0)
    row_sums = probabilities.sum(axis=1)
    if np.any(~np.isfinite(row_sums)) or np.any(row_sums <= 0):
        raise ValueError("model returned invalid conditional probabilities")
    probabilities /= row_sums[:, None]
    log_probability = np.log(np.maximum(probabilities, np.finfo(float).tiny))
    log_likelihood = np.sum(observed * log_probability, axis=1)
    valid_positions = valid.sum(axis=1)
    null_log_likelihood = -totals * np.log(valid_positions)

    with np.errstate(divide="ignore", invalid="ignore"):
        saturated = np.where(
            observed > 0,
            observed * np.log(observed / totals[:, None]),
            0.0,
        )
    deviance = 2.0 * np.sum(saturated - observed * log_probability, axis=1)

    expected = probabilities * totals[:, None]
    aggregate_observed = observed.sum(axis=0)
    aggregate_expected = expected.sum(axis=0)
    aggregate_observed /= max(float(aggregate_observed.sum()), 1.0)
    aggregate_expected /= max(float(aggregate_expected.sum()), 1.0)
    midpoint = 0.5 * (aggregate_observed + aggregate_expected)
    jsd = 0.5 * (
        float(np.sum(rel_entr(aggregate_observed, midpoint)))
        + float(np.sum(rel_entr(aggregate_expected, midpoint)))
    )

    positive = expected > 0
    ratios = np.divide(observed, expected, out=np.zeros_like(observed), where=positive)
    probability_values = probabilities[positive]
    ratio_values = ratios[positive]
    if len(probability_values) >= 10:
        quantiles = np.unique(
            np.quantile(probability_values, np.linspace(0.0, 1.0, 11))
        )
        groups = np.digitize(probability_values, quantiles[1:-1])
        calibration = float(
            np.mean(
                [
                    abs(float(np.mean(ratio_values[groups == group])) - 1.0)
                    for group in np.unique(groups)
                ]
            )
        )
    else:
        calibration = np.nan
    return ConditionalControlScores(
        log_likelihood=log_likelihood,
        null_log_likelihood=null_log_likelihood,
        deviance=deviance,
        totals=totals,
        calibration_error=calibration,
        aggregate_jsd=jsd,
    )


def paired_block_bootstrap_gain(
    candidate_log_likelihood: np.ndarray,
    reference_log_likelihood: np.ndarray,
    totals: np.ndarray,
    blocks: Iterable[object],
    *,
    bootstraps: int = 2000,
    seed: int = 2026,
) -> dict[str, float | int | bool]:
    """Bootstrap the paired per-cut likelihood gain over independent blocks."""

    candidate = np.asarray(candidate_log_likelihood, dtype=np.float64)
    reference = np.asarray(reference_log_likelihood, dtype=np.float64)
    weights = np.asarray(totals, dtype=np.float64)
    block_values = np.asarray(list(blocks), dtype=str)
    if candidate.shape != reference.shape or candidate.shape != weights.shape:
        raise ValueError("paired likelihood and total arrays must have equal shapes")
    if block_values.shape != candidate.shape:
        raise ValueError("one block identifier is required per window")
    if bootstraps < 20 or np.any(weights <= 0):
        raise ValueError("bootstraps must be >=20 and all totals must be positive")
    identities = np.unique(block_values)
    if len(identities) < 2:
        raise ValueError("paired block bootstrap requires at least two blocks")
    numerators = np.asarray(
        [np.sum((candidate - reference)[block_values == block]) for block in identities]
    )
    denominators = np.asarray(
        [np.sum(weights[block_values == block]) for block in identities]
    )
    point = float(np.sum(numerators) / np.sum(denominators))
    rng = np.random.default_rng(seed)
    values = np.empty(bootstraps, dtype=np.float64)
    for index in range(bootstraps):
        selected = rng.integers(0, len(identities), size=len(identities))
        values[index] = np.sum(numerators[selected]) / np.sum(denominators[selected])
    lower, upper = np.quantile(values, [0.025, 0.975])
    return {
        "blocks": int(len(identities)),
        "windows": int(len(candidate)),
        "paired_log_likelihood_gain_per_cut": point,
        "paired_gain_lower_95": float(lower),
        "paired_gain_upper_95": float(upper),
        "paired_gain_interval_excludes_zero": bool(lower > 0 or upper < 0),
    }


def retain_control_candidates(
    metrics: pd.DataFrame,
    *,
    maximum_retained: int = 2,
    maximum_model_size_mb: float = 25.0,
) -> pd.DataFrame:
    """Apply the frozen one-SE and significant-larger-model retention rule."""

    required = {
        "candidate_id",
        "context_length",
        "mean_conditional_nll",
        "standard_error_conditional_nll",
        "minimum_library_nll_gain",
        "model_size_mb",
    }
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError("candidate metrics are missing: " + ", ".join(sorted(missing)))
    output = metrics.copy()
    output["within_one_standard_error"] = False
    output["retained"] = False
    output["retention_reason"] = ""
    eligible = output[
        (output["minimum_library_nll_gain"].astype(float) > 0)
        & (output["model_size_mb"].astype(float) <= maximum_model_size_mb)
    ].copy()
    if eligible.empty:
        return output
    best = eligible.sort_values(
        ["mean_conditional_nll", "context_length", "candidate_id"],
        kind="mergesort",
    ).iloc[0]
    threshold = float(best["mean_conditional_nll"]) + float(
        best["standard_error_conditional_nll"]
    )
    within = eligible["mean_conditional_nll"].astype(float) <= threshold
    output.loc[eligible.index[within], "within_one_standard_error"] = True
    smallest_index = (
        eligible.loc[within]
        .sort_values(
            ["context_length", "model_size_mb", "mean_conditional_nll", "candidate_id"],
            kind="mergesort",
        )
        .index[0]
    )
    output.loc[smallest_index, ["retained", "retention_reason"]] = [
        True,
        "smallest model within one SE of best control likelihood",
    ]
    if maximum_retained <= 1:
        return output

    larger = eligible[
        (
            eligible["context_length"].astype(int)
            > int(output.at[smallest_index, "context_length"])
        )
    ].copy()
    if "gain_over_smallest_lower_95" in larger.columns:
        larger = larger[larger["gain_over_smallest_lower_95"].astype(float) > 0]
    else:
        larger = larger.iloc[0:0]
    if not larger.empty:
        larger_index = larger.sort_values(
            ["mean_conditional_nll", "context_length", "candidate_id"],
            kind="mergesort",
        ).index[0]
        output.loc[larger_index, ["retained", "retention_reason"]] = [
            True,
            "larger model with significant paired likelihood gain",
        ]
    return output


def motif_residual_effect(
    observed: np.ndarray,
    expected: np.ndarray,
    positions: np.ndarray,
    *,
    bootstraps: int = 1000,
    seed: int = 2026,
    threshold: float = 0.25,
    pseudocount: float = 0.5,
) -> dict[str, float | int | bool]:
    """Test motif-centered observed-minus-predicted log-ratio structure."""

    y = np.asarray(observed, dtype=np.float64)
    mu = np.asarray(expected, dtype=np.float64)
    x = np.asarray(positions, dtype=np.float64)
    if y.shape != mu.shape or y.ndim != 2 or y.shape[1] != len(x):
        raise ValueError("observed and expected motif profiles must match positions")
    if len(y) < 2 or np.any(y < 0) or np.any(mu < 0):
        raise ValueError("at least two non-negative profiles are required")
    center = np.abs(x) <= 15
    flanks = (np.abs(x) >= 40) & (np.abs(x) <= 80)
    if not center.any() or not flanks.any():
        raise ValueError("positions must cover the center and 40-80 bp flanks")
    residual = np.log((y + pseudocount) / (mu + pseudocount))
    effects = residual[:, flanks].mean(axis=1) - residual[:, center].mean(axis=1)
    point = float(np.mean(effects))
    rng = np.random.default_rng(seed)
    values = np.empty(bootstraps, dtype=np.float64)
    for index in range(bootstraps):
        selected = rng.integers(0, len(effects), size=len(effects))
        values[index] = np.mean(effects[selected])
    lower, upper = np.quantile(values, [0.025, 0.975])
    flag = bool(abs(point) >= threshold and (lower > 0 or upper < 0))
    return {
        "sites": int(len(y)),
        "observed_minus_predicted_center_flank_effect": point,
        "effect_lower_95": float(lower),
        "effect_upper_95": float(upper),
        "motif_residual_flag": flag,
    }


def wilson_interval(
    positive: int,
    total: int,
    *,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Return the Wilson binomial confidence interval without continuity correction."""

    if total < 1 or positive < 0 or positive > total or not 0 < confidence < 1:
        raise ValueError("invalid binomial counts or confidence")
    z = float(norm.ppf(0.5 + confidence / 2.0))
    p = positive / total
    denominator = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / denominator
    half = (
        z * np.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / denominator
    )
    lower = 0.0 if positive == 0 else max(0.0, center - half)
    upper = 1.0 if positive == total else min(1.0, center + half)
    return float(lower), float(upper)


def negative_control_safety(
    candidate_calls: np.ndarray,
    reference_calls: np.ndarray,
    *,
    maximum_rate: float = 0.05,
    maximum_upper: float = 0.05,
    maximum_increase: float = 0.01,
) -> dict[str, float | int | bool]:
    """Apply the frozen independent-null false-positive gates."""

    candidate = np.asarray(candidate_calls, dtype=bool)
    reference = np.asarray(reference_calls, dtype=bool)
    if candidate.shape != reference.shape or candidate.ndim != 1 or len(candidate) == 0:
        raise ValueError("candidate and reference calls must be equal nonempty vectors")
    candidate_positive = int(np.sum(candidate))
    reference_positive = int(np.sum(reference))
    candidate_rate = candidate_positive / len(candidate)
    reference_rate = reference_positive / len(reference)
    _lower, upper = wilson_interval(candidate_positive, len(candidate))
    return {
        "total": int(len(candidate)),
        "candidate_positive": candidate_positive,
        "reference_positive": reference_positive,
        "candidate_false_positive_rate": float(candidate_rate),
        "reference_false_positive_rate": float(reference_rate),
        "candidate_wilson_upper_95": upper,
        "false_positive_rate_increase": float(candidate_rate - reference_rate),
        "passed_negative_control_safety": bool(
            candidate_rate <= maximum_rate
            and upper <= maximum_upper
            and candidate_rate - reference_rate <= maximum_increase
        ),
    }
