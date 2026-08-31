#!/usr/bin/env python
"""Fast parametric sequence-bias models and calibrated cut-site residuals.

The production :mod:`fp_tools.tools.atacorrect` path continues to use the
established PWM/DWM implementation.  This module provides the opt-in research
primitives used to compare a conditional log-linear sequence model with those
baselines without introducing a deep-learning dependency.

The model conditions on the total cuts in each window.  It therefore learns
where Tn5 cuts within a window rather than using total accessibility as a
training target::

    cuts_w | total_w ~ Multinomial(total_w, softmax(sequence_features @ beta))

Reverse-strand contexts must be reverse-complemented before fitting/scoring so
one coefficient set describes both orientations.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.special import gammaln, logsumexp


MODEL_SCHEMA = "fp-tools-conditional-bias-v1"
BASES = "ACGT"
BASE_TO_CODE = np.full(256, 4, dtype=np.uint8)
for _index, _base in enumerate(BASES):
    BASE_TO_CODE[ord(_base)] = _index
    BASE_TO_CODE[ord(_base.lower())] = _index
COMPLEMENT = np.asarray([3, 2, 1, 0, 4], dtype=np.uint8)


@dataclass(frozen=True)
class BiasFeatureSpec:
    """Sequence context and interaction geometry for a bias model."""

    name: str
    context_length: int
    pair_distances: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.context_length < 2:
            raise ValueError("context_length must be at least 2")
        distances = tuple(sorted(set(int(value) for value in self.pair_distances)))
        if any(value < 1 or value >= self.context_length for value in distances):
            raise ValueError("pair distances must be between 1 and context_length - 1")
        object.__setattr__(self, "pair_distances", distances)

    @classmethod
    def selma10(cls) -> "BiasFeatureSpec":
        """Compact 10-mer model with adjacent and Tn5-dimer interactions."""

        return cls("selma10", 10, (1, 9))

    @classmethod
    def loglinear81(cls) -> "BiasFeatureSpec":
        """ChromBPNet-inspired 81-bp parametric receptive field."""

        return cls("loglinear81", 81, (1, 5, 9))


def encode_sequence(sequence: str) -> np.ndarray:
    """Encode A/C/G/T as 0/1/2/3 and every other symbol as 4."""

    raw = np.frombuffer(sequence.encode("ascii", errors="replace"), dtype=np.uint8)
    return BASE_TO_CODE[raw]


def reverse_complement_contexts(contexts: np.ndarray) -> np.ndarray:
    """Reverse-complement one context or an arbitrary stack of contexts."""

    values = np.asarray(contexts, dtype=np.uint8)
    if values.ndim < 1:
        raise ValueError("contexts must have at least one dimension")
    if np.any(values > 4):
        raise ValueError("encoded contexts must contain only values 0..4")
    return COMPLEMENT[values[..., ::-1]]


def contexts_from_sequence(
    sequence: str | np.ndarray,
    positions: Iterable[int],
    context_length: int,
    strands: Iterable[str | int] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract oriented sequence contexts centered on genomic cut positions.

    Returns ``(contexts, valid)``.  Invalid boundary contexts are filled with
    the unknown-base code and marked false.  A reverse strand may be encoded as
    ``"-"``, ``-1``, ``"reverse"``, or ``True``.
    """

    codes = encode_sequence(sequence) if isinstance(sequence, str) else np.asarray(sequence, dtype=np.uint8)
    positions_array = np.asarray(list(positions), dtype=np.int64)
    if strands is None:
        strands_array = np.zeros(len(positions_array), dtype=bool)
    else:
        raw_strands = list(strands)
        if len(raw_strands) != len(positions_array):
            raise ValueError("strands must contain one value per position")
        strands_array = np.asarray(
            [value in ("-", -1, "reverse", True) for value in raw_strands],
            dtype=bool,
        )

    left = context_length // 2
    right = context_length - left
    contexts = np.full((len(positions_array), context_length), 4, dtype=np.uint8)
    valid = (positions_array >= left) & (positions_array + right <= len(codes))
    for row in np.flatnonzero(valid):
        position = int(positions_array[row])
        contexts[row] = codes[position - left:position + right]
    reverse_rows = valid & strands_array
    contexts[reverse_rows] = reverse_complement_contexts(contexts[reverse_rows])
    valid &= np.all(contexts < 4, axis=1)
    return contexts, valid


def _double_center(array: np.ndarray) -> np.ndarray:
    """Apply sum-to-zero identifiability constraints to a 4x4 interaction."""

    return (
        array
        - array.mean(axis=-1, keepdims=True)
        - array.mean(axis=-2, keepdims=True)
        + array.mean(axis=(-2, -1), keepdims=True)
    )


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class ConditionalSequenceBiasModel:
    """Conditional log-linear sequence model fitted with deterministic Adam."""

    def __init__(self, feature_spec: BiasFeatureSpec):
        self.feature_spec = feature_spec
        length = feature_spec.context_length
        self.main = np.zeros((length, 4), dtype=np.float64)
        self.pairs = {
            distance: np.zeros((length - distance, 4, 4), dtype=np.float64)
            for distance in feature_spec.pair_distances
        }
        self.training_history: list[dict[str, float | int]] = []
        self.metadata: dict[str, Any] = {}

    def _validate_contexts(self, contexts: np.ndarray) -> np.ndarray:
        values = np.asarray(contexts, dtype=np.uint8)
        if values.ndim != 3:
            raise ValueError("contexts must have shape (windows, positions, context_length)")
        if values.shape[2] != self.feature_spec.context_length:
            raise ValueError(
                f"expected context length {self.feature_spec.context_length}, got {values.shape[2]}"
            )
        if np.any(values > 4):
            raise ValueError("encoded contexts must contain only values 0..4")
        return values

    def log_scores(self, contexts: np.ndarray) -> np.ndarray:
        """Return unnormalized log propensity for each candidate cut position."""

        values = self._validate_contexts(contexts)
        flat = values.reshape(-1, values.shape[-1])
        valid = np.all(flat < 4, axis=1)
        scores = np.zeros(len(flat), dtype=np.float64)
        valid_flat = flat[valid]
        if len(valid_flat):
            row_index = np.arange(len(valid_flat))
            for position in range(self.feature_spec.context_length):
                scores[valid] += self.main[position, valid_flat[:, position]]
            for distance, coefficients in self.pairs.items():
                for position in range(self.feature_spec.context_length - distance):
                    scores[valid] += coefficients[
                        position,
                        valid_flat[:, position],
                        valid_flat[:, position + distance],
                    ]
        scores[~valid] = -np.inf
        return scores.reshape(values.shape[:2])

    def probabilities(self, contexts: np.ndarray) -> np.ndarray:
        """Return conditional cut probabilities within each input window."""

        scores = self.log_scores(contexts)
        normalizer = logsumexp(scores, axis=1, keepdims=True)
        probabilities = np.exp(scores - normalizer)
        probabilities[~np.isfinite(probabilities)] = 0.0
        return probabilities

    def conditional_nll(self, contexts: np.ndarray, counts: np.ndarray) -> float:
        """Mean conditional negative log likelihood per observed cut."""

        values = self._validate_contexts(contexts)
        observed = np.asarray(counts, dtype=np.float64)
        if observed.shape != values.shape[:2]:
            raise ValueError("counts must match the first two context dimensions")
        if np.any(observed < 0) or np.any(~np.isfinite(observed)):
            raise ValueError("counts must be finite and non-negative")
        scores = self.log_scores(values)
        valid = np.isfinite(scores)
        if np.any((observed > 0) & ~valid):
            raise ValueError("positive counts cannot have invalid sequence contexts")
        totals = observed.sum(axis=1)
        keep = (totals > 0) & np.any(valid, axis=1)
        if not np.any(keep):
            raise ValueError("at least one window must contain observed cuts")
        log_probabilities = scores[keep] - logsumexp(scores[keep], axis=1, keepdims=True)
        contribution = np.where(observed[keep] > 0, observed[keep] * log_probabilities, 0.0)
        return float(-np.sum(contribution[np.isfinite(contribution)]) / np.sum(totals[keep]))

    def _gradient(
        self,
        contexts: np.ndarray,
        counts: np.ndarray,
    ) -> tuple[np.ndarray, dict[int, np.ndarray], float]:
        probabilities = self.probabilities(contexts)
        totals = counts.sum(axis=1, keepdims=True)
        error = probabilities * totals - counts
        flat_contexts = contexts.reshape(-1, contexts.shape[-1])
        flat_error = error.reshape(-1)
        valid_rows = np.all(flat_contexts < 4, axis=1)
        flat_error = np.where(valid_rows, flat_error, 0.0)
        scale = max(float(np.sum(totals)), 1.0)

        main_gradient = np.zeros_like(self.main)
        for position in range(self.feature_spec.context_length):
            codes = flat_contexts[:, position]
            valid = codes < 4
            main_gradient[position] = np.bincount(
                codes[valid], weights=flat_error[valid], minlength=4
            )[:4]
        main_gradient /= scale

        pair_gradients: dict[int, np.ndarray] = {}
        for distance, coefficients in self.pairs.items():
            gradient = np.zeros_like(coefficients)
            for position in range(self.feature_spec.context_length - distance):
                left = flat_contexts[:, position]
                right = flat_contexts[:, position + distance]
                valid = (left < 4) & (right < 4)
                code = left[valid].astype(np.int64) * 4 + right[valid].astype(np.int64)
                gradient[position] = np.bincount(
                    code, weights=flat_error[valid], minlength=16
                ).reshape(4, 4)
            pair_gradients[distance] = gradient / scale
        return main_gradient, pair_gradients, scale

    def fit(
        self,
        contexts: np.ndarray,
        counts: np.ndarray,
        *,
        epochs: int = 200,
        batch_windows: int = 64,
        learning_rate: float = 0.03,
        l2: float = 1e-3,
        prior: "ConditionalSequenceBiasModel | None" = None,
        prior_strength: float = 0.0,
        seed: int = 2026,
        tolerance: float = 1e-7,
        patience: int = 12,
    ) -> "ConditionalSequenceBiasModel":
        """Fit the model and retain a compact convergence history.

        ``prior`` supplies pooled coefficients for sample-specific adaptation.
        ``prior_strength`` penalizes deviations from those coefficients.
        """

        values = self._validate_contexts(contexts)
        observed = np.asarray(counts, dtype=np.float64)
        if observed.shape != values.shape[:2]:
            raise ValueError("counts must match the first two context dimensions")
        if np.any(observed < 0) or np.any(~np.isfinite(observed)):
            raise ValueError("counts must be finite and non-negative")
        valid_positions = np.all(values < 4, axis=2)
        if np.any((observed > 0) & ~valid_positions):
            raise ValueError("positive counts cannot have invalid sequence contexts")
        keep = (observed.sum(axis=1) > 0) & valid_positions.any(axis=1)
        values = values[keep]
        observed = observed[keep]
        if len(values) == 0:
            raise ValueError("at least one valid window with cuts is required")
        if epochs < 1 or batch_windows < 1:
            raise ValueError("epochs and batch_windows must be positive")
        if prior is not None and prior.feature_spec != self.feature_spec:
            raise ValueError("prior and fitted model must use the same feature specification")
        if prior_strength < 0 or l2 < 0:
            raise ValueError("regularization strengths must be non-negative")

        prior_main = prior.main if prior is not None else np.zeros_like(self.main)
        prior_pairs = {
            distance: prior.pairs[distance] if prior is not None else np.zeros_like(coefficients)
            for distance, coefficients in self.pairs.items()
        }
        first_main = np.zeros_like(self.main)
        second_main = np.zeros_like(self.main)
        first_pairs = {distance: np.zeros_like(value) for distance, value in self.pairs.items()}
        second_pairs = {distance: np.zeros_like(value) for distance, value in self.pairs.items()}
        rng = np.random.default_rng(seed)
        beta1, beta2, epsilon = 0.9, 0.999, 1e-8
        step = 0
        previous = np.inf
        stale = 0
        self.training_history = []

        for epoch in range(1, epochs + 1):
            order = rng.permutation(len(values))
            for start in range(0, len(order), batch_windows):
                selected = order[start:start + batch_windows]
                gradient_main, gradient_pairs, _ = self._gradient(values[selected], observed[selected])
                gradient_main += l2 * self.main + prior_strength * (self.main - prior_main)
                for distance in gradient_pairs:
                    gradient_pairs[distance] += (
                        l2 * self.pairs[distance]
                        + prior_strength * (self.pairs[distance] - prior_pairs[distance])
                    )

                step += 1
                first_main = beta1 * first_main + (1.0 - beta1) * gradient_main
                second_main = beta2 * second_main + (1.0 - beta2) * np.square(gradient_main)
                main_update = (first_main / (1.0 - beta1**step)) / (
                    np.sqrt(second_main / (1.0 - beta2**step)) + epsilon
                )
                self.main -= learning_rate * main_update
                self.main -= self.main.mean(axis=1, keepdims=True)

                for distance, gradient in gradient_pairs.items():
                    first_pairs[distance] = beta1 * first_pairs[distance] + (1.0 - beta1) * gradient
                    second_pairs[distance] = beta2 * second_pairs[distance] + (1.0 - beta2) * np.square(gradient)
                    update = (first_pairs[distance] / (1.0 - beta1**step)) / (
                        np.sqrt(second_pairs[distance] / (1.0 - beta2**step)) + epsilon
                    )
                    self.pairs[distance] -= learning_rate * update
                    self.pairs[distance] = _double_center(self.pairs[distance])

            nll = self.conditional_nll(values, observed)
            penalty = 0.5 * l2 * (
                float(np.sum(np.square(self.main)))
                + sum(float(np.sum(np.square(value))) for value in self.pairs.values())
            )
            objective = nll + penalty
            self.training_history.append({"epoch": epoch, "nll": nll, "objective": objective})
            improvement = previous - objective
            if improvement >= 0 and improvement < tolerance:
                stale += 1
            elif improvement < 0 and abs(improvement) < tolerance:
                stale += 1
            else:
                stale = 0
            previous = objective
            if stale >= patience:
                break

        self.metadata.update(
            {
                "seed": int(seed),
                "epochs_completed": int(len(self.training_history)),
                "learning_rate": float(learning_rate),
                "l2": float(l2),
                "prior_strength": float(prior_strength),
                "training_windows": int(len(values)),
                "training_cuts": float(np.sum(observed)),
            }
        )
        return self

    def save(self, path: str | Path, metadata: dict[str, Any] | None = None) -> tuple[Path, Path]:
        """Serialize coefficients to NPZ and provenance/integrity data to JSON."""

        npz_path = Path(path)
        if npz_path.suffix != ".npz":
            npz_path = npz_path.with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        arrays: dict[str, np.ndarray] = {"main": self.main}
        arrays.update({f"pair_{distance}": value for distance, value in self.pairs.items()})
        np.savez_compressed(npz_path, **arrays)
        document = {
            "schema": MODEL_SCHEMA,
            "feature_spec": asdict(self.feature_spec),
            "npz_sha256": _sha256_file(npz_path),
            "training_history": self.training_history,
            "metadata": {**self.metadata, **dict(metadata or {})},
        }
        json_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "ConditionalSequenceBiasModel":
        """Load a versioned model after validating its checksum and shapes."""

        npz_path = Path(path)
        if npz_path.suffix != ".npz":
            npz_path = npz_path.with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        document = json.loads(json_path.read_text(encoding="utf-8"))
        if document.get("schema") != MODEL_SCHEMA:
            raise ValueError(f"unsupported bias model schema: {document.get('schema')}")
        if document.get("npz_sha256") != _sha256_file(npz_path):
            raise ValueError("bias model checksum does not match its JSON metadata")
        raw_spec = document["feature_spec"]
        spec = BiasFeatureSpec(
            str(raw_spec["name"]),
            int(raw_spec["context_length"]),
            tuple(int(value) for value in raw_spec["pair_distances"]),
        )
        model = cls(spec)
        with np.load(npz_path, allow_pickle=False) as arrays:
            if arrays["main"].shape != model.main.shape:
                raise ValueError("bias model main-effect shape is incompatible with its metadata")
            model.main = np.asarray(arrays["main"], dtype=np.float64)
            for distance in spec.pair_distances:
                key = f"pair_{distance}"
                if key not in arrays or arrays[key].shape != model.pairs[distance].shape:
                    raise ValueError(f"bias model interaction {distance} is missing or malformed")
                model.pairs[distance] = np.asarray(arrays[key], dtype=np.float64)
        model.training_history = list(document.get("training_history", []))
        model.metadata = dict(document.get("metadata", {}))
        return model


def expected_from_log_bias(
    observed: np.ndarray,
    log_bias: np.ndarray,
    window: int = 100,
    epsilon: float = 1e-12,
) -> np.ndarray:
    """Scale sequence propensities toward rolling local observed cut totals."""

    cuts = np.asarray(observed, dtype=np.float64)
    scores = np.asarray(log_bias, dtype=np.float64)
    if cuts.shape != scores.shape or cuts.ndim != 1:
        raise ValueError("observed and log_bias must be one-dimensional arrays of equal length")
    if window < 1 or np.any(cuts < 0):
        raise ValueError("window must be positive and observed cuts non-negative")
    finite = np.isfinite(scores)
    if not finite.any():
        return np.zeros_like(cuts)
    centered = np.where(finite, scores - np.max(scores[finite]), -np.inf)
    propensity = np.exp(centered)
    kernel = np.ones(int(window), dtype=np.float64)
    local_cuts = np.convolve(cuts, kernel, mode="same")
    local_propensity = np.convolve(propensity, kernel, mode="same")
    expected = propensity * local_cuts / np.maximum(local_propensity, epsilon)
    expected[~np.isfinite(expected)] = 0.0
    expected_sum = float(expected.sum())
    observed_sum = float(cuts.sum())
    if expected_sum > 0:
        expected *= observed_sum / expected_sum
    return expected


def estimate_nb_dispersion(observed: np.ndarray, expected: np.ndarray) -> float:
    """Estimate NB2 dispersion by a non-negative method-of-moments equation."""

    y = np.asarray(observed, dtype=np.float64)
    mu = np.asarray(expected, dtype=np.float64)
    finite = np.isfinite(y) & np.isfinite(mu) & (y >= 0) & (mu > 0)
    if not finite.any():
        return 0.0
    numerator = np.sum(np.square(y[finite] - mu[finite]) - y[finite])
    denominator = np.sum(np.square(mu[finite]))
    return float(max(0.0, numerator / denominator)) if denominator > 0 else 0.0


def calibrated_residuals(
    observed: np.ndarray,
    expected: np.ndarray,
    mode: str,
    *,
    dispersion: float = 0.0,
    pseudocount: float = 0.5,
) -> np.ndarray:
    """Calculate raw, Pearson, deviance, or log-ratio cut-site residuals."""

    y = np.asarray(observed, dtype=np.float64)
    mu = np.asarray(expected, dtype=np.float64)
    if y.shape != mu.shape:
        raise ValueError("observed and expected arrays must have equal shapes")
    if np.any(y < 0) or np.any(mu < 0):
        raise ValueError("observed and expected signals must be non-negative")
    if dispersion < 0 or pseudocount <= 0:
        raise ValueError("dispersion must be non-negative and pseudocount positive")
    name = mode.lower().replace("_", "-")
    if name in {"raw", "difference"}:
        return y - mu
    if name == "pearson":
        variance = mu + dispersion * np.square(mu)
        return (y - mu) / np.sqrt(np.maximum(variance, np.finfo(float).eps))
    if name in {"log-ratio", "logratio"}:
        return np.log((y + pseudocount) / (mu + pseudocount))
    if name not in {"deviance", "signed-deviance"}:
        raise ValueError(f"unknown residual mode: {mode}")

    safe_mu = np.maximum(mu, np.finfo(float).tiny)
    with np.errstate(divide="ignore", invalid="ignore"):
        first = np.where(y > 0, y * np.log(y / safe_mu), 0.0)
        if dispersion == 0:
            deviance = 2.0 * (first - (y - safe_mu))
        else:
            size = 1.0 / dispersion
            second = (y + size) * np.log((y + size) / (safe_mu + size))
            deviance = 2.0 * (first - second)
    deviance = np.maximum(np.nan_to_num(deviance, nan=0.0, posinf=0.0, neginf=0.0), 0.0)
    return np.sign(y - mu) * np.sqrt(deviance)


def negative_binomial_log_likelihood(
    observed: np.ndarray,
    expected: np.ndarray,
    dispersion: float = 0.0,
) -> np.ndarray:
    """Elementwise Poisson/NB2 log likelihood, excluding no constants."""

    y = np.asarray(observed, dtype=np.float64)
    mu = np.maximum(np.asarray(expected, dtype=np.float64), np.finfo(float).tiny)
    if y.shape != mu.shape or np.any(y < 0) or dispersion < 0:
        raise ValueError("invalid observed/expected arrays or dispersion")
    if dispersion == 0:
        return y * np.log(mu) - mu - gammaln(y + 1.0)
    size = 1.0 / dispersion
    return (
        gammaln(y + size)
        - gammaln(size)
        - gammaln(y + 1.0)
        + size * np.log(size / (size + mu))
        + y * np.log(mu / (size + mu))
    )


def center_flank_likelihood_score(
    observed_profiles: np.ndarray,
    expected_profiles: np.ndarray,
    *,
    center_width: int = 15,
    flank_width: int = 30,
    gap: int = 5,
    dispersion: float = 0.0,
) -> np.ndarray:
    """Signed NB likelihood-ratio score for central protection versus expectation."""

    observed = np.asarray(observed_profiles, dtype=np.float64)
    expected = np.asarray(expected_profiles, dtype=np.float64)
    if observed.shape != expected.shape or observed.ndim != 2:
        raise ValueError("profiles must be equal two-dimensional arrays")
    width = observed.shape[1]
    midpoint = width // 2
    half_center = max(1, center_width // 2)
    center = np.zeros(width, dtype=bool)
    center[max(0, midpoint - half_center):min(width, midpoint + half_center + 1)] = True
    flanks = np.zeros(width, dtype=bool)
    left_end = max(0, midpoint - half_center - gap)
    left_start = max(0, left_end - flank_width)
    right_start = min(width, midpoint + half_center + gap + 1)
    right_end = min(width, right_start + flank_width)
    flanks[left_start:left_end] = True
    flanks[right_start:right_end] = True
    if not center.any() or not flanks.any():
        raise ValueError("profile is too short for the requested center/flank geometry")

    epsilon = np.finfo(float).eps
    observed_center = observed[:, center].sum(axis=1)
    expected_center = expected[:, center].sum(axis=1)
    observed_flank = observed[:, flanks].sum(axis=1)
    expected_flank = expected[:, flanks].sum(axis=1)
    flank_scale = (observed_flank + 0.5) / (expected_flank + 0.5)
    null = expected * flank_scale[:, None]
    center_scale = (observed_center + 0.5) / (expected_center + 0.5)
    alternative = null.copy()
    alternative[:, center] = expected[:, center] * center_scale[:, None]
    null_ll = negative_binomial_log_likelihood(observed, np.maximum(null, epsilon), dispersion).sum(axis=1)
    alternative_ll = negative_binomial_log_likelihood(
        observed, np.maximum(alternative, epsilon), dispersion
    ).sum(axis=1)
    statistic = np.maximum(0.0, 2.0 * (alternative_ll - null_ll))
    direction = np.sign(expected_center * flank_scale - observed_center)
    return direction * np.sqrt(statistic)
