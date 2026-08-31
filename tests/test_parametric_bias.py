from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from fp_tools.tools.parametric_bias import (
    BiasFeatureSpec,
    ConditionalSequenceBiasModel,
    calibrated_residuals,
    center_flank_likelihood_score,
    contexts_from_sequence,
    encode_sequence,
    estimate_nb_dispersion,
    expected_from_log_bias,
    reverse_complement_contexts,
)


def test_sequence_encoding_and_orientation() -> None:
    assert encode_sequence("ACGTNacgt").tolist() == [0, 1, 2, 3, 4, 0, 1, 2, 3]
    contexts = np.asarray([[0, 1, 2, 3, 4]], dtype=np.uint8)
    assert reverse_complement_contexts(contexts).tolist() == [[4, 0, 1, 2, 3]]

    extracted, valid = contexts_from_sequence(
        "AACCGGTTAACC",
        [4, 4, 0],
        5,
        strands=["+", "-", "+"],
    )
    assert valid.tolist() == [True, True, False]
    assert extracted[1].tolist() == reverse_complement_contexts(extracted[[0]])[0].tolist()


def _synthetic_conditional_data(seed: int = 7) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    contexts = rng.integers(0, 4, size=(180, 14, 5), dtype=np.uint8)
    center = contexts[:, :, 2]
    neighbor = contexts[:, :, 3]
    scores = 1.4 * (center == 0) - 0.9 * (center == 3) + 0.7 * (
        (center == 1) & (neighbor == 2)
    )
    probabilities = np.exp(scores - scores.max(axis=1, keepdims=True))
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    counts = np.vstack([rng.multinomial(35, row) for row in probabilities]).astype(float)
    return contexts, counts


def test_conditional_model_recovers_sequence_signal() -> None:
    contexts, counts = _synthetic_conditional_data()
    model = ConditionalSequenceBiasModel(BiasFeatureSpec("test", 5, (1,)))
    null_nll = model.conditional_nll(contexts, counts)
    model.fit(
        contexts,
        counts,
        epochs=100,
        batch_windows=45,
        learning_rate=0.04,
        l2=2e-3,
        seed=11,
        patience=20,
    )
    fitted_nll = model.conditional_nll(contexts, counts)
    assert fitted_nll < null_nll - 0.08
    assert model.main.shape == (5, 4)
    assert np.allclose(model.main.mean(axis=1), 0.0, atol=1e-10)
    assert model.training_history[-1]["nll"] == pytest.approx(fitted_nll)


def test_pooled_prior_adaptation_and_safe_roundtrip(tmp_path: Path) -> None:
    contexts, counts = _synthetic_conditional_data()
    pooled = ConditionalSequenceBiasModel(BiasFeatureSpec("test", 5, (1,))).fit(
        contexts[:120], counts[:120], epochs=25, batch_windows=40, seed=2
    )
    adapted = ConditionalSequenceBiasModel(pooled.feature_spec).fit(
        contexts[120:],
        counts[120:],
        epochs=20,
        batch_windows=30,
        prior=pooled,
        prior_strength=0.2,
        seed=3,
    )
    npz_path, json_path = adapted.save(
        tmp_path / "bias_model.npz",
        metadata={"read_shift": [4, -4], "training_source": "synthetic"},
    )
    assert npz_path.is_file() and json_path.is_file()
    loaded = ConditionalSequenceBiasModel.load(npz_path)
    assert np.allclose(loaded.main, adapted.main)
    assert np.allclose(loaded.log_scores(contexts[:2]), adapted.log_scores(contexts[:2]))
    assert loaded.metadata["read_shift"] == [4, -4]

    with npz_path.open("ab") as handle:
        handle.write(b"corrupt")
    with pytest.raises(ValueError, match="checksum"):
        ConditionalSequenceBiasModel.load(npz_path)


def test_calibrated_residuals_and_expected_signal() -> None:
    observed = np.asarray([0.0, 2.0, 8.0, 1.0, 0.0])
    expected = np.asarray([0.5, 2.0, 3.0, 2.0, 0.5])
    assert np.allclose(calibrated_residuals(observed, expected, "raw"), observed - expected)
    pearson = calibrated_residuals(observed, expected, "pearson", dispersion=0.1)
    deviance = calibrated_residuals(observed, expected, "deviance", dispersion=0.1)
    log_ratio = calibrated_residuals(observed, expected, "log-ratio")
    assert np.all(np.isfinite(pearson))
    assert np.all(np.isfinite(deviance))
    assert np.all(np.isfinite(log_ratio))
    assert deviance[2] > 0 and deviance[0] < 0
    assert estimate_nb_dispersion(observed, expected) >= 0

    predicted = expected_from_log_bias(
        np.asarray([0.0, 1.0, 4.0, 1.0, 0.0]),
        np.log(np.asarray([1.0, 1.0, 5.0, 1.0, 1.0])),
        window=5,
    )
    assert predicted[2] > predicted[1]
    assert predicted.sum() == pytest.approx(6.0)


def test_center_flank_likelihood_detects_protection() -> None:
    expected = np.full((2, 101), 5.0)
    observed = expected.copy()
    observed[0, 45:56] = 1.0
    observed[1, 45:56] = 9.0
    scores = center_flank_likelihood_score(
        observed,
        expected,
        center_width=11,
        flank_width=20,
        gap=5,
    )
    assert scores[0] > 0
    assert scores[1] < 0

