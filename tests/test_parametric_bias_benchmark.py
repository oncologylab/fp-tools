from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_parametric_bias import (  # noqa: E402
    ControlWindowDataset,
    conditional_metrics,
    cut_position,
    gc_matched_indices,
    select_bias_configurations,
    split_mitochondrial_dataset,
    stable_u64,
    summarize_bias_depth_stability,
    thin_counts,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    BiasFeatureSpec,
    ConditionalSequenceBiasModel,
)


def _dataset(n_windows: int = 12, seed: int = 7) -> ControlWindowDataset:
    rng = np.random.default_rng(seed)
    window_size = 20
    margin = 41
    sequences = np.asarray(
        ["".join(rng.choice(list("ACGT"), size=window_size + 2 * margin)) for _ in range(n_windows)]
    )
    forward = rng.poisson(0.5, size=(n_windows, window_size)).astype(np.uint16)
    reverse = rng.poisson(0.5, size=(n_windows, window_size)).astype(np.uint16)
    return ControlWindowDataset(
        sample="sample",
        split="train",
        source="synthetic",
        shift=(4, -5),
        window_size=window_size,
        margin=margin,
        chromosomes=np.asarray(["chrM"] * n_windows),
        starts=np.arange(n_windows) * window_size + margin,
        sequences=sequences,
        forward_counts=forward,
        reverse_counts=reverse,
        gc_fraction=np.full(n_windows, 0.5),
    )


def test_cut_position_preserves_one_read_soft_clip_convention() -> None:
    forward = SimpleNamespace(
        query_length=50,
        query_alignment_start=3,
        query_alignment_end=50,
        reference_start=100,
        reference_end=147,
        is_reverse=False,
        infer_query_length=lambda: 50,
    )
    reverse = SimpleNamespace(**{**vars(forward), "is_reverse": True})
    assert cut_position(forward, (4, -5)) == 101
    assert cut_position(reverse, (4, -5)) == 142
    assert cut_position(reverse, (4, -4)) == 143


def test_control_dataset_builds_oriented_model_arrays_and_roundtrips(tmp_path: Path) -> None:
    dataset = _dataset()
    contexts, counts = dataset.model_arrays(BiasFeatureSpec.selma10())
    assert contexts.shape[1:] == (20, 10)
    assert counts.shape == contexts.shape[:2]
    assert len(contexts) <= len(dataset.starts) * 2
    npz, metadata = dataset.save(tmp_path / "windows")
    loaded = ControlWindowDataset.load(npz)
    assert metadata.is_file()
    assert np.array_equal(loaded.forward_counts, dataset.forward_counts)
    assert np.array_equal(loaded.sequences, dataset.sequences)

    with npz.open("ab") as handle:
        handle.write(b"broken")
    with pytest.raises(ValueError, match="checksum"):
        ControlWindowDataset.load(npz)


def test_gc_matching_moves_candidates_toward_target() -> None:
    candidate = np.linspace(0.0, 1.0, 1000)
    target = np.linspace(0.72, 0.88, 200)
    selected = gc_matched_indices(candidate, target, 150, seed=11)
    random_selected = gc_matched_indices(candidate, np.asarray([]), 150, seed=11)
    assert len(selected) == 150
    assert abs(np.mean(candidate[selected]) - np.mean(target)) < abs(
        np.mean(candidate[random_selected]) - np.mean(target)
    )


def test_mitochondrial_hash_split_and_thinning_are_deterministic() -> None:
    dataset = _dataset(20)
    train_a, validation_a = split_mitochondrial_dataset(dataset, 14, 6, seed=5)
    train_b, validation_b = split_mitochondrial_dataset(dataset, 14, 6, seed=5)
    assert np.array_equal(train_a.starts, train_b.starts)
    assert np.array_equal(validation_a.starts, validation_b.starts)
    assert not set(train_a.starts).intersection(validation_a.starts)

    counts = np.full((20, 10), 5.0)
    first = thin_counts(counts, 250, seed=9)
    second = thin_counts(counts, 250, seed=9)
    assert np.array_equal(first, second)
    assert first.sum() < counts.sum()
    assert stable_u64("x", 1) == stable_u64("x", 1)


def test_conditional_metrics_reward_fitted_sequence_model() -> None:
    rng = np.random.default_rng(13)
    contexts = rng.integers(0, 4, size=(80, 20, 5), dtype=np.uint8)
    center = contexts[:, :, 2]
    scores = 1.5 * (center == 0) - 1.0 * (center == 3)
    probabilities = np.exp(scores - scores.max(axis=1, keepdims=True))
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    counts = np.vstack([rng.multinomial(40, row) for row in probabilities]).astype(float)
    model = ConditionalSequenceBiasModel(BiasFeatureSpec("test", 5, (1,))).fit(
        contexts,
        counts,
        epochs=50,
        batch_windows=40,
        seed=3,
    )
    metrics = conditional_metrics(model, contexts, counts)
    assert metrics["nll_gain"] > 0.05
    assert metrics["conditional_nll"] < metrics["null_nll"]
    assert metrics["multinomial_deviance_per_cut"] >= 0
    assert np.isfinite(metrics["aggregate_jsd"])


def test_selection_uses_validation_only_and_keeps_two() -> None:
    rows = []
    for model, nll in (("selma10", 4.0), ("loglinear81", 3.8), ("other", 4.2)):
        for sample in ("A", "B"):
            for split in ("train", "validation"):
                rows.append(
                    {
                        "source": "nonpeak",
                        "sample": sample,
                        "split": split,
                        "shift_forward": 4,
                        "shift_reverse": -4,
                        "model": model,
                        "configuration": "pooled",
                        "l2": 0.001,
                        "training_depth": "full",
                        "seed": 1,
                        "conditional_nll": nll if split == "validation" else 0.1,
                        "nll_gain": 0.2,
                        "multinomial_deviance_per_cut": 1.0,
                        "calibration_error": 0.1,
                        "runtime_seconds": 2.0,
                        "model_size_mb": 1.0,
                    }
                )
    selection = select_bias_configurations(pd.DataFrame(rows))
    retained = selection[selection["retained_for_functional_screen"]]
    assert retained["model"].tolist() == ["loglinear81", "selma10"]
    assert retained["rank"].tolist() == [1, 2]


def test_selection_rejects_models_that_do_not_beat_uniform_control() -> None:
    frame = pd.DataFrame(
        [
            {
                "source": "nonpeak",
                "sample": "A",
                "split": "validation",
                "shift_forward": 4,
                "shift_reverse": -5,
                "model": "selma10",
                "configuration": "sample",
                "l2": 0.001,
                "training_depth": "full",
                "seed": 1,
                "conditional_nll": 5.4,
                "nll_gain": -0.1,
                "multinomial_deviance_per_cut": 2.0,
                "calibration_error": 0.2,
                "runtime_seconds": 1.0,
                "model_size_mb": 1.0,
            }
        ]
    )
    selected = select_bias_configurations(frame)
    assert not selected.loc[0, "passed_control_likelihood"]
    assert not selected.loc[0, "retained_for_functional_screen"]


def test_depth_recommendation_aggregates_seeds_and_uses_smallest_stable_depth() -> None:
    rows = []
    means = {"1000": 4.20, "5000": 4.002, "25000": 4.00, "full": 4.01}
    for depth, mean in means.items():
        for seed, offset in enumerate((-0.01, -0.005, 0.0, 0.005, 0.01), start=1):
            for sample in ("A", "B"):
                rows.append(
                    {
                        "source": "mitochondrial",
                        "sample": sample,
                        "split": "validation",
                        "shift_forward": 4,
                        "shift_reverse": -4,
                        "model": "selma10",
                        "configuration": "pooled",
                        "l2": 0.001,
                        "training_depth": depth,
                        "seed": seed,
                        "conditional_nll": mean + offset,
                        "nll_gain": 4.6 - mean - offset,
                        "calibration_error": 0.1,
                        "runtime_seconds": 1.0,
                        "model_size_mb": 1.0,
                    }
                )
    stability, recommendations = summarize_bias_depth_stability(pd.DataFrame(rows))
    assert len(recommendations) == 1
    assert str(recommendations.loc[0, "training_depth"]) == "5000"
    selected = stability[stability["recommended_minimum_depth"]]
    assert selected["seed_count"].tolist() == [5]
    assert selected["passed_all_seed_control_likelihood"].all()
