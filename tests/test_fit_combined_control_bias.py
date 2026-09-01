from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from fp_tools.tools.parametric_bias import BiasFeatureSpec


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "fit_combined_control_bias.py"
spec = importlib.util.spec_from_file_location("fit_combined_control_bias", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def dataset(name: str, split: str, shift: tuple[int, int], seed: int):
    rng = np.random.default_rng(seed)
    width = 20
    margin = 41
    count = 8
    sequence_length = width + 2 * margin
    sequences = np.asarray(
        ["".join(rng.choice(list("ACGT"), sequence_length)) for _ in range(count)]
    )
    forward = rng.poisson(1.0, size=(count, width))
    reverse = rng.poisson(1.0, size=(count, width))
    return module.ControlWindowDataset(
        sample=name,
        split=split,
        source="fixture",
        shift=shift,
        window_size=width,
        margin=margin,
        chromosomes=np.repeat("chr1", count),
        starts=np.arange(count) * width,
        sequences=sequences,
        forward_counts=forward,
        reverse_counts=reverse,
        gc_fraction=np.full(count, 0.5),
    )


def test_combined_grid_fits_all_sources_and_builds_safe_ensemble(tmp_path) -> None:
    shift = (4, -4)
    training = {
        (shift, "naked"): dataset("naked", "train", shift, 1),
        (shift, "mito"): dataset("mito", "train", shift, 2),
    }
    validation = {
        (shift, "naked"): dataset("naked", "validation", shift, 3),
        (shift, "mito"): dataset("mito", "validation", shift, 4),
    }
    artifacts, ensembles, metrics = module.fit_combined_grid(
        training,
        validation,
        tmp_path,
        models=["selma10"],
        l2_values=[1e-3],
        seeds=[2026, 2027],
        epochs=2,
        batch_windows=8,
        jobs=2,
    )
    assert len(artifacts) == 2
    assert len(ensembles) == 1
    assert len(metrics) == 2
    assert set(metrics["library"]) == {"naked", "mito"}
    model = module.ConditionalSequenceBiasModel.load(ensembles.loc[0, "model_npz"])
    assert model.feature_spec == BiasFeatureSpec.selma10()
    with np.load(ensembles.loc[0, "model_npz"], allow_pickle=False) as arrays:
        assert all(arrays[key].dtype != object for key in arrays.files)

    member_mtimes = {
        Path(path): Path(path).stat().st_mtime_ns for path in artifacts["model_npz"]
    }
    ensemble_mtime = Path(ensembles.loc[0, "model_npz"]).stat().st_mtime_ns
    resumed_artifacts, resumed_ensembles, resumed_metrics = module.fit_combined_grid(
        training,
        validation,
        tmp_path,
        models=["selma10"],
        l2_values=[1e-3],
        seeds=[2026, 2027],
        epochs=2,
        batch_windows=8,
        jobs=2,
    )
    assert resumed_artifacts["resumed"].all()
    assert resumed_ensembles["resumed"].all()
    assert len(resumed_metrics) == len(metrics)
    assert all(path.stat().st_mtime_ns == mtime for path, mtime in member_mtimes.items())
    assert Path(ensembles.loc[0, "model_npz"]).stat().st_mtime_ns == ensemble_mtime


def test_combined_grid_rejects_nonpositive_jobs(tmp_path) -> None:
    shift = (4, -4)
    training = {
        (shift, "naked"): dataset("naked", "train", shift, 1),
        (shift, "mito"): dataset("mito", "train", shift, 2),
    }
    validation = {
        (shift, "naked"): dataset("naked", "validation", shift, 3),
        (shift, "mito"): dataset("mito", "validation", shift, 4),
    }
    try:
        module.fit_combined_grid(
            training,
            validation,
            tmp_path,
            models=["selma10"],
            l2_values=[1e-3],
            seeds=[2026],
            epochs=1,
            batch_windows=8,
            jobs=0,
        )
    except ValueError as exc:
        assert str(exc) == "jobs must be positive"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected a ValueError")
