from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "fit_tobias_dwm_control_reference.py"
spec = importlib.util.spec_from_file_location(
    "fit_tobias_dwm_control_reference", SCRIPT
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def make_dataset(name: str, split: str, seed: int):
    rng = np.random.default_rng(seed)
    width = 20
    margin = 41
    count = 10
    sequences = np.asarray(
        ["".join(rng.choice(list("ACGT"), width + 2 * margin)) for _ in range(count)]
    )
    return module.ControlWindowDataset(
        sample=name,
        split=split,
        source="fixture",
        shift=(4, -4),
        window_size=width,
        margin=margin,
        chromosomes=np.repeat("chr16", count),
        starts=np.arange(count) * width,
        sequences=sequences,
        forward_counts=rng.poisson(2, size=(count, width)),
        reverse_counts=rng.poisson(2, size=(count, width)),
        gc_fraction=np.full(count, 0.5),
    )


def test_fit_reference_scores_each_validation_library(tmp_path) -> None:
    shift = (4, -4)
    training = {
        (shift, "naked"): make_dataset("naked", "train", 1),
        (shift, "mito"): make_dataset("mito", "train", 2),
    }
    validation = {
        (shift, "naked"): make_dataset("naked", "validation", 3),
        (shift, "mito"): make_dataset("mito", "validation", 4),
    }
    artifacts, metrics = module.fit_references(training, validation, tmp_path)
    assert len(artifacts) == 1
    assert set(metrics["library"]) == {"naked", "mito"}
    model = module.TobiasDwmReferenceModel.load(artifacts.loc[0, "model_npz"])
    assert model.metadata["configuration"] == "conventional_tobias_style_dwm"
