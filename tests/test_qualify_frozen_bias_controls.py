from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from fp_tools.tools.parametric_bias import BiasFeatureSpec, ConditionalSequenceBiasModel


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "qualify_frozen_bias_controls.py"
spec = importlib.util.spec_from_file_location("qualify_frozen_bias_controls", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def make_dataset(path: Path, name: str, seed: int) -> Path:
    rng = np.random.default_rng(seed)
    width = 20
    margin = 41
    count = 12
    sequences = np.asarray(
        ["".join(rng.choice(list("ACGT"), width + 2 * margin)) for _ in range(count)]
    )
    forward = rng.poisson(2, size=(count, width))
    reverse = rng.poisson(2, size=(count, width))
    dataset = module.ControlWindowDataset(
        sample=name,
        split="validation",
        source="fixture",
        shift=(4, -4),
        window_size=width,
        margin=margin,
        chromosomes=np.repeat(["chr16", "chr17", "chr18"], 4),
        starts=np.arange(count) * width,
        sequences=sequences,
        forward_counts=forward,
        reverse_counts=reverse,
        gc_fraction=np.full(count, 0.5),
    )
    dataset.save(path)
    return path


def save_model(path: Path, strength: float) -> Path:
    model = ConditionalSequenceBiasModel(BiasFeatureSpec.selma10())
    model.main[4, 0] = strength
    model.main[4] -= model.main[4].mean()
    npz, _json = model.save(path, metadata={"read_shift": [4, -4]})
    return npz


def test_control_qualification_pairs_windows_and_is_deterministic(tmp_path) -> None:
    first = make_dataset(tmp_path / "first.npz", "first", 1)
    second = make_dataset(tmp_path / "second.npz", "second", 2)
    reference = save_model(tmp_path / "reference", 0.0)
    candidate = save_model(tmp_path / "candidate", 0.2)
    result = module.qualify(
        [("reference", reference), ("candidate", candidate)],
        [("first", first), ("second", second)],
        reference_id="reference",
        bootstraps=100,
        seed=2026,
        maximum_model_size_mb=25,
    )
    windows, libraries, paired, selection = result
    assert set(windows["candidate_id"]) == {"reference", "candidate"}
    assert len(libraries) == 4
    assert len(paired) == 1
    assert paired.loc[0, "paired_support_fraction"] == 1.0
    assert len(selection) == 2
    repeated = module.qualify(
        [("reference", reference), ("candidate", candidate)],
        [("first", first), ("second", second)],
        reference_id="reference",
        bootstraps=100,
        seed=2026,
        maximum_model_size_mb=25,
    )
    assert paired.equals(repeated[2])


def test_model_arrays_keep_stable_strand_window_keys(tmp_path) -> None:
    dataset_path = make_dataset(tmp_path / "data.npz", "data", 3)
    dataset = module.ControlWindowDataset.load(dataset_path)
    model = ConditionalSequenceBiasModel(BiasFeatureSpec.selma10())
    contexts, counts, keys, blocks = module.model_arrays_with_keys(
        dataset, model, "library"
    )
    assert len(contexts) == len(counts) == len(keys) == len(blocks)
    assert len(keys) == len(set(keys))
    assert any("|forward|" in value for value in keys)
    assert any("|reverse|" in value for value in keys)


def test_model_arrays_exclude_unresolved_windows_on_common_support(tmp_path) -> None:
    dataset_path = make_dataset(tmp_path / "data.npz", "data", 4)
    dataset = module.ControlWindowDataset.load(dataset_path)
    dataset.sequences[3] = "N" + str(dataset.sequences[3])[1:]
    model = ConditionalSequenceBiasModel(BiasFeatureSpec.selma10())
    _contexts, _counts, keys, _blocks = module.model_arrays_with_keys(
        dataset, model, "library"
    )
    assert not any(value.endswith("|3") for value in keys)


def test_qualifier_loads_safe_dwm_reference(tmp_path) -> None:
    dataset_path = make_dataset(tmp_path / "data.npz", "data", 5)
    dataset = module.ControlWindowDataset.load(dataset_path)
    dwm = module.TobiasDwmReferenceModel()
    contexts, counts = dataset.model_arrays(dwm.feature_spec)
    dwm.fit(contexts, counts)
    model_path, _ = dwm.save(tmp_path / "dwm", metadata={"read_shift": [4, -4]})
    restored = module.load_control_model(model_path)
    assert isinstance(restored, module.TobiasDwmReferenceModel)
