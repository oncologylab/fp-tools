from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "evaluate_shifted_atac_null_policy.py"
spec = importlib.util.spec_from_file_location("evaluate_shifted_atac_null_policy", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_cyclic_shift_preserves_rows_and_moves_center() -> None:
    values = np.arange(15, dtype=float).reshape(3, 5)
    shifted = module.cyclic_shift_profiles({"observed": values}, 2)["observed"]
    assert np.array_equal(shifted, np.roll(values, 2, axis=1))
    assert np.array_equal(shifted.sum(axis=1), values.sum(axis=1))


def test_cyclic_shift_rejects_zero_and_bad_shapes() -> None:
    with pytest.raises(ValueError, match="zero"):
        module.cyclic_shift_profiles({"observed": np.ones((2, 3))}, 0)
    with pytest.raises(ValueError, match="two-dimensional"):
        module.cyclic_shift_profiles({"observed": np.ones(3)}, 1)


def test_hash_folds_are_deterministic_and_complete() -> None:
    hashes = np.array([1, 2, 3, 4, 5], dtype=np.uint64)
    folds = np.mod(hashes, 4).astype(int)
    assert folds.tolist() == [1, 2, 3, 0, 1]
    assert sum(int(np.sum(folds == fold)) for fold in range(4)) == len(hashes)
