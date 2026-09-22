from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from crossfit_unit_floor_geometry import normalize_unit_floor
from fixed_shape_normalization import detrend


def test_floor_never_amplifies_and_keeps_zero_profiles_zero():
    rng = np.random.default_rng(2026)
    values = rng.normal(size=(5, 201)) * np.array([0, .001, .1, 1, 20])[:, None]
    residual, scale = detrend(values, np.arange(-100, 101))
    actual = normalize_unit_floor(values)
    np.testing.assert_allclose(actual, np.clip(residual / np.maximum(scale, 1)[:, None], -12, 12))
    assert (abs(actual) <= abs(residual) + 1e-12).all()
    assert not actual[0].any()


def test_floor_does_not_depend_on_batch_and_retains_sparse_amplitude():
    values = np.zeros((2, 201))
    values[0, 0] = 1
    values[1] = values[0] * .001
    result = normalize_unit_floor(values)
    np.testing.assert_allclose(result[1], result[0] * .001, atol=1e-15)
    for i in range(2):
        np.testing.assert_array_equal(normalize_unit_floor(values[i:i+1])[0], result[i])
