from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from diagnose_low_count_normalization import components, count_stratum, POSITIONS
from fixed_shape_normalization import normalize_fixed


def test_components_reproduce_normalization_and_zero_fallback():
    values = np.zeros((3, 201))
    values[1, 0] = 1
    values[2] = values[1] * 1e-4
    c = components(values)
    assert c['fallback'].tolist() == [True, False, False]
    assert c['denominator'][0] == 1
    np.testing.assert_allclose(c['amplification'][2], c['amplification'][1] * 1e4)
    np.testing.assert_allclose(normalize_fixed(values[1:2], POSITIONS, 1), normalize_fixed(values[2:3], POSITIONS, 1), atol=1e-12)


def test_count_bins_include_zero_and_exact_boundaries():
    assert count_stratum([0, 1, 2, 4, 5, 9, 10, 19, 20, 49, 50, 99, 100, 101]).tolist() == [0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7]
    with pytest.raises(ValueError):
        count_stratum([-1])
