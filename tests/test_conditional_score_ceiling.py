from itertools import product
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from audit_conditional_score_ceiling import score_ceiling
from conditional_protection_score import depletion_score


def test_analytic_ceiling_matches_zero_center_and_bounds_all_counts():
    for n, m, p, q in product([0, 1, 2, 5, 10], [0, 1, 4], [0., .2, .8, 1.], [.1, .7, 1.]):
        ceiling = float(score_ceiling(n, m, p, q))
        np.testing.assert_allclose(ceiling, np.round(depletion_score(n, m, 0, p, q), 10), atol=1.01e-10, rtol=0)
        for center in range(n+m+1):
            assert np.round(depletion_score(n, m, center, p, q), 10) <= ceiling + 1.01e-10


def test_zero_times_log_zero_is_zero_and_ceiling_is_capped():
    np.testing.assert_array_equal(score_ceiling([0, 1], [0, 0], [1., 1.], [1., 1.]), [0., 300.])
    with pytest.raises(ValueError):
        score_ceiling(.5, 1, .2, .3)
    with pytest.raises(ValueError):
        score_ceiling(1, 1, np.nan, .3)
