from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from audit_rep1_region_calibration import region_probabilities, predictive_interval


def test_uniform_categories_partition_and_logits_are_shift_invariant():
    x = np.zeros((2, 201)); x[1] += 1000
    p, valid = region_probabilities(x)
    np.testing.assert_allclose(p, np.tile(np.array([15, 29, 29])/73, (2, 1)))
    np.testing.assert_allclose(p.sum(1), 1)
    assert valid.all()
    x[1, 100] = np.nan
    p, valid = region_probabilities(x)
    assert valid.tolist() == [True, False] and np.isnan(p[1]).all()


def test_predictive_interval_matches_binomial_and_retains_zero_count_case():
    assert predictive_interval([0, 0], [.1, .9]) == (0, 0)
    assert predictive_interval([3, 7], [.2, .2]) == tuple(binom.ppf([.025, .975], 10, .2).astype(int))
    assert predictive_interval([3, 7], [1., 0.]) == (3, 3)
    with pytest.raises(ValueError):
        predictive_interval([1.5], [.2])
