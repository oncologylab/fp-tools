from itertools import product
from math import comb
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from conditional_protection_score import depletion_score, profile_scores, region_components


def mass(n, p):
    return np.array([comb(n, k) * p**k * (1-p)**(n-k) for k in range(n+1)])


def test_exact_enumeration_tail_and_superuniform_null():
    for n, m, p, q in product(range(6), range(6), [0., .01, .2, .5, .8, .99, 1.], [.1, .5, .9]):
        probabilities = np.convolve(mass(n, p), mass(m, q))
        tails = probabilities.cumsum()
        scores = np.array([depletion_score(n, m, k, p, q) for k in range(n+m+1)])
        np.testing.assert_allclose(10**(-scores), tails, rtol=1e-11, atol=1e-14)
        assert np.all(np.diff(scores) <= 1e-12)
        for alpha in [.01, .05, .1, .5]:
            assert probabilities[scores > -np.log10(alpha) + 1e-12].sum() <= alpha + 1e-12
            assert probabilities[np.round(scores, 10) > np.round(-np.log10(alpha), 10)].sum() <= alpha + 1e-12


def test_strand_exchange_and_zero_counts():
    assert depletion_score(0, 0, 0, .1, .9) == 0
    assert depletion_score(4, 7, 2, .2, .6) == depletion_score(7, 4, 2, .6, .2)
    assert depletion_score(10, 10, 20, .5, .5) == 0
    with pytest.raises(ValueError):
        depletion_score(1, 2, 4, .2, .3)
    with pytest.raises(ValueError):
        depletion_score(1.5, 2, 1, .2, .3)


def test_uniform_propensity_and_region_counts_exclude_gap_and_outer_flanks():
    plus = np.zeros((1, 201)); minus = plus.copy(); logits = plus.copy()
    plus[0, 100] = 1; plus[0, 108] = 9; plus[0, 115] = 2; minus[0, 0] = 7
    c = region_components(plus, minus, logits, logits, (7, 12, 40))
    assert c['n_plus'][0] == 3 and c['n_minus'][0] == 0 and c['central'][0] == 1
    np.testing.assert_allclose(c['p_plus'], 15/73)
    score, _ = profile_scores(plus, minus, logits, logits, (7, 12, 40))
    np.testing.assert_allclose(score[0], depletion_score(3, 0, 1, 15/73, 15/73), atol=1e-10)
    logits[0, 100] = np.nan
    score, c = profile_scores(plus, minus, logits, logits, (7, 12, 40))
    assert not c['valid'][0] and np.isnan(score[0])
