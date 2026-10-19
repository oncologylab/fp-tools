"""Focused tests for the calibration-round model, scorer wrapper and call rule."""
from itertools import product
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.special import expit, logit
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import (apply_recalibration, binomial_loglik, call_status, count_tail,  # noqa: E402
                  fit_recalibration)
from protocol import THRESHOLD  # noqa: E402


def simulate(a, b, n_sites=40000, seed=7):
    rng = np.random.default_rng(seed)
    p = rng.uniform(.01, .5, n_sites)
    total = rng.integers(0, 12, n_sites)
    q = expit(a + b * logit(p))
    return rng.binomial(total, q), total, p


def test_recovers_intercept_and_slope():
    central, total, p = simulate(.3, .7)
    fit = fit_recalibration(central, total, p)
    assert abs(fit['a'] - .3) < .06 and abs(fit['b'] - .7) < .04
    assert fit['rows'] == int((total > 0).sum())


def test_identity_null_is_near_zero_one():
    central, total, p = simulate(0., 1.)
    fit = fit_recalibration(central, total, p)
    assert abs(fit['a']) < .05 and abs(fit['b'] - 1) < .03


def test_intercept_only_closed_form_for_constant_propensity():
    central = np.array([1, 0, 3, 2]); total = np.array([5, 4, 6, 5]); p = np.full(4, .2)
    fit = fit_recalibration(central, total, p, slope=False)
    assert fit['b'] == 1.
    assert fit['a'] == pytest.approx(logit(6 / 20) - logit(.2), abs=1e-10)


def test_fit_is_likelihood_maximum():
    central, total, p = simulate(-.2, 1.2, n_sites=3000, seed=3)
    fit = fit_recalibration(central, total, p)
    best = binomial_loglik(central, total, apply_recalibration(p, fit['a'], fit['b']))
    for da, db in product([-.01, .01], [-.01, .01]):
        assert binomial_loglik(central, total, apply_recalibration(p, fit['a'] + da, fit['b'] + db)) < best


def test_zero_information_is_rejected():
    with pytest.raises(ValueError):
        fit_recalibration([0, 0], [0, 0], [.1, .2])
    with pytest.raises(ValueError):
        fit_recalibration([0, 0], [3, 2], [.1, .2])
    with pytest.raises(ValueError):
        fit_recalibration([4, 0], [3, 2], [.1, .2])


def brute_tail(n_plus, n_minus, central, p_plus, p_minus):
    total = 0.
    for i in range(n_plus + 1):
        for j in range(n_minus + 1):
            if i + j <= central:
                total += binom.pmf(i, n_plus, p_plus) * binom.pmf(j, n_minus, p_minus)
    return min(300., max(0., -np.log10(total)))


@pytest.mark.parametrize('n_plus,n_minus,p_plus,p_minus', [(0, 3, .2, .1), (4, 2, .03, .05), (7, 9, .25, .15)])
def test_count_tail_matches_brute_force_and_ceiling(n_plus, n_minus, p_plus, p_minus):
    for central in range(n_plus + n_minus + 1):
        score, ceiling = count_tail(n_plus, n_minus, central, p_plus, p_minus)
        assert score == pytest.approx(round(brute_tail(n_plus, n_minus, central, p_plus, p_minus), 10), abs=1e-9)
        assert ceiling >= score
        assert ceiling == pytest.approx(-np.log10((1 - p_plus) ** n_plus * (1 - p_minus) ** n_minus), abs=1e-9)


def test_call_status_separates_abstention_from_calls():
    assert call_status(0, 0., 0.) == 'NO_INFORMATION'
    assert call_status(3, .2, .5) == 'ABSTAIN_UNREACHABLE'
    assert call_status(3, THRESHOLD, THRESHOLD) == 'ABSTAIN_UNREACHABLE'
    assert call_status(9, THRESHOLD, 2.) == 'NO_CALL'
    assert call_status(9, THRESHOLD + 1e-9, 2.) == 'CALL'
    assert call_status(9, np.nan, np.nan) == 'UNAVAILABLE'


def test_recalibration_is_monotone_for_positive_slope():
    p = np.linspace(.01, .9, 50)
    q = apply_recalibration(p, .4, .6)
    assert np.all(np.diff(q) > 0)
    assert np.allclose(apply_recalibration(p, 0., 1.), p)


# ---------------------------------------------------------------- pipeline helpers
def synthetic_sites(n=600, seed=11):
    import pandas as pd
    from protocol import CHROMOSOMES
    rng = np.random.default_rng(seed)
    frame = pd.DataFrame(dict(site_id=[f's{i}' for i in range(n)], cell='HepG2', tf='CTCF',
                              TFBS_chr=rng.choice(CHROMOSOMES, n), chip_label=rng.integers(0, 2, n),
                              physical_id=[f'p{i}' for i in range(n)], deduplicated_fit_representative=True,
                              eligible=True, LOG21_plus_valid=True, LOG21_minus_valid=True,
                              DWM_plus_valid=True, DWM_minus_valid=True))
    for strand in ['plus', 'minus']:
        frame[f'LOG21_{strand}_p'] = rng.uniform(.02, .4, n)
        frame[f'DWM_{strand}_p'] = rng.uniform(.02, .4, n)
        frame[f'{strand}_total'] = rng.integers(0, 15, n)
        frame[f'{strand}_center'] = rng.binomial(frame[f'{strand}_total'], expit(.2 + .8 * logit(frame[f'LOG21_{strand}_p'])))
    frame['total'] = frame.plus_total + frame.minus_total
    frame['count_stratum'] = '5-9'
    frame['library'] = 'HepG2_rep1'; frame['depth'] = 'full'
    return frame


def test_fold_fits_exclude_held_chromosome_and_score_all_arms():
    from concurrent.futures import ThreadPoolExecutor
    from run import ARMS, coefficient_lookup, fit_folds, score_frame
    frame = synthetic_sites()
    records = fit_folds(frame)
    assert {r['fold'] for r in records} >= {'chr1', 'deployable_all_chr1_15'}
    for r in records:
        if r['fold'].startswith('chr'):
            assert r['fold'] not in r['fitting_chromosomes']
    coefficients = coefficient_lookup(records)
    with ThreadPoolExecutor(2) as pool:
        scored = score_frame(frame, coefficients, pool, ['site_id', 'TFBS_chr', 'total'])
    assert set(scored.arm) == set(ARMS) and len(scored) == len(frame) * len(ARMS)
    reference = scored[scored.arm.eq('LOG21_count_tail')].reset_index(drop=True)
    for i in [0, 5, 77]:
        r = frame.iloc[i]
        expected = count_tail(r.plus_total, r.minus_total, r.plus_center + r.minus_center, r.LOG21_plus_p, r.LOG21_minus_p)
        assert reference.score[i] == expected[0] and reference.ceiling[i] == expected[1]
    primary = scored[scored.arm.eq('LOG21NC_count_tail')].reset_index(drop=True)
    a, b = coefficients[('LOG21NC_count_tail', frame.TFBS_chr[3])]
    assert primary.q_plus[3] == pytest.approx(apply_recalibration(frame.LOG21_plus_p[3], a, b))
    assert (scored.status.eq('CALL') == scored.call).all()


def test_stratified_auc_only_compares_within_strata():
    from run import stratified_auc
    labels = np.array([1, 0, 1, 0]); scores = np.array([10., 9., 1., 0.]); strata = np.array(['a', 'a', 'b', 'b'])
    assert stratified_auc(labels, scores, strata) == (1., 2.)
    scores = np.array([10., 9., 0., 1.])
    assert stratified_auc(labels, scores, strata) == (.5, 2.)


def test_screen_pairs_keeps_missing_tasks_unavailable():
    import pandas as pd
    from protocol import CHROMOSOMES
    from run import screen_pairs
    ledger = pd.DataFrame(dict(cell=['HepG2', 'HepG2'], tf=['CTCF', 'HMGA1']))
    rows = []
    for rep in (1, 2, 3):
        for fold in CHROMOSOMES:
            for arm, auc in [('x', .6), ('y', .5)]:
                rows.append(dict(depth='full', cell='HepG2', tf='CTCF', library=f'HepG2_rep{rep}', fold=fold, arm=arm,
                                 auroc=auc, auprc=auc))
    detail, summary = screen_pairs(pd.DataFrame(rows), ledger, 'full', 'x', 'y')
    assert summary['available_tasks'] == 1 and summary['original_tasks'] == 21
    assert summary['mean_auroc_gain'] == pytest.approx(.1)
    assert summary['mean_relative_auprc_gain'] == pytest.approx(.2)
    assert [r['status'] for r in detail] == ['DEVELOPMENT_DIAGNOSTIC', 'UNAVAILABLE']
