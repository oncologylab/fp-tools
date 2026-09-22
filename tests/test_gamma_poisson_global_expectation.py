from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import gamma, nbinom, poisson

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import gamma_poisson_global_expectation as model


def test_posterior_parameters_match_scalar_conjugacy_and_zero_background():
    counts = np.zeros(401)
    bias = np.zeros(401)
    result = model.posterior_track(counts, bias)
    assert result['valid'].sum() == 201
    assert result['shape'][200] == .5 and result['prior_fraction'][200] == 1
    assert result['mean'][200] == pytest.approx(.5/102)
    assert result['mean_variance'][200] == pytest.approx(.5/102**2)
    counts[100] = 3
    other = model.posterior_track(counts, bias)
    assert other['shape'][200] == 3.5 and other['mean'][200] == pytest.approx(3.5/102)


def test_predictive_mass_matches_independent_gamma_poisson_integration():
    for shape, ratio in ((.5, .1), (2.5, .3), (10.5, 1.)):
        for count in (0, 1, 4):
            integrated = quad(lambda mu: poisson.pmf(count, mu)*gamma.pdf(mu, shape, scale=ratio),
                              0, np.inf, epsabs=1e-9)[0]
            assert nbinom.pmf(count, shape, 1/(1+ratio)) == pytest.approx(integrated, rel=1e-7, abs=1e-9)
        low, high = model.predictive_interval(shape, ratio)
        mass = nbinom.cdf(high, shape, 1/(1+ratio))-nbinom.cdf(low-1, shape, 1/(1+ratio))
        assert mass >= .95


def test_chunk_halo_and_log_offset_invariance_for_complete_score():
    rng = np.random.default_rng(2026)
    plus, minus = rng.poisson(.05, (2, 900))
    forward, reverse = rng.normal(0, .5, (2, 900))
    full = model.score_tracks(plus, minus, forward, reverse)
    local = model.score_tracks(plus[60:840], minus[60:840], forward[60:840]+100, reverse[60:840]-100)
    for field in full:
        np.testing.assert_allclose(local[field][140:-140], full[field][200:700], rtol=1e-10, atol=1e-10, equal_nan=True)


def test_zero_reads_are_ties_despite_nonconstant_sequence_propensity():
    bias = np.random.default_rng(2).normal(size=401)
    result = model.score_tracks(np.zeros(401), np.zeros(401), bias, -bias)
    valid = result['footprint_valid']
    assert valid.any() and not result['footprints'][valid].any()
    assert not result['informative'].any()
    assert (result['predictive_variance'][result['expectation_valid']] > 0).all()


def test_missing_sequence_is_not_repaired_by_the_prior():
    counts, bias = np.zeros(401), np.zeros(401)
    bias[200] = np.nan
    result = model.posterior_track(counts, bias)
    assert not result['valid'].any()
    assert np.isnan(result['mean']).all()


def test_invalid_counts_and_posterior_parameters_reject():
    with pytest.raises(ValueError): model.posterior_track([.5], [0])
    with pytest.raises(ValueError): model.posterior_track([-1], [0])
    with pytest.raises(ValueError): model.predictive_interval(0, 1)
    with pytest.raises(ValueError): model.predictive_interval(1, np.nan)
