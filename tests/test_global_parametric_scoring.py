from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import global_parametric_scoring as scorer
from corrected_shape_candidates import geometry_weights


def test_conventional_weights_match_frozen_detector():
    np.testing.assert_array_equal(scorer.conventional_weights(), geometry_weights((7, 12, 40))[60:141])


def test_signed_poisson_deviance_and_scalar_geometry_without_rescaling():
    counts = np.ones(401)
    expected = np.full(401, 2.)
    counts[193:208] = 0
    residual, score, valid = scorer.residual_and_footprint(counts, expected, np.ones(401, bool), dispersion=0)
    assert residual[200] == pytest.approx(-2.)
    assert residual[100] == pytest.approx(-np.sqrt(2*(np.log(.5)-1+2)))
    assert score[200] == pytest.approx(residual[160:241] @ scorer.conventional_weights())
    assert valid[200] and score[200] > 0
    # Per-window total normalization would change this expectation and value.
    assert residual[100] != 0


@pytest.mark.parametrize('dispersion', [0., 1.6144318781609228])
def test_count_to_score_halo_parity_includes_negative_binomial(dispersion):
    rng = np.random.default_rng(2026)
    plus, minus = rng.poisson(2, (2, 1000))
    forward, reverse = rng.normal(0, .3, (2, 1000))
    full = scorer.score_genomic_counts(plus, minus, forward, reverse, dispersion=dispersion)
    for start, end in ((140, 400), (350, 650), (600, 860)):
        left, right = start-scorer.SCORING_HALO, end+scorer.SCORING_HALO
        local = scorer.score_genomic_counts(plus[left:right], minus[left:right], forward[left:right], reverse[left:right], dispersion=dispersion)
        for field in full:
            np.testing.assert_allclose(local[field][140:-140], full[field][start:end], rtol=1e-10, atol=1e-10, equal_nan=True)


def test_unavailable_background_cannot_create_a_finite_footprint():
    result = scorer.score_genomic_counts(np.zeros(401), np.zeros(401), np.zeros(401), np.zeros(401), dispersion=0)
    assert not result['expectation_valid'].any() and not result['footprint_valid'].any()
    assert np.isnan(result['footprints']).all()
    counts, expected, valid = np.ones(401), np.ones(401), np.ones(401, bool)
    expected[200] = np.nan
    _, score, mask = scorer.residual_and_footprint(counts, expected, valid, dispersion=0)
    assert not mask[200] and np.isnan(score[200])


def test_invalid_dispersion_and_misaligned_strands_reject():
    with pytest.raises(ValueError):
        scorer.residual_and_footprint(np.ones(10), np.ones(10), np.ones(10, bool), dispersion=np.nan)
    with pytest.raises(ValueError):
        scorer.score_genomic_counts(np.ones(10), np.ones(11), np.zeros(10), np.zeros(11), dispersion=0)
