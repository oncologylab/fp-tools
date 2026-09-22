from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from global_parametric_expectation import expected_track


def test_matches_independent_scalar_annular_reference():
    rng = np.random.default_rng(2026)
    observed = rng.poisson(3, 501)
    log_bias = rng.normal(0, .3, 501)
    expected, valid = expected_track(observed, log_bias)
    assert valid.sum() == 301
    for i in np.flatnonzero(valid):
        ix = np.r_[np.arange(i-100, i-49), np.arange(i+50, i+101)]
        reference = np.exp(log_bias[i]) * observed[ix].sum() / np.exp(log_bias[ix]).sum()
        assert expected[i] == pytest.approx(reference, rel=1e-12)


def test_halo_chunking_and_overlapping_windows_give_identical_values():
    rng = np.random.default_rng(7)
    counts = rng.poisson(2, 1000)
    bias = rng.normal(0, 1, 1000)
    full, mask = expected_track(counts, bias)
    for start, end in ((100, 400), (300, 600), (600, 900)):
        local, valid = expected_track(counts[start-100:end+100], bias[start-100:end+100])
        np.testing.assert_array_equal(valid[100:-100], mask[start:end])
        np.testing.assert_allclose(local[100:-100], full[start:end], rtol=1e-12, atol=1e-12)


def test_reverse_genomic_direction_and_log_offset_are_invariant():
    rng = np.random.default_rng(9)
    counts = rng.poisson(2, 501)
    bias = rng.normal(size=501)
    full, valid = expected_track(counts, bias)
    reverse, reverse_valid = expected_track(counts[::-1], bias[::-1])
    shifted, _ = expected_track(counts, bias+500)
    np.testing.assert_array_equal(valid, reverse_valid[::-1])
    np.testing.assert_allclose(full, reverse[::-1], rtol=1e-12, equal_nan=True)
    np.testing.assert_allclose(full, shifted, rtol=1e-12, equal_nan=True)


def test_center_counts_do_not_set_own_background_and_missing_support_is_explicit():
    counts = np.ones(401)
    bias = np.zeros(401)
    before, _ = expected_track(counts, bias)
    counts[200] = 1000
    after, _ = expected_track(counts, bias)
    assert before[200] == after[200] == 1
    bias[150] = np.nan
    missing, valid = expected_track(counts, bias)
    assert not valid[200] and np.isnan(missing[200])
    empty, valid = expected_track(np.zeros(401), np.zeros(401))
    assert not valid.any() and np.isnan(empty).all()


def test_invalid_arrays_and_underflow_fail_closed():
    for counts, bias in (([1, -1], [0, 0]), ([1], [np.inf]), ([], []), ([1], [0, 1])):
        with pytest.raises(ValueError):
            expected_track(counts, bias)
    with pytest.raises(ValueError, match='underflows'):
        expected_track(np.ones(401), np.r_[0., np.full(400, -1000.)])
    with pytest.raises(ValueError):
        expected_track(np.ones(401), np.zeros(401), inner=0)
