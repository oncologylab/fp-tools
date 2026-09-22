from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import screen_fixed_shape_operators as screen


def test_calibration_excludes_missing_expectations_and_zero_counts():
    arrays = dict(valid=np.array([True, True, True, True, False, True]),
                  observed=np.ones((6, 201)), expected=np.ones((6, 201)))
    dwm = np.ones((6, 201))
    dwm[1] = 0
    arrays["expected"][2] = 0
    arrays["observed"][3] = 0
    np.testing.assert_array_equal(screen.support(arrays, dwm, np.array([1, 1, 1, 1, 1, 0], bool)),
                                  [True, False, False, False, False, False])


def summarize(n=200, candidate=None, comparable=None, eligible=True):
    return screen.rate_summary(np.zeros(n) if candidate is None else candidate, np.zeros(n), np.ones(n),
        np.ones(n, bool), np.ones(n, bool) if comparable is None else comparable, 1., 1., eligible)


def test_zero_calls_require_enough_support_for_wilson_bound():
    assert summarize()[0]["passes_fpr_screen"]
    assert not summarize(20)[0]["passes_fpr_screen"]


def test_missing_baseline_cannot_pass_relative_screen():
    row, _, _ = summarize(comparable=np.zeros(200, bool))
    assert row["false_positive_rate"] == 0
    assert row["baseline_comparable_sites"] == 0
    assert not row["passes_fpr_screen"]


def test_absolute_calls_not_hidden_by_baseline_coverage():
    scores = np.zeros(200)
    scores[:20] = 2
    common = np.ones(200, bool)
    common[:20] = False
    row, calls, _ = summarize(candidate=scores, comparable=common)
    assert calls.sum() == row["false_positive_calls"] == 20
    assert row["paired_candidate_calls"] == 0
    assert not row["passes_fpr_screen"]


def test_zero_observed_sites_never_called_and_do_not_inflate_informative_support():
    totals = np.zeros(200)
    totals[:20] = 1
    row, calls, _ = screen.rate_summary(np.full(200, 2.), np.zeros(200), totals,
        np.ones(200, bool), np.ones(200, bool), 1., 1., True)
    assert row["candidate_valid_sites"] == 200
    assert row["informative_sites"] == calls.sum() == 20
    assert not row["passes_fpr_screen"]


def test_paired_increase_and_ineligible_thresholds_fail_closed():
    values = np.zeros(1000)
    values[:20] = 2
    row, _, _ = summarize(1000, candidate=values)
    assert row["wilson_upper_bound"] < .05
    assert row["paired_increase_over_DWM"] == .02
    assert not row["passes_fpr_screen"]
    unavailable = summarize(eligible=False)[0]
    assert not unavailable["passes_fpr_screen"]
    assert not unavailable["threshold_evaluated"]
    assert np.isnan(unavailable["false_positive_rate"])


def test_nonfinite_scores_are_explicitly_unavailable():
    row, calls, _ = summarize(candidate=np.full(200, np.nan))
    assert row["candidate_valid_sites"] == 0
    assert not calls.any() and not row["passes_fpr_screen"]
