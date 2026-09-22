from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from audit_stable_geometry_direction import direction_summary


def test_direction_retains_missing_and_tied_chromosomes():
    result = direction_summary([.1, -.2, 0, np.nan], [True, False, True, True])
    assert result['evaluable'] == 3 and result['positive'] == result['negative'] == result['tied'] == 1
    assert result['fixed_protection'] == 2 and result['positive_and_protective'] == 1
    np.testing.assert_allclose([result['omitted_mean_min'], result['omitted_mean_max']], [-.1, .05])


def test_unavailable_task_is_not_zero_gain():
    result = direction_summary([np.nan] * 15, [False] * 15)
    assert result['evaluable'] == result['positive'] == 0
    assert np.isnan(result['mean']) and np.isnan(result['omitted_mean_min'])
