from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from audit_control_track_scope import window_support


def test_absent_track_is_not_recorded_zero_or_positive_expectation():
    v=np.full(201,np.nan)
    assert window_support(v)==dict(finite_bases=0,positive_bases=0,category='absent_track_entries')
    v[0]=0
    assert window_support(v)==dict(finite_bases=1,positive_bases=0,category='recorded_zero_only')
    v[1]=2
    assert window_support(v)==dict(finite_bases=2,positive_bases=1,category='positive_expectation')


def test_invalid_values_and_shapes_reject():
    for v in [np.zeros(200),np.full(201,-1.),np.full(201,np.inf)]:
        with pytest.raises(ValueError):window_support(v)
