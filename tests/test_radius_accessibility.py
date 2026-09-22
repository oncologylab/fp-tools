from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from audit_radius_accessibility import count_features,association


def test_flank_features_have_exact_boundaries_and_no_central_overlap():
    plus=np.ones((2,1201));minus=np.zeros_like(plus)
    values=count_features(plus,minus)
    for key,count in [('central_201',201),('outer_50_100',102),('outer_50_200',302),('outer_50_500',902),('outer_101_500',800)]:
        np.testing.assert_array_equal(values[key],[count,count])
    plus[:,500:701]=5
    np.testing.assert_array_equal(count_features(plus,minus)['outer_101_500'],[800,800])
    with pytest.raises(ValueError): count_features(plus[:,:100],minus[:,:100])


def test_associations_report_both_directions_without_selecting_one():
    result=association(np.array([0,1,2,3]),np.array([0,0,1,1]))
    assert result['auroc_increasing']==1 and result['auroc_decreasing']==0 and result['smd']>0
    tied=association(np.ones(4),np.array([0,0,1,1]))
    assert tied['smd']==0 and tied['auroc_increasing']==tied['auroc_decreasing']==.5
    with pytest.raises(ValueError): association(np.ones(4),np.ones(4))
