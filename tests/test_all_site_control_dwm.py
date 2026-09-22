from pathlib import Path
import sys
import pandas as pd
import pytest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from build_all_site_control_dwm import site_regions
from review_all_site_control_dwm import compare_expectations


def test_union_preserves_coordinates_without_duplicate_prediction_regions():
    a=pd.DataFrame(dict(TFBS_chr=['chr2','chr1'],TFBS_start=[10,20],TFBS_end=[25,35]))
    result=site_regions([a,a.iloc[:1]])
    assert result.values.tolist()==[['chr1',20,35],['chr2',10,25]]


def test_control_region_builder_rejects_labels_and_invalid_coordinates():
    a=pd.DataFrame(dict(TFBS_chr=['chr1'],TFBS_start=[10],TFBS_end=[25],chip_label=[1]))
    with pytest.raises(ValueError,match='label-free'):site_regions([a])
    a=a.drop(columns='chip_label');a['TFBS_end']=5
    with pytest.raises(ValueError,match='invalid'):site_regions([a])


def test_overlap_allows_global_scaling_but_rejects_shape_drift_and_support_loss():
    old=np.array([[1,2,3],[0,0,0],[1,2,3],[1,2,3]],float)
    new=np.array([[2,4,6],[1,2,3],[3,2,1],[0,0,0]],float)
    positive,match,error=compare_expectations(old,new)
    np.testing.assert_array_equal(positive,[True,False,True,True])
    np.testing.assert_array_equal(match,[True,False,False,False])
    assert error[0]==0 and error[2]>0
