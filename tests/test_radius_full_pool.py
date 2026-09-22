from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from prepare_radius_full_pool import match


def test_three_covariate_matching_recovers_exact_pairs_without_subsampling():
    x=np.array([[0,0,0],[1,1,1],[2,2,2],[0,0,0],[1,1,1],[2,2,2],[10,10,10]],float)
    p,n,d,smd,total=match(x,np.array([1,1,1,0,0,0,1]))
    assert len(p)==total==3
    np.testing.assert_array_equal(x[p],x[n]);np.testing.assert_array_equal(smd,[0,0,0])
    assert (d==0).all() and len(set(p)&set(n))==0


def test_unmatchable_class_and_invalid_features_are_explicit():
    x=np.r_[np.zeros((3,3)),np.ones((3,3))*10]
    p,n,d,smd,total=match(x,np.array([1,1,1,0,0,0]))
    assert len(p)==total==0 and np.isnan(smd).all()
    assert len(match(x,np.ones(6))[0])==0
    with pytest.raises(ValueError): match(x[:,:2],np.ones(6))
