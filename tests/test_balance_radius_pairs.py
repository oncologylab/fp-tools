from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from balance_radius_training_pairs import select_pairs


def test_already_balanced_pairs_are_all_preserved():
    p=np.array([[1,2,3],[2,3,4],[3,4,5]],float)
    ix,smd,distance=select_pairs(p,p.copy(),np.array([3,1,2]))
    np.testing.assert_array_equal(ix,[1,2,0]);np.testing.assert_array_equal(smd,[0,0,0])


def test_outlier_pair_is_removed_without_splitting_pairs_or_using_scores():
    p=np.array([[1,1,1],[2,2,2],[3,3,30]],float)
    n=np.array([[1,1,1],[2,2,2],[3,3,3]],float)
    ix,smd,_=select_pairs(p,n,np.arange(3))
    np.testing.assert_array_equal(ix,[0,1]);np.testing.assert_array_equal(smd,[0,0,0])


def test_no_balanced_prefix_and_invalid_pairs_fail_explicitly():
    ix,smd,_=select_pairs(np.ones((3,3)),np.zeros((3,3)),np.arange(3))
    assert len(ix)==0 and np.isnan(smd).all()
    with pytest.raises(ValueError): select_pairs(np.ones((3,3)),np.zeros((3,3)),np.zeros(3))
    with pytest.raises(ValueError): select_pairs(np.ones((3,2)),np.zeros((3,2)),np.arange(3))
