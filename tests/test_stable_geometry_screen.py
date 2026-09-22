from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from screen_stable_geometry_data_root import available,baseline_scores


def test_missing_wide_controls_and_models_remain_unavailable():
    assert not available({'status':'selected','arm':'radius_500'})
    assert not available({'status':'unavailable','arm':'DWM_expected'})
    assert available({'status':'underpowered_fallback','arm':'DWM_expected'})


def test_baseline_is_unnormalized_conventional_and_zero_read_is_tied():
    p=np.zeros((3,201));p[:,93:108]=-2;p[:,60:89]=1;p[:,112:141]=1
    result=baseline_scores({'DWM_expected':p},np.array([1,0,2]))
    np.testing.assert_array_equal(result,[3,0,3])


def test_relative_data_references_use_frozen_original_root(monkeypatch):
    import screen_stable_geometry_data_root as screen
    calls=[]
    def load(config,root,cell,control):
        calls.append((root,cell,control));return 'profiles'
    monkeypatch.setattr(screen.prior,'load_cell',load)
    config={'protocol':screen.PROTOCOL,'prior':{}}
    assert screen.load_cell(config,'HepG2',False)=='profiles'
    assert screen.load_cell(config,'K562',True)=='profiles'
    assert calls==[(screen.DATA_ROOT,'HepG2',False),(screen.DATA_ROOT,'K562',True)]
    assert screen.DATA_ROOT!=screen.ROOT
