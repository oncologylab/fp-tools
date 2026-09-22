from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from evaluate_radius_full_pool import evaluate_panels,ARMS


def test_independent_normalization_agrees_for_variable_and_flat_profiles():
    from review_radius_full_pool import normalized
    from fixed_shape_normalization import normalize_fixed
    rng=np.random.default_rng(2026);values=np.vstack([rng.normal(size=(10,201)),np.ones((2,201)),np.zeros((2,201))])
    np.testing.assert_allclose(normalized(values),normalize_fixed(values,np.arange(-100,101),1.),atol=1e-12)


def test_full_denominator_empty_task_zero_ties_and_recomputed_shapes():
    rows=[];diag=[]
    for i in range(21):
        diag.append(dict(cell='HepG2',tf=f'TF{i}',matched_pairs=2 if i<20 else 0,status='underpowered'))
        if i==20: continue
        for j,label in enumerate([0,0,1,1]):
            rows.append(dict(cell='HepG2',tf=f'TF{i}',TFBS_chr='chr1',TFBS_start=i*100+j*10,
                TFBS_end=i*100+j*10+5,TFBS_strand='+',chip_label=label))
    sites=pd.DataFrame(rows);panels={arm:np.zeros((80,201)) for arm in ARMS}
    metrics,curves,scores=evaluate_panels(sites,panels,pd.DataFrame(diag))
    assert len(metrics)==147 and metrics.auroc.notna().sum()==140
    assert metrics.auroc.dropna().eq(.5).all() and not metrics.canonical_protection.dropna().any()
    assert len(curves)==20*7*201 and scores.score.eq(0).all()
    panels['radius_500'][0,0]=np.nan
    with pytest.raises(ValueError,match='finite'): evaluate_panels(sites,panels,pd.DataFrame(diag))
