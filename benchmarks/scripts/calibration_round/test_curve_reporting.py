import numpy as np
import pandas as pd
from curve_reporting import fold_curve, curves, realised_power


def test_ties_not_precision_envelope():
    v = fold_curve([1,0,1,0], [1,1,0,0], np.array([0,.5,1]))
    assert v['auroc'] == .5
    assert v['ap'] == .5
    np.testing.assert_allclose(v['pr'], [1,.5,.5])
    w = fold_curve([0,1,1,0], [4,3,2,1], np.array([.5,1]))
    np.testing.assert_allclose(w['pr'], [.5,2/3])


def test_no_silent_nan_deletion():
    assert fold_curve([0,1], [0,np.nan], [0,1]) is None


def test_equal_fold_then_library_and_missing_task():
    rows=[]
    for lib, chrom, repeats, good in [('a','chr1',1,True),('a','chr2',10,False),('b','chr1',1,True)]:
        for i in range(repeats):
            for label in [0,1]:
                rows.append(dict(cell='K562',tf='CTCF',depth='full',library=lib,TFBS_chr=chrom,
                    site_id=f'{chrom}_{i}_{label}',chip_label=label,score=label if good else 1-label,arm='ref'))
    s=pd.DataFrame(rows); c=s[['site_id','library','score','depth']]
    ledger=pd.DataFrame([['K562','CTCF'],['HepG2','HMGA1']],columns=['cell','tf'])
    _,values,_,metrics=curves(s,c,ledger,['ref','covariates_only'],'ref')
    assert values['K562','CTCF','ref']['auroc'] == .75
    assert values['HepG2','HMGA1','ref'] is None
    assert len(metrics[metrics.status.eq('UNAVAILABLE')]) == 2


def test_power_uses_realised_effect_not_weight():
    p=pd.DataFrame(dict(block=[0,1],site_id=['x','y'],r=[.5,.5]))
    e=p.assign(realised_marginal_center_ratio=[.7,np.nan])
    result=realised_power(p,e)
    assert np.isclose(result.realised_depletion.iloc[0],.3)
    assert result.effect_status.iloc[1]=='UNDEFINED_ZERO_SOURCE_CENTER'
