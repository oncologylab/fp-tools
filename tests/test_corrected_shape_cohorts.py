import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import prepare_corrected_shape_cohorts as m


def fixture():
    sites = pd.DataFrame(dict(cell=["HepG2"]*8, tf=["A"]*8, TFBS_chr=["chr1"]*4+["chr16"]*4,
        chip_label=[1,1,0,0]*2, motif_score=[1,2,1,2]*2))
    arrays = [np.ones((8,201)) for _ in range(3)]
    protocol = json.loads(m.PROTOCOL.read_text())
    protocol["minimum_sites_per_class"] = 2
    return sites, arrays, protocol


def test_pairs_are_balanced_and_reproducible_without_crossing_split():
    sites, a, p = fixture()
    args = (sites,*a,np.ones(8,dtype=bool),np.ones(8,dtype=bool),p)
    pairs, rows, support = m.match_cohorts(*args)
    pd.testing.assert_frame_equal(pairs, m.match_cohorts(*args)[0])
    assert rows.status.eq("eligible").all() and support.all()
    for _, group in pairs.groupby(["split","pair_id"]):
        assert sorted(group.chip_label)==[0,1]
        assert len(sites.iloc[group.site_index].TFBS_chr.unique())==1


def test_zero_counts_retained_but_missing_expected_support_excluded():
    sites,a,p=fixture()
    a[0][[0,2]]=0
    a[1][[0,2]]=0
    a[2][1]=0
    _,rows,support=m.match_cohorts(sites,*a,np.ones(8,dtype=bool),np.ones(8,dtype=bool),p)
    assert support[[0,2]].all() and not support[1]
    assert len(rows)==2 and rows.iloc[0].status=="underpowered"


def test_constant_but_shifted_features_are_not_balanced():
    values=np.array([[1.,0.],[1.,0.],[2.,0.],[2.,0.]])
    result=m.smd(values,np.array([1,1,0,0]))
    assert np.isneginf(result[0]) and result[1]==0


def test_absent_tuning_class_and_rejected_pairs_remain_in_ledger():
    sites,a,p=fixture()
    sites.loc[sites.TFBS_chr.eq("chr16"),"chip_label"]=1
    sites.loc[sites.chip_label.eq(0),"motif_score"]=100
    pairs,rows,_=m.match_cohorts(sites,*a,np.ones(8,dtype=bool),np.ones(8,dtype=bool),p)
    assert len(rows)==2 and pairs.empty
    assert rows.status.tolist()==["underpowered","missing_class"]


def test_test_chromosomes_rejected():
    sites,a,p=fixture()
    sites.loc[0,"TFBS_chr"]="chr19"
    with pytest.raises(ValueError,match="development chromosomes"):
        m.match_cohorts(sites,*a,np.ones(8,dtype=bool),np.ones(8,dtype=bool),p)
