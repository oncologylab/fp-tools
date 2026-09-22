import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import balance_corrected_shape_cohorts as m


def test_largest_balanced_prefix_preserves_exact_pairs():
    a=np.array([[0.,0.],[1.,1.],[2.,2.],[100.,100.]])
    b=np.array([[0.,0.],[1.,1.],[2.,2.],[3.,3.]])
    ix,smd=m.balanced_prefix(a,b,[0,0,0,.5],[0,1,2,3],.1)
    assert ix.tolist()==[0,1,2]
    np.testing.assert_array_equal(smd,[0,0])


def test_ties_have_deterministic_pair_id_order():
    a=np.array([[2.,2.],[0.,0.],[1.,1.]])
    ix,_=m.balanced_prefix(a,a,[0,0,0],[2,0,1],.1)
    assert ix.tolist()==[1,2,0]


def test_no_balanced_subset_does_not_relax_gate():
    a=np.ones((4,2));b=np.zeros((4,2))
    assert not len(m.balanced_prefix(a,b,[.1]*4,range(4),.1)[0])
    with pytest.raises(ValueError):m.balanced_prefix(a,b,[.1]*4,range(4),.2)


def fixture():
    sites=pd.DataFrame(dict(cell=["HepG2"]*6,tf=["A"]*6,TFBS_chr=["chr1"]*6,
        chip_label=[1,1,1,0,0,0],motif_score=[0,1,10,0,1,2]))
    pairs=pd.DataFrame(dict(cell=["HepG2"]*6,tf=["A"]*6,split=["train"]*6,
        pair_id=[0,1,2,0,1,2],site_index=range(6),chip_label=sites.chip_label,
        match_distance=[0,0,.4,0,0,.4],status=["unbalanced"]*6))
    diagnostics=pd.DataFrame([dict(cell="HepG2",tf="A",split="train",matched_pairs=3,status="unbalanced"),
                              dict(cell="HepG2",tf="A",split="tune",matched_pairs=0,status="missing_class")])
    protocol=json.loads(m.parent.base.PROTOCOL.read_text())
    protocol["minimum_sites_per_class"]=2
    return sites,pairs,diagnostics,protocol


def test_balance_retains_empty_tasks_and_power_gate():
    sites,pairs,diagnostics,p=fixture()
    out,rows=m.balance(sites,np.ones(6),pairs,diagnostics,p)
    assert rows.status.tolist()==["eligible","missing_class"]
    assert sorted(out.site_index)==[0,1,3,4]
    assert out.groupby('pair_id').chip_label.sum().tolist()==[1,1]
    p["minimum_sites_per_class"]=200
    _,rows=m.balance(sites,np.ones(6),pairs,diagnostics,p)
    assert rows.iloc[0].status=="underpowered"


def test_broken_pair_and_changed_identity_fail():
    sites,pairs,diagnostics,p=fixture()
    with pytest.raises(ValueError,match="both classes"):
        m.balance(sites,np.ones(6),pairs.iloc[:-1],diagnostics,p)
    sites.loc[0,"TFBS_chr"]="chr19"
    with pytest.raises(ValueError,match="identity or split"):
        m.balance(sites,np.ones(6),pairs,diagnostics,p)
