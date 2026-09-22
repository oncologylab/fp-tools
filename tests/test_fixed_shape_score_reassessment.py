from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks/scripts"))
import reassess_fixed_shape_scores as m


def fixture():
    sites=pd.DataFrame(dict(cell=["K562"]*4,tf=["CTCF"]*4,TFBS_chr=["chr16"]*4,
        TFBS_start=[100,200,300,400],TFBS_end=[110,210,310,410],TFBS_strand=["+"]*4,chip_label=[0,0,1,1]))
    frames=[]
    for d in ("conventional","trained_geometry","hierarchical_spline"):
        f=sites.copy();f["correction"]="raw";f["detector"]=d;f["score"]=[0,1,2,3];frames.append(f)
    metrics=pd.DataFrame([dict(cell="K562",tf="CTCF",correction="raw",detector=d,auroc=1.,auprc=1.,cohort_status="underpowered") for d in ("conventional","trained_geometry","hierarchical_spline")])
    return sites,metrics,pd.concat(frames,ignore_index=True)


def test_only_spline_metrics_and_scores_change():
    s,t,v=fixture();before=t.copy()
    audit=m.replace_scores(t,v,s,dict(conventional=np.arange(4),trained_geometry=np.arange(4),hierarchical_spline=np.array([3,2,1,0])),"K562","CTCF","raw")
    assert len(t)==3 and len(v)==12 and len(audit)==1
    assert t.iloc[2].auroc==0 and t.iloc[:2].equals(before.iloc[:2])
    assert t.cohort_status.eq("underpowered").all()


def test_changed_site_order_is_rejected():
    s,t,v=fixture()
    with pytest.raises(AssertionError):
        m.replace_scores(t,v,s.iloc[::-1],{"hierarchical_spline":np.arange(4)},"K562","CTCF","raw")


def test_geometry_change_or_new_detector_is_rejected():
    s,t,v=fixture()
    with pytest.raises(AssertionError):m.replace_scores(t,v,s,{"trained_geometry":np.zeros(4)},"K562","CTCF","raw")
    with pytest.raises(ValueError):m.replace_scores(t,v,s,{"invented_detector":np.zeros(4)},"K562","CTCF","raw")
