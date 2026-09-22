from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks/scripts"))
import render_corrected_shape_tuning as m


def tables():
    metrics=pd.DataFrame([dict(cell="K562",tf="TF",correction=c,detector=d,matched_pairs=4,cohort_status="underpowered")
                          for c in m.CORRECTIONS for d in ("conventional","trained_geometry","hierarchical_spline")])
    curves=pd.concat([pd.DataFrame(dict(cell="K562",tf="TF",correction=c,position=np.arange(-100,101),
        positive_mean=np.zeros(201),negative_mean=np.ones(201),positive_minus_negative=-np.ones(201))) for c in m.CORRECTIONS],ignore_index=True)
    return metrics,curves


def test_incomplete_or_inconsistent_curves_rejected():
    metrics,curves=tables()
    with pytest.raises(ValueError):m.validate_tables(metrics,curves.iloc[1:])
    curves.loc[0,"positive_minus_negative"]=0
    with pytest.raises(ValueError):m.validate_tables(metrics,curves)


def test_task_denominator_and_factorial_are_complete():
    metrics,curves=tables()
    with pytest.raises(ValueError):m.validate_tables(metrics.iloc[1:],curves)
    with pytest.raises(ValueError):m.validate_tables(metrics,pd.concat([curves,curves.iloc[:1]]))


def test_underpowered_task_is_rendered_and_no_overwrite(tmp_path):
    metrics,curves=tables()
    out=tmp_path/"plots"
    paths=m.render(metrics,curves,out)
    assert [p.suffix for p in paths]==[".pdf",".png"]
    assert all(p.stat().st_size>1000 for p in paths)
    with pytest.raises(ValueError):m.render(metrics,curves,out)
