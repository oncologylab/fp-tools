from pathlib import Path
import sys
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from bootstrap_radius_training import weighted_ranks,weighted_shape


def test_weighted_tie_metrics_match_explicit_block_replication():
    score=np.array([0,0,1,2,2,3.]);labels=np.array([0,1,1,0,1,0]);blocks=np.array([0,0,1,1,2,2])
    draws=np.array([[1,1,1],[2,0,1],[0,3,0]])
    roc,ap=weighted_ranks(score,labels,blocks,draws)
    for i,draw in enumerate(draws):
        ix=np.repeat(np.arange(6),draw[blocks]);np.testing.assert_allclose([roc[i],ap[i]],
            [roc_auc_score(labels[ix],score[ix]),average_precision_score(labels[ix],score[ix])],atol=1e-12)


def test_missing_classes_are_retained_and_shape_uses_same_weights():
    y=np.array([0,0,1,1]);b=np.array([0,0,1,1]);draws=np.array([[2,0],[0,2],[1,1]])
    roc,ap=weighted_ranks(np.arange(4),y,b,draws)
    assert np.isnan(roc[:2]).all() and np.isnan(ap[:2]).all() and roc[2]==ap[2]==1
    values=np.arange(12).reshape(4,3);result=weighted_shape(values,y,b,draws)
    assert np.isnan(result[:2]).all()
    np.testing.assert_allclose(result[2],values[2:].mean(0)-values[:2].mean(0))
