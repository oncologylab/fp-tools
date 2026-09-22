from pathlib import Path
import sys
import json
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from freeze_stable_geometry_models import score_model


def test_serialized_operator_preserves_sparse_scores_and_zero_rule():
    model=json.loads(json.dumps(dict(status='selected',geometry=[3,8,60],normalization_fallback=1.,decimals=10)))
    p=np.random.default_rng(2026).poisson(.1,(10,201)).astype(float);zero=np.zeros(10,bool);zero[0]=True
    scores=score_model(model,p,zero)
    assert scores[0]==0
    np.testing.assert_array_equal(scores,[score_model(model,p[i:i+1],zero[i:i+1])[0] for i in range(10)])
    model['decimals']=8
    with pytest.raises(ValueError,match='identity'):score_model(model,p,zero)


def test_unavailable_model_cannot_silently_score():
    with pytest.raises(ValueError,match='unavailable'):score_model({'status':'unavailable'},np.zeros((1,201)),[False])
