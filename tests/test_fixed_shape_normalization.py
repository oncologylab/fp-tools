from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"benchmarks/scripts"))
import corrected_shape_candidates as candidates
from fixed_shape_normalization import normalize_fixed, training_fallback


def profiles():
    rng=np.random.default_rng(2026)
    values=rng.normal(0,.2,(20,201))
    values[-1]=0
    values[-1,100]=-1
    return values


def test_legacy_batch_dependence_is_reproduced():
    x=profiles();p=candidates.POSITIONS
    full=candidates.normalize_functional_profiles(x,p)[-1]
    singleton=candidates.normalize_functional_profiles(x[-1:],p)[0]
    assert not np.allclose(full,singleton)


def test_fixed_fallback_matches_training_and_is_batch_independent():
    x=profiles();p=candidates.POSITIONS
    fallback=training_fallback(x,p)
    full=normalize_fixed(x,p,fallback)
    np.testing.assert_allclose(full,candidates.normalize_functional_profiles(x,p),rtol=1e-12,atol=1e-12)
    singleton=np.concatenate([normalize_fixed(row[None,:],p,fallback) for row in x])
    np.testing.assert_allclose(full,singleton,rtol=1e-12,atol=1e-12)
    unrelated=np.random.default_rng(7).normal(0,10,(10,201))
    together=normalize_fixed(np.concatenate([x,unrelated]),p,fallback)
    np.testing.assert_allclose(full,together[:len(x)],rtol=1e-12,atol=1e-12)


def test_zero_variation_training_has_fixed_unit_fallback():
    x=np.zeros((3,201));p=candidates.POSITIONS
    assert training_fallback(x,p)==1
    assert not normalize_fixed(x,p,1).any()


def test_invalid_data_and_fallback_are_rejected():
    x=profiles();p=candidates.POSITIONS
    with pytest.raises(ValueError):normalize_fixed(x,p,0)
    with pytest.raises(ValueError):training_fallback(x[:0],p)
    x[0,0]=np.nan
    with pytest.raises(ValueError):normalize_fixed(x,p,1)


def test_artifact_copy_preserves_weights_and_source_bytes(tmp_path):
    import freeze_shape_normalization as m
    source=tmp_path/"source.npz";target=tmp_path/"new.npz"
    np.savez(source,geometry=np.array([7,12,40]),discriminant=np.ones(201),midpoint=np.array(.2),score_scale=np.array(2.))
    original=source.read_bytes();x=profiles()
    row=m.copy_model(source,target,x)
    assert source.read_bytes()==original and row['learned_weights_unchanged']
    assert row['legacy_batch_dependent_sites']>0 and row['fixed_batch_invariant']
    with np.load(source) as a,np.load(target) as b:
        for k in a.files:np.testing.assert_array_equal(a[k],b[k])
    full=m.scores(target,x)['hierarchical_spline']
    single=np.array([m.scores(target,r[None,:])['hierarchical_spline'][0] for r in x])
    np.testing.assert_allclose(full,single,rtol=1e-12,atol=1e-12)


def test_fixed_scorer_rejects_legacy_spline_and_model_overwrite(tmp_path):
    import freeze_shape_normalization as m
    source=tmp_path/"legacy.npz"
    np.savez(source,geometry=np.array([7,12,40]),discriminant=np.ones(201),midpoint=np.array(0.),score_scale=np.array(1.))
    with pytest.raises(ValueError,match="training-frozen"):
        m.scores(source,profiles())
    with pytest.raises(ValueError,match="overwrite"):
        m.copy_model(source,source,profiles())
