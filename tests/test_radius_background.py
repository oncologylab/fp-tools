from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from radius_background import posterior_mean,score_tracks
from gamma_poisson_global_expectation import posterior_track


@pytest.mark.parametrize('radius',[100,200,500])
def test_radius_matches_independent_scalar_and_halo(radius):
    rng=np.random.default_rng(2026)
    counts=rng.poisson(.2,2401);bias=rng.normal(size=len(counts))
    mean,n,valid=posterior_mean(counts,bias,radius)
    for center in (radius,1200,len(counts)-radius-1):
        indices=np.r_[np.arange(center-radius,center-49),np.arange(center+50,center+radius+1)]
        expected=(counts[indices].sum()+.5)*np.exp(bias[center])/np.exp(bias[indices]).sum()
        assert mean[center]==pytest.approx(expected,rel=1e-11)
        assert n[center]==counts[indices].sum() and valid[center]
        local=posterior_mean(counts[center-radius:center+radius+1],bias[center-radius:center+radius+1],radius)
        assert local[0][radius]==pytest.approx(mean[center],rel=1e-11)
    shifted=posterior_mean(counts,bias+100,radius)
    np.testing.assert_allclose(mean,shifted[0],rtol=1e-11,atol=1e-11,equal_nan=True)
    score=score_tracks(counts,counts,bias,bias,radius,1.6)
    halo=radius+40;center=1200
    chunk=score_tracks(counts[center-halo:center+halo+1],counts[center-halo:center+halo+1],
                       bias[center-halo:center+halo+1],bias[center-halo:center+halo+1],radius,1.6)
    assert score['footprints'][center]==pytest.approx(chunk['footprints'][halo],abs=1e-11)


def test_radius100_matches_previous_mean_and_zero_read_ties():
    rng=np.random.default_rng(9);counts=rng.poisson(.1,1201);bias=rng.normal(size=len(counts))
    new=posterior_mean(counts,bias,100);old=posterior_track(counts,bias)
    np.testing.assert_array_equal(new[2],old['valid'])
    np.testing.assert_allclose(new[0],old['mean'],rtol=1e-11,atol=1e-11,equal_nan=True)
    zero=np.zeros_like(counts)
    scored=score_tracks(zero,zero,bias,bias,500,0)
    assert scored['footprints'][600]==0 and scored['expected'][600]>0
    bias[600]=np.nan
    assert not posterior_mean(counts,bias,500)[2][600]
    with pytest.raises(ValueError): posterior_mean(counts,bias,300)
    with pytest.raises(ValueError): posterior_mean(counts+.1,bias,100)
