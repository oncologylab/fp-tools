from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
from audit_wide_background import prefix_diagnostic, summary
from build_wide_training_profiles import background_diagnostic


@pytest.mark.parametrize('radius',[100,200,500])
def test_independent_prefix_and_logspace_calculation(radius):
    rng=np.random.default_rng(2026)
    counts=rng.poisson(.1,(12,1201));bias=rng.normal(size=counts.shape)
    bias[0,100]=np.nan;bias[1,600]=np.nan
    first=prefix_diagnostic(counts,bias,600,radius)
    second=background_diagnostic(counts,bias,600,radius)
    for key in first:
        np.testing.assert_allclose(first[key],second[key],rtol=1e-10,atol=1e-10,equal_nan=True)
    shifted=prefix_diagnostic(counts,bias+1000,600,radius)
    np.testing.assert_allclose(first['disagreement'],shifted['disagreement'],atol=1e-10,equal_nan=True)


def test_common_support_does_not_silently_fill_missing_sequence():
    counts=np.zeros((2,1201));bias=np.zeros_like(counts);bias[0,100]=np.nan
    small=prefix_diagnostic(counts,bias,600,100)
    large=prefix_diagnostic(counts,bias,600,500)
    common=small['available']&large['available']
    assert small['available'].sum()==2 and common.sum()==1
    for diagnostic in (small,large):
        result=summary(diagnostic,common)
        assert result['supported']==result['zero_background']==1
        assert result['median_absolute_log2_rate_disagreement']==0
    with pytest.raises(ValueError): prefix_diagnostic(counts,bias,600,700)
