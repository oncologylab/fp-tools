"""Research-only genomic Gamma background with a prespecified flank radius.

This changes only the background radius. The half-count prior and fixed
deviance dispersion are inherited development choices, not calibrated
biological uncertainty. No labels are used by this numerical primitive.
"""
import numpy as np
from scipy.ndimage import correlate1d
from fp_tools.tools.parametric_bias import calibrated_residuals
from global_parametric_scoring import conventional_weights

RADII=(100,200,500)


def posterior_mean(counts,log_bias,radius):
    counts,log_bias=np.asarray(counts,float),np.asarray(log_bias,float)
    if (radius not in RADII or counts.ndim!=1 or counts.shape!=log_bias.shape
            or not counts.size or not np.isfinite(counts).all() or (counts<0).any()
            or (counts!=np.floor(counts)).any() or np.isinf(log_bias).any()):
        raise ValueError('aligned integer counts, finite-or-missing bias and frozen radius required')
    mean=np.full(counts.shape,np.nan);background=np.full_like(mean,np.nan)
    valid=np.zeros(counts.shape,bool);finite=np.isfinite(log_bias)
    if len(counts)<2*radius+1 or not finite.any(): return mean,background,valid
    bias=np.zeros(counts.shape);bias[finite]=np.exp(log_bias[finite]-log_bias[finite].max())
    if (bias[finite]==0).any(): raise ValueError('bias underflow')
    centers=np.arange(radius,len(counts)-radius)
    def cumulative(x): return np.r_[0,np.cumsum(x)]
    def flank_sum(x):
        p=cumulative(x)
        return p[centers-49]-p[centers-radius]+p[centers+radius+1]-p[centers+50]
    missing=cumulative((~finite).astype(int))
    available=(missing[centers+radius+1]-missing[centers-radius])==0
    n=flank_sum(counts);exposure=flank_sum(bias)
    if (exposure[available]<=0).any() or (n[available]<0).any():
        raise ValueError('numerically invalid background sum')
    indices=centers[available]
    predicted=(n[available]+.5)*bias[indices]/exposure[available]
    if not np.isfinite(predicted).all() or (predicted<=0).any():
        raise ValueError('invalid posterior mean')
    mean[indices],background[indices],valid[indices]=predicted,n[available],True
    return mean,background,valid


def score_tracks(plus,minus,plus_bias,minus_bias,radius,dispersion):
    plus,minus=np.asarray(plus,float),np.asarray(minus,float)
    if plus.shape!=minus.shape: raise ValueError('strand shapes differ')
    first,_,a=posterior_mean(plus,plus_bias,radius)
    second,_,b=posterior_mean(minus,minus_bias,radius)
    observed=plus+minus;expected=first+second;valid=a&b
    residual=np.full(observed.shape,np.nan)
    residual[valid]=calibrated_residuals(observed[valid],expected[valid],'deviance',dispersion=dispersion)
    weights=conventional_weights();required=(weights!=0).astype(float)
    supported=correlate1d(valid.astype(float),required,mode='constant',cval=0)==required.sum()
    score=correlate1d(np.where(valid,residual,0),weights,mode='constant',cval=0)
    informative=correlate1d(observed,np.ones(201),mode='constant',cval=0)>0
    score[~informative]=0;score[~supported]=np.nan
    return dict(observed=observed,expected=expected,corrected=residual,footprints=score,
                expectation_valid=valid,footprint_valid=supported,informative=informative)
