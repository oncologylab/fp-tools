"""Research-only uncertain local background; no calibration claim is implied.

For background count N and bias exposure B, the scale-invariant Poisson
Jeffreys prior gives lambda | N ~ Gamma(N+1/2, rate=B). At center exposure b,
the count prediction is negative binomial, with mean (N+1/2)b/B and variance
mean + (N+1/2)(b/B)^2. This is a model-based prediction, not a binding
probability or a calibrated simultaneous functional band.
"""

import numpy as np
from scipy.ndimage import correlate1d
from scipy.stats import nbinom

from global_parametric_scoring import conventional_weights

PRIOR_SHAPE = .5


def posterior_track(observed, log_bias):
    observed, log_bias = np.asarray(observed, float), np.asarray(log_bias, float)
    if observed.ndim != 1 or observed.shape != log_bias.shape or not observed.size:
        raise ValueError("aligned nonempty one-dimensional tracks required")
    if (not np.isfinite(observed).all() or np.any(observed < 0) or np.any(observed != np.floor(observed))
            or np.isinf(log_bias).any()):
        raise ValueError("integer counts and finite-or-NaN log bias required")
    finite = np.isfinite(log_bias)
    result = {key: np.full(observed.shape, np.nan) for key in
              ("mean", "mean_variance", "predictive_variance", "shape", "exposure_ratio", "background_count", "prior_fraction")}
    result["valid"] = np.zeros(observed.shape, bool)
    if not finite.any():
        return result
    bias = np.zeros(observed.shape)
    bias[finite] = np.exp(log_bias[finite] - log_bias[finite].max())
    if np.any(bias[finite] == 0):
        raise ValueError("bias underflow; do not silently clip")
    kernel = (np.abs(np.arange(-100, 101)) >= 50).astype(float)
    counts = correlate1d(observed, kernel, mode="constant", cval=0.)
    exposure = correlate1d(bias, kernel, mode="constant", cval=0.)
    support = correlate1d(finite.astype(float), np.ones(201), mode="constant", cval=0.)
    valid = (support == 201) & (exposure > 0)
    ratio, shape = bias[valid] / exposure[valid], counts[valid] + PRIOR_SHAPE
    mean, mean_variance = shape * ratio, shape * ratio**2
    if np.any(ratio <= 0) or np.any(mean <= 0):
        raise ValueError("posterior expectation underflows; do not silently replace it")
    for key, values in (("mean", mean), ("mean_variance", mean_variance),
                        ("predictive_variance", mean + mean_variance), ("shape", shape),
                        ("exposure_ratio", ratio), ("background_count", counts[valid]),
                        ("prior_fraction", PRIOR_SHAPE / shape)):
        if not np.isfinite(values).all():
            raise ValueError("nonfinite posterior calculation")
        result[key][valid] = values
    result["valid"] = valid
    return result


def predictive_interval(shape, exposure_ratio, level=.95):
    shape, ratio = np.broadcast_arrays(np.asarray(shape, float), np.asarray(exposure_ratio, float))
    if (not 0 < level < 1 or not np.isfinite(shape).all() or not np.isfinite(ratio).all()
            or np.any(shape <= 0) or np.any(ratio <= 0)):
        raise ValueError("positive finite posterior parameters and interval level required")
    tail = (1-level)/2
    probability = 1 / (1 + ratio)
    return nbinom.ppf(tail, shape, probability), nbinom.ppf(1-tail, shape, probability)


def score_tracks(plus, minus, plus_log_bias, minus_log_bias):
    """Moment-standardized combined counts and conventional geometry.

    Strand variances add under the working independence assumption. Their
    sum is NOT asserted to be negative binomial when exposure ratios differ.
    Zero observed 201-position neighborhoods are retained as zero-score ties.
    """
    plus, minus = np.asarray(plus, float), np.asarray(minus, float)
    if plus.shape != minus.shape:
        raise ValueError("strand shapes differ")
    first, second = posterior_track(plus, plus_log_bias), posterior_track(minus, minus_log_bias)
    observed = plus + minus
    mean = first["mean"] + second["mean"]
    variance = first["predictive_variance"] + second["predictive_variance"]
    valid = first["valid"] & second["valid"]
    residual = np.full(observed.shape, np.nan)
    residual[valid] = (observed[valid]-mean[valid]) / np.sqrt(variance[valid])
    weights = conventional_weights()
    required = (weights != 0).astype(float)
    footprint_valid = correlate1d(valid.astype(float), required, mode="constant", cval=0.) == required.sum()
    footprints = correlate1d(np.where(valid, residual, 0.), weights, mode="constant", cval=0.)
    informative = correlate1d(observed, np.ones(201), mode="constant", cval=0.) > 0
    footprints[~informative] = 0.
    footprints[~footprint_valid] = np.nan
    return dict(observed=observed, expected=mean, predictive_variance=variance,
                corrected=residual, expectation_valid=valid, footprints=footprints,
                footprint_valid=footprint_valid, informative=informative)
