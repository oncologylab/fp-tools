"""Research prototype: one strand-specific expectation per genomic position.

No fitting or candidate selection occurs here. Real-data execution must use a
separate frozen configuration. Earlier conditional-profile evidence does not
qualify this numerically different operator.
"""

import numpy as np
from scipy.ndimage import correlate1d


def expected_track(observed, log_bias, *, inner=50, outer=100):
    """Estimate background from symmetric outer flanks, excluding the center.

    E[i] = exp(log_bias[i]) * sum(observed[j]) / sum(exp(log_bias[j])),
    where inner <= abs(j-i) <= outer. Both sums use the same strand.
    A common exponential offset cancels. A full valid 2*outer+1 neighborhood
    is required; no padding is treated as data. A zero background remains
    explicitly unavailable rather than claiming a supported zero expectation.

    Chunked callers must supply an `outer` halo and retain only the interior.
    Inputs are genomic-oriented, not motif-oriented. Invalid expectation is
    NaN, and the returned mask must be propagated into downstream coverage.
    """
    observed, log_bias = np.asarray(observed, float), np.asarray(log_bias, float)
    if observed.ndim != 1 or observed.shape != log_bias.shape or observed.size == 0:
        raise ValueError("aligned nonempty one-dimensional tracks required")
    if (not isinstance(inner, (int, np.integer)) or not isinstance(outer, (int, np.integer))
            or not 1 <= inner <= outer):
        raise ValueError("require integer 1 <= inner <= outer")
    if not np.isfinite(observed).all() or (observed < 0).any() or np.isinf(log_bias).any():
        raise ValueError("finite nonnegative counts and finite-or-NaN log bias required")
    finite = np.isfinite(log_bias)
    expected = np.full(observed.shape, np.nan)
    if not finite.any():
        return expected, np.zeros(observed.shape, bool)
    bias = np.zeros(observed.shape)
    bias[finite] = np.exp(log_bias[finite] - log_bias[finite].max())
    if np.any(bias[finite] == 0):
        raise ValueError("bias dynamic range underflows; do not silently clip")
    offsets = np.arange(-outer, outer + 1)
    kernel = (np.abs(offsets) >= inner).astype(float)
    count_sum = correlate1d(observed, kernel, mode="constant", cval=0.)
    bias_sum = correlate1d(bias, kernel, mode="constant", cval=0.)
    support = correlate1d(finite.astype(float), np.ones(len(offsets)), mode="constant", cval=0.)
    valid = (support == len(offsets)) & (count_sum > 0) & (bias_sum > 0)
    # Form the bounded bias ratio before multiplying by counts.
    expected[valid] = (bias[valid] / bias_sum[valid]) * count_sum[valid]
    if not np.isfinite(expected[valid]).all():
        raise ValueError("nonfinite expected counts")
    return expected, valid
