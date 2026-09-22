"""Research-only genomic residual and conventional footprint track primitives.

No per-window total rescaling or learned detector is applied. These functions
do not fit, calibrate, export a pipeline, or establish biological performance.
"""

from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import correlate1d

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from fp_tools.tools.parametric_bias import calibrated_residuals
from global_parametric_expectation import expected_track

EXPECTATION_HALO = 100
SCORING_HALO = 140


def conventional_weights():
    positions = np.arange(-40, 41)
    center = np.abs(positions) <= 7
    shoulders = np.abs(positions) >= 12
    return shoulders / shoulders.sum() - center / center.sum()


def residual_and_footprint(observed, expected, valid, *, dispersion):
    observed, expected = np.asarray(observed, float), np.asarray(expected, float)
    valid = np.asarray(valid, bool)
    if observed.ndim != 1 or observed.shape != expected.shape or valid.shape != observed.shape or not observed.size:
        raise ValueError("aligned nonempty one-dimensional arrays required")
    if (not np.isfinite(observed).all() or np.any(observed < 0) or np.isinf(expected).any()
            or np.any(expected < 0) or not np.isfinite(dispersion) or dispersion < 0):
        raise ValueError("invalid counts, expectation or dispersion")
    usable = valid & np.isfinite(expected) & (expected > 0)
    residual = np.full(observed.shape, np.nan)
    residual[usable] = calibrated_residuals(observed[usable], expected[usable], "deviance", dispersion=dispersion)
    weights = conventional_weights()
    required = (weights != 0).astype(float)
    support = correlate1d(usable.astype(float), required, mode="constant", cval=0.)
    score_valid = support == required.sum()
    # Zeroes are internal convolution padding only. Every used position must
    # pass coverage before its score can be returned as finite evidence.
    score = correlate1d(np.where(usable, residual, 0.), weights, mode="constant", cval=0.)
    score[~score_valid] = np.nan
    return residual, score, score_valid


def score_genomic_counts(plus, minus, plus_log_bias, minus_log_bias, *, dispersion):
    plus, minus = np.asarray(plus, float), np.asarray(minus, float)
    if plus.shape != minus.shape:
        raise ValueError("strand count shapes differ")
    first, first_valid = expected_track(plus, plus_log_bias)
    second, second_valid = expected_track(minus, minus_log_bias)
    expected = first + second
    valid = first_valid & second_valid
    expected[~valid] = np.nan
    observed = plus + minus
    residual, scores, score_valid = residual_and_footprint(observed, expected, valid, dispersion=dispersion)
    return dict(observed=observed, expected=expected, expectation_valid=valid,
                corrected=residual, footprints=scores, footprint_valid=score_valid)
