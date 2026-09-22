"""Training-frozen normalization for research spline model inference.

The legacy normalizer uses the current batch's median outer-flank scale for
profiles with no outer-flank variation. Freeze that fallback on training rows
so prediction is independent of batching. Existing model artifacts are not
modified by this module.
"""

import numpy as np


def detrend(profiles, positions):
    values, x = np.asarray(profiles, dtype=float), np.asarray(positions, dtype=float)
    if (values.ndim != 2 or x.shape != (values.shape[1],) or not np.isfinite(values).all()
            or not np.isfinite(x).all() or len(x) < 5):
        raise ValueError("normalization requires finite profiles and matching positions")
    outer = np.abs(x) >= min(60.0, .6 * float(np.max(np.abs(x))))
    if outer.sum() < 4:
        raise ValueError("insufficient outer-flank positions")
    design = np.column_stack([np.ones(outer.sum()), x[outer]])
    coefficients = values[:, outer] @ np.linalg.pinv(design).T
    baseline = coefficients[:, [0]] + coefficients[:, [1]] * x[None, :]
    residual = values - baseline
    scale = np.sqrt(np.mean(residual[:, outer] ** 2, axis=1))
    return residual, scale


def training_fallback(profiles, positions):
    if not len(profiles):
        raise ValueError("training profiles are required to freeze normalization")
    _, scale = detrend(profiles, positions)
    positive = scale[scale > 1e-8]
    return float(np.median(positive)) if len(positive) else 1.0


def normalize_fixed(profiles, positions, fallback):
    if not np.isfinite(fallback) or fallback <= 0:
        raise ValueError("training-frozen fallback must be finite and positive")
    residual, scale = detrend(profiles, positions)
    scale = np.where(scale > 1e-8, scale, fallback)
    return np.clip(residual / scale[:, None], -12.0, 12.0)
