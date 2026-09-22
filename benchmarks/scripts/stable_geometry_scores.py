"""Fixed numerical precision for research geometry ranks, not probabilities."""
import numpy as np
from scipy.stats import rankdata

DECIMALS = 10
SHAPE_MARGIN = 1e-10


def stable_scores(values, weights, zero_read):
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    zero = np.asarray(zero_read, bool)
    if (values.ndim != 2 or weights.shape != (values.shape[1],)
            or zero.shape != (len(values),) or not np.isfinite(values).all()
            or not np.isfinite(weights).all()):
        raise ValueError('finite aligned profiles, weights and zero-read mask required')
    scores = np.round(values @ weights, decimals=DECIMALS)
    scores[zero | (scores == 0)] = 0.
    return scores


def rank_auc(labels, scores):
    y, score = np.asarray(labels), np.asarray(scores, float)
    if y.shape != score.shape or y.ndim != 1 or not np.isin(y, [0, 1]).all() or not np.isfinite(score).all():
        raise ValueError('finite aligned scores and binary labels required')
    n = int((y == 1).sum())
    if not 0 < n < len(y):
        raise ValueError('both classes required')
    # Half-integer ranks give an exact integer numerator, so mathematically
    # tied candidate AUROCs cannot be separated by integration roundoff.
    numerator = int((2 * rankdata(score)[y == 1]).sum()) - n * (n + 1)
    return numerator / (2 * n * (len(y) - n))
