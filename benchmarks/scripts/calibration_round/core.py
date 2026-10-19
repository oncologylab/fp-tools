"""Label-free null-mean recalibration of LOG21 central shares and shared scoring.

The recalibration changes only the probability that a conditioned cut falls in
the fixed centre (positions -7..7) rather than a shoulder (12..40 either side).
It is fitted on protein-free naked-DNA windows and never sees occupancy labels.
"""
from pathlib import Path
import sys

import numpy as np
from scipy.special import expit, logit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from conditional_protection_score import depletion_score  # noqa: E402
from protocol import CAP, DECIMALS, THRESHOLD  # noqa: E402

P_CLIP = (1e-6, 1 - 1e-6)


def clipped_logit(p):
    return logit(np.clip(np.asarray(p, float), *P_CLIP))


def binomial_loglik(central, total, q):
    central = np.asarray(central, float)
    total = np.asarray(total, float)
    q = np.clip(np.asarray(q, float), 1e-300, 1 - 1e-16)
    return float(np.sum(central * np.log(q) + (total - central) * np.log1p(-q)))


def fit_recalibration(central, total, p, slope=True, tol=1e-12, max_iter=100):
    """Maximum-likelihood (a, b) for central ~ Binomial(total, expit(a + b*logit(p))).

    With slope=False, b is fixed at 1 and only the intercept is estimated.
    Rows with zero total carry no information and are removed.
    """
    central = np.asarray(central, float)
    total = np.asarray(total, float)
    x = clipped_logit(p)
    keep = total > 0
    if (central < 0).any() or (central > total).any() or not np.isfinite(x[keep]).all():
        raise ValueError('invalid recalibration inputs')
    central, total, x = central[keep], total[keep], x[keep]
    if not len(total) or central.sum() in (0, total.sum()):
        raise ValueError('recalibration is not identified')
    design = np.column_stack([np.ones_like(x), x]) if slope else np.ones((len(x), 1))
    offset = np.zeros_like(x) if slope else x
    theta = np.zeros(design.shape[1])
    if slope:
        theta[1] = 1.
    for iteration in range(1, max_iter + 1):
        q = expit(offset + design @ theta)
        gradient = design.T @ (central - total * q)
        hessian = (design * (total * q * (1 - q))[:, None]).T @ design
        step = np.linalg.solve(hessian, gradient)
        theta = theta + step
        if np.max(np.abs(step)) < tol:
            break
    else:
        raise ArithmeticError('recalibration did not converge')
    a, b = (float(theta[0]), float(theta[1])) if slope else (float(theta[0]), 1.)
    q = expit(a + b * x)
    return dict(a=a, b=b, iterations=iteration, rows=int(len(x)),
                central=int(central.sum()), total=int(total.sum()),
                loglik=binomial_loglik(central, total, q))


def apply_recalibration(p, a, b):
    return expit(a + b * clipped_logit(p))


def count_tail(n_plus, n_minus, central, q_plus, q_minus):
    """Frozen count-tail score and attainable ceiling, rounded as historically."""
    if not np.isfinite([q_plus, q_minus]).all():
        return np.nan, np.nan
    score = round(depletion_score(n_plus, n_minus, central, q_plus, q_minus), DECIMALS)
    ceiling = round(depletion_score(n_plus, n_minus, 0, q_plus, q_minus), DECIMALS)
    return min(score, CAP), min(ceiling, CAP)


def call_status(total, score, ceiling):
    """Separate calls from honest abstention; abstention is never a negative call."""
    if not np.isfinite(score):
        return 'UNAVAILABLE'
    if total == 0:
        return 'NO_INFORMATION'
    if ceiling <= THRESHOLD:
        return 'ABSTAIN_UNREACHABLE'
    return 'CALL' if score > THRESHOLD else 'NO_CALL'
