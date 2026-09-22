"""Research central-depletion score under a strand-conditioned binomial null.

The null assumes independent cuts within and between strands and correct fixed
sequence propensities. It does not guarantee calibration for paired fragments,
overdispersed biological data, or an imperfect bias model. Scores are not
binding probabilities. Geometry must be fixed without the evaluated labels.
"""
import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom

POSITIONS = np.arange(-100, 101)
MAX_SCORE = 300.


def depletion_score(n_plus, n_minus, central, p_plus, p_minus):
    """Return capped -log10 P(Bin(n_plus,p_plus)+Bin(n_minus,p_minus)<=central)."""
    counts = np.asarray([n_plus, n_minus, central], float)
    probabilities = np.asarray([p_plus, p_minus], float)
    if (not np.isfinite(counts).all() or (counts < 0).any() or (counts != np.floor(counts)).any()
            or not np.isfinite(probabilities).all() or (probabilities < 0).any() or (probabilities > 1).any()):
        raise ValueError('nonnegative integer counts and finite probabilities in [0,1] required')
    n_plus, n_minus, central = map(int, counts)
    if central > n_plus + n_minus:
        raise ValueError('central counts exceed conditioned counts')
    if central == n_plus + n_minus:
        return 0.
    # Sum over the smaller strand for bounded practical work on sparse windows.
    if n_plus > n_minus:
        n_plus, n_minus, p_plus, p_minus = n_minus, n_plus, p_minus, p_plus
    possible = np.arange(min(n_plus, central) + 1)
    terms = binom.logpmf(possible, n_plus, p_plus) + binom.logcdf(central - possible, n_minus, p_minus)
    log_p = min(0., float(logsumexp(terms)))
    return min(MAX_SCORE, max(0., -log_p / np.log(10.)))


def region_components(plus, minus, plus_log_bias, minus_log_bias, geometry):
    arrays = [np.asarray(a, float) for a in [plus, minus, plus_log_bias, minus_log_bias]]
    plus, minus, plus_log_bias, minus_log_bias = arrays
    if plus.ndim != 2 or plus.shape[1] != 201 or any(a.shape != plus.shape for a in arrays):
        raise ValueError('aligned 201-position strand arrays required')
    if any(not np.isfinite(a).all() or (a < 0).any() or (a != np.floor(a)).any() for a in [plus, minus]):
        raise ValueError('finite nonnegative integer strand counts required')
    center, left, right = geometry
    if not 0 < center < left < right <= 100:
        raise ValueError('invalid geometry')
    inner = abs(POSITIONS) <= center
    region = inner | ((abs(POSITIONS) >= left) & (abs(POSITIONS) <= right))
    probabilities, validities = [], []
    for logits in [plus_log_bias, minus_log_bias]:
        valid = np.isfinite(logits[:, region]).all(1)
        safe = np.where(np.isfinite(logits), logits, -np.inf)
        denominator = logsumexp(safe[:, region], axis=1)
        numerator = logsumexp(safe[:, inner], axis=1)
        probability = np.zeros(len(plus))
        probability[valid] = np.exp(numerator[valid] - denominator[valid])
        probabilities.append(np.clip(probability, 0., 1.))
        validities.append(valid)
    return dict(n_plus=plus[:, region].sum(1).astype(np.int64),
                n_minus=minus[:, region].sum(1).astype(np.int64),
                central=(plus[:, inner].sum(1) + minus[:, inner].sum(1)).astype(np.int64),
                p_plus=probabilities[0], p_minus=probabilities[1],
                valid=validities[0] & validities[1])


def profile_scores(plus, minus, plus_log_bias, minus_log_bias, geometry):
    components = region_components(plus, minus, plus_log_bias, minus_log_bias, geometry)
    scores = np.full(len(plus), np.nan)
    for i in np.flatnonzero(components['valid']):
        scores[i] = depletion_score(*(components[k][i] for k in ['n_plus', 'n_minus', 'central', 'p_plus', 'p_minus']))
    # Explicit numerical rank precision; the unrounded null probabilities are
    # tested separately. This does not claim exactly calibrated rounded tails.
    return np.round(scores, 10), components
