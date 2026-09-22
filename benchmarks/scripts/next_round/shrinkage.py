"""Research-only signed depletion under a fixed zero-centred normal prior.

The two strand factors are a working likelihood, not independent biological
evidence. Posterior intervals are conditional model intervals, not binding
probabilities or calibrated functional uncertainty.
"""
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import betaln, expit, log_expit, logit
from scipy.stats import norm


def validate(c, n, p, rho):
    c, n, p = [np.asarray(v, float) for v in (c, n, p)]
    if (c.shape != (2,) or n.shape != c.shape or p.shape != c.shape
            or not np.isfinite(np.r_[c, n, p, rho]).all()
            or (c < 0).any() or (n < c).any()
            or (c != np.floor(c)).any() or (n != np.floor(n)).any()
            or (p < 0).any() or (p > 1).any() or not 0 <= rho <= .95):
        raise ValueError('two finite integer strand counts and valid probabilities required')
    if ((p == 0) & (c != 0)).any() or ((p == 1) & (c != n)).any():
        raise ValueError('observed count outside exact probability support')
    return c, n, p


def log_likelihood(theta, c, n, p, rho):
    """Drop theta-constant combinatorial terms; exact p=0/1 contributions are 0."""
    result = np.zeros_like(np.asarray(theta, float))
    for ci, ni, pi in zip(c, n, p):
        if ni == 0 or pi in (0, 1):
            continue
        eta = logit(pi) - theta
        if rho == 0:
            result = result + ci * log_expit(eta) + (ni-ci) * log_expit(-eta)
        else:
            concentration = 1/rho - 1
            a = expit(eta) * concentration
            b = expit(-eta) * concentration
            result = result + betaln(ci+a, ni-ci+b) - betaln(a, b)
    return result


@lru_cache(maxsize=200000)
def posterior(c0, c1, n0, n1, p0, p1, rho=0., prior_sd=1., tolerance=2e-11):
    c, n, p = validate([c0, c1], [n0, n1], [p0, p1], rho)
    if prior_sd <= 0 or not np.isfinite(prior_sd):
        raise ValueError('positive prior standard deviation required')
    informative = bool(np.any((n > 0) & (p > 0) & (p < 1)))
    if not informative:
        return dict(mean=0., median=0., sd=prior_sd,
                    low=float(norm.ppf(.025)*prior_sd), high=float(norm.ppf(.975)*prior_sd),
                    probability_positive=.5, no_information=True, prior_only=True,
                    variance_reduction=0., integration_error=0.)
    def log_density(t):
        return float(log_likelihood(t, c, n, p, rho) - .5*(t/prior_sd)**2)
    fit = minimize_scalar(lambda t: -log_density(t), bounds=(-40*prior_sd, 40*prior_sd),
                          method='bounded', options={'xatol': 1e-10})
    if not fit.success or abs(fit.x) >= 39*prior_sd:
        raise ArithmeticError('posterior mode not qualified')
    mode = float(fit.x)
    peak = log_density(mode)
    # Integration is split at the mode so adaptive quadrature cannot miss a
    # narrow posterior far from the prior centre. Twelve prior SDs beyond the
    # mode bracket the tails; independently verify by a different quadrature.
    lo, hi = min(0., mode)-12*prior_sd, max(0., mode)+12*prior_sd
    def density(t):
        return np.exp(log_density(t)-peak)
    def integral(fun, a, b):
        v, e = quad(fun, a, b, epsabs=tolerance, epsrel=tolerance, limit=200)
        if not np.isfinite(v) or not np.isfinite(e) or e > max(1e-9, abs(v)*1e-8):
            raise ArithmeticError('quadrature did not meet error budget')
        return v, e
    left, el = integral(density, lo, mode)
    right, er = integral(density, mode, hi)
    z = left+right
    def moment(k):
        f = lambda t: (t-mode)**k*density(t)
        return (integral(f, lo, mode)[0]+integral(f, mode, hi)[0])/z
    shift = moment(1)
    mean = mode+shift
    variance = max(0., moment(2)-shift**2)
    def cdf(t):
        if t <= mode:
            return integral(density, lo, t)[0]/z
        return (left+integral(density, mode, t)[0])/z
    quantiles = [brentq(lambda t: cdf(t)-q, lo, hi, xtol=2e-9) for q in (.025, .5, .975)]
    return dict(mean=mean, median=quantiles[1], sd=float(np.sqrt(variance)),
                low=quantiles[0], high=quantiles[2], probability_positive=1-cdf(0.),
                no_information=False, prior_only=False,
                variance_reduction=1-variance/prior_sd**2, integration_error=(el+er)/z)


def fit_dispersion(c, n, p):
    """Null theta=0 constrained composite likelihood; profile is diagnostic."""
    c, n, p = (np.asarray(a, float) for a in (c, n, p))
    if c.ndim != 2 or c.shape[1] != 2 or c.shape != n.shape or c.shape != p.shape:
        raise ValueError('aligned site by strand arrays required')
    for row in zip(c, n, p):
        validate(*row, 0.)
    def objective(rho):
        if rho == 0:
            from scipy.special import xlogy, xlog1py
            return -float((xlogy(c, p)+xlog1py(n-c, -p)).sum())
        good = (p > 0) & (p < 1) & (n > 0)
        a, b = p[good]*(1/rho-1), (1-p[good])*(1/rho-1)
        return -float((betaln(c[good]+a, n[good]-c[good]+b)-betaln(a, b)).sum())
    fit = minimize_scalar(objective, bounds=(1e-7, .95), method='bounded',
                          options={'xatol': 1e-9})
    candidates = [(0., objective(0.)), (float(fit.x), float(fit.fun)), (.95, objective(.95))]
    rho, minimum = min(candidates, key=lambda x: x[1])
    # Dense predetermined support grid and exact boundary roots. The usual
    # likelihood-ratio cutoff is an identification diagnostic, not a valid CI
    # under arbitrary dependence or a boundary parameter.
    cutoff = 1.920729410347062
    def delta(r):
        return objective(r)-minimum-cutoff
    low = brentq(delta, 0., rho) if rho > 0 and delta(0.) > 0 else 0.
    high = brentq(delta, rho, .95) if rho < .95 and delta(.95) > 0 else .95
    return dict(rho=rho, profile_low=low, profile_high=high,
                n_unique_sites=len(c), informative_sites=int((n.sum(1)>0).sum()),
                identified=bool(fit.success and high-low <= .25 and high < .95),
                profile_is_confidence_interval=False,
                nll=minimum, binomial_nll=objective(0.))


def threshold(scores):
    scores = np.asarray(scores, float)
    if not np.isfinite(scores).all():
        raise ValueError('calibration scores must be finite')
    if len(scores) < 200:
        return np.nan
    k = int(np.ceil(.95*(len(scores)+1)))
    return float(np.sort(scores)[k-1]) if k <= len(scores) else np.nan
