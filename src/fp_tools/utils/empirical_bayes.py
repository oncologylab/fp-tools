"""Empirical-Bayes variance moderation for replicate-level motif scores.

The implementation follows the scaled-inverse-chi-square variance model used
for moderated t statistics. Biological samples are the units of replication;
motif sites are summarized within each sample before this module is called.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import digamma, polygamma
from scipy.stats import t as student_t


@dataclass(frozen=True)
class VariancePrior:
    """Scaled-inverse-chi-square prior for residual variances."""

    degrees_of_freedom: float
    variance: float


def benjamini_hochberg(pvalues: Sequence[float]) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted values, preserving missing entries."""

    values = np.asarray(pvalues, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(values))
    if not len(finite):
        return result
    clipped = np.clip(values[finite], 0.0, 1.0)
    order = np.argsort(clipped, kind="stable")
    ranked = clipped[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result[finite[order]] = np.minimum(adjusted, 1.0)
    return result


def estimate_variance_prior(
    residual_variances: Sequence[float],
    residual_degrees_of_freedom: int,
) -> VariancePrior:
    """Estimate a scaled-inverse-chi-square prior by log-moment matching.

    If the observed log-variance dispersion is no greater than the sampling
    dispersion expected at ``residual_degrees_of_freedom``, the prior degrees
    of freedom is infinite and every posterior variance equals the common
    prior variance.
    """

    if residual_degrees_of_freedom <= 0:
        raise ValueError("residual_degrees_of_freedom must be positive")
    variances = np.asarray(residual_variances, dtype=float)
    variances = variances[np.isfinite(variances) & (variances >= 0)]
    if len(variances) < 2:
        raise ValueError("At least two finite residual variances are required")

    positive = variances[variances > 0]
    if not len(positive):
        raise ValueError("At least one positive residual variance is required")
    variance_floor = max(float(np.median(positive)) * np.finfo(float).eps, np.finfo(float).tiny)
    log_variances = np.log(np.maximum(variances, variance_floor))
    mean_log_variance = float(np.mean(log_variances))
    observed_log_variance = float(np.var(log_variances, ddof=1))
    sampling_log_variance = float(polygamma(1, residual_degrees_of_freedom / 2.0))
    target = observed_log_variance - sampling_log_variance

    if target <= 1e-12:
        prior_df = np.inf
        prior_log_adjustment = 0.0
    else:
        prior_half_df = brentq(
            lambda value: float(polygamma(1, value)) - target,
            1e-8,
            1e8,
        )
        prior_df = float(2.0 * prior_half_df)
        prior_log_adjustment = float(digamma(prior_half_df) - np.log(prior_half_df))

    residual_half_df = residual_degrees_of_freedom / 2.0
    log_prior_variance = (
        mean_log_variance
        - float(digamma(residual_half_df))
        + float(np.log(residual_half_df))
        + prior_log_adjustment
    )
    return VariancePrior(prior_df, float(np.exp(log_prior_variance)))


def fit_moderated_contrast(
    score_matrix: pd.DataFrame,
    sample_conditions: Mapping[str, str],
    condition_1: str,
    condition_2: str,
    *,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Fit a two-condition moderated t model to per-sample motif scores.

    Parameters
    ----------
    score_matrix
        Rows are motifs and columns are biological samples. Values must already
        be motif-level summaries, such as the mean footprint score across the
        common motif-site universe for each sample.
    sample_conditions
        Mapping from every score-matrix column used in the model to its
        biological condition.
    condition_1, condition_2
        Contrast direction, reported as ``condition_1 - condition_2``.
    alpha
        False-discovery-rate threshold for the convenience significance column.
    """

    samples_1 = [sample for sample, condition in sample_conditions.items() if condition == condition_1]
    samples_2 = [sample for sample, condition in sample_conditions.items() if condition == condition_2]
    if len(samples_1) < 2 or len(samples_2) < 2:
        raise ValueError("Empirical-Bayes testing requires at least two biological replicates per condition")
    missing = [sample for sample in samples_1 + samples_2 if sample not in score_matrix.columns]
    if missing:
        raise ValueError(f"Score matrix is missing sample column(s): {', '.join(missing)}")

    values_1 = score_matrix[samples_1].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    values_2 = score_matrix[samples_2].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    complete = np.isfinite(values_1).all(axis=1) & np.isfinite(values_2).all(axis=1)
    residual_df = len(samples_1) + len(samples_2) - 2
    output = pd.DataFrame(index=score_matrix.index)
    for column in (
        "condition_1_mean",
        "condition_2_mean",
        "effect",
        "residual_variance",
        "prior_variance",
        "prior_df",
        "posterior_variance",
        "moderated_se",
        "moderated_t",
        "moderated_df",
        "pvalue",
        "qvalue_bh",
        "ci_lower",
        "ci_upper",
    ):
        output[column] = np.nan
    output["significant_fdr05"] = False
    if np.count_nonzero(complete) < 2:
        raise ValueError("At least two motifs with complete replicate scores are required")

    mean_1 = np.mean(values_1[complete], axis=1)
    mean_2 = np.mean(values_2[complete], axis=1)
    effect = mean_1 - mean_2
    residual_sum_squares = (
        np.sum((values_1[complete] - mean_1[:, None]) ** 2, axis=1)
        + np.sum((values_2[complete] - mean_2[:, None]) ** 2, axis=1)
    )
    residual_variance = residual_sum_squares / residual_df
    prior = estimate_variance_prior(residual_variance, residual_df)
    if np.isinf(prior.degrees_of_freedom):
        posterior_variance = np.full_like(residual_variance, prior.variance)
        moderated_df = np.full_like(residual_variance, np.inf)
    else:
        posterior_variance = (
            prior.degrees_of_freedom * prior.variance + residual_df * residual_variance
        ) / (prior.degrees_of_freedom + residual_df)
        moderated_df = np.full_like(residual_variance, prior.degrees_of_freedom + residual_df)

    unscaled_se = np.sqrt(1.0 / len(samples_1) + 1.0 / len(samples_2))
    moderated_se = np.sqrt(posterior_variance) * unscaled_se
    moderated_t = np.divide(
        effect,
        moderated_se,
        out=np.full_like(effect, np.nan),
        where=moderated_se > 0,
    )
    pvalue = 2.0 * student_t.sf(np.abs(moderated_t), moderated_df)
    qvalue = benjamini_hochberg(pvalue)
    critical = student_t.ppf(1.0 - alpha / 2.0, moderated_df)

    valid_index = output.index[complete]
    output.loc[valid_index, "condition_1_mean"] = mean_1
    output.loc[valid_index, "condition_2_mean"] = mean_2
    output.loc[valid_index, "effect"] = effect
    output.loc[valid_index, "residual_variance"] = residual_variance
    output.loc[valid_index, "prior_variance"] = prior.variance
    output.loc[valid_index, "prior_df"] = prior.degrees_of_freedom
    output.loc[valid_index, "posterior_variance"] = posterior_variance
    output.loc[valid_index, "moderated_se"] = moderated_se
    output.loc[valid_index, "moderated_t"] = moderated_t
    output.loc[valid_index, "moderated_df"] = moderated_df
    output.loc[valid_index, "pvalue"] = pvalue
    output.loc[valid_index, "qvalue_bh"] = qvalue
    output.loc[valid_index, "ci_lower"] = effect - critical * moderated_se
    output.loc[valid_index, "ci_upper"] = effect + critical * moderated_se
    output.loc[valid_index, "significant_fdr05"] = qvalue <= alpha
    return output
