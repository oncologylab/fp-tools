#!/usr/bin/env python3
"""Audit label-free TF detector selection and the balanced bias/shift factorial.

The detector fits consumed no ChIP labels.  Development labels are opened here
only to evaluate frozen selection rules and to measure the gap to a development
oracle.  Locked GM12878 and IMR-90 labels are never read by this program.
"""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd


ARMS = (
    "MT_SELMA10_4m4",
    "MT_SELMA10_4m5",
    "MT_LOG81_4m4",
    "MT_LOG81_4m5",
)
TASK_KEYS = ("cell", "tf", "motif_family")
METRICS = ("auroc", "auprc", "brier")
FIXED_COMMON_CANDIDATES = (
    "count_gp.bg_none.window_30",
    "count_gp.bg_gp-long.window_30",
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*values: object, seed: int = 2026) -> int:
    digest = blake2b(digest_size=8)
    for value in (seed, *values):
        digest.update(str(value).encode())
        digest.update(b"\0")
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def _passing_evaluation(metrics: pd.DataFrame) -> pd.DataFrame:
    required = {
        *TASK_KEYS,
        "bias_configuration",
        "candidate_id",
        "status",
        "converged",
        "selection_score",
        *METRICS,
    }
    missing = sorted(required.difference(metrics.columns))
    if missing:
        raise ValueError(f"label-free evaluation metrics lack columns: {missing}")
    output = metrics[
        metrics["status"].eq("ok")
        & metrics["converged"].eq(True)
        & metrics["bias_configuration"].isin(ARMS)
    ].copy()
    duplicate = list(TASK_KEYS) + ["bias_configuration", "candidate_id"]
    if output.duplicated(duplicate).any():
        raise ValueError("evaluation metrics contain duplicate task/arm/candidate rows")
    return output


def select_development_oracle_common(metrics: pd.DataFrame) -> pd.DataFrame:
    """Choose one common detector across arms using development labels."""

    passing = _passing_evaluation(metrics)
    wide = passing.pivot(
        index=list(TASK_KEYS) + ["candidate_id"],
        columns="bias_configuration",
        values="selection_score",
    ).dropna(subset=list(ARMS))
    if wide.empty:
        return pd.DataFrame()
    complete = wide.reset_index()
    complete["selection_policy"] = "development_oracle_common"
    complete["selection_labels_used"] = True
    complete["selection_value"] = complete[list(ARMS)].mean(axis=1)
    complete["minimum_arm_value"] = complete[list(ARMS)].min(axis=1)
    return (
        complete.sort_values(
            list(TASK_KEYS)
            + ["selection_value", "minimum_arm_value", "candidate_id"],
            ascending=[True] * len(TASK_KEYS) + [False, False, True],
            kind="mergesort",
        )
        .groupby(list(TASK_KEYS), as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )


def _passing_unlabeled_cv(cv_metrics: pd.DataFrame) -> pd.DataFrame:
    required = {
        *TASK_KEYS,
        "bias_configuration",
        "candidate_id",
        "status",
        "converged",
        "profile_plausible",
        "heldout_gain_per_site",
        "heldout_gain_per_site_position",
        "heldout_posterior_entropy",
        "profile_depletion",
    }
    missing = sorted(required.difference(cv_metrics.columns))
    if missing:
        raise ValueError(f"unlabeled CV metrics lack columns: {missing}")
    output = cv_metrics[
        cv_metrics["status"].eq("ok")
        & cv_metrics["converged"].eq(True)
        & cv_metrics["profile_plausible"].eq(True)
        & cv_metrics["bias_configuration"].isin(ARMS)
    ].copy()
    duplicate = list(TASK_KEYS) + ["bias_configuration", "candidate_id"]
    if output.duplicated(duplicate).any():
        raise ValueError("unlabeled CV metrics contain duplicate task/arm/candidate rows")
    return output


def select_unlabeled_common(cv_metrics: pd.DataFrame) -> pd.DataFrame:
    """Choose a common detector by mean held-out unlabeled likelihood."""

    passing = _passing_unlabeled_cv(cv_metrics)
    wide = passing.pivot(
        index=list(TASK_KEYS) + ["candidate_id"],
        columns="bias_configuration",
        values="heldout_gain_per_site_position",
    ).dropna(subset=list(ARMS))
    if wide.empty:
        return pd.DataFrame()
    complete = wide.reset_index()
    complete["selection_policy"] = "unlabeled_likelihood_common"
    complete["selection_labels_used"] = False
    complete["selection_value"] = complete[list(ARMS)].mean(axis=1)
    complete["minimum_arm_value"] = complete[list(ARMS)].min(axis=1)
    return (
        complete.sort_values(
            list(TASK_KEYS)
            + ["selection_value", "minimum_arm_value", "candidate_id"],
            ascending=[True] * len(TASK_KEYS) + [False, False, True],
            kind="mergesort",
        )
        .groupby(list(TASK_KEYS), as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )


def fixed_common_settings(metrics: pd.DataFrame, candidate_id: str) -> pd.DataFrame:
    passing = _passing_evaluation(metrics)
    selected = passing[passing["candidate_id"].eq(candidate_id)].copy()
    counts = (
        selected.groupby(list(TASK_KEYS))["bias_configuration"]
        .nunique()
        .rename("arm_count")
        .reset_index()
    )
    complete = counts[counts["arm_count"].eq(len(ARMS))].drop(columns="arm_count")
    complete["candidate_id"] = candidate_id
    complete["selection_policy"] = f"prespecified_{candidate_id}"
    complete["selection_labels_used"] = False
    complete["selection_value"] = np.nan
    complete["minimum_arm_value"] = np.nan
    return complete


def factorial_effects(metrics: pd.DataFrame, settings: pd.DataFrame) -> pd.DataFrame:
    """Calculate balanced model, shift, interaction, and diagonal effects."""

    if settings.empty:
        return pd.DataFrame()
    passing = _passing_evaluation(metrics)
    join_keys = list(TASK_KEYS) + ["candidate_id"]
    policy_keys = join_keys + ["selection_policy"]
    if settings.duplicated(policy_keys).any():
        raise ValueError("factorial settings contain duplicate task/candidate/policy rows")
    selected = passing.merge(
        settings[
            join_keys
            + ["selection_policy", "selection_labels_used", "selection_value"]
        ],
        on=join_keys,
        how="inner",
        validate="many_to_many",
    )
    rows: list[dict] = []
    group_keys = join_keys + [
        "selection_policy",
        "selection_labels_used",
        "selection_value",
    ]
    for keys, group in selected.groupby(group_keys, sort=True, dropna=False):
        identity = dict(zip(group_keys, keys))
        indexed = group.set_index("bias_configuration")
        if not set(ARMS).issubset(indexed.index):
            continue
        for metric in METRICS:
            values = {arm: float(indexed.loc[arm, metric]) for arm in ARMS}
            if not all(np.isfinite(list(values.values()))):
                continue
            direction = -1.0 if metric == "brier" else 1.0
            oriented = {arm: direction * value for arm, value in values.items()}
            shift = 0.5 * (
                oriented["MT_SELMA10_4m4"]
                - oriented["MT_SELMA10_4m5"]
                + oriented["MT_LOG81_4m4"]
                - oriented["MT_LOG81_4m5"]
            )
            bias = 0.5 * (
                oriented["MT_SELMA10_4m4"]
                - oriented["MT_LOG81_4m4"]
                + oriented["MT_SELMA10_4m5"]
                - oriented["MT_LOG81_4m5"]
            )
            interaction = (
                oriented["MT_SELMA10_4m4"]
                - oriented["MT_SELMA10_4m5"]
                - oriented["MT_LOG81_4m4"]
                + oriented["MT_LOG81_4m5"]
            )
            total = oriented["MT_SELMA10_4m4"] - oriented["MT_LOG81_4m5"]
            if not np.isclose(total, shift + bias, atol=1e-12):
                raise AssertionError("balanced factorial contributions do not sum")
            best_arm = max(oriented, key=lambda arm: (oriented[arm], arm))
            rows.append(
                {
                    **identity,
                    "metric": metric,
                    "total_improvement": total,
                    "shift_contribution": shift,
                    "bias_model_contribution": bias,
                    "model_by_shift_interaction": interaction,
                    "best_arm": best_arm,
                    **{f"value_{arm.lower()}": values[arm] for arm in ARMS},
                }
            )
    return pd.DataFrame(rows)


def bootstrap_mean_interval(
    values: np.ndarray, *, iterations: int, seed: int
) -> tuple[float, float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return np.nan, np.nan, np.nan
    observed = float(np.mean(finite))
    if len(finite) < 2:
        return observed, observed, observed
    rng = np.random.default_rng(seed)
    draws = finite[rng.integers(0, len(finite), size=(iterations, len(finite)))]
    lower, upper = np.quantile(np.mean(draws, axis=1), [0.025, 0.975])
    return observed, float(lower), float(upper)


def summarize_factorial(
    effects: pd.DataFrame, *, bootstraps: int, seed: int
) -> pd.DataFrame:
    rows: list[dict] = []
    components = (
        "total_improvement",
        "shift_contribution",
        "bias_model_contribution",
        "model_by_shift_interaction",
    )
    for (policy, labels, metric), group in effects.groupby(
        ["selection_policy", "selection_labels_used", "metric"], sort=True
    ):
        for component in components:
            mean, lower, upper = bootstrap_mean_interval(
                group[component].to_numpy(dtype=float),
                iterations=bootstraps,
                seed=stable_seed(policy, metric, component, seed=seed),
            )
            rows.append(
                {
                    "selection_policy": policy,
                    "selection_labels_used": bool(labels),
                    "metric": metric,
                    "component": component,
                    "cell_tf_pairs": int(len(group)),
                    "mean": mean,
                    "bootstrap_lower": lower,
                    "bootstrap_upper": upper,
                    "median": float(group[component].median()),
                    "positive_fraction": float(np.mean(group[component] > 0)),
                }
            )
    return pd.DataFrame(rows)


def merge_cv_with_evaluation(
    cv_metrics: pd.DataFrame, evaluation_metrics: pd.DataFrame
) -> pd.DataFrame:
    cv = _passing_unlabeled_cv(cv_metrics)
    evaluation = _passing_evaluation(evaluation_metrics)
    columns = [
        *TASK_KEYS,
        "bias_configuration",
        "candidate_id",
        "family",
        "auroc",
        "auprc",
        "brier",
        "prevalence",
        "selection_score",
    ]
    available = [column for column in columns if column in evaluation]
    return cv.merge(
        evaluation[available],
        on=list(TASK_KEYS) + ["bias_configuration", "candidate_id"],
        how="inner",
        validate="one_to_one",
        suffixes=("_cv", "_evaluation"),
    )


def rank_correlations(merged: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    pairs = (
        ("heldout_gain_per_site_position", "auroc"),
        ("heldout_gain_per_site_position", "auprc"),
        ("heldout_gain_per_site_position", "selection_score"),
        ("profile_depletion", "auroc"),
        ("heldout_posterior_entropy", "auroc"),
    )
    for keys, group in merged.groupby(list(TASK_KEYS), sort=True):
        row = {**dict(zip(TASK_KEYS, keys)), "candidate_arm_rows": int(len(group))}
        for left, right in pairs:
            row[f"spearman_{left}_vs_{right}"] = float(
                group[left].corr(group[right], method="spearman")
            )
        rows.append(row)
    return pd.DataFrame(rows)


def _rank_within_task(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame.groupby(list(TASK_KEYS))[column].transform(
        lambda values: values.rank(method="average", pct=True)
    )


def _top_by_score(frame: pd.DataFrame, score: str, rule: str) -> pd.DataFrame:
    output = (
        frame.sort_values(
            list(TASK_KEYS) + [score, "bias_configuration", "candidate_id"],
            ascending=[True] * len(TASK_KEYS) + [False, True, True],
            kind="mergesort",
        )
        .groupby(list(TASK_KEYS), as_index=False, sort=True)
        .head(1)
        .copy()
    )
    output["selection_rule"] = rule
    output["selection_rule_score"] = output[score]
    output["selection_labels_used"] = False
    return output


def selection_rule_choices(merged: pd.DataFrame) -> pd.DataFrame:
    """Apply fixed label-free selection rules; labels only travel as outcomes."""

    values = merged.copy()
    values["rank_gain_position"] = _rank_within_task(
        values, "heldout_gain_per_site_position"
    )
    values["rank_depletion"] = _rank_within_task(values, "profile_depletion")
    values["rank_entropy"] = _rank_within_task(
        values, "heldout_posterior_entropy"
    )
    values["score_gain_depletion"] = (
        values["rank_gain_position"] + 0.5 * values["rank_depletion"]
    )
    values["score_depletion_low_entropy"] = (
        values["rank_depletion"] - 0.5 * values["rank_entropy"]
    )
    values["score_gain_depletion_entropy"] = (
        values["rank_gain_position"]
        + 0.5 * values["rank_depletion"]
        - 0.5 * values["rank_entropy"]
    )
    choices = [
        _top_by_score(
            values,
            "heldout_gain_per_site_position",
            "likelihood_per_position",
        ),
        _top_by_score(values, "heldout_gain_per_site", "likelihood_per_site"),
        _top_by_score(values, "profile_depletion", "profile_depletion"),
        _top_by_score(
            values, "score_gain_depletion", "likelihood_plus_depletion"
        ),
        _top_by_score(
            values,
            "score_depletion_low_entropy",
            "depletion_minus_entropy",
        ),
        _top_by_score(
            values,
            "score_gain_depletion_entropy",
            "likelihood_depletion_entropy",
        ),
    ]

    candidate_stability = (
        values.groupby(list(TASK_KEYS) + ["candidate_id"], as_index=False)
        .agg(
            arm_coverage=("bias_configuration", "nunique"),
            arm_mean_gain=("heldout_gain_per_site_position", "mean"),
            arm_sd_gain=("heldout_gain_per_site_position", "std"),
        )
    )
    maximum_coverage = candidate_stability.groupby(list(TASK_KEYS))[
        "arm_coverage"
    ].transform("max")
    candidate_stability = candidate_stability[
        candidate_stability["arm_coverage"].eq(maximum_coverage)
    ].copy()
    candidate_stability["cross_arm_mean_score"] = candidate_stability[
        "arm_mean_gain"
    ]
    candidate_stability["cross_arm_robust_score"] = candidate_stability[
        "arm_mean_gain"
    ] - 0.5 * candidate_stability["arm_sd_gain"].fillna(0.0)
    for score, rule in (
        ("cross_arm_mean_score", "cross_arm_mean_likelihood"),
        ("cross_arm_robust_score", "cross_arm_robust_likelihood"),
    ):
        stable = (
            candidate_stability.sort_values(
                list(TASK_KEYS) + [score, "candidate_id"],
                ascending=[True] * len(TASK_KEYS) + [False, True],
                kind="mergesort",
            )
            .groupby(list(TASK_KEYS), as_index=False, sort=True)
            .head(1)
        )
        selected = values.merge(
            stable[list(TASK_KEYS) + ["candidate_id", score, "arm_coverage"]],
            on=list(TASK_KEYS) + ["candidate_id"],
            how="inner",
            validate="many_to_one",
        )
        selected = _top_by_score(selected, "heldout_gain_per_site_position", rule)
        selected["selection_rule_score"] = selected[score]
        choices.append(selected)

    arm_stability = (
        values.groupby(list(TASK_KEYS) + ["bias_configuration"], as_index=False)
        .agg(arm_mean_candidate_gain=("heldout_gain_per_site_position", "mean"))
    )
    selected_arms = (
        arm_stability.sort_values(
            list(TASK_KEYS) + ["arm_mean_candidate_gain", "bias_configuration"],
            ascending=[True] * len(TASK_KEYS) + [False, True],
            kind="mergesort",
        )
        .groupby(list(TASK_KEYS), as_index=False, sort=True)
        .head(1)
    )
    arm_rows = values.merge(
        selected_arms[list(TASK_KEYS) + ["bias_configuration"]],
        on=list(TASK_KEYS) + ["bias_configuration"],
        how="inner",
        validate="many_to_one",
    )
    choices.append(
        _top_by_score(arm_rows, "heldout_gain_per_site_position", "arm_consensus")
    )
    return pd.concat(choices, ignore_index=True)


def fixed_dwm_rows(
    candidate_matrix: pd.DataFrame, fixed_candidate_id: str
) -> pd.DataFrame:
    required = {*TASK_KEYS, "bias_configuration", "candidate_id", *METRICS}
    missing = sorted(required.difference(candidate_matrix.columns))
    if missing:
        raise ValueError(f"candidate matrix lacks columns: {missing}")
    rows = candidate_matrix[
        candidate_matrix["bias_configuration"].eq("DWM")
        & candidate_matrix["candidate_id"].eq(fixed_candidate_id)
    ][list(TASK_KEYS) + list(METRICS)].copy()
    if rows.duplicated(list(TASK_KEYS)).any():
        raise ValueError("fixed DWM reference contains duplicate tasks")
    return rows.rename(columns={metric: f"dwm_{metric}" for metric in METRICS})


def score_rule_choices(
    choices: pd.DataFrame, dwm: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scored = choices.merge(dwm, on=list(TASK_KEYS), how="inner", validate="many_to_one")
    scored["auroc_gain_vs_dwm"] = scored["auroc"] - scored["dwm_auroc"]
    scored["auprc_gain_vs_dwm"] = scored["auprc"] - scored["dwm_auprc"]
    scored["relative_auprc_gain_vs_dwm"] = (
        scored["auprc"] / np.maximum(scored["dwm_auprc"], 1e-8) - 1.0
    )
    summary = (
        scored.groupby("selection_rule", sort=True)
        .agg(
            cell_tf_pairs=("tf", "size"),
            mean_auroc=("auroc", "mean"),
            mean_auprc=("auprc", "mean"),
            mean_brier=("brier", "mean"),
            mean_dwm_auroc=("dwm_auroc", "mean"),
            mean_dwm_auprc=("dwm_auprc", "mean"),
            mean_auroc_gain=("auroc_gain_vs_dwm", "mean"),
            mean_auprc_gain=("auprc_gain_vs_dwm", "mean"),
            relative_mean_auprc_gain=(
                "relative_auprc_gain_vs_dwm",
                "mean",
            ),
            maximum_auroc_loss=(
                "auroc_gain_vs_dwm",
                lambda values: float(np.maximum(-values, 0.0).max()),
            ),
            fraction_auroc_above_half=(
                "auroc",
                lambda values: float(np.mean(values > 0.5)),
            ),
            fraction_improving_dwm=(
                "auroc_gain_vs_dwm",
                lambda values: float(np.mean(values > 0.0)),
            ),
        )
        .reset_index()
    )
    summary["passes_development_effect_gates"] = (
        (summary["mean_auroc_gain"] >= 0.03)
        & (summary["relative_mean_auprc_gain"] >= 0.10)
        & (summary["maximum_auroc_loss"] <= 0.02)
    )
    return scored, summary


def development_oracle_rows(evaluation_metrics: pd.DataFrame) -> pd.DataFrame:
    passing = _passing_evaluation(evaluation_metrics)
    count = passing[passing["family"].eq("count")].copy()
    return (
        count.sort_values(
            list(TASK_KEYS)
            + ["selection_score", "auprc", "auroc", "bias_configuration", "candidate_id"],
            ascending=[True] * len(TASK_KEYS) + [False, False, False, True, True],
            kind="mergesort",
        )
        .groupby(list(TASK_KEYS), as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )


def selector_regret(
    rule_rows: pd.DataFrame, evaluation_metrics: pd.DataFrame
) -> pd.DataFrame:
    oracle = development_oracle_rows(evaluation_metrics)[
        list(TASK_KEYS)
        + ["bias_configuration", "candidate_id", "auroc", "auprc", "brier"]
    ].rename(
        columns={
            "bias_configuration": "oracle_bias_configuration",
            "candidate_id": "oracle_candidate_id",
            "auroc": "oracle_auroc",
            "auprc": "oracle_auprc",
            "brier": "oracle_brier",
        }
    )
    output = rule_rows.merge(oracle, on=list(TASK_KEYS), how="left", validate="many_to_one")
    output["auroc_regret_to_count_oracle"] = output["oracle_auroc"] - output["auroc"]
    output["auprc_regret_to_count_oracle"] = output["oracle_auprc"] - output["auprc"]
    output["selection_diagnosis"] = np.select(
        [
            (output["auroc"] < 0.5)
            & (output["oracle_auroc"] >= 0.6),
            output["auroc_regret_to_count_oracle"] >= 0.05,
            output["auroc_gain_vs_dwm"] < -0.02,
            (output["auroc_gain_vs_dwm"] >= 0.03)
            & (output["relative_auprc_gain_vs_dwm"] >= 0.10),
        ],
        [
            "unlabeled_criterion_or_orientation_failure",
            "unlabeled_selection_regret",
            "dwm_nonregression_failure",
            "promising_development_effect",
        ],
        default="small_or_mixed_effect",
    )
    return output


def fixed_route_summary(merged: pd.DataFrame, dwm: pd.DataFrame) -> pd.DataFrame:
    scored = merged.merge(dwm, on=list(TASK_KEYS), how="inner", validate="many_to_one")
    scored["auroc_gain_vs_dwm"] = scored["auroc"] - scored["dwm_auroc"]
    scored["auprc_gain_vs_dwm"] = scored["auprc"] - scored["dwm_auprc"]
    return (
        scored.groupby(["bias_configuration", "candidate_id"], sort=True)
        .agg(
            cell_tf_pairs=("tf", "size"),
            mean_auroc=("auroc", "mean"),
            mean_auprc=("auprc", "mean"),
            mean_auroc_gain=("auroc_gain_vs_dwm", "mean"),
            mean_auprc_gain=("auprc_gain_vs_dwm", "mean"),
            minimum_auroc_gain=("auroc_gain_vs_dwm", "min"),
        )
        .reset_index()
        .sort_values(
            ["cell_tf_pairs", "mean_auroc_gain", "mean_auprc_gain"],
            ascending=[False, False, False],
            kind="mergesort",
        )
    )


def render_audit(
    effects: pd.DataFrame,
    rule_summary: pd.DataFrame,
    rule_rows: pd.DataFrame,
    output: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output) as pdf:
        for policy in sorted(effects["selection_policy"].unique()):
            subset = effects[
                effects["selection_policy"].eq(policy)
                & effects["metric"].eq("auroc")
            ].copy()
            if subset.empty:
                continue
            subset["cell_tf"] = subset["cell"] + " — " + subset["tf"]
            subset = subset.sort_values("total_improvement", ascending=False)
            columns = [
                "bias_model_contribution",
                "shift_contribution",
                "total_improvement",
            ]
            matrix = subset[columns].to_numpy(dtype=float)
            limit = max(float(np.nanquantile(np.abs(matrix), 0.98)), 0.01)
            figure, axis = plt.subplots(
                figsize=(9.5, max(4.5, 0.32 * len(subset) + 1.8))
            )
            image = axis.imshow(
                matrix, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit
            )
            axis.set_xticks(
                range(3), labels=["Bias model", "Cut shift", "Diagonal total"]
            )
            axis.set_yticks(range(len(subset)), labels=subset["cell_tf"])
            axis.set_title(f"Balanced label-free-fit AUROC factorial — {policy}")
            for row in range(len(subset)):
                for column in range(3):
                    axis.text(
                        column,
                        row,
                        f"{matrix[row, column]:+.3f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                    )
            figure.colorbar(image, ax=axis, label="Positive is better")
            figure.tight_layout()
            pdf.savefig(figure)
            plt.close(figure)

        ordered = rule_summary.sort_values("mean_auroc_gain", ascending=True)
        figure, axis = plt.subplots(figsize=(10.0, 5.5))
        axis.barh(ordered["selection_rule"], ordered["mean_auroc_gain"], color="#4472C4")
        axis.axvline(0.0, color="black", linewidth=0.8)
        axis.axvline(0.03, color="#B22222", linewidth=0.8, linestyle="--")
        axis.set_xlabel("Mean AUROC gain versus fixed DWM")
        axis.set_title("Unlabeled detector-selection rules on development labels")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        best_rule = rule_summary.sort_values(
            ["passes_development_effect_gates", "mean_auroc_gain", "mean_auprc_gain"],
            ascending=[False, False, False],
            kind="mergesort",
        ).iloc[0]["selection_rule"]
        subset = rule_rows[rule_rows["selection_rule"].eq(best_rule)].copy()
        subset["cell_tf"] = subset["cell"] + " — " + subset["tf"]
        subset = subset.sort_values("auroc_gain_vs_dwm", ascending=True)
        figure, axis = plt.subplots(
            figsize=(9.5, max(4.5, 0.32 * len(subset) + 1.8))
        )
        colors = np.where(subset["auroc_gain_vs_dwm"] >= 0, "#2E8B57", "#B22222")
        axis.barh(subset["cell_tf"], subset["auroc_gain_vs_dwm"], color=colors)
        axis.axvline(0.0, color="black", linewidth=0.8)
        axis.set_xlabel("AUROC gain versus fixed DWM")
        axis.set_title(f"Per-TF development effects — {best_rule}")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-metrics", type=Path, required=True)
    parser.add_argument("--unlabeled-cv-metrics", type=Path, required=True)
    parser.add_argument("--candidate-matrix", type=Path, required=True)
    parser.add_argument("--policy-manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bootstraps", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    if args.bootstraps < 2:
        raise SystemExit("--bootstraps must be at least two")
    evaluation = pd.read_csv(args.evaluation_metrics, sep="\t")
    cv_metrics = pd.read_csv(args.unlabeled_cv_metrics, sep="\t")
    candidate_matrix = pd.read_csv(args.candidate_matrix, sep="\t")
    policy_document = json.loads(args.policy_manifest.read_text(encoding="utf-8"))
    fixed_dwm_candidate = policy_document["global_dwm_fallback"]["candidate_id"]

    settings_frames = [
        select_development_oracle_common(evaluation),
        select_unlabeled_common(cv_metrics),
    ]
    settings_frames.extend(
        fixed_common_settings(evaluation, candidate_id)
        for candidate_id in FIXED_COMMON_CANDIDATES
    )
    settings = pd.concat(
        [frame for frame in settings_frames if not frame.empty], ignore_index=True
    )
    effects = factorial_effects(evaluation, settings)
    factorial_summary = summarize_factorial(
        effects, bootstraps=args.bootstraps, seed=args.seed
    )
    merged = merge_cv_with_evaluation(cv_metrics, evaluation)
    correlations = rank_correlations(merged)
    dwm = fixed_dwm_rows(candidate_matrix, fixed_dwm_candidate)
    choices = selection_rule_choices(merged)
    scored_choices, rule_summary = score_rule_choices(choices, dwm)
    regret = selector_regret(scored_choices, evaluation)
    routes = fixed_route_summary(merged, dwm)

    args.outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "common_settings": args.outdir / "label_free_factorial_common_settings.tsv",
        "factorial_effects": args.outdir / "label_free_factorial_effects.tsv",
        "factorial_summary": args.outdir / "label_free_factorial_summary.tsv",
        "criterion_correlations": args.outdir / "unlabeled_criterion_correlations.tsv",
        "selection_rule_rows": args.outdir / "unlabeled_selection_rule_rows.tsv.gz",
        "selection_rule_summary": args.outdir / "unlabeled_selection_rule_summary.tsv",
        "selector_regret": args.outdir / "unlabeled_selector_regret.tsv.gz",
        "fixed_route_summary": args.outdir / "unlabeled_fixed_route_summary.tsv",
        "audit_pdf": args.outdir / "label_free_factorial_selection_audit.pdf",
    }
    settings.to_csv(outputs["common_settings"], sep="\t", index=False)
    effects.to_csv(outputs["factorial_effects"], sep="\t", index=False)
    factorial_summary.to_csv(outputs["factorial_summary"], sep="\t", index=False)
    correlations.to_csv(outputs["criterion_correlations"], sep="\t", index=False)
    scored_choices.to_csv(outputs["selection_rule_rows"], sep="\t", index=False)
    rule_summary.to_csv(outputs["selection_rule_summary"], sep="\t", index=False)
    regret.to_csv(outputs["selector_regret"], sep="\t", index=False)
    routes.to_csv(outputs["fixed_route_summary"], sep="\t", index=False)
    render_audit(effects, rule_summary, scored_choices, outputs["audit_pdf"])
    manifest = {
        "schema": "fp-tools-label-free-factorial-selection-audit-v1",
        "locked_holdout_labels_read": False,
        "model_training_labels_used": False,
        "development_labels_used_for_evaluation": True,
        "development_oracle_is_not_deployable": True,
        "selection_rules_use_labels": False,
        "arms": list(ARMS),
        "fixed_common_candidates": list(FIXED_COMMON_CANDIDATES),
        "fixed_dwm_candidate": fixed_dwm_candidate,
        "inputs": {
            str(path): file_sha256(path)
            for path in (
                args.evaluation_metrics,
                args.unlabeled_cv_metrics,
                args.candidate_matrix,
                args.policy_manifest,
            )
        },
        "bootstraps": args.bootstraps,
        "seed": args.seed,
        "outputs": {
            name: {"path": str(path), "sha256": file_sha256(path)}
            for name, path in outputs.items()
        },
    }
    manifest_path = args.outdir / "label_free_factorial_selection_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(factorial_summary.to_string(index=False))
    print("\nUnlabeled selection rules")
    print(rule_summary.sort_values("mean_auroc_gain", ascending=False).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
