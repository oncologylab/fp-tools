#!/usr/bin/env python3
"""Integrate bias, residual, FDA, GP, and aggregate evidence per TF.

All model fitting and route selection inputs are label-free unless a column is
explicitly named ``development_best`` or ``supervised``.  Those columns are
development diagnostics, never deployable selections.  Locked holdout labels
are not accepted or read.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd


TASK_KEYS = ("cell", "tf", "motif_family")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def add_selection_score(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    if "selection_score" not in output.columns:
        prevalence = output["prevalence"].to_numpy(dtype=float)
        output["selection_score"] = output["auroc"].to_numpy(dtype=float) + (
            output["auprc"].to_numpy(dtype=float) - prevalence
        ) / np.maximum(1.0 - prevalence, 1e-8)
    return output


def best_rows(
    frame: pd.DataFrame,
    *,
    group_keys: tuple[str, ...] = TASK_KEYS,
    score: str = "selection_score",
) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    return (
        frame.sort_values(
            list(group_keys) + [score, "auprc", "auroc"],
            ascending=[True] * len(group_keys) + [False, False, False],
            kind="mergesort",
        )
        .groupby(list(group_keys), as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )


def current_dwm(
    candidate_matrix: pd.DataFrame, fixed_candidate_id: str
) -> pd.DataFrame:
    selected = candidate_matrix[
        candidate_matrix["bias_configuration"].eq("DWM")
        & candidate_matrix["candidate_id"].eq(fixed_candidate_id)
    ].copy()
    if selected.duplicated(list(TASK_KEYS)).any():
        raise ValueError("fixed DWM reference contains duplicate tasks")
    rename = {
        "candidate_id": "current_candidate_id",
        "auroc": "current_auroc",
        "auprc": "current_auprc",
        "brier": "current_brier",
        "prevalence": "prevalence",
        "n_sites": "evaluation_sites",
        "positive_sites": "positive_sites",
        "negative_sites": "negative_sites",
        "depletion": "current_profile_depletion",
        "width": "current_profile_width",
        "asymmetry": "current_profile_asymmetry",
        "periodicity": "current_profile_periodicity",
    }
    columns = list(TASK_KEYS) + [column for column in rename if column in selected]
    return selected[columns].rename(columns=rename)


def supervised_ceiling(pilot_metrics: pd.DataFrame) -> pd.DataFrame:
    metrics = add_selection_score(pilot_metrics)
    baseline = best_rows(metrics[metrics["method"].eq("supervised_baseline")]).copy()
    ceiling = best_rows(metrics[metrics["method"].eq("supervised_fpca")]).copy()
    baseline = baseline[
        list(TASK_KEYS) + ["auroc", "auprc", "brier", "correction"]
    ].rename(
        columns={
            "auroc": "supervised_baseline_auroc",
            "auprc": "supervised_baseline_auprc",
            "brier": "supervised_baseline_brier",
            "correction": "supervised_baseline_correction",
        }
    )
    ceiling = ceiling[
        list(TASK_KEYS) + ["auroc", "auprc", "brier", "correction"]
    ].rename(
        columns={
            "auroc": "supervised_ceiling_auroc",
            "auprc": "supervised_ceiling_auprc",
            "brier": "supervised_ceiling_brier",
            "correction": "supervised_ceiling_correction",
        }
    )
    output = baseline.merge(ceiling, on=list(TASK_KEYS), how="outer", validate="one_to_one")
    output["supervised_functional_relative_auprc_gain"] = (
        output["supervised_ceiling_auprc"]
        / np.maximum(output["supervised_baseline_auprc"], 1e-8)
        - 1.0
    )
    return output


def residual_screen(pilot_metrics: pd.DataFrame) -> pd.DataFrame:
    metrics = add_selection_score(pilot_metrics)
    residuals = metrics[metrics["method"].astype(str).str.startswith("residual_")].copy()
    dwm = residuals[residuals["correction"].eq("DWM")]
    difference = best_rows(dwm[dwm["method"].eq("residual_difference")]).copy()
    best_dwm = best_rows(dwm).copy()
    best_parametric = best_rows(residuals[~residuals["correction"].eq("DWM")]).copy()

    def choose(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
        columns = list(TASK_KEYS) + ["method", "correction", "auroc", "auprc", "brier"]
        return frame[columns].rename(
            columns={
                "method": f"{prefix}_method",
                "correction": f"{prefix}_correction",
                "auroc": f"{prefix}_auroc",
                "auprc": f"{prefix}_auprc",
                "brier": f"{prefix}_brier",
            }
        )

    output = choose(difference, "difference_residual").merge(
        choose(best_dwm, "best_dwm_residual"),
        on=list(TASK_KEYS),
        how="outer",
        validate="one_to_one",
    ).merge(
        choose(best_parametric, "best_legacy_parametric_residual"),
        on=list(TASK_KEYS),
        how="outer",
        validate="one_to_one",
    )
    output["same_bias_residual_auroc_gain"] = (
        output["best_dwm_residual_auroc"] - output["difference_residual_auroc"]
    )
    output["same_bias_residual_auprc_gain"] = (
        output["best_dwm_residual_auprc"] - output["difference_residual_auprc"]
    )
    output["legacy_parametric_residual_auroc_gain_vs_dwm"] = (
        output["best_legacy_parametric_residual_auroc"]
        - output["best_dwm_residual_auroc"]
    )
    return output


def development_models(evaluation_metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    passing = add_selection_score(
        evaluation_metrics[
            evaluation_metrics["status"].eq("ok")
            & evaluation_metrics["converged"].eq(True)
        ]
    )
    overall = best_rows(passing).copy()
    overall_columns = list(TASK_KEYS) + [
        "bias_configuration",
        "candidate_id",
        "family",
        "smoother",
        "background",
        "window",
        "channel",
        "training_pool",
        "auroc",
        "auprc",
        "brier",
        "functional_separation",
        "depletion",
        "width",
        "shoulder_distance",
        "asymmetry",
        "periodicity",
        "fit_seconds",
    ]
    overall = overall[[column for column in overall_columns if column in overall]].rename(
        columns={
            "bias_configuration": "development_best_bias_configuration",
            "candidate_id": "development_best_candidate_id",
            "family": "development_best_family",
            "smoother": "development_best_smoother",
            "background": "development_best_background",
            "window": "development_best_window",
            "channel": "development_best_channel",
            "training_pool": "development_best_training_pool",
            "auroc": "development_best_auroc",
            "auprc": "development_best_auprc",
            "brier": "development_best_brier",
            "functional_separation": "development_best_functional_separation",
            "depletion": "development_best_profile_depletion",
            "width": "development_best_profile_width",
            "shoulder_distance": "development_best_shoulder_distance",
            "asymmetry": "development_best_asymmetry",
            "periodicity": "development_best_periodicity",
            "fit_seconds": "development_best_fit_seconds",
        }
    )

    family_rows = []
    comparisons = []
    for family_label, selector in (
        ("count_gp", passing[passing["family"].eq("count") & passing["smoother"].eq("gp")]),
        (
            "count_spline",
            passing[passing["family"].eq("count") & passing["smoother"].eq("spline")],
        ),
        ("fda", passing[passing["family"].eq("fda")]),
        ("hybrid", passing[passing["family"].eq("hybrid")]),
    ):
        selected = best_rows(selector).copy()
        selected["model_class"] = family_label
        family_rows.append(selected)
    family_best = pd.concat(family_rows, ignore_index=True)
    for keys, group in family_best.groupby(list(TASK_KEYS), sort=True):
        identity = dict(zip(TASK_KEYS, keys))
        row: dict[str, object] = dict(identity)
        for model_class, selected in group.groupby("model_class", sort=True):
            best = selected.iloc[0]
            prefix = str(model_class)
            row[f"{prefix}_bias_configuration"] = best["bias_configuration"]
            row[f"{prefix}_candidate_id"] = best["candidate_id"]
            row[f"{prefix}_auroc"] = best["auroc"]
            row[f"{prefix}_auprc"] = best["auprc"]
            row[f"{prefix}_brier"] = best["brier"]
            row[f"{prefix}_fit_seconds"] = best["fit_seconds"]
        comparisons.append(row)
    comparison = pd.DataFrame(comparisons)
    comparison["gp_relative_auprc_gain_over_spline"] = (
        comparison["count_gp_auprc"]
        / np.maximum(comparison["count_spline_auprc"], 1e-8)
        - 1.0
    )
    comparison["gp_auroc_gain_over_spline"] = (
        comparison["count_gp_auroc"] - comparison["count_spline_auroc"]
    )
    comparison["gp_runtime_ratio_over_spline"] = (
        comparison["count_gp_fit_seconds"]
        / np.maximum(comparison["count_spline_fit_seconds"], 1e-8)
    )
    comparison["gp_passes_relative_auprc_gate"] = (
        comparison["gp_relative_auprc_gain_over_spline"] >= 0.05
    )
    return overall, comparison


def unlabeled_selected(rule_rows: pd.DataFrame, rule: str) -> pd.DataFrame:
    selected = rule_rows[rule_rows["selection_rule"].eq(rule)].copy()
    if selected.duplicated(list(TASK_KEYS)).any():
        raise ValueError("unlabeled rule contains duplicate task selections")
    columns = list(TASK_KEYS) + [
        "bias_configuration",
        "candidate_id",
        "family_evaluation",
        "smoother",
        "background",
        "window",
        "channel",
        "heldout_gain_per_site_position",
        "profile_depletion",
        "profile_width",
        "heldout_posterior_entropy",
        "auroc",
        "auprc",
        "brier",
        "auroc_gain_vs_dwm",
        "auprc_gain_vs_dwm",
        "relative_auprc_gain_vs_dwm",
    ]
    selected = selected[[column for column in columns if column in selected]].rename(
        columns={
            "bias_configuration": "unlabeled_bias_configuration",
            "candidate_id": "unlabeled_candidate_id",
            "family_evaluation": "unlabeled_family",
            "smoother": "unlabeled_smoother",
            "background": "unlabeled_background",
            "window": "unlabeled_window",
            "channel": "unlabeled_channel",
            "heldout_gain_per_site_position": "unlabeled_likelihood_gain_per_position",
            "profile_depletion": "unlabeled_profile_depletion",
            "profile_width": "unlabeled_profile_width",
            "heldout_posterior_entropy": "unlabeled_posterior_entropy",
            "auroc": "unlabeled_auroc",
            "auprc": "unlabeled_auprc",
            "brier": "unlabeled_brier",
            "auroc_gain_vs_dwm": "unlabeled_auroc_gain_vs_dwm",
            "auprc_gain_vs_dwm": "unlabeled_auprc_gain_vs_dwm",
            "relative_auprc_gain_vs_dwm": "unlabeled_relative_auprc_gain_vs_dwm",
        }
    )
    selected["unlabeled_selection_rule"] = rule
    return selected


def per_tf_factorial(factorial_effects: pd.DataFrame, policy: str) -> pd.DataFrame:
    selected = factorial_effects[factorial_effects["selection_policy"].eq(policy)].copy()
    if selected.empty:
        raise ValueError(f"factorial effects lack policy {policy}")
    values = selected.pivot(
        index=list(TASK_KEYS),
        columns="metric",
        values=[
            "bias_model_contribution",
            "shift_contribution",
            "total_improvement",
            "model_by_shift_interaction",
        ],
    )
    values.columns = [f"factorial_{component}_{metric}" for component, metric in values.columns]
    return values.reset_index()


def aggregate_evidence(summary: pd.DataFrame, unlabeled_method: str) -> pd.DataFrame:
    method_map = {
        "Current DWM": "current",
        unlabeled_method: "unlabeled",
        "Development best (diagnostic)": "development_best",
    }
    rows = []
    for keys, group in summary.groupby(list(TASK_KEYS), sort=True):
        row: dict[str, object] = dict(zip(TASK_KEYS, keys))
        for method, prefix in method_map.items():
            selected = group[group["method"].eq(method)]
            if selected.empty:
                continue
            item = selected.iloc[0]
            for column in (
                "visual_rms_difference",
                "center_mean_difference",
                "shoulder_mean_difference",
                "protection_shape_contrast",
                "band_exclusion_fraction_within_50bp",
            ):
                row[f"{prefix}_aggregate_{column}"] = item[column]
        rows.append(row)
    output = pd.DataFrame(rows)
    output["unlabeled_protection_contrast_gain"] = (
        output["unlabeled_aggregate_protection_shape_contrast"]
        - output["current_aggregate_protection_shape_contrast"]
    )
    output["development_best_protection_contrast_gain"] = (
        output["development_best_aggregate_protection_shape_contrast"]
        - output["current_aggregate_protection_shape_contrast"]
    )
    return output


def bias_control_table(selection: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for item in selection.itertuples(index=False):
        model = "SELMA10" if str(item.model) == "selma10" else "LOG81"
        label = f"MT_{model}_4m{abs(int(item.shift_reverse))}"
        rows.append(
            {
                "bias_configuration": label,
                "control_conditional_nll": float(item.mean_conditional_nll),
                "control_nll_gain": float(item.mean_nll_gain),
                "control_passed_likelihood": bool(item.passed_control_likelihood),
                "control_model_size_mb": float(item.model_size_mb),
                "control_runtime_seconds": float(item.runtime_seconds),
            }
        )
    return pd.DataFrame(rows).drop_duplicates("bias_configuration")


def motif_response_flags(summary: pd.DataFrame) -> pd.DataFrame:
    columns = [
        *TASK_KEYS,
        "model_label",
        "center_flank_log_bias_effect",
        "potential_motif_response_requires_review",
    ]
    return summary[columns].rename(
        columns={
            "model_label": "bias_configuration",
            "center_flank_log_bias_effect": "motif_bias_response_effect",
            "potential_motif_response_requires_review": "motif_bias_response_requires_review",
        }
    )


def classify_diagnosis(row: pd.Series) -> tuple[str, str, str]:
    positives = float(row.get("positive_sites", np.nan))
    current_auroc = float(row.get("current_auroc", np.nan))
    current_ap = float(row.get("current_auprc", np.nan))
    prevalence = float(row.get("prevalence", np.nan))
    development_auroc = float(row.get("development_best_auroc", np.nan))
    unlabeled_auroc = float(row.get("unlabeled_auroc", np.nan))
    supervised_gain = float(row.get("supervised_functional_relative_auprc_gain", np.nan))
    residual_gain = float(row.get("same_bias_residual_auroc_gain", np.nan))
    bias_effect = float(row.get("factorial_bias_model_contribution_auroc", np.nan))
    current_adjusted_ap = (current_ap - prevalence) / max(1.0 - prevalence, 1e-8)
    development_gain = development_auroc - current_auroc
    selection_regret = development_auroc - unlabeled_auroc
    flags = []
    if np.isfinite(bias_effect) and bias_effect >= 0.03:
        flags.append("functional_bias_sensitive")
    if np.isfinite(residual_gain) and residual_gain >= 0.03:
        flags.append("residual_sensitive")
    if np.isfinite(supervised_gain) and supervised_gain < 0.10:
        flags.append("low_supervised_shape_increment")
    if np.isfinite(selection_regret) and selection_regret >= 0.05:
        flags.append("unlabeled_selection_regret")

    if np.isfinite(positives) and positives < 25:
        classification = "power_limited"
        recommendation = "Acquire more positive sites or a richer occupancy benchmark before model conclusions."
    elif current_auroc >= 0.65 and current_adjusted_ap >= 0.10:
        classification = "detectable_current"
        recommendation = "Keep the DWM route until external non-regression and naked-DNA gates pass."
    elif development_auroc >= 0.60 and development_gain >= 0.03:
        if unlabeled_auroc >= 0.58 and development_auroc - unlabeled_auroc < 0.05:
            classification = "detectable_with_label_free_route"
            recommendation = "Externally validate the frozen label-free route; do not promote before safety gates."
        else:
            classification = "shape_model_selection_limited"
            recommendation = "The ATAC shape is informative, but the unlabeled selector must be made TF/family stable."
    elif np.isfinite(supervised_gain) and supervised_gain < 0.10:
        classification = "assay_limited_or_motif_ambiguous"
        recommendation = "Test depth and orthogonal occupancy; ATAC shape adds little over motif plus accessibility."
    elif np.isfinite(residual_gain) and residual_gain >= 0.03:
        classification = "residual_formulation_limited"
        recommendation = "Prioritize residual calibration with the same bias model before changing sequence bias."
    elif np.isfinite(bias_effect) and bias_effect >= 0.03:
        classification = "bias_model_sensitive_unconfirmed"
        recommendation = "Bias choice changes performance, but require superior control likelihood and naked-DNA safety."
    else:
        classification = "currently_undetectable_or_low_information"
        recommendation = "Treat as unresolved; test depth, replicate stability, and orthogonal occupancy evidence."
    return classification, ",".join(flags), recommendation


def assemble_diagnosis(
    current: pd.DataFrame,
    supervised: pd.DataFrame,
    residual: pd.DataFrame,
    development: pd.DataFrame,
    comparison: pd.DataFrame,
    unlabeled: pd.DataFrame,
    factorial: pd.DataFrame,
    aggregates: pd.DataFrame,
    control: pd.DataFrame,
    motif_response: pd.DataFrame | None,
) -> pd.DataFrame:
    output = current.copy()
    for frame in (supervised, residual, development, comparison, unlabeled, factorial, aggregates):
        output = output.merge(frame, on=list(TASK_KEYS), how="left", validate="one_to_one")
    output["development_best_auroc_gain_vs_dwm"] = (
        output["development_best_auroc"] - output["current_auroc"]
    )
    output["development_best_auprc_gain_vs_dwm"] = (
        output["development_best_auprc"] - output["current_auprc"]
    )
    output["development_best_relative_auprc_gain_vs_dwm"] = (
        output["development_best_auprc"] / np.maximum(output["current_auprc"], 1e-8) - 1.0
    )
    output["unlabeled_selection_auroc_regret"] = (
        output["development_best_auroc"] - output["unlabeled_auroc"]
    )

    control_by_config = control.set_index("bias_configuration")
    for route in ("development_best", "unlabeled"):
        config = output[f"{route}_bias_configuration"]
        for column in control.columns:
            if column == "bias_configuration":
                continue
            output[f"{route}_{column}"] = config.map(control_by_config[column])
    if motif_response is not None and not motif_response.empty:
        for route in ("development_best", "unlabeled"):
            routes = output[list(TASK_KEYS) + [f"{route}_bias_configuration"]].rename(
                columns={f"{route}_bias_configuration": "bias_configuration"}
            )
            joined = routes.merge(
                motif_response,
                on=list(TASK_KEYS) + ["bias_configuration"],
                how="left",
                validate="one_to_one",
            )
            output = output.merge(
                joined[
                    list(TASK_KEYS)
                    + ["motif_bias_response_effect", "motif_bias_response_requires_review"]
                ].rename(
                    columns={
                        "motif_bias_response_effect": f"{route}_motif_bias_response_effect",
                        "motif_bias_response_requires_review": f"{route}_motif_bias_response_requires_review",
                    }
                ),
                on=list(TASK_KEYS),
                how="left",
                validate="one_to_one",
            )

    diagnoses = output.apply(classify_diagnosis, axis=1)
    output["failure_classification"] = [item[0] for item in diagnoses]
    output["evidence_flags"] = [item[1] for item in diagnoses]
    output["recommended_next_experiment"] = [item[2] for item in diagnoses]
    candidate_control = output.get(
        "development_best_control_conditional_nll", pd.Series(np.nan, index=output.index)
    )
    reference_nll = float(
        control.loc[
            control["bias_configuration"].eq("MT_LOG81_4m5"),
            "control_conditional_nll",
        ].iloc[0]
    )
    output["bias_limited_confirmed"] = (
        (output["factorial_bias_model_contribution_auroc"] >= 0.03)
        & (candidate_control < reference_nll)
        & output.get(
            "development_best_motif_bias_response_requires_review",
            pd.Series(False, index=output.index),
        ).eq(False)
    )
    output["motif_ambiguity_status"] = "not_identifiable_from_current_single_motif_labels"
    output["depth_status"] = "pending_full_depth_matrix"
    return output.sort_values(list(TASK_KEYS)).reset_index(drop=True)


def summarize_tf_consistency(diagnosis: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (tf, family), group in diagnosis.groupby(["tf", "motif_family"], sort=True):
        current_hard = group["current_auroc"] < 0.60
        development_rescue = (
            (group["development_best_auroc"] >= 0.60)
            & (group["development_best_auroc_gain_vs_dwm"] >= 0.03)
        )
        unlabeled_rescue = (
            (group["unlabeled_auroc"] >= 0.58)
            & (group["unlabeled_auroc_gain_vs_dwm"] >= 0.0)
        )
        rows.append(
            {
                "tf": tf,
                "motif_family": family,
                "contexts": int(len(group)),
                "cells": ",".join(sorted(group["cell"].astype(str))),
                "mean_current_auroc": float(group["current_auroc"].mean()),
                "minimum_current_auroc": float(group["current_auroc"].min()),
                "mean_development_best_auroc": float(group["development_best_auroc"].mean()),
                "mean_development_auroc_gain": float(
                    group["development_best_auroc_gain_vs_dwm"].mean()
                ),
                "mean_unlabeled_auroc": float(group["unlabeled_auroc"].mean()),
                "mean_unlabeled_auroc_gain": float(
                    group["unlabeled_auroc_gain_vs_dwm"].mean()
                ),
                "consistently_hard_current": bool(current_hard.all()),
                "development_rescue_in_all_contexts": bool(development_rescue.all()),
                "unlabeled_rescue_in_all_contexts": bool(unlabeled_rescue.all()),
                "classifications": ",".join(
                    sorted(group["failure_classification"].dropna().unique())
                ),
                "recommended_bias_configurations": ",".join(
                    sorted(group["development_best_bias_configuration"].dropna().unique())
                ),
                "recommended_candidate_ids": ",".join(
                    sorted(group["development_best_candidate_id"].dropna().unique())
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["consistently_hard_current", "mean_development_auroc_gain"],
        ascending=[False, False],
        kind="mergesort",
    )


def render_diagnosis(diagnosis: pd.DataFrame, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output.parent.mkdir(parents=True, exist_ok=True)
    frame = diagnosis.copy()
    frame["cell_tf"] = frame["cell"] + " — " + frame["tf"]
    with PdfPages(output) as pdf:
        columns = ["current_auroc", "unlabeled_auroc", "development_best_auroc"]
        matrix = frame[columns].to_numpy(dtype=float)
        figure, axis = plt.subplots(figsize=(9.5, max(5.0, 0.34 * len(frame) + 2.0)))
        image = axis.imshow(matrix, aspect="auto", cmap="viridis", vmin=0.3, vmax=0.85)
        axis.set_xticks(range(3), labels=["Current DWM", "Unlabeled route", "Development best"])
        axis.set_yticks(range(len(frame)), labels=frame["cell_tf"])
        axis.set_title("Per-TF AUROC: deployable evidence versus development potential")
        for row in range(len(frame)):
            for column in range(3):
                value = matrix[row, column]
                axis.text(column, row, "NA" if not np.isfinite(value) else f"{value:.3f}", ha="center", va="center", color="white" if np.isfinite(value) and value < 0.6 else "black", fontsize=7)
        figure.colorbar(image, ax=axis, label="AUROC")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        effect_columns = [
            "factorial_bias_model_contribution_auroc",
            "factorial_shift_contribution_auroc",
            "same_bias_residual_auroc_gain",
            "gp_auroc_gain_over_spline",
            "unlabeled_auroc_gain_vs_dwm",
        ]
        effect_labels = ["Bias model", "Cut shift", "Residual", "GP vs spline", "Unlabeled route"]
        matrix = frame[effect_columns].to_numpy(dtype=float)
        limit = max(float(np.nanquantile(np.abs(matrix), 0.98)), 0.03)
        figure, axis = plt.subplots(figsize=(11.0, max(5.0, 0.34 * len(frame) + 2.0)))
        image = axis.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
        axis.set_xticks(range(len(effect_labels)), labels=effect_labels)
        axis.set_yticks(range(len(frame)), labels=frame["cell_tf"])
        axis.set_title("Mechanistic AUROC effects (positive is better)")
        for row in range(len(frame)):
            for column in range(len(effect_columns)):
                value = matrix[row, column]
                axis.text(column, row, "NA" if not np.isfinite(value) else f"{value:+.3f}", ha="center", va="center", fontsize=7)
        figure.colorbar(image, ax=axis, label="AUROC change")
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-matrix", type=Path, required=True)
    parser.add_argument("--policy-manifest", type=Path, required=True)
    parser.add_argument("--pilot-metrics", type=Path, required=True)
    parser.add_argument("--evaluation-metrics", type=Path, required=True)
    parser.add_argument("--selection-rule-rows", type=Path, required=True)
    parser.add_argument("--selection-rule", default="likelihood_plus_depletion")
    parser.add_argument("--factorial-effects", type=Path, required=True)
    parser.add_argument(
        "--factorial-policy",
        default="prespecified_count_gp.bg_gp-long.window_30",
    )
    parser.add_argument("--aggregate-summary", type=Path, required=True)
    parser.add_argument("--bias-control-selection", type=Path, required=True)
    parser.add_argument("--motif-response", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    candidate_matrix = pd.read_csv(args.candidate_matrix, sep="\t")
    policy = json.loads(args.policy_manifest.read_text(encoding="utf-8"))
    fixed_candidate = policy["global_dwm_fallback"]["candidate_id"]
    pilot = pd.read_csv(args.pilot_metrics, sep="\t")
    evaluation = pd.read_csv(args.evaluation_metrics, sep="\t")
    rule_rows = pd.read_csv(args.selection_rule_rows, sep="\t")
    effects = pd.read_csv(args.factorial_effects, sep="\t")
    aggregate = pd.read_csv(args.aggregate_summary, sep="\t")
    control_selection = pd.read_csv(args.bias_control_selection, sep="\t")
    response = (
        motif_response_flags(pd.read_csv(args.motif_response, sep="\t"))
        if args.motif_response
        else None
    )

    development, comparison = development_models(evaluation)
    unlabeled = unlabeled_selected(rule_rows, args.selection_rule)
    unlabeled_method = f"Unlabeled: {args.selection_rule}"
    diagnosis = assemble_diagnosis(
        current_dwm(candidate_matrix, fixed_candidate),
        supervised_ceiling(pilot),
        residual_screen(pilot),
        development,
        comparison,
        unlabeled,
        per_tf_factorial(effects, args.factorial_policy),
        aggregate_evidence(aggregate, unlabeled_method),
        bias_control_table(control_selection),
        response,
    )
    consistency = summarize_tf_consistency(diagnosis)
    classification_counts = (
        diagnosis.groupby("failure_classification", sort=True)
        .agg(cell_tf_pairs=("tf", "size"), tfs=("tf", "nunique"))
        .reset_index()
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "diagnosis": args.outdir / "per_tf_detectability_diagnosis.tsv",
        "consistency": args.outdir / "tf_consistency_summary.tsv",
        "classification_counts": args.outdir / "failure_classification_counts.tsv",
        "pdf": args.outdir / "per_tf_detectability_diagnosis.pdf",
    }
    diagnosis.to_csv(paths["diagnosis"], sep="\t", index=False)
    consistency.to_csv(paths["consistency"], sep="\t", index=False)
    classification_counts.to_csv(paths["classification_counts"], sep="\t", index=False)
    render_diagnosis(diagnosis, paths["pdf"])
    input_paths = [
        args.candidate_matrix,
        args.policy_manifest,
        args.pilot_metrics,
        args.evaluation_metrics,
        args.selection_rule_rows,
        args.factorial_effects,
        args.aggregate_summary,
        args.bias_control_selection,
    ]
    if args.motif_response:
        input_paths.append(args.motif_response)
    manifest = {
        "schema": "fp-tools-per-tf-detectability-diagnosis-v1",
        "locked_holdout_labels_read": False,
        "development_labels_used_for_diagnostic_oracle": True,
        "development_oracle_is_not_deployable": True,
        "unlabeled_selection_rule": args.selection_rule,
        "factorial_policy": args.factorial_policy,
        "depth_status": "pending_full_depth_matrix",
        "motif_ambiguity_status": "not_identifiable_from_current_single_motif_labels",
        "inputs": {str(path): file_sha256(path) for path in input_paths},
        "outputs": {
            name: {"path": str(path), "sha256": file_sha256(path)}
            for name, path in paths.items()
        },
    }
    (args.outdir / "per_tf_detectability_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    columns = [
        "cell",
        "tf",
        "failure_classification",
        "current_auroc",
        "unlabeled_auroc",
        "development_best_auroc",
        "factorial_bias_model_contribution_auroc",
        "factorial_shift_contribution_auroc",
        "gp_relative_auprc_gain_over_spline",
    ]
    print(diagnosis[columns].to_string(index=False))
    print("\nTF consistency")
    print(consistency.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
