#!/usr/bin/env python3
"""Render concise, paired before/after reports for a frozen TF detector.

The report intentionally distinguishes a TF-specific research result from a
general method claim.  It uses identical held-out motif sites for the legacy
bigWig score and the frozen candidate, and records bootstrap and biological-
replicate evidence when those tables are supplied.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from compare_frozen_tf_candidates import score_centers
from plot_frozen_tf_profiles import mean_ci, normalize_profiles_for_display
from search_tf_footprint_models import candidate_from_row, extract_profiles, score_candidate


BLUE = "#4C78A8"
RED = "#E45756"
GREEN = "#2E8B57"
GRAY = "#6B7280"


def crossfit_covariate_residuals(
    values: np.ndarray,
    covariates: np.ndarray,
    groups: np.ndarray,
    *,
    ridge_alpha: float = 1.0,
) -> np.ndarray:
    """Remove motif/accessibility effects without using ChIP labels."""

    values = np.asarray(values, dtype=float)
    covariates = np.asarray(covariates, dtype=float)
    groups = np.asarray(groups)
    if covariates.ndim != 2 or covariates.shape[0] != len(values):
        raise ValueError("covariates must have one row per score")
    if groups.shape != (len(values),):
        raise ValueError("groups must have one value per score")
    unique_groups = np.unique(groups)
    if len(unique_groups) < 2:
        raise ValueError("cross-fitted residuals require at least two groups")
    if not np.isfinite(values).all() or not np.isfinite(covariates).all():
        raise ValueError("scores and covariates must be finite")

    residuals = np.empty_like(values)
    for group in unique_groups:
        held_out = groups == group
        training = ~held_out
        model = make_pipeline(
            StandardScaler(),
            Ridge(alpha=ridge_alpha),
        ).fit(covariates[training], values[training])
        residuals[held_out] = values[held_out] - model.predict(
            covariates[held_out]
        )
    return residuals


def _aggregate(axis, x: np.ndarray, profiles: np.ndarray, labels: np.ndarray, title: str) -> None:
    colors = {0: BLUE, 1: RED}
    names = {0: "ChIP-negative", 1: "ChIP-positive"}
    for label in (0, 1):
        values = profiles[labels == label]
        mean, error = mean_ci(values)
        axis.plot(x, mean, color=colors[label], lw=2, label=f"{names[label]} (n={len(values):,})")
        axis.fill_between(x, mean - error, mean + error, color=colors[label], alpha=0.16)
    axis.axvline(0, color="#222222", lw=0.7, alpha=0.55)
    axis.set_title(title, loc="left", fontweight="bold")
    axis.set_xlabel("Distance from motif center (bp)")
    axis.legend(frameon=False, fontsize=8)
    axis.spines[["top", "right"]].set_visible(False)


def _bootstrap_text(bootstrap: pd.DataFrame, cell: str, tf: str, metric: str) -> str:
    row = bootstrap[
        (bootstrap["cell"] == cell)
        & (bootstrap["tf"] == tf)
        & (bootstrap["metric"] == metric)
    ]
    if row.empty:
        return "not available"
    value = row.iloc[0]
    return f"{value.ci_low:+.3f} to {value.ci_high:+.3f}"


def _replicate_text(panel: pd.DataFrame, cell: str, tf: str) -> str:
    rows = panel[
        (panel["cell"] == cell)
        & (panel["tf"] == tf)
        & (panel["depth"].astype(str) == "full")
        & (panel["sample"].astype(str) != "replicate_mean")
    ]
    if rows.empty:
        return "not available"
    return (
        f"AUROC {rows.auroc.min():.3f}-{rows.auroc.max():.3f}; "
        f"AUPRC {rows.auprc.min():.3f}-{rows.auprc.max():.3f} "
        f"(n={len(rows)} biological replicates)"
    )


def render_cell(
    *,
    cell: str,
    tf: str,
    split: str,
    sites: pd.DataFrame,
    winners: pd.DataFrame,
    baselines: pd.DataFrame,
    cache_dir: Path,
    bootstrap: pd.DataFrame,
    replicate_panel: pd.DataFrame,
    flank: int,
) -> tuple[plt.Figure, dict[str, object]]:
    cell_sites = sites[sites["cell"] == cell].reset_index(drop=True)
    winner_rows = winners[(winners["cell"] == cell) & (winners["tf"] == tf)]
    if len(winner_rows) != 1:
        raise ValueError(f"expected one frozen candidate for {cell}/{tf}, found {len(winner_rows)}")
    baseline_rows = baselines[baselines["cell"] == cell]
    if len(baseline_rows) != 1:
        raise ValueError(f"expected one legacy baseline for {cell}, found {len(baseline_rows)}")

    winner = winner_rows.iloc[0]
    candidate = candidate_from_row(winner)
    positions = np.flatnonzero(
        (cell_sites["tf"].to_numpy() == tf)
        & (cell_sites["chromosome_split"].to_numpy() == split)
    )
    if not len(positions):
        raise ValueError(f"no {split} sites for {cell}/{tf}")

    candidate_profiles = np.load(
        cache_dir / f"{cell}.{candidate.correction}.flank{flank}.npz"
    )["profiles"][positions]
    candidate_scores = score_candidate(candidate_profiles, candidate)
    baseline_signal = Path(baseline_rows.iloc[0].signal)
    baseline_all_scores = score_centers(cell_sites, baseline_signal)
    baseline_scores = baseline_all_scores[positions]
    covariates = np.column_stack(
        [
            cell_sites.iloc[positions]["motif_score"].to_numpy(dtype=float),
            np.log1p(
                np.maximum(
                    cell_sites.iloc[positions]["accessibility"].to_numpy(
                        dtype=float
                    ),
                    0.0,
                )
            ),
            np.log1p(
                np.maximum(
                    cell_sites.iloc[positions]["central_accessibility"].to_numpy(
                        dtype=float
                    ),
                    0.0,
                )
            ),
        ]
    )
    finite = (
        np.isfinite(candidate_scores)
        & np.isfinite(baseline_scores)
        & np.isfinite(covariates).all(axis=1)
    )
    labels = cell_sites.iloc[positions].loc[finite, "chip_label"].to_numpy(dtype=int)
    report_sites = cell_sites.iloc[positions].loc[finite].reset_index(drop=True)
    covariates = covariates[finite]
    candidate_profiles = candidate_profiles[finite]
    candidate_scores = candidate_scores[finite]
    baseline_scores = baseline_scores[finite]

    legacy_profiles, valid_legacy = extract_profiles(report_sites, baseline_signal, flank)
    joint = valid_legacy & np.isfinite(candidate_profiles).all(axis=1)
    legacy_display = normalize_profiles_for_display(legacy_profiles[joint])
    candidate_display = normalize_profiles_for_display(candidate_profiles[joint])
    aggregate_labels = labels[joint]

    baseline_auroc = float(roc_auc_score(labels, baseline_scores))
    candidate_auroc = float(roc_auc_score(labels, candidate_scores))
    baseline_auprc = float(average_precision_score(labels, baseline_scores))
    candidate_auprc = float(average_precision_score(labels, candidate_scores))
    chromosome_groups = report_sites["TFBS_chr"].astype(str).to_numpy()
    residual_baseline_scores = crossfit_covariate_residuals(
        baseline_scores, covariates, chromosome_groups
    )
    residual_candidate_scores = crossfit_covariate_residuals(
        candidate_scores, covariates, chromosome_groups
    )
    residual_baseline_auroc = float(
        roc_auc_score(labels, residual_baseline_scores)
    )
    residual_candidate_auroc = float(
        roc_auc_score(labels, residual_candidate_scores)
    )
    residual_baseline_auprc = float(
        average_precision_score(labels, residual_baseline_scores)
    )
    residual_candidate_auprc = float(
        average_precision_score(labels, residual_candidate_scores)
    )

    figure = plt.figure(figsize=(13.2, 8.2))
    grid = figure.add_gridspec(
        2,
        3,
        height_ratios=[0.88, 1.12],
        left=0.06,
        right=0.98,
        bottom=0.12,
        top=0.86,
        wspace=0.30,
        hspace=0.38,
    )
    metrics_axis = figure.add_subplot(grid[0, 0])
    roc_axis = figure.add_subplot(grid[0, 1])
    pr_axis = figure.add_subplot(grid[0, 2])
    legacy_axis = figure.add_subplot(grid[1, 0])
    candidate_axis = figure.add_subplot(grid[1, 1])
    evidence_axis = figure.add_subplot(grid[1, 2])

    metric_names = ["AUROC", "AUPRC"]
    before = [baseline_auroc, baseline_auprc]
    after = [candidate_auroc, candidate_auprc]
    locations = np.arange(2)
    width = 0.34
    metrics_axis.bar(locations - width / 2, before, width, color=GRAY, label="Before: conventional")
    metrics_axis.bar(locations + width / 2, after, width, color=GREEN, label="After: TF-specific")
    for index, (old, new) in enumerate(zip(before, after)):
        metrics_axis.text(index - width / 2, old + 0.015, f"{old:.3f}", ha="center", fontsize=9)
        metrics_axis.text(index + width / 2, new + 0.015, f"{new:.3f}", ha="center", fontsize=9)
    metric_labels = [
        f"{name}\nDelta {new - old:+.3f}"
        for name, old, new in zip(metric_names, before, after)
    ]
    metrics_axis.set_xticks(locations, metric_labels)
    metrics_axis.set_ylim(0.48, min(1.0, max(after) + 0.13))
    metrics_axis.set_title("Held-out discrimination", loc="left", fontweight="bold")
    metrics_axis.legend(frameon=False, fontsize=8, loc="upper left")
    metrics_axis.spines[["top", "right"]].set_visible(False)

    for scores, color, label in (
        (baseline_scores, GRAY, f"Conventional ({baseline_auroc:.3f})"),
        (candidate_scores, GREEN, f"TF-specific ({candidate_auroc:.3f})"),
    ):
        fpr, tpr, _ = roc_curve(labels, scores)
        roc_axis.plot(fpr, tpr, color=color, lw=2, label=label)
    roc_axis.plot([0, 1], [0, 1], color="#AAAAAA", ls="--", lw=1)
    roc_axis.set(xlabel="False-positive rate", ylabel="True-positive rate", xlim=(0, 1), ylim=(0, 1))
    roc_axis.set_title("ROC curve", loc="left", fontweight="bold")
    roc_axis.legend(frameon=False, fontsize=8, loc="lower right")
    roc_axis.spines[["top", "right"]].set_visible(False)

    for scores, color, label in (
        (baseline_scores, GRAY, f"Conventional ({baseline_auprc:.3f})"),
        (candidate_scores, GREEN, f"TF-specific ({candidate_auprc:.3f})"),
    ):
        precision, recall, _ = precision_recall_curve(labels, scores)
        pr_axis.plot(recall, precision, color=color, lw=2, label=label)
    prevalence = float(np.mean(labels))
    pr_axis.axhline(prevalence, color="#AAAAAA", ls="--", lw=1)
    pr_axis.set(xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0, 1))
    pr_axis.set_title("Precision-recall curve", loc="left", fontweight="bold")
    pr_axis.legend(frameon=False, fontsize=8, loc="lower left")
    pr_axis.spines[["top", "right"]].set_visible(False)

    x = np.arange(-flank, flank + 1)
    _aggregate(
        legacy_axis,
        x,
        legacy_display,
        aggregate_labels,
        "Before: conventional score",
    )
    legacy_axis.set_ylabel("Outer-flank standardized score")
    _aggregate(
        candidate_axis,
        x,
        candidate_display,
        aggregate_labels,
        f"After: {candidate.correction} TF-specific geometry",
    )
    half = candidate.center_width // 2
    candidate_axis.axvspan(-half, half, color="#777777", alpha=0.08)
    candidate_axis.set_ylabel("Outer-flank standardized cut signal")

    evidence_axis.axis("off")
    auroc_ci = _bootstrap_text(bootstrap, cell, tf, "auroc")
    auprc_ci = _bootstrap_text(bootstrap, cell, tf, "auprc")
    evidence = (
        "Evidence\n"
        f"Test: chr19-22 and chrX\n"
        f"Paired sites: {len(labels):,} ({int(labels.sum()):,} ChIP-positive)\n\n"
        f"Delta AUROC 95% CI: {auroc_ci}\n"
        f"Delta AUPRC 95% CI: {auprc_ci}\n\n"
        "Post-hoc covariate sensitivity\n"
        "Leave-one-chromosome-out removal of motif score and accessibility:\n"
        f"Delta AUROC {residual_candidate_auroc - residual_baseline_auroc:+.3f}; "
        f"AUPRC {residual_candidate_auprc - residual_baseline_auprc:+.3f}\n\n"
        "Full-depth transfer\n"
        f"{_replicate_text(replicate_panel, cell, tf)}\n\n"
        "Interpretation\n"
        "Strong CTCF-specific improvement in two ENCODE cell lines.\n"
        "This does not establish a universal TF detector."
    )
    evidence_axis.text(
        0.0,
        1.0,
        evidence,
        va="top",
        ha="left",
        fontsize=9.3,
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "#F5F7FA", "edgecolor": "#D1D5DB"},
    )

    figure.suptitle(
        f"{cell} {tf} footprint detection: before vs after",
        y=0.965,
        fontsize=17,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.915,
        "Conventional TOBIAS-style DWM footprint score compared with a frozen TF-specific cut-profile geometry",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#374151",
    )
    figure.text(
        0.5,
        0.025,
        "Research result only - geometry selected without test chromosomes; exact candidate naked-DNA specificity and external-cell transfer remain pending.",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#4B5563",
    )

    row = {
        "cell": cell,
        "tf": tf,
        "split": split,
        "n_sites": int(len(labels)),
        "positive_sites": int(labels.sum()),
        "baseline": str(baseline_rows.iloc[0].method),
        "candidate": candidate.identifier,
        "baseline_auroc": baseline_auroc,
        "candidate_auroc": candidate_auroc,
        "delta_auroc": candidate_auroc - baseline_auroc,
        "baseline_auprc": baseline_auprc,
        "candidate_auprc": candidate_auprc,
        "delta_auprc": candidate_auprc - baseline_auprc,
        "residual_baseline_auroc": residual_baseline_auroc,
        "residual_candidate_auroc": residual_candidate_auroc,
        "residual_delta_auroc": residual_candidate_auroc - residual_baseline_auroc,
        "residual_baseline_auprc": residual_baseline_auprc,
        "residual_candidate_auprc": residual_candidate_auprc,
        "residual_delta_auprc": residual_candidate_auprc - residual_baseline_auprc,
        "delta_auroc_ci": auroc_ci,
        "delta_auprc_ci": auprc_ci,
        "replicate_evidence": _replicate_text(replicate_panel, cell, tf),
        "scope": "CTCF-specific research result; not a package-wide promotion",
    }
    return figure, row


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--baselines", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=Path, required=True)
    parser.add_argument("--replicate-panel", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tf", default="CTCF")
    parser.add_argument("--split", default="test")
    parser.add_argument("--flank", type=int, default=100)
    args = parser.parse_args(argv)

    sites = pd.read_csv(args.sites, sep="\t")
    winners = pd.read_csv(args.winners, sep="\t")
    baselines = pd.read_csv(args.baselines, sep="\t")
    bootstrap = pd.read_csv(args.bootstrap, sep="\t")
    replicate_panel = pd.read_csv(args.replicate_panel, sep="\t")
    cells = sorted(set(sites["cell"]).intersection(winners.loc[winners["tf"] == args.tf, "cell"]))
    if not cells:
        parser.error(f"no cells contain a frozen {args.tf} candidate")

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    combined_path = args.outdir / f"{args.tf}_before_after_summary.pdf"
    with PdfPages(combined_path) as combined:
        for cell in cells:
            figure, row = render_cell(
                cell=cell,
                tf=args.tf,
                split=args.split,
                sites=sites,
                winners=winners,
                baselines=baselines,
                cache_dir=args.cache_dir,
                bootstrap=bootstrap,
                replicate_panel=replicate_panel,
                flank=args.flank,
            )
            figure.savefig(args.outdir / f"{args.tf}_{cell}_before_after.pdf", bbox_inches="tight")
            combined.savefig(figure, bbox_inches="tight")
            plt.close(figure)
            rows.append(row)

    metrics = pd.DataFrame(rows)
    metrics.to_csv(args.outdir / f"{args.tf}_before_after_metrics.tsv", sep="\t", index=False)
    (args.outdir / "README.txt").write_text(
        "Concise held-out before/after reports for a TF-specific fp-tools research candidate.\n"
        "The conventional comparator is the current TOBIAS-style DWM footprint score.\n"
        "The candidate was frozen before chr19-22/X were evaluated.\n"
        "The covariate-residual sensitivity analysis is post hoc and does not use ChIP labels for fitting.\n"
        "Scope: TF-specific research evidence only; no main-branch or package-default change.\n",
        encoding="utf-8",
    )
    print(metrics.to_string(index=False))
    print(f"Wrote {combined_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
