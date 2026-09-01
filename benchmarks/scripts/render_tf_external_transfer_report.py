#!/usr/bin/env python3
"""Render a concise before/after report for frozen external TF transfer."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_tf_geometry_external_transfer import (  # noqa: E402
    conventional_profile_scores,
)
from fp_tools.utils.signals import footprint_score_array_fast  # noqa: E402
from plot_frozen_tf_profiles import (  # noqa: E402
    mean_ci,
    normalize_profiles_for_display,
)


BLUE = "#4C78A8"
RED = "#E45756"
GREEN = "#2E8B57"
GRAY = "#6B7280"


def parse_cell_path(value: str) -> tuple[str, Path]:
    fields = value.split(",", 1)
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError("profile must use CELL,PATH")
    return fields[0], Path(fields[1])


def conventional_score_profiles(profiles: np.ndarray) -> np.ndarray:
    values = np.asarray(profiles, dtype=float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("profiles must be a finite matrix")
    return np.stack(
        [footprint_score_array_fast(row, 10, 30, 20, 50) for row in values]
    )


def plot_aggregate(
    axis: plt.Axes,
    profiles: np.ndarray,
    labels: np.ndarray,
    *,
    title: str,
) -> None:
    x = np.arange(-(profiles.shape[1] // 2), profiles.shape[1] // 2 + 1)
    for label, color, name in (
        (0, BLUE, "ChIP-negative"),
        (1, RED, "ChIP-positive"),
    ):
        subset = profiles[labels == label]
        mean, error = mean_ci(subset)
        axis.plot(x, mean, color=color, lw=1.8, label=f"{name} (n={len(subset):,})")
        axis.fill_between(x, mean - error, mean + error, color=color, alpha=0.15)
    axis.axvline(0, color="#222222", lw=0.7, alpha=0.55)
    axis.set_title(title, loc="left", fontweight="bold")
    axis.set_xlabel("Distance from motif center (bp)")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, fontsize=7.5)


def evidence_text(row: pd.Series) -> str:
    return (
        f"Sites: {int(row.sites):,} ({int(row.positive_sites):,} ChIP-positive)\n\n"
        "AUROC\n"
        f"{row.conventional_auroc:.3f} -> {row.candidate_auroc:.3f}  "
        f"(Delta {row.delta_auroc:+.3f})\n"
        f"95% CI {row.delta_auroc_ci_low:+.3f} to {row.delta_auroc_ci_high:+.3f}\n\n"
        "AUPRC\n"
        f"{row.conventional_auprc:.3f} -> {row.candidate_auprc:.3f}  "
        f"(Delta {row.delta_auprc:+.3f})\n"
        f"95% CI {row.delta_auprc_ci_low:+.3f} to {row.delta_auprc_ci_high:+.3f}\n\n"
        "Label-free covariate residual\n"
        f"Delta AUROC {row.residual_delta_auroc:+.3f}; "
        f"AUPRC {row.residual_delta_auprc:+.3f}\n\n"
        "Both biological replicates improve\n"
        f"Minimum Delta AUROC {row.replicate_min_delta_auroc:+.3f}; "
        f"AUPRC {row.replicate_min_delta_auprc:+.3f}"
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument(
        "--profile",
        action="append",
        type=parse_cell_path,
        required=True,
        metavar="CELL,NPZ",
    )
    parser.add_argument("--tf", default="CTCF")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    summary = pd.read_csv(args.summary, sep="\t")
    profiles = dict(args.profile)
    cells = list(profiles)
    if len(cells) != 2:
        parser.error("the concise report requires exactly two external cells")

    figure, axes = plt.subplots(2, 3, figsize=(14.8, 8.2))
    figure.subplots_adjust(
        left=0.055,
        right=0.985,
        bottom=0.11,
        top=0.86,
        hspace=0.48,
        wspace=0.30,
    )
    for row_index, cell in enumerate(cells):
        rows = summary[
            summary["cell"].astype(str).eq(cell)
            & summary["tf"].astype(str).eq(args.tf)
        ]
        if len(rows) != 1:
            raise ValueError(f"expected one summary row for {cell}/{args.tf}")
        data = np.load(profiles[cell])
        corrected = np.asarray(data["profiles"], dtype=float)
        labels = np.asarray(data["labels"], dtype=int)
        conventional = conventional_score_profiles(corrected)
        center = corrected.shape[1] // 2
        expected = conventional_profile_scores(corrected, center)
        if not np.allclose(conventional[:, center], expected):
            raise AssertionError("conventional profile center mismatch")
        plot_aggregate(
            axes[row_index, 0],
            normalize_profiles_for_display(conventional),
            labels,
            title=f"{cell}: conventional score",
        )
        axes[row_index, 0].set_ylabel("Outer-flank standardized score")
        plot_aggregate(
            axes[row_index, 1],
            normalize_profiles_for_display(corrected),
            labels,
            title=f"{cell}: frozen TF geometry input",
        )
        axes[row_index, 1].axvspan(-20, 20, color="#777777", alpha=0.08)
        axes[row_index, 1].set_ylabel("Outer-flank standardized cut signal")
        axes[row_index, 2].axis("off")
        axes[row_index, 2].text(
            0,
            1,
            evidence_text(rows.iloc[0]),
            va="top",
            ha="left",
            fontsize=9.1,
            linespacing=1.25,
            bbox={
                "boxstyle": "round,pad=0.55",
                "facecolor": "#F5F7FA",
                "edgecolor": "#D1D5DB",
            },
        )

    figure.suptitle(
        f"{args.tf} footprint detection: external transfer before vs after",
        y=0.965,
        fontsize=17,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.915,
        "K562-frozen DWM cut-profile geometry applied without retuning to two external ENCODE cell types",
        ha="center",
        fontsize=10.5,
        color="#374151",
    )
    figure.text(
        0.5,
        0.025,
        "Post-hoc no-retuning transfer evidence; the candidate was selected without external labels. CTCF-specific research result, not package-wide promotion.",
        ha="center",
        fontsize=8.5,
        color="#4B5563",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
