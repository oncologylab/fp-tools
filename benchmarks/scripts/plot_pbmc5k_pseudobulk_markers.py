#!/usr/bin/env python
"""Plot PBMC5k pseudobulk marker footprint summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import matplotlib.patheffects as patheffects
import numpy as np
import pandas as pd


MARKERS = {
    "B_cell": "PAX5",
    "Monocyte": "CEBPB",
    "T_NK_cell": "TCF7",
}

COLORS = {
    "B_cell": "#3B82F6",
    "Monocyte": "#D97706",
    "T_NK_cell": "#059669",
    "background": "#D1D5DB",
    "up": "#B91C1C",
    "down": "#1D4ED8",
}


def label_groups(ax: plt.Axes, annotations: pd.DataFrame, column: str, fontsize: int = 8) -> None:
    for label, group_df in annotations.groupby(column, sort=True):
        text = ax.text(
            group_df["umap_1"].median(),
            group_df["umap_2"].median(),
            str(label),
            fontsize=fontsize,
            weight="bold",
            ha="center",
            va="center",
            color="#111827",
        )
        text.set_path_effects([patheffects.withStroke(linewidth=2.2, foreground="white")])


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def plot_celltype_umap(annotations: pd.DataFrame, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for cell_type in ["B_cell", "Monocyte", "T_NK_cell"]:
        subset = annotations[annotations["cell_type"] == cell_type]
        ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=2.0,
            alpha=0.9,
            linewidths=0,
            color=COLORS[cell_type],
            label=cell_type,
            rasterized=True,
        )
    label_groups(ax, annotations, "cell_type", fontsize=9)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("PBMC5k broad cell types")
    ax.legend(frameon=False, markerscale=3, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_marker_umap(annotations: pd.DataFrame, aggregate_screen: pd.DataFrame, output_prefix: Path) -> None:
    aggregate_screen = aggregate_screen.set_index("tf")
    groups = ["B_cell", "Monocyte", "T_NK_cell"]
    values: dict[str, dict[str, float]] = {}
    for group, tf in MARKERS.items():
        row = aggregate_screen.loc[tf]
        values[tf] = {}
        for cell_type in groups:
            col = f"{cell_type}_center_protection"
            # More negative center footprint score is easier to read as stronger protection.
            values[tf][cell_type] = -float(row[col])

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.8), sharex=True, sharey=True)
    for ax, (cell_type, tf) in zip(axes, MARKERS.items(), strict=True):
        score = annotations["cell_type"].map(values[tf]).astype(float)
        vmax = float(np.nanmax(np.abs(score))) if np.isfinite(score).any() else 1.0
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        sc = ax.scatter(
            annotations["umap_1"],
            annotations["umap_2"],
            c=score,
            cmap="RdBu_r",
            norm=norm,
            s=2.0,
            alpha=0.92,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{tf} footprint signature")
        ax.set_xlabel("UMAP 1")
        ax.spines[["top", "right"]].set_visible(False)
        label_groups(ax, annotations, "cell_type", fontsize=7)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label("- center footprint score", fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    axes[0].set_ylabel("UMAP 2")
    fig.suptitle("PBMC5k marker footprint signatures projected onto broad-cell UMAP", y=1.03)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    rows = []
    for tf, group_values in values.items():
        for group, value in group_values.items():
            rows.append({"tf": tf, "cell_type": group, "projected_signature": value})
    pd.DataFrame(rows).to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)


def plot_volcano(results: pd.DataFrame, comparison: str, first: str, second: str, output_prefix: Path, markers: list[str]) -> None:
    change_col = f"{comparison}_change"
    p_col = f"{comparison}_pvalue"
    q_col = f"{comparison}_qvalue_bh"
    plot_df = results.copy()
    plot_df["change"] = coerce_numeric(plot_df[change_col])
    plot_df["pvalue"] = coerce_numeric(plot_df[p_col]).clip(lower=np.nextafter(0, 1))
    plot_df["qvalue"] = coerce_numeric(plot_df[q_col])
    plot_df["neg_log10_p"] = -np.log10(plot_df["pvalue"])
    plot_df["status"] = "not significant"
    plot_df.loc[(plot_df["qvalue"] <= 0.05) & (plot_df["change"] > 0), "status"] = "higher in first"
    plot_df.loc[(plot_df["qvalue"] <= 0.05) & (plot_df["change"] < 0), "status"] = "higher in second"

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for status, color, size, alpha in [
        ("not significant", COLORS["background"], 12, 0.45),
        ("higher in second", COLORS["down"], 18, 0.75),
        ("higher in first", COLORS["up"], 18, 0.75),
    ]:
        subset = plot_df[plot_df["status"] == status]
        ax.scatter(subset["change"], subset["neg_log10_p"], s=size, color=color, alpha=alpha, linewidths=0, label=status)

    for marker in markers:
        marker_rows = plot_df[plot_df["name"].astype(str).str.upper() == marker.upper()]
        for _, row in marker_rows.iterrows():
            ax.scatter(row["change"], row["neg_log10_p"], s=60, color="#111827", edgecolor="white", linewidth=0.7, zorder=4)
            ax.annotate(
                row["name"],
                (row["change"], row["neg_log10_p"]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                weight="bold",
            )

    ax.axvline(0, color="#374151", linewidth=0.8)
    ax.axhline(-np.log10(0.05), color="#6B7280", linewidth=0.8, linestyle="--")
    ax.set_xlabel(f"BINDetect change ({first} vs {second})")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{first} vs {second}")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), dpi=220)
    fig.savefig(output_prefix.with_suffix(".pdf"))
    plt.close(fig)

    plot_df[["output_prefix", "name", "total_tfbs", "change", "pvalue", "qvalue", "status"]].to_csv(
        output_prefix.with_suffix(".tsv"),
        sep="\t",
        index=False,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--aggregate-screen", required=True)
    parser.add_argument("--bindetect-results", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--markers", default="PAX5,CEBPB,TCF7")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    annotations = read_table(Path(args.annotations))
    aggregate_screen = read_table(Path(args.aggregate_screen))
    results = read_table(Path(args.bindetect_results))
    markers = [marker.strip() for marker in args.markers.split(",") if marker.strip()]

    plot_celltype_umap(annotations, outdir / "pbmc5k_umap_broad_celltypes_scprinter_style")
    plot_marker_umap(annotations, aggregate_screen, outdir / "pbmc5k_marker_footprint_umap")
    for comparison, first, second in [
        ("B_cell_Monocyte", "B_cell", "Monocyte"),
        ("B_cell_T_NK_cell", "B_cell", "T_NK_cell"),
        ("Monocyte_T_NK_cell", "Monocyte", "T_NK_cell"),
    ]:
        plot_volcano(results, comparison, first, second, outdir / f"pbmc5k_volcano_{comparison}", markers)
    print(f"Wrote PBMC5k marker UMAP and volcano plots to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
