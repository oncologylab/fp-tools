#!/usr/bin/env python
"""Plot PBMC5k pseudobulk marker footprint summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

from adjustText import adjust_text
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import matplotlib.patheffects as patheffects
import matplotlib.text as mtext
import numpy as np
import pandas as pd


plt.rcParams.update(
    {
        "font.size": 8,
        "axes.titlesize": 8.8,
        "axes.labelsize": 8.4,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.2,
        "font.family": "Arial",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "svg.fonttype": "none",
    }
)

MARKERS = {
    "B_cell": "PAX5",
    "Monocyte": "CEBPB",
    "T_NK_cell": "TCF7",
}
PAIRWISE_COMPARISONS = [
    ("B_cell_Monocyte", "B_cell", "Monocyte"),
    ("B_cell_T_NK_cell", "B_cell", "T_NK_cell"),
    ("Monocyte_T_NK_cell", "Monocyte", "T_NK_cell"),
]
CELL_TYPES = ("B_cell", "Monocyte", "T_NK_cell")
PSEUDOBULK_HEATMAP_ROW_HEIGHT_IN = 0.065
PSEUDOBULK_LABELED_HEATMAP_ROW_HEIGHT_IN = 0.14
PSEUDOBULK_HEATMAP_MIN_BODY_HEIGHT_IN = 3.2
PSEUDOBULK_HEATMAP_EXTRA_HEIGHT_IN = 1.7
MARKER_GROUPS = {
    "B_cell": ("PAX5", "EBF1", "POU2F2", "POU2AF1", "BCL6", "SPIB"),
    "Monocyte": ("CEBPB", "CEBPA"),
    "T_NK_cell": ("TCF7", "LEF1", "ZBTB7B", "RUNX3", "GATA3"),
}

COLORS = {
    "B_cell": "#3B82F6",
    "Monocyte": "#D97706",
    "T_NK_cell": "#059669",
    "background": "#D1D5DB",
    "up": "#B91C1C",
    "down": "#1D4ED8",
}


def save_illustrator_svg(fig: plt.Figure, output_prefix: Path) -> None:
    for text in fig.findobj(match=mtext.Text):
        text.set_fontfamily("Arial")
        text.set_fontsize(9)
        text.set_fontweight("bold")
    fig.savefig(output_prefix.with_suffix(".svg"), bbox_inches="tight")


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


def marker_to_group(marker_groups: dict[str, tuple[str, ...]]) -> dict[str, str]:
    return {marker.upper(): group for group, markers in marker_groups.items() for marker in markers}


def prepare_volcano_df(results: pd.DataFrame, comparison: str) -> pd.DataFrame:
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
    return plot_df


def directional_marker_rows(
    plot_df: pd.DataFrame,
    first: str,
    second: str,
    marker_groups: dict[str, tuple[str, ...]] = MARKER_GROUPS,
    qvalue_threshold: float = 0.05,
) -> pd.DataFrame:
    group_by_marker = marker_to_group(marker_groups)
    labels = []
    for _, row in plot_df.iterrows():
        marker = str(row["name"]).upper()
        marker_group = group_by_marker.get(marker)
        if marker_group not in {first, second}:
            continue
        change = float(row["change"])
        qvalue = float(row["qvalue"])
        if not np.isfinite(change) or not np.isfinite(qvalue) or qvalue > qvalue_threshold:
            continue
        if (marker_group == first and change > 0) or (marker_group == second and change < 0):
            out = row.copy()
            out["marker_group"] = marker_group
            labels.append(out)
    if not labels:
        return plot_df.iloc[0:0].assign(marker_group=pd.Series(dtype=str))
    return pd.DataFrame(labels)


def annotate_marker_labels(ax: plt.Axes, labels: pd.DataFrame, fontsize: int) -> list[dict[str, object]]:
    source_rows = []
    texts = []
    for _, row in labels.iterrows():
        ax.scatter(row["change"], row["neg_log10_p"], s=52, color="#111827", edgecolor="white", linewidth=0.7, zorder=4)
        text = ax.text(
            float(row["change"]),
            float(row["neg_log10_p"]),
            row["name"],
            fontsize=fontsize,
            weight="bold",
            ha="center",
            va="bottom",
            clip_on=False,
            zorder=5,
        )
        text.set_path_effects([patheffects.withStroke(linewidth=2.0, foreground="white")])
        texts.append(text)
        source_rows.append(row.to_dict())
    if texts:
        adjust_text(
            texts,
            ax=ax,
            x=labels["change"].to_numpy(dtype=float),
            y=labels["neg_log10_p"].to_numpy(dtype=float),
            only_move={"text": "xy", "static": "xy", "explode": "xy", "pull": "xy"},
            force_text=(0.28, 0.42),
            force_static=(0.12, 0.18),
            force_pull=(0.03, 0.05),
            expand=(1.12, 1.25),
            max_move=(12, 12),
            min_arrow_len=6,
            ensure_inside_axes=True,
            expand_axes=False,
            arrowprops={"arrowstyle": "-", "color": "#4B5563", "linewidth": 0.45, "alpha": 0.7},
        )
    return source_rows


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
    save_illustrator_svg(fig, output_prefix)
    plt.close(fig)


def plot_marker_umap(annotations: pd.DataFrame, aggregate_screen: pd.DataFrame, output_prefix: Path) -> None:
    aggregate_screen = aggregate_screen.set_index("tf")
    values: dict[str, dict[str, float]] = {}
    for group, tf in MARKERS.items():
        row = aggregate_screen.loc[tf]
        values[tf] = {}
        for cell_type in CELL_TYPES:
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
    save_illustrator_svg(fig, output_prefix)
    plt.close(fig)

    rows = []
    for tf, group_values in values.items():
        for group, value in group_values.items():
            rows.append({"tf": tf, "cell_type": group, "projected_signature": value})
    pd.DataFrame(rows).to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)


def prepare_all_signature_heatmap(results: pd.DataFrame) -> pd.DataFrame:
    score_cols = [f"{cell_type}_mean_score" for cell_type in CELL_TYPES]
    change_cols = [f"{comparison}_change" for comparison, _, _ in PAIRWISE_COMPARISONS]
    required = {"output_prefix", "name", "motif_id", "total_tfbs", *score_cols, *change_cols}
    missing = required.difference(results.columns)
    if missing:
        raise SystemExit(f"BINDetect results table is missing required columns: {', '.join(sorted(missing))}")

    table = results[["output_prefix", "name", "motif_id", "total_tfbs", *score_cols, *change_cols]].copy()
    for col in ["total_tfbs", *score_cols, *change_cols]:
        table[col] = coerce_numeric(table[col])
    table = table.dropna(subset=change_cols).reset_index(drop=True)

    # Convert pairwise differential footprint scores into relative cell-type scores.
    # The fitted values have sum zero per motif and preserve the pairwise directions.
    design = np.array(
        [
            [1.0, -1.0, 0.0],
            [1.0, 0.0, -1.0],
            [0.0, 1.0, -1.0],
            [1.0, 1.0, 1.0],
        ]
    )
    pairwise_values = np.column_stack([table[col].to_numpy(dtype=float) for col in change_cols])
    fitted = []
    for row in pairwise_values:
        effect, *_ = np.linalg.lstsq(design, np.concatenate([row, [0.0]]), rcond=None)
        fitted.append(effect)
    raw_values = np.vstack(fitted)
    row_mean = np.nanmean(raw_values, axis=1, keepdims=True)
    row_sd = np.nanstd(raw_values, axis=1, keepdims=True)
    z_values = np.divide(raw_values - row_mean, row_sd, out=np.zeros_like(raw_values), where=row_sd > 0)

    table["mean_differential_score"] = row_mean[:, 0]
    table["dynamic_range"] = np.nanmax(raw_values, axis=1) - np.nanmin(raw_values, axis=1)
    table["dominant_cell_type"] = [CELL_TYPES[index] for index in np.nanargmax(z_values, axis=1)]
    for index, cell_type in enumerate(CELL_TYPES):
        table[f"{cell_type}_differential_score"] = raw_values[:, index]
        table[f"{cell_type}_signature_z"] = z_values[:, index]

    group_order = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    table["_group_order"] = table["dominant_cell_type"].map(group_order)
    table = table.sort_values(["_group_order", "dynamic_range", "name", "output_prefix"], ascending=[True, False, True, True])
    table = table.drop(columns="_group_order").reset_index(drop=True)
    table.insert(0, "rank", np.arange(1, len(table) + 1))
    table.insert(1, "signature_id", table["output_prefix"].astype(str))
    return table


def plot_signature_heatmap(
    table: pd.DataFrame,
    output_prefix: Path,
    title: str,
    *,
    top_n: int | None = None,
    label_rows: bool = False,
) -> None:
    plot_table = table.head(top_n).copy() if top_n is not None else table.copy()
    z_cols = [f"{cell_type}_signature_z" for cell_type in CELL_TYPES]
    values = plot_table[z_cols].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = float(np.nanpercentile(np.abs(finite), 99)) if finite.size else 1.0
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    n_rows = len(plot_table)
    row_height = PSEUDOBULK_LABELED_HEATMAP_ROW_HEIGHT_IN if label_rows else PSEUDOBULK_HEATMAP_ROW_HEIGHT_IN
    height = max(PSEUDOBULK_HEATMAP_MIN_BODY_HEIGHT_IN, row_height * n_rows) + PSEUDOBULK_HEATMAP_EXTRA_HEIGHT_IN
    width = 4.8 if not label_rows else 6.0
    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest", rasterized=True)
    ax.set_xticks(range(len(CELL_TYPES)), labels=CELL_TYPES)
    ax.tick_params(axis="x", labelrotation=0, labelsize=8.2, length=0)
    if label_rows:
        labels = plot_table["name"].astype(str)
        duplicate_names = labels.duplicated(keep=False)
        row_labels = labels.where(~duplicate_names, plot_table["signature_id"].astype(str))
        ax.set_yticks(range(n_rows), labels=row_labels)
        ax.tick_params(axis="y", labelsize=6.2, length=0)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{n_rows:,} motif signatures")

    boundaries = np.flatnonzero(plot_table["dominant_cell_type"].to_numpy()[1:] != plot_table["dominant_cell_type"].to_numpy()[:-1]) + 0.5
    for boundary in boundaries:
        ax.axhline(boundary, color="black", linewidth=0.45, alpha=0.55)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=9.6, fontweight="bold", pad=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("row-standardized differential score", fontsize=7.4)
    cbar.ax.tick_params(labelsize=7.0, length=2.2, width=0.6)
    fig.tight_layout()
    save_illustrator_svg(fig, output_prefix)
    plt.close(fig)


def balanced_top_signatures(table: pd.DataFrame, top_n: int) -> pd.DataFrame:
    per_group = max(1, int(np.ceil(top_n / len(CELL_TYPES))))
    groups = []
    for cell_type in CELL_TYPES:
        groups.append(table[table["dominant_cell_type"] == cell_type].head(per_group))
    top = pd.concat(groups, ignore_index=True).head(top_n)
    return top


def plot_all_signature_heatmaps(results: pd.DataFrame, output_prefix: Path, top_n: int) -> None:
    table = prepare_all_signature_heatmap(results)
    top_table = balanced_top_signatures(table, top_n)
    table.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)
    plot_signature_heatmap(
        table,
        output_prefix,
        "PBMC5k all pseudobulk footprint signatures",
    )
    plot_signature_heatmap(
        top_table,
        output_prefix.with_name(output_prefix.name.replace("_all_", "_top_")),
        f"PBMC5k top pseudobulk footprint signatures ({len(top_table)} total)",
        label_rows=True,
    )


def plot_volcano(results: pd.DataFrame, comparison: str, first: str, second: str, output_prefix: Path, markers: list[str]) -> None:
    plot_df = prepare_volcano_df(results, comparison)

    fig, ax = plt.subplots(figsize=(3.15, 3.15))
    for status, color, size, alpha in [
        ("not significant", COLORS["background"], 12, 0.45),
        ("higher in second", COLORS["down"], 18, 0.75),
        ("higher in first", COLORS["up"], 18, 0.75),
    ]:
        subset = plot_df[plot_df["status"] == status]
        ax.scatter(subset["change"], subset["neg_log10_p"], s=size, color=color, alpha=alpha, linewidths=0, label=status)

    ax.axvline(0, color="#374151", linewidth=0.8)
    ax.axhline(-np.log10(0.05), color="#6B7280", linewidth=0.8, linestyle="--")
    labels = directional_marker_rows(plot_df, first, second)
    if markers:
        requested = {marker.upper() for marker in markers}
        labels = labels[labels["name"].astype(str).str.upper().isin(requested)]
    annotate_marker_labels(ax, labels, fontsize=7.6)
    ax.set_xlabel("Differential footprint score")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{first} vs {second}", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_box_aspect(1)
    fig.tight_layout()
    save_illustrator_svg(fig, output_prefix)
    plt.close(fig)

    plot_df[["output_prefix", "name", "total_tfbs", "change", "pvalue", "qvalue", "status"]].to_csv(
        output_prefix.with_suffix(".tsv"),
        sep="\t",
        index=False,
    )


def plot_directional_pairwise_volcano(results: pd.DataFrame, output_prefix: Path) -> None:
    fig, axes = plt.subplots(
        1,
        len(PAIRWISE_COMPARISONS),
        figsize=(7.7, 3.25),
        sharey=True,
        gridspec_kw={"wspace": 0.28},
    )
    source_rows = []
    for ax, (comparison, first, second) in zip(axes, PAIRWISE_COMPARISONS, strict=True):
        plot_df = prepare_volcano_df(results, comparison)
        labels = directional_marker_rows(plot_df, first, second)
        for status, color, size, alpha in [
            ("not significant", COLORS["background"], 10, 0.42),
            ("higher in second", COLORS["down"], 14, 0.72),
            ("higher in first", COLORS["up"], 14, 0.72),
        ]:
            subset = plot_df[plot_df["status"] == status]
            ax.scatter(subset["change"], subset["neg_log10_p"], s=size, color=color, alpha=alpha, linewidths=0, label=status)

        ax.axvline(0, color="#374151", linewidth=0.8)
        ax.axhline(-np.log10(0.05), color="#6B7280", linewidth=0.8, linestyle="--")
        for row in annotate_marker_labels(ax, labels, fontsize=7.5):
            source_rows.append(
                {
                    "comparison": comparison,
                    "first": first,
                    "second": second,
                    "output_prefix": row["output_prefix"],
                    "name": row["name"],
                    "marker_group": row["marker_group"],
                    "total_tfbs": row["total_tfbs"],
                    "change": row["change"],
                    "pvalue": row["pvalue"],
                    "qvalue": row["qvalue"],
                    "status": row["status"],
                }
            )

        ax.set_xlabel("Differential footprint score")
        ax.set_title(f"{first} vs {second}", fontsize=8.8, fontweight="bold", pad=4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=7.3, length=2.5, width=0.7)
        ax.set_box_aspect(1)
    axes[0].set_ylabel("-log10(p-value)")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=7.2, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle("PBMC5k pseudobulk differential footprint scores", y=1.14, fontsize=9.8, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_illustrator_svg(fig, output_prefix)
    plt.close(fig)
    pd.DataFrame(source_rows).to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--aggregate-screen", required=True)
    parser.add_argument("--bindetect-results", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--markers", default="PAX5,CEBPB,TCF7")
    parser.add_argument("--top-heatmap-n", type=int, default=60, help="Number of high-contrast signatures to label in the companion heatmap.")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    annotations = read_table(Path(args.annotations))
    aggregate_screen = read_table(Path(args.aggregate_screen))
    results = read_table(Path(args.bindetect_results))
    markers = [marker.strip() for marker in args.markers.split(",") if marker.strip()]

    plot_celltype_umap(annotations, outdir / "pbmc5k_umap_broad_celltypes_scprinter_style")
    plot_marker_umap(annotations, aggregate_screen, outdir / "pbmc5k_marker_footprint_umap")
    for comparison, first, second in PAIRWISE_COMPARISONS:
        plot_volcano(results, comparison, first, second, outdir / f"pbmc5k_volcano_{comparison}", markers)
    plot_directional_pairwise_volcano(results, outdir / "pbmc5k_volcano_pairwise_directional_markers")
    plot_all_signature_heatmaps(results, outdir / "pbmc5k_all_footprint_signature_heatmap", args.top_heatmap_n)
    print(f"Wrote PBMC5k marker UMAP and volcano plots to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
