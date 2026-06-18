#!/usr/bin/env python
"""Plot per-cell PBMC5k marker footprint and motif-activity signatures."""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
import gzip
from pathlib import Path

import anndata as ad
import matplotlib.patheffects as patheffects
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import pysam
from sklearn.neighbors import NearestNeighbors


plt.rcParams.update(
    {
        "font.size": 8.5,
        "axes.titlesize": 9.2,
        "axes.labelsize": 8.8,
        "xtick.labelsize": 7.6,
        "ytick.labelsize": 7.6,
        "legend.fontsize": 7.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

MARKERS = ("PAX5", "CEBPB", "TCF7", "CEBPA", "SPIB", "ZBTB7B", "POU2F2")
SUMMARY_MARKER_ORDER = ("PAX5", "POU2F2", "SPIB", "CEBPB", "CEBPA", "TCF7", "ZBTB7B")
CELL_TYPES = ("B_cell", "Monocyte", "T_NK_cell")
MARKER_GROUPS = {
    "PAX5": "B_cell",
    "CEBPA": "Monocyte",
    "CEBPB": "Monocyte",
    "POU2F2": "B_cell",
    "SPIB": "B_cell",
    "TCF7": "T_NK_cell",
    "ZBTB7B": "T_NK_cell",
}
CELL_TYPE_COLORS = {
    "B_cell": "#3B82F6",
    "Monocyte": "#D97706",
    "T_NK_cell": "#059669",
}
REFERENCE_LABEL_POSITIONS = {
    "B_cell": (-5.15, 5.15),
    "Monocyte": (12.8, 8.4),
    "T_NK_cell": (0.2, 11.0),
}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_annotations(path: Path) -> pd.DataFrame:
    annotations = pd.read_csv(path, sep="\t")
    required = {"barcode", "cell_type", "umap_1", "umap_2"}
    missing = required.difference(annotations.columns)
    if missing:
        raise SystemExit(f"Annotation table is missing required columns: {', '.join(sorted(missing))}")
    return annotations


def site_path(site_dir: Path, tf: str) -> Path:
    for suffix in ("motif_hits", "motif_peaks"):
        path = site_dir / f"{tf}.{suffix}.bed"
        if path.exists():
            return path
    raise FileNotFoundError(f"No motif BED found for {tf} in {site_dir}")


def read_sites(site_dir: Path, markers: list[str], max_sites: int | None) -> dict[str, list[tuple[str, int]]]:
    sites: dict[str, list[tuple[str, int]]] = {}
    for tf in markers:
        rows = []
        with site_path(site_dir, tf).open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, start, end, *_ = line.rstrip("\n").split("\t")
                rows.append((chrom, (int(start) + int(end)) // 2))
                if max_sites is not None and len(rows) >= max_sites:
                    break
        if not rows:
            raise SystemExit(f"No sites found for {tf}")
        sites[tf] = rows
    return sites


def smooth_profile(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values
    if window % 2 == 0:
        window += 1
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(values, kernel, mode="same")


def protection_score(profile: np.ndarray, center_half_width: int, flank_inner: int, flank_outer: int) -> float:
    arr = smooth_profile(profile.astype(float), center_half_width * 2 + 1)
    mid = arr.shape[0] // 2
    center = arr[max(0, mid - center_half_width) : mid + center_half_width + 1]
    left = arr[max(0, mid - flank_outer) : max(0, mid - flank_inner)]
    right = arr[min(arr.shape[0], mid + flank_inner + 1) : min(arr.shape[0], mid + flank_outer + 1)]
    flanks = [x for x in (left, right) if x.size]
    if not flanks or center.size == 0:
        return float("nan")
    return float(np.nanmean(np.concatenate(flanks)) - np.nanmean(center))


def build_site_lookup(sites: dict[str, list[tuple[str, int]]]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for tf_index, tf in enumerate(sites):
        for chrom, center in sites[tf]:
            by_chrom.setdefault(chrom, []).append((center, tf_index))
    lookup = {}
    for chrom, rows in by_chrom.items():
        rows = sorted(rows)
        lookup[chrom] = (
            np.asarray([row[0] for row in rows], dtype=np.int64),
            np.asarray([row[1] for row in rows], dtype=np.int16),
        )
    return lookup


def merged_site_windows(sites: dict[str, list[tuple[str, int]]], flank: int) -> dict[str, list[tuple[int, int]]]:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for tf_sites in sites.values():
        for chrom, center in tf_sites:
            by_chrom.setdefault(chrom, []).append((max(0, center - flank), center + flank))
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in by_chrom.items():
        chrom_merged = []
        for start, end in sorted(intervals):
            if not chrom_merged or start > chrom_merged[-1][1]:
                chrom_merged.append([start, end])
            else:
                chrom_merged[-1][1] = max(chrom_merged[-1][1], end)
        merged[chrom] = [(int(start), int(end)) for start, end in chrom_merged]
    return merged


def add_cut(profiles: np.ndarray, cell_index: int, chrom: str, position: int, multiplicity: int, lookup: dict[str, tuple[np.ndarray, np.ndarray]], flank: int) -> int:
    if chrom not in lookup:
        return 0
    centers, tf_indices = lookup[chrom]
    left = bisect_left(centers, position - flank + 1)
    right = bisect_right(centers, position + flank)
    if right <= left:
        return 0
    offsets = position - centers[left:right] + flank
    tfs = tf_indices[left:right]
    for tf_index in np.unique(tfs):
        tf_offsets = offsets[tfs == tf_index]
        np.add.at(profiles[cell_index, tf_index], tf_offsets, multiplicity)
    return right - left


def ensure_tabix_index(fragments: Path, create_index: bool) -> bool:
    index_path = Path(str(fragments) + ".tbi")
    if index_path.exists():
        return True
    if not create_index:
        return False
    pysam.tabix_index(str(fragments), preset="bed", force=True, keep_original=True)
    return index_path.exists()


def iter_fragment_rows(fragments: Path, sites: dict[str, list[tuple[str, int]]], flank: int, create_index: bool):
    if ensure_tabix_index(fragments, create_index):
        tabix = pysam.TabixFile(str(fragments))
        try:
            windows = merged_site_windows(sites, flank)
            for chrom, intervals in windows.items():
                if chrom not in tabix.contigs:
                    continue
                for start, end in intervals:
                    yield from tabix.fetch(chrom, start, end)
        finally:
            tabix.close()
        return

    with open_text(fragments) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            yield line.rstrip("\n")


def count_fragment_profiles(
    fragments: Path,
    annotations: pd.DataFrame,
    sites: dict[str, list[tuple[str, int]]],
    flank: int,
    create_index: bool,
) -> np.ndarray:
    barcode_to_index = {barcode: idx for idx, barcode in enumerate(annotations["barcode"].astype(str))}
    profiles = np.zeros((annotations.shape[0], len(sites), 2 * flank), dtype=np.float32)
    lookup = build_site_lookup(sites)
    matched_cuts = 0
    for line in iter_fragment_rows(fragments, sites, flank, create_index):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        cell_index = barcode_to_index.get(fields[3])
        if cell_index is None:
            continue
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        if end <= start:
            continue
        multiplicity = int(fields[4]) if len(fields) > 4 and fields[4].isdigit() else 1
        matched_cuts += add_cut(profiles, cell_index, chrom, start, multiplicity, lookup, flank)
        matched_cuts += add_cut(profiles, cell_index, chrom, end - 1, multiplicity, lookup, flank)
    if matched_cuts == 0:
        raise SystemExit("No fragment cuts overlapped the requested motif windows.")
    return profiles


def knn_indices(annotations: pd.DataFrame, h5ad_path: Path, k: int) -> np.ndarray:
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        if "X_spectral" in adata.obsm:
            embedding = np.asarray(adata.obsm["X_spectral"])
            indexer = pd.Index(adata.obs_names.astype(str)).get_indexer(annotations["barcode"].astype(str))
            if (indexer < 0).any():
                raise SystemExit("Some annotation barcodes are missing from the h5ad obs_names.")
            embedding = embedding[indexer]
        else:
            embedding = annotations[["umap_1", "umap_2"]].to_numpy(dtype=float)
    finally:
        adata.file.close()
    n_neighbors = min(k, annotations.shape[0])
    model = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    model.fit(embedding)
    return model.kneighbors(return_distance=False)


def score_knn_profiles(
    profiles: np.ndarray,
    neighbors: np.ndarray,
    markers: list[str],
    center_half_width: int,
    flank_inner: int,
    flank_outer: int,
) -> pd.DataFrame:
    rows = []
    for cell_index in range(profiles.shape[0]):
        neighborhood = profiles[neighbors[cell_index]].sum(axis=0)
        for tf_index, tf in enumerate(markers):
            profile = neighborhood[tf_index]
            total = float(profile.sum())
            score = float("nan") if total <= 0 else protection_score(profile / total * 10_000.0, center_half_width, flank_inner, flank_outer)
            rows.append({"cell_index": cell_index, "tf": tf, "knn_footprint_score": score, "neighborhood_cut_count": total})
    return pd.DataFrame(rows)


def parse_marker_groups(text: str) -> dict[str, str]:
    groups = dict(MARKER_GROUPS)
    if not text:
        return groups
    for item in text.split(","):
        if not item.strip():
            continue
        if ":" not in item:
            raise SystemExit("--marker-groups entries must be TF:cell_type pairs.")
        tf, group = item.split(":", 1)
        groups[tf.strip()] = group.strip()
    return groups


def add_oriented_knn_score(scores: pd.DataFrame, marker_groups: dict[str, str]) -> pd.DataFrame:
    scores = scores.copy()
    scores["knn_footprint_oriented_z"] = np.nan
    rows = []
    for tf, subset in scores.groupby("tf", sort=False):
        expected_group = marker_groups.get(str(tf))
        raw = pd.to_numeric(subset["knn_footprint_score"], errors="coerce")
        sign = 1.0
        if expected_group:
            expected = raw[subset["cell_type"] == expected_group]
            other = raw[subset["cell_type"] != expected_group]
            if expected.notna().any() and other.notna().any() and float(expected.median()) < float(other.median()):
                sign = -1.0
        oriented = raw * sign
        finite = oriented[np.isfinite(oriented)]
        if finite.empty or float(finite.std(ddof=0)) == 0.0:
            z = oriented * np.nan
        else:
            z = (oriented - float(finite.mean())) / float(finite.std(ddof=0))
        scores.loc[subset.index, "knn_footprint_oriented_z"] = z
        rows.append(
            {
                "tf": tf,
                "expected_cell_type": expected_group or "",
                "orientation_sign": sign,
                "expected_median_raw": float(raw[subset["cell_type"] == expected_group].median()) if expected_group else np.nan,
                "other_median_raw": float(raw[subset["cell_type"] != expected_group].median()) if expected_group else np.nan,
            }
        )
    scores.attrs["orientation_summary"] = pd.DataFrame(rows)
    return scores


def bin_name_for_center(chrom: str, center: int, bin_size: int) -> str:
    start = (center // bin_size) * bin_size
    return f"{chrom}:{start}-{start + bin_size}"


def chromvar_like_scores(h5ad_path: Path, annotations: pd.DataFrame, sites: dict[str, list[tuple[str, int]]], markers: list[str], bin_size: int) -> pd.DataFrame:
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        obs_index = pd.Index(adata.obs_names.astype(str)).get_indexer(annotations["barcode"].astype(str))
        if (obs_index < 0).any():
            raise SystemExit("Some annotation barcodes are missing from the h5ad obs_names.")
        selected = np.asarray(adata.var["selected"], dtype=bool)
        selected_indices = np.flatnonzero(selected)
        selected_depth = np.asarray(adata.X[:, selected_indices].sum(axis=1)).ravel()[obs_index]
        rows = []
        for tf in markers:
            names = sorted({bin_name_for_center(chrom, center, bin_size) for chrom, center in sites[tf]})
            tf_indices = adata.var_names.get_indexer(names)
            tf_indices = np.asarray(sorted({int(idx) for idx in tf_indices if idx >= 0 and selected[idx]}), dtype=int)
            if tf_indices.size == 0:
                activity = np.full(annotations.shape[0], np.nan)
                counts = np.zeros(annotations.shape[0], dtype=float)
            else:
                counts = np.asarray(adata.X[:, tf_indices].sum(axis=1)).ravel()[obs_index]
                expected = selected_depth * (tf_indices.size / max(1, selected_indices.size))
                activity = (counts - expected) / np.sqrt(expected + 1.0)
                finite = np.isfinite(activity)
                if finite.any() and float(np.nanstd(activity[finite])) > 0:
                    activity = (activity - np.nanmean(activity[finite])) / np.nanstd(activity[finite])
            for cell_index, value in enumerate(activity):
                rows.append(
                    {
                        "cell_index": cell_index,
                        "tf": tf,
                        "chromvar_like_activity_z": float(value),
                        "motif_bin_count": int(tf_indices.size),
                        "motif_accessibility_count": float(counts[cell_index]),
                    }
                )
    finally:
        adata.file.close()
    return pd.DataFrame(rows)


def attach_annotations(scores: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    annot = annotations.reset_index(names="cell_index")
    return scores.merge(annot[["cell_index", "barcode", "cell_type", "snap_cell_type", "umap_1", "umap_2"]], on="cell_index", how="left")


def label_groups(
    ax: plt.Axes,
    annotations: pd.DataFrame,
    fontsize: int = 7,
    offsets: dict[str, tuple[float, float]] | None = None,
    positions: dict[str, tuple[float, float]] | None = None,
    arrows: bool = False,
) -> None:
    for label, group in annotations.groupby("cell_type", sort=True):
        target = (float(group["umap_1"].median()), float(group["umap_2"].median()))
        if positions and str(label) in positions:
            xytext = positions[str(label)]
        else:
            dx, dy = (offsets or {}).get(str(label), (0.0, 0.0))
            xytext = (target[0] + dx, target[1] + dy)
        text = ax.annotate(
            str(label),
            xy=target,
            xytext=xytext,
            fontsize=fontsize,
            weight="bold",
            ha="center",
            va="center",
            color="#111827",
            arrowprops=(
                {
                    "arrowstyle": "-",
                    "color": "#4B5563",
                    "linewidth": 0.6,
                    "shrinkA": 2,
                    "shrinkB": 2,
                }
                if arrows
                else None
            ),
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
            clip_on=False,
        )
        text.set_path_effects([patheffects.withStroke(linewidth=2.1, foreground="white")])


def robust_norm(values: pd.Series) -> TwoSlopeNorm:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    vmax = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 1.0
    vmax = max(vmax, 1e-6)
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def plot_score_grid(
    annotations: pd.DataFrame,
    score_tables: list[tuple[str, pd.DataFrame, str, str]],
    markers: list[str],
    output_prefix: Path,
) -> None:
    fig, axes = plt.subplots(len(score_tables), len(markers), figsize=(4.2 * len(markers), 3.55 * len(score_tables)), sharex=True, sharey=True)
    if len(score_tables) == 1:
        axes = np.asarray([axes])
    for row_index, (row_label, scores, value_column, colorbar_label) in enumerate(score_tables):
        for col_index, tf in enumerate(markers):
            ax = axes[row_index, col_index]
            subset = scores[scores["tf"] == tf].copy()
            subset = annotations[["barcode", "cell_type", "umap_1", "umap_2"]].merge(
                subset[["barcode", value_column]], on="barcode", how="left"
            )
            sc = ax.scatter(
                subset["umap_1"],
                subset["umap_2"],
                c=subset[value_column],
                cmap="RdBu_r",
                norm=robust_norm(subset[value_column]),
                s=2.0,
                alpha=0.92,
                linewidths=0,
                rasterized=True,
            )
            label_groups(ax, annotations, fontsize=7)
            if row_index == 0:
                ax.set_title(tf)
            if col_index == 0:
                ax.set_ylabel(f"{row_label}\nUMAP 2")
            ax.set_xlabel("UMAP 1")
            ax.spines[["top", "right"]].set_visible(False)
            cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
            cbar.set_label(colorbar_label, fontsize=7)
            cbar.ax.tick_params(labelsize=6)
    fig.suptitle("PBMC5k per-cell marker signatures", y=1.01)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), dpi=240, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_knn_with_reference(
    annotations: pd.DataFrame,
    knn_scores: pd.DataFrame,
    markers: list[str],
    output_prefix: Path,
) -> None:
    panels = [None] + markers
    ncols = 3 if len(panels) > 4 else len(panels)
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.8 * ncols, 2.55 * nrows), sharex=True, sharey=True)
    axes_flat = np.asarray(axes).reshape(-1)
    reference_ax = axes_flat[0]
    for cell_type in CELL_TYPES:
        subset = annotations[annotations["cell_type"] == cell_type]
        reference_ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=2.8,
            alpha=0.9,
            linewidths=0,
            color=CELL_TYPE_COLORS.get(cell_type),
            label=cell_type,
            rasterized=True,
        )
    label_groups(reference_ax, annotations, fontsize=8.6, positions=REFERENCE_LABEL_POSITIONS)
    reference_ax.set_title("Broad cell types", fontweight="bold", pad=4)
    reference_ax.set_xlabel("UMAP 1")
    reference_ax.set_ylabel("UMAP 2")
    reference_ax.legend(frameon=False, markerscale=2.8, fontsize=7.0, loc="center left", bbox_to_anchor=(1.02, 0.5))
    reference_ax.spines[["top", "right"]].set_visible(False)

    marker_axes = axes_flat[1 : 1 + len(markers)]
    for panel_index, (ax, tf) in enumerate(zip(marker_axes, markers, strict=True), start=1):
        subset = knn_scores[knn_scores["tf"] == tf].copy()
        subset = annotations[["barcode", "cell_type", "umap_1", "umap_2"]].merge(
            subset[["barcode", "knn_footprint_oriented_z"]],
            on="barcode",
            how="left",
        )
        sc = ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            c=subset["knn_footprint_oriented_z"],
            cmap="RdBu_r",
            norm=robust_norm(subset["knn_footprint_oriented_z"]),
            s=2.8,
            alpha=0.92,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{tf} footprint signature", fontweight="bold", pad=4)
        ax.set_xlabel("UMAP 1")
        if panel_index % ncols == 0:
            ax.set_ylabel("UMAP 2")
        ax.spines[["top", "right"]].set_visible(False)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label("footprint signature z-score", fontsize=7.2)
        cbar.ax.tick_params(labelsize=6.8, length=2.2, width=0.6)

    for ax in axes_flat[len(panels) :]:
        ax.set_visible(False)
    fig.suptitle("PBMC5k per-cell KNN footprint signature scores", y=1.02, fontsize=10.8, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(output_prefix.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=450, bbox_inches="tight")
    plt.close(fig)


def plot_single_cell_footprinting_summary(
    annotations: pd.DataFrame,
    knn_scores: pd.DataFrame,
    markers: list[str],
    representative_markers: list[str],
    output_prefix: Path,
) -> None:
    marker_by_cell = (
        knn_scores.groupby(["tf", "cell_type"], sort=False)["knn_footprint_oriented_z"]
        .mean()
        .unstack("cell_type")
        .reindex(index=markers, columns=CELL_TYPES)
    )
    heatmap = marker_by_cell.T
    fig = plt.figure(figsize=(10.8, 6.6))

    heat_left = 0.22
    heat_bottom = 0.55
    heat_width = 0.52
    heat_height = heat_width * fig.get_figwidth() * len(CELL_TYPES) / (len(markers) * fig.get_figheight())
    heat_ax = fig.add_axes([heat_left, heat_bottom, heat_width, heat_height])
    values = heatmap.to_numpy(dtype=float)
    vmax = float(np.nanpercentile(np.abs(values[np.isfinite(values)]), 98)) if np.isfinite(values).any() else 1.0
    vmax = max(vmax, 0.75)
    im = heat_ax.imshow(values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
    heat_ax.set_xticks(range(len(markers)), labels=markers)
    heat_ax.set_yticks(range(len(CELL_TYPES)), labels=CELL_TYPES)
    heat_ax.tick_params(axis="x", labelrotation=28, labelsize=8.0, length=0, pad=3)
    heat_ax.tick_params(axis="y", labelsize=8.2, length=0, pad=3)
    heat_ax.set_title("Marker footprint signatures by broad cell type", fontsize=10.0, fontweight="bold", pad=7)
    for label in heat_ax.get_xticklabels():
        label.set_ha("right")
        label.set_rotation_mode("anchor")
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            value = values[row, col]
            if np.isfinite(value):
                heat_ax.text(col, row, f"{value:.1f}", ha="center", va="center", fontsize=6.9, color="#111827")
    heat_ax.set_xticks(np.arange(-0.5, len(markers), 1), minor=True)
    heat_ax.set_yticks(np.arange(-0.5, len(CELL_TYPES), 1), minor=True)
    heat_ax.grid(which="minor", color="white", linewidth=1.0)
    heat_ax.tick_params(which="minor", bottom=False, left=False)
    for spine in heat_ax.spines.values():
        spine.set_visible(False)
    cax = fig.add_axes([heat_left + heat_width + 0.012, heat_bottom, 0.010, heat_height])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("mean signature z-score", fontsize=7.2)
    cbar.ax.tick_params(labelsize=6.8, length=2.2, width=0.6)
    fig.text(heat_left - 0.03, heat_bottom + heat_height + 0.012, "A", fontsize=13, fontweight="bold", va="top")

    umap_bottom = 0.08
    umap_width = 0.17
    umap_height = 0.34
    umap_lefts = [0.06, 0.29, 0.52, 0.75]
    umap_axes = [fig.add_axes([left, umap_bottom, umap_width, umap_height]) for left in umap_lefts]
    x_padding = float((annotations["umap_1"].max() - annotations["umap_1"].min()) * 0.03)
    y_padding = float((annotations["umap_2"].max() - annotations["umap_2"].min()) * 0.03)
    xlim = (float(annotations["umap_1"].min() - x_padding), float(annotations["umap_1"].max() + x_padding))
    ylim = (float(annotations["umap_2"].min() - y_padding), float(annotations["umap_2"].max() + y_padding))
    reference_ax = umap_axes[0]
    for cell_type in CELL_TYPES:
        subset = annotations[annotations["cell_type"] == cell_type]
        reference_ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=2.8,
            alpha=0.9,
            linewidths=0,
            color=CELL_TYPE_COLORS.get(cell_type),
            rasterized=True,
        )
    label_groups(reference_ax, annotations, fontsize=8.4, positions=REFERENCE_LABEL_POSITIONS)
    reference_ax.set_title("Cell type\nannotation", fontsize=9.0, fontweight="bold", pad=5)
    reference_ax.set_xlabel("UMAP 1")
    reference_ax.set_ylabel("UMAP 2")
    reference_ax.set_xlim(xlim)
    reference_ax.set_ylim(ylim)
    reference_ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.032, umap_bottom + umap_height + 0.025, "B", fontsize=13, fontweight="bold", va="top")

    marker_values = knn_scores[knn_scores["tf"].isin(representative_markers)]["knn_footprint_oriented_z"]
    marker_arr = pd.to_numeric(marker_values, errors="coerce").to_numpy(dtype=float)
    finite = marker_arr[np.isfinite(marker_arr)]
    marker_vmax = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 2.0
    marker_vmax = max(marker_vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-marker_vmax, vcenter=0.0, vmax=marker_vmax)
    last_scatter = None
    for ax, tf in zip(umap_axes[1:], representative_markers, strict=True):
        subset = knn_scores[knn_scores["tf"] == tf].copy()
        subset = annotations[["barcode", "cell_type", "umap_1", "umap_2"]].merge(
            subset[["barcode", "knn_footprint_oriented_z"]],
            on="barcode",
            how="left",
        )
        last_scatter = ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            c=subset["knn_footprint_oriented_z"],
            cmap="RdBu_r",
            norm=norm,
            s=2.8,
            alpha=0.92,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{tf}\nfootprint signature", fontsize=9.0, fontweight="bold", pad=5)
        ax.set_xlabel("UMAP 1")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelleft=False)

    if last_scatter is not None:
        cax = fig.add_axes([0.935, umap_bottom, 0.014, umap_height])
        cbar = fig.colorbar(last_scatter, cax=cax)
        cbar.set_label("footprint signature z-score", fontsize=7.2)
        cbar.ax.tick_params(labelsize=6.8, length=2.2, width=0.6)

    for ax in umap_axes:
        ax.tick_params(axis="both", labelsize=8.0)

    fig.suptitle("Single-cell footprinting", y=0.995, fontsize=13.2, fontweight="bold")
    fig.savefig(output_prefix.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=450, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--fragments", required=True)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--tf-site-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--markers", default=",".join(MARKERS))
    parser.add_argument("--max-sites-per-tf", type=int, default=1500)
    parser.add_argument("--knn", type=int, default=75)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--center-half-width", type=int, default=10)
    parser.add_argument("--flank-inner", type=int, default=25)
    parser.add_argument("--flank-outer", type=int, default=100)
    parser.add_argument("--bin-size", type=int, default=500)
    parser.add_argument(
        "--marker-groups",
        default="PAX5:B_cell,CEBPA:Monocyte,CEBPB:Monocyte,POU2F2:B_cell,SPIB:B_cell,TCF7:T_NK_cell,ZBTB7B:T_NK_cell",
        help="Comma-separated TF:cell_type pairs used to orient KNN marker scores for UMAP review.",
    )
    parser.add_argument("--no-create-fragment-index", action="store_true", help="Do not create a tabix index for the fragment file when it is missing.")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    markers = [marker.strip() for marker in args.markers.split(",") if marker.strip()]
    annotations = read_annotations(Path(args.annotations))
    sites = read_sites(Path(args.tf_site_dir), markers, args.max_sites_per_tf)

    profiles = count_fragment_profiles(Path(args.fragments), annotations, sites, args.flank, create_index=not args.no_create_fragment_index)
    neighbors = knn_indices(annotations, Path(args.h5ad), args.knn)
    knn_scores = attach_annotations(
        score_knn_profiles(profiles, neighbors, markers, args.center_half_width, args.flank_inner, args.flank_outer),
        annotations,
    )
    knn_scores = add_oriented_knn_score(knn_scores, parse_marker_groups(args.marker_groups))
    chromvar_scores = attach_annotations(chromvar_like_scores(Path(args.h5ad), annotations, sites, markers, args.bin_size), annotations)

    knn_scores.to_csv(outdir / "knn_footprint_signature_scores.tsv", sep="\t", index=False)
    if "orientation_summary" in knn_scores.attrs:
        knn_scores.attrs["orientation_summary"].to_csv(outdir / "knn_footprint_orientation_summary.tsv", sep="\t", index=False)
    chromvar_scores.to_csv(outdir / "chromvar_like_motif_activity_scores.tsv", sep="\t", index=False)

    plot_knn_with_reference(annotations, knn_scores, markers, outdir / "pbmc5k_knn_footprint_signature_umap")
    summary_markers = [marker for marker in SUMMARY_MARKER_ORDER if marker in markers]
    summary_markers.extend(marker for marker in markers if marker not in summary_markers)
    representative_markers = [marker for marker in ("PAX5", "CEBPB", "TCF7") if marker in markers]
    plot_single_cell_footprinting_summary(
        annotations,
        knn_scores,
        summary_markers,
        representative_markers,
        outdir / "pbmc5k_single_cell_footprinting_summary",
    )
    plot_score_grid(
        annotations,
        [("ChromVAR-like", chromvar_scores, "chromvar_like_activity_z", "activity z")],
        markers,
        outdir / "pbmc5k_chromvar_like_motif_activity_umap",
    )
    plot_score_grid(
        annotations,
        [
            ("KNN footprint", knn_scores, "knn_footprint_score", "protection score"),
            ("KNN oriented", knn_scores, "knn_footprint_oriented_z", "oriented z"),
            ("ChromVAR-like", chromvar_scores, "chromvar_like_activity_z", "activity z"),
        ],
        markers,
        outdir / "pbmc5k_per_cell_signature_review",
    )
    print(f"Wrote per-cell PBMC5k signature plots and score tables to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
