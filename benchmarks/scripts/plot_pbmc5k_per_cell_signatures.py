#!/usr/bin/env python
"""Plot per-cell PBMC5k marker footprint and motif-activity signatures."""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
import gzip
from pathlib import Path

import anndata as ad
import matplotlib.patches as mpatches
import matplotlib.patheffects as patheffects
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
from matplotlib.colors import ListedColormap, TwoSlopeNorm
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
        "font.family": "Arial",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "svg.fonttype": "none",
    }
)

MARKERS = ("PAX5", "CEBPB", "TCF7", "CEBPA", "SPIB", "ZBTB7B", "POU2F2")
SUMMARY_MARKER_ORDER = ("PAX5", "POU2F2", "SPIB", "CEBPB", "CEBPA", "TCF7", "ZBTB7B")
UMAP_MARKER_LABELS = ("PAX5", "POU2F2", "SPIB", "CEBPB", "CEBPA", "TCF7", "ZBTB7B")
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
EXPECTED_MARKER_GROUPS = {
    "B_cell": ("PAX5", "POU2F2", "SPIB"),
    "Monocyte": ("CEBPB", "CEBPA"),
    "T_NK_cell": ("TCF7", "ZBTB7B"),
}
CELL_TYPE_COLORS = {
    "B_cell": "#3B82F6",
    "Monocyte": "#D97706",
    "T_NK_cell": "#059669",
}
MARKER_HEATMAP_ROW_HEIGHT_IN = 0.34
MOTIF_HEATMAP_ROW_HEIGHT_IN = 0.035
HEATMAP_MIN_BODY_HEIGHT_IN = 3.2
HEATMAP_TOP_MARGIN_IN = 1.65
HEATMAP_BOTTOM_MARGIN_IN = 0.55
REFERENCE_LABEL_POSITIONS = {
    "B_cell": (-5.15, 5.15),
    "Monocyte": (12.8, 8.4),
    "T_NK_cell": (0.2, 11.0),
}


def save_illustrator_svg(fig: plt.Figure, output_prefix: Path, *, tight: bool = True) -> None:
    for text in fig.findobj(match=mtext.Text):
        text.set_fontfamily("Arial")
        text.set_fontsize(9)
        text.set_fontweight("bold")
    save_kwargs = {"bbox_inches": "tight"} if tight else {}
    fig.savefig(output_prefix.with_suffix(".svg"), **save_kwargs)


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


def batched(items: list[str], batch_size: int):
    for start in range(0, len(items), batch_size):
        yield items[start : start + batch_size]


def read_bindetect_motif_table(results_path: Path, bindetect_dir: Path, max_motifs: int | None) -> pd.DataFrame:
    results = pd.read_csv(results_path, sep="\t")
    required = {"output_prefix", "name", "motif_id", "total_tfbs"}
    missing = required.difference(results.columns)
    if missing:
        raise SystemExit(f"BINDetect results table is missing required columns: {', '.join(sorted(missing))}")
    rows = []
    for _, row in results.iterrows():
        motif_id = str(row["output_prefix"])
        bed = bindetect_dir / motif_id / "beds" / f"{motif_id}_all.bed"
        if not bed.exists():
            continue
        rows.append(
            {
                "motif_id": motif_id,
                "tf_name": str(row["name"]),
                "jaspar_id": str(row["motif_id"]),
                "total_tfbs": int(row["total_tfbs"]),
                "bed_path": str(bed),
            }
        )
        if max_motifs is not None and len(rows) >= max_motifs:
            break
    if not rows:
        raise SystemExit(f"No BINDetect motif BED files found under {bindetect_dir}")
    return pd.DataFrame(rows)


def read_motif_bed_sites(path: Path, max_sites: int | None) -> list[tuple[str, int]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[:3]
            score = float(fields[4]) if len(fields) > 4 else 0.0
            rows.append((chrom, (int(start) + int(end)) // 2, score))
    if max_sites is not None and max_sites > 0 and len(rows) > max_sites:
        rows = sorted(rows, key=lambda row: row[2], reverse=True)[:max_sites]
    return [(chrom, center) for chrom, center, _ in rows]


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
    save_illustrator_svg(fig, output_prefix)
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
    save_illustrator_svg(fig, output_prefix)
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
    save_illustrator_svg(fig, output_prefix, tight=False)
    plt.close(fig)


def ordered_cell_annotations(annotations: pd.DataFrame) -> pd.DataFrame:
    ordered = annotations.copy()
    ordered["cell_type_order"] = pd.Categorical(ordered["cell_type"], categories=CELL_TYPES, ordered=True)
    return ordered.sort_values(["cell_type_order", "umap_1", "umap_2", "barcode"]).drop(columns="cell_type_order")


def plot_per_cell_signature_heatmap(
    annotations: pd.DataFrame,
    knn_scores: pd.DataFrame,
    markers: list[str],
    output_prefix: Path,
) -> None:
    ordered = ordered_cell_annotations(annotations)
    cell_order = ordered["barcode"].astype(str).tolist()
    matrix = (
        knn_scores.pivot_table(index="tf", columns="barcode", values="knn_footprint_oriented_z", aggfunc="mean")
        .reindex(index=markers, columns=cell_order)
    )
    matrix.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index_label="tf")

    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = float(np.nanpercentile(np.abs(finite), 99)) if finite.size else 1.0
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    heat_body_height = max(HEATMAP_MIN_BODY_HEIGHT_IN, MARKER_HEATMAP_ROW_HEIGHT_IN * len(markers))
    fig = plt.figure(figsize=(10.8, heat_body_height + HEATMAP_TOP_MARGIN_IN + HEATMAP_BOTTOM_MARGIN_IN))
    grid = fig.add_gridspec(2, 2, height_ratios=[0.15, 1.0], width_ratios=[1.0, 0.025], hspace=0.05, wspace=0.025)
    strip_ax = fig.add_subplot(grid[0, 0])
    heat_ax = fig.add_subplot(grid[1, 0])
    cax = fig.add_subplot(grid[1, 1])

    cell_type_to_code = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    strip = ordered["cell_type"].map(cell_type_to_code).to_numpy(dtype=int)[np.newaxis, :]
    strip_cmap = ListedColormap([CELL_TYPE_COLORS[cell_type] for cell_type in CELL_TYPES])
    strip_ax.imshow(strip, aspect="auto", interpolation="nearest", cmap=strip_cmap, vmin=-0.5, vmax=len(CELL_TYPES) - 0.5)
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])
    strip_ax.set_title("Cells ordered by broad cell type and UMAP position", fontsize=8.8, fontweight="bold", pad=3)
    for spine in strip_ax.spines.values():
        spine.set_visible(False)

    im = heat_ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest", rasterized=True)
    heat_ax.set_yticks(range(len(markers)), labels=markers)
    heat_ax.set_xticks([])
    heat_ax.set_xlabel(f"{len(cell_order):,} cells")
    heat_ax.set_ylabel("Footprint signature")
    heat_ax.tick_params(axis="y", labelsize=8.2, length=0)
    boundaries = np.flatnonzero(ordered["cell_type"].to_numpy()[1:] != ordered["cell_type"].to_numpy()[:-1]) + 0.5
    for boundary in boundaries:
        strip_ax.axvline(boundary, color="white", linewidth=1.0)
        heat_ax.axvline(boundary, color="black", linewidth=0.35, alpha=0.5)
    for spine in heat_ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("footprint signature z-score", fontsize=7.4)
    cbar.ax.tick_params(labelsize=7.0, length=2.2, width=0.6)

    handles = [mpatches.Patch(color=CELL_TYPE_COLORS[cell_type], label=cell_type) for cell_type in CELL_TYPES]
    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=len(CELL_TYPES), bbox_to_anchor=(0.5, 1.02), fontsize=7.6)
    fig.suptitle("PBMC5k per-cell KNN footprint signature heatmap", y=1.08, fontsize=11.0, fontweight="bold")
    save_illustrator_svg(fig, output_prefix)
    plt.close(fig)


def row_zscore(values: np.ndarray) -> np.ndarray:
    means = np.nanmean(values, axis=1, keepdims=True)
    sds = np.nanstd(values, axis=1, keepdims=True)
    return np.divide(values - means, sds, out=np.zeros_like(values), where=sds > 0)


def add_motif_specificity(metadata: pd.DataFrame) -> pd.DataFrame:
    metadata = metadata.copy()
    mean_columns = {cell_type: f"{cell_type}_mean_z" for cell_type in CELL_TYPES}
    specificity_scores = []
    for _, row in metadata.iterrows():
        group_means = {cell_type: float(row[column]) for cell_type, column in mean_columns.items()}
        dominant = str(row.get("dominant_cell_type") or max(group_means, key=group_means.get))
        other_means = [value for cell_type, value in group_means.items() if cell_type != dominant]
        specificity_scores.append(float(group_means[dominant] - max(other_means)))
    metadata["cell_type_specificity"] = specificity_scores
    return metadata


def add_repelled_row_labels(
    ax: plt.Axes,
    metadata: pd.DataFrame,
    label_tfs: list[str] | tuple[str, ...],
    *,
    min_gap_rows: float = 4.0,
) -> None:
    label_set = {str(tf).upper() for tf in label_tfs}
    label_rows = [
        (row_index, str(row["tf_name"]))
        for row_index, row in metadata.reset_index(drop=True).iterrows()
        if str(row.get("tf_name", "")).upper() in label_set
    ]
    if not label_rows:
        return

    adjusted = []
    previous_y = -float("inf")
    for row_index, label in label_rows:
        label_y = max(float(row_index), previous_y + min_gap_rows)
        adjusted.append((row_index, label_y, label))
        previous_y = label_y

    overflow = adjusted[-1][1] - (len(metadata) - 1)
    if overflow > 0:
        adjusted = [(row_index, label_y - overflow, label) for row_index, label_y, label in adjusted]
        for idx in range(len(adjusted) - 2, -1, -1):
            next_y = adjusted[idx + 1][1]
            row_index, label_y, label = adjusted[idx]
            adjusted[idx] = (row_index, min(label_y, next_y - min_gap_rows), label)

    transform = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    for row_index, label_y, label in adjusted:
        ax.annotate(
            label,
            xy=(-0.002, row_index),
            xytext=(-0.018, label_y),
            xycoords=transform,
            textcoords=transform,
            ha="right",
            va="center",
            fontsize=9,
            fontweight="bold",
            fontfamily="Arial",
            arrowprops={"arrowstyle": "-", "linewidth": 0.35, "color": "#374151", "shrinkA": 0, "shrinkB": 0},
            clip_on=False,
        )


def select_top_motif_signatures(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    per_cell_type: int,
    min_specificity: float,
    marker_tfs: list[str] | tuple[str, ...] = UMAP_MARKER_LABELS,
    correlation_threshold: float = 0.96,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = add_motif_specificity(metadata)
    marker_set = {str(tf).upper() for tf in marker_tfs}
    matrix_values = matrix.to_numpy(dtype=float)
    row_lookup = {motif_id: idx for idx, motif_id in enumerate(matrix.index.astype(str))}

    def row_corr(first: np.ndarray, second: np.ndarray) -> float:
        mask = np.isfinite(first) & np.isfinite(second)
        if mask.sum() < 3:
            return 0.0
        x = first[mask]
        y = second[mask]
        x_sd = float(np.std(x))
        y_sd = float(np.std(y))
        if x_sd == 0.0 or y_sd == 0.0:
            return 0.0
        return float(np.corrcoef(x, y)[0, 1])

    selected = []
    selected_vectors: list[np.ndarray] = []
    for cell_type in CELL_TYPES:
        subset = metadata[
            (metadata["dominant_cell_type"] == cell_type)
            & (metadata["cell_type_specificity"] >= float(min_specificity))
        ].sort_values(["cell_type_specificity", "dynamic_range", "motif_id"], ascending=[False, False, True])
        marker_subset = metadata[
            (metadata["dominant_cell_type"] == cell_type)
            & (metadata["tf_name"].astype(str).str.upper().isin(marker_set))
        ].sort_values(["cell_type_specificity", "dynamic_range", "motif_id"], ascending=[False, False, True])
        candidates = pd.concat([marker_subset, subset], axis=0).drop_duplicates("motif_id", keep="first")
        kept_rows = []
        kept_vectors = list(selected_vectors)
        for _, row in candidates.iterrows():
            motif_id = str(row["motif_id"])
            tf_name = str(row.get("tf_name", "")).upper()
            vector = matrix_values[row_lookup[motif_id]]
            force_keep = tf_name in marker_set
            redundant = False
            if not force_keep:
                redundant = any(abs(row_corr(vector, existing)) >= correlation_threshold for existing in kept_vectors)
            if redundant:
                continue
            kept_rows.append(row)
            kept_vectors.append(vector)
            if len(kept_rows) >= per_cell_type:
                break
        selected_vectors.extend(matrix_values[row_lookup[str(row["motif_id"])]] for row in kept_rows)
        selected.append(pd.DataFrame(kept_rows).head(per_cell_type))
    top_metadata = pd.concat(selected, axis=0) if selected else metadata.head(0)
    if top_metadata.empty:
        raise SystemExit("No top motif signatures passed the specificity threshold.")
    group_order = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    top_metadata = top_metadata.copy()
    top_metadata["_group_order"] = top_metadata["dominant_cell_type"].map(group_order)
    top_metadata = top_metadata.sort_values(
        ["_group_order", "cell_type_specificity", "dynamic_range", "motif_id"],
        ascending=[True, False, False, True],
    ).drop(columns="_group_order")
    return matrix.loc[top_metadata["motif_id"]], top_metadata.reset_index(drop=True)


def plot_all_motif_per_cell_heatmap(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    ordered_annotations: pd.DataFrame,
    output_prefix: Path,
    *,
    title: str = "PBMC5k all-motif per-cell KNN footprint signatures",
    label_tfs: list[str] | tuple[str, ...] | None = None,
    show_axis_counts: bool = True,
) -> None:
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = float(np.nanpercentile(np.abs(finite), 99)) if finite.size else 1.0
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    heat_body_height = max(HEATMAP_MIN_BODY_HEIGHT_IN, MOTIF_HEATMAP_ROW_HEIGHT_IN * matrix.shape[0])
    fig = plt.figure(figsize=(12.2, heat_body_height + HEATMAP_TOP_MARGIN_IN + HEATMAP_BOTTOM_MARGIN_IN))
    strip_height = 0.26
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[strip_height, heat_body_height],
        width_ratios=[1.0, 0.022],
        hspace=0.025,
        wspace=0.025,
    )
    strip_ax = fig.add_subplot(grid[0, 0])
    heat_ax = fig.add_subplot(grid[1, 0])
    cax = fig.add_subplot(grid[1, 1])

    cell_type_to_code = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    strip = ordered_annotations["cell_type"].map(cell_type_to_code).to_numpy(dtype=int)[np.newaxis, :]
    strip_cmap = ListedColormap([CELL_TYPE_COLORS[cell_type] for cell_type in CELL_TYPES])
    strip_ax.imshow(strip, aspect="auto", interpolation="nearest", cmap=strip_cmap, vmin=-0.5, vmax=len(CELL_TYPES) - 0.5)
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])
    strip_ax.set_title("Cells ordered by broad cell type and UMAP position", fontsize=9.0, fontweight="bold", pad=4)
    for spine in strip_ax.spines.values():
        spine.set_visible(False)

    im = heat_ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest", rasterized=True)
    heat_ax.set_xticks([])
    if label_tfs:
        heat_ax.set_yticks([])
        add_repelled_row_labels(heat_ax, metadata, label_tfs)
    else:
        heat_ax.set_yticks([])
    heat_ax.set_xlabel(f"{matrix.shape[1]:,} cells" if show_axis_counts else "")
    heat_ax.set_ylabel(f"{matrix.shape[0]:,} motif signatures" if show_axis_counts else "")
    cell_boundaries = np.flatnonzero(ordered_annotations["cell_type"].to_numpy()[1:] != ordered_annotations["cell_type"].to_numpy()[:-1]) + 0.5
    for boundary in cell_boundaries:
        strip_ax.axvline(boundary, color="white", linewidth=1.0)
        heat_ax.axvline(boundary, color="black", linewidth=0.35, alpha=0.55)
    motif_boundaries = np.flatnonzero(metadata["dominant_cell_type"].to_numpy()[1:] != metadata["dominant_cell_type"].to_numpy()[:-1]) + 0.5
    for boundary in motif_boundaries:
        heat_ax.axhline(boundary, color="black", linewidth=0.45, alpha=0.6)
    for spine in heat_ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("row-standardized KNN footprint score", fontsize=7.6)
    cbar.ax.tick_params(labelsize=7.0, length=2.2, width=0.6)
    handles = [mpatches.Patch(color=CELL_TYPE_COLORS[cell_type], label=cell_type) for cell_type in CELL_TYPES]
    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=len(CELL_TYPES), bbox_to_anchor=(0.5, 1.01), fontsize=7.8)
    fig.suptitle(title, y=1.045, fontsize=11.0, fontweight="bold")
    save_illustrator_svg(fig, output_prefix)
    plt.close(fig)


def read_all_motif_score_table(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = pd.read_csv(path, sep="\t")
    required = {"motif_id", "dominant_cell_type", "dynamic_range", *[f"{cell_type}_mean_z" for cell_type in CELL_TYPES]}
    missing = required.difference(table.columns)
    if missing:
        raise SystemExit(f"All-motif score table is missing required columns: {', '.join(sorted(missing))}")
    first_cell_col = next((idx for idx, column in enumerate(table.columns) if str(column).endswith("-1")), None)
    if first_cell_col is None:
        raise SystemExit("Could not find per-cell columns in all-motif score table.")
    metadata = table.iloc[:, :first_cell_col].copy()
    matrix = table.iloc[:, first_cell_col:].copy()
    matrix.index = metadata["motif_id"].astype(str)
    return matrix, metadata


def expected_group_for_marker(tf: str) -> str:
    marker = str(tf).upper()
    for cell_type, markers in EXPECTED_MARKER_GROUPS.items():
        if marker in {item.upper() for item in markers}:
            return cell_type
    raise SystemExit(f"No expected cell type configured for marker TF {tf}.")


def marker_signature_rows(
    marker_scores: pd.DataFrame,
    ordered: pd.DataFrame,
    marker_tfs: list[str] | tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"tf", "barcode", "knn_footprint_oriented_z"}
    missing = required.difference(marker_scores.columns)
    if missing:
        raise SystemExit(f"Marker score table is missing required columns: {', '.join(sorted(missing))}")

    marker_set = [str(tf) for tf in marker_tfs]
    score_tfs = {str(tf).upper() for tf in marker_scores["tf"].dropna().unique()}
    missing_markers = [tf for tf in marker_set if tf.upper() not in score_tfs]
    if missing_markers:
        raise SystemExit(f"Marker score table is missing selected TFs: {', '.join(missing_markers)}")

    cell_order = ordered["barcode"].astype(str).tolist()
    matrix = (
        marker_scores[marker_scores["tf"].astype(str).str.upper().isin({tf.upper() for tf in marker_set})]
        .pivot_table(index="tf", columns="barcode", values="knn_footprint_oriented_z", aggfunc="mean")
        .reindex(index=marker_set, columns=cell_order)
    )
    metadata_rows = []
    for tf in marker_set:
        values = matrix.loc[tf].to_numpy(dtype=float)
        group_means = {
            cell_type: float(np.nanmean(values[ordered["cell_type"].to_numpy() == cell_type]))
            for cell_type in CELL_TYPES
        }
        expected_group = expected_group_for_marker(tf)
        other_means = [value for cell_type, value in group_means.items() if cell_type != expected_group]
        metadata_rows.append(
            {
                "motif_id": f"{tf}_marker_oriented_knn",
                "tf_name": tf,
                "name": tf,
                "dominant_cell_type": expected_group,
                "dynamic_range": float(np.nanmax(values) - np.nanmin(values)),
                "cell_type_specificity": float(group_means[expected_group] - max(other_means)),
                "source": "marker_oriented_knn",
                **{f"{cell_type}_mean_z": group_means[cell_type] for cell_type in CELL_TYPES},
            }
        )
    matrix.index = [row["motif_id"] for row in metadata_rows]
    return matrix, pd.DataFrame(metadata_rows)


def replace_selected_marker_rows(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    ordered: pd.DataFrame,
    marker_scores: pd.DataFrame | None,
    marker_tfs: list[str] | tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if marker_scores is None:
        return matrix, metadata
    marker_set = {str(tf).upper() for tf in marker_tfs}
    keep_mask = ~metadata["tf_name"].astype(str).str.upper().isin(marker_set)
    base_metadata = metadata.loc[keep_mask].copy()
    base_matrix = matrix.loc[base_metadata["motif_id"].astype(str)]
    marker_matrix, marker_metadata = marker_signature_rows(marker_scores, ordered, marker_tfs)
    merged_matrix = pd.concat([base_matrix, marker_matrix], axis=0)
    merged_metadata = pd.concat([base_metadata, marker_metadata], axis=0, ignore_index=True)
    return merged_matrix, merged_metadata


def top_motif_signature_table(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    ordered: pd.DataFrame,
    per_cell_type: int,
    min_specificity: float,
    marker_tfs: list[str] | tuple[str, ...] = UMAP_MARKER_LABELS,
    correlation_threshold: float = 0.96,
    marker_scores: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix, metadata = replace_selected_marker_rows(matrix, metadata, ordered, marker_scores, marker_tfs)
    return select_top_motif_signatures(
        matrix,
        metadata,
        per_cell_type,
        min_specificity,
        marker_tfs=marker_tfs,
        correlation_threshold=correlation_threshold,
    )


def write_top_motif_heatmap(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    ordered: pd.DataFrame,
    output_prefix: Path,
    per_cell_type: int,
    min_specificity: float,
    marker_tfs: list[str] | tuple[str, ...] = UMAP_MARKER_LABELS,
    correlation_threshold: float = 0.96,
    marker_scores: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    top_matrix, top_metadata = top_motif_signature_table(
        matrix,
        metadata,
        ordered,
        per_cell_type,
        min_specificity,
        marker_tfs=marker_tfs,
        correlation_threshold=correlation_threshold,
        marker_scores=marker_scores,
    )
    top_output = pd.concat([top_metadata.reset_index(drop=True), top_matrix.reset_index(drop=True)], axis=1)
    top_output.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)
    plot_all_motif_per_cell_heatmap(
        top_matrix,
        top_metadata,
        ordered,
        output_prefix,
        title=f"PBMC5k top cell-type-specific per-cell KNN footprint signatures",
        label_tfs=marker_tfs,
        show_axis_counts=False,
    )
    return top_matrix, top_metadata


def draw_celltype_strip(ax: plt.Axes, ordered_annotations: pd.DataFrame) -> None:
    cell_type_to_code = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    strip = ordered_annotations["cell_type"].map(cell_type_to_code).to_numpy(dtype=int)[np.newaxis, :]
    strip_cmap = ListedColormap([CELL_TYPE_COLORS[cell_type] for cell_type in CELL_TYPES])
    ax.imshow(strip, aspect="auto", interpolation="nearest", cmap=strip_cmap, vmin=-0.5, vmax=len(CELL_TYPES) - 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Cells ordered by broad cell type and UMAP position", fontsize=9, fontweight="bold", pad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_top_motif_heatmap_panel(
    fig: plt.Figure,
    strip_ax: plt.Axes,
    heat_ax: plt.Axes,
    cax: plt.Axes,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    ordered_annotations: pd.DataFrame,
    label_tfs: list[str] | tuple[str, ...],
) -> None:
    values = matrix.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    vmax = float(np.nanpercentile(np.abs(finite), 99)) if finite.size else 1.0
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    draw_celltype_strip(strip_ax, ordered_annotations)
    im = heat_ax.imshow(values, cmap="RdBu_r", norm=norm, aspect="auto", interpolation="nearest", rasterized=True)
    heat_ax.set_xticks([])
    heat_ax.set_yticks([])
    add_repelled_row_labels(heat_ax, metadata, label_tfs)
    cell_boundaries = np.flatnonzero(ordered_annotations["cell_type"].to_numpy()[1:] != ordered_annotations["cell_type"].to_numpy()[:-1]) + 0.5
    for boundary in cell_boundaries:
        strip_ax.axvline(boundary, color="white", linewidth=1.0)
        heat_ax.axvline(boundary, color="black", linewidth=0.35, alpha=0.55)
    motif_boundaries = np.flatnonzero(metadata["dominant_cell_type"].to_numpy()[1:] != metadata["dominant_cell_type"].to_numpy()[:-1]) + 0.5
    for boundary in motif_boundaries:
        heat_ax.axhline(boundary, color="black", linewidth=0.45, alpha=0.6)
    for spine in heat_ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("row-standardized KNN footprint score", fontsize=9)
    cbar.ax.tick_params(labelsize=9, length=2.2, width=0.6)


def draw_knn_umap_panel(
    fig: plt.Figure,
    axes: list[plt.Axes],
    annotations: pd.DataFrame,
    knn_scores: pd.DataFrame,
    markers: list[str],
    cax: plt.Axes | None = None,
) -> None:
    x_padding = float((annotations["umap_1"].max() - annotations["umap_1"].min()) * 0.03)
    y_padding = float((annotations["umap_2"].max() - annotations["umap_2"].min()) * 0.03)
    xlim = (float(annotations["umap_1"].min() - x_padding), float(annotations["umap_1"].max() + x_padding))
    ylim = (float(annotations["umap_2"].min() - y_padding), float(annotations["umap_2"].max() + y_padding))

    reference_ax = axes[0]
    for cell_type in CELL_TYPES:
        subset = annotations[annotations["cell_type"] == cell_type]
        reference_ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=1.8,
            alpha=0.9,
            linewidths=0,
            color=CELL_TYPE_COLORS.get(cell_type),
            label=cell_type,
            rasterized=True,
        )
    label_groups(reference_ax, annotations, fontsize=9, positions=REFERENCE_LABEL_POSITIONS)
    reference_ax.set_title("Broad cell types", fontweight="bold", pad=4)
    reference_ax.set_xlim(xlim)
    reference_ax.set_ylim(ylim)
    reference_ax.set_xticks([])
    reference_ax.set_yticks([])
    reference_ax.spines[["top", "right"]].set_visible(False)

    marker_values = knn_scores[knn_scores["tf"].isin(markers)]["knn_footprint_oriented_z"]
    marker_arr = pd.to_numeric(marker_values, errors="coerce").to_numpy(dtype=float)
    finite = marker_arr[np.isfinite(marker_arr)]
    marker_vmax = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 2.0
    marker_vmax = max(marker_vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-marker_vmax, vcenter=0.0, vmax=marker_vmax)
    last_scatter = None
    for ax, tf in zip(axes[1:], markers, strict=True):
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
            norm=norm,
            s=1.8,
            alpha=0.92,
            linewidths=0,
            rasterized=True,
        )
        last_scatter = sc
        ax.set_title(tf, fontweight="bold", pad=4)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[["top", "right"]].set_visible(False)

    if last_scatter is not None:
        if cax is not None:
            cbar = fig.colorbar(last_scatter, cax=cax)
        else:
            cbar = fig.colorbar(last_scatter, ax=axes[1 : 1 + len(markers)], fraction=0.022, pad=0.018)
        cbar.ax.set_title("Footprint\nsignature", fontsize=8.5, fontweight="bold", pad=5)
        cbar.ax.tick_params(labelsize=9, length=2.2, width=0.6)


def plot_all_tf_signature_review_pdfs(
    annotations: pd.DataFrame,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    output_prefix: Path,
    *,
    panels_per_page: int = 12,
) -> None:
    ordered_meta = add_motif_specificity(metadata).copy()
    ordered_meta["_tf_sort"] = ordered_meta["tf_name"].astype(str)
    ordered_meta = ordered_meta.sort_values(
        ["dominant_cell_type", "cell_type_specificity", "dynamic_range", "_tf_sort", "motif_id"],
        ascending=[True, False, False, True, True],
    ).drop(columns="_tf_sort")

    x_padding = float((annotations["umap_1"].max() - annotations["umap_1"].min()) * 0.03)
    y_padding = float((annotations["umap_2"].max() - annotations["umap_2"].min()) * 0.03)
    xlim = (float(annotations["umap_1"].min() - x_padding), float(annotations["umap_1"].max() + x_padding))
    ylim = (float(annotations["umap_2"].min() - y_padding), float(annotations["umap_2"].max() + y_padding))

    all_values = matrix.to_numpy(dtype=float)
    finite = all_values[np.isfinite(all_values)]
    vmax = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 2.0
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    for cell_type in CELL_TYPES:
        subset_meta = ordered_meta[ordered_meta["dominant_cell_type"] == cell_type].reset_index(drop=True)
        output_pdf = output_prefix.with_name(f"{output_prefix.name}_{cell_type}.pdf")
        with PdfPages(output_pdf) as pdf:
            for page_start in range(0, subset_meta.shape[0], panels_per_page):
                page_meta = subset_meta.iloc[page_start : page_start + panels_per_page]
                fig, axes = plt.subplots(3, 4, figsize=(11.0, 8.5), sharex=True, sharey=True)
                fig.subplots_adjust(left=0.035, right=0.86, bottom=0.055, top=0.90, wspace=0.06, hspace=0.14)
                axes_flat = axes.reshape(-1)
                last_scatter = None
                for ax_index, ax in enumerate(axes_flat):
                    if ax_index >= page_meta.shape[0]:
                        ax.set_visible(False)
                        continue
                    row = page_meta.iloc[ax_index]
                    motif_id = str(row["motif_id"])
                    values = pd.to_numeric(matrix.loc[motif_id], errors="coerce")
                    value_table = values.rename("signature_z").reset_index()
                    value_table.columns = ["barcode", "signature_z"]
                    plot_df = annotations[["barcode", "umap_1", "umap_2"]].merge(
                        value_table,
                        on="barcode",
                        how="left",
                    )
                    last_scatter = ax.scatter(
                        plot_df["umap_1"],
                        plot_df["umap_2"],
                        c=plot_df["signature_z"],
                        cmap="RdBu_r",
                        norm=norm,
                        s=1.8,
                        alpha=0.92,
                        linewidths=0,
                        rasterized=True,
                    )
                    title = f"{row['tf_name']} ({motif_id})"
                    if len(title) > 42:
                        title = title[:39] + "..."
                    ax.set_title(title, fontsize=8.6, fontweight="bold", pad=3)
                    ax.set_xlim(xlim)
                    ax.set_ylim(ylim)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.spines[["top", "right"]].set_visible(False)
                    ax.spines[["left", "bottom"]].set_linewidth(0.5)
                if last_scatter is not None:
                    cax = fig.add_axes([0.895, 0.25, 0.018, 0.48])
                    cbar = fig.colorbar(last_scatter, cax=cax)
                    cbar.set_label("footprint signature z-score", fontsize=9)
                    cbar.ax.tick_params(labelsize=8, length=2.2, width=0.6)
                page_number = page_start // panels_per_page + 1
                total_pages = int(np.ceil(subset_meta.shape[0] / panels_per_page))
                fig.suptitle(
                    f"PBMC5k {cell_type} dominant footprint signatures ({page_number}/{total_pages})",
                    fontsize=11.0,
                    fontweight="bold",
                    y=0.98,
                )
                pdf.savefig(fig)
                plt.close(fig)
        print(f"Wrote {output_pdf}")


def plot_fig4_single_cell_footprinting(
    annotations: pd.DataFrame,
    knn_scores: pd.DataFrame,
    top_matrix: pd.DataFrame,
    top_metadata: pd.DataFrame,
    output_prefix: Path,
    markers: list[str],
) -> None:
    fig = plt.figure(figsize=(8.5, 11.0))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[0.36, 0.64],
        hspace=0.20,
        left=0.11,
        right=0.95,
        top=0.955,
        bottom=0.045,
    )
    top_grid = outer[0].subgridspec(
        2,
        2,
        height_ratios=[0.16, 1.0],
        width_ratios=[1.0, 0.026],
        hspace=0.04,
        wspace=0.03,
    )
    strip_ax = fig.add_subplot(top_grid[0, 0])
    heat_ax = fig.add_subplot(top_grid[1, 0])
    cax = fig.add_subplot(top_grid[1, 1])
    draw_top_motif_heatmap_panel(fig, strip_ax, heat_ax, cax, top_matrix, top_metadata, ordered_cell_annotations(annotations), UMAP_MARKER_LABELS)
    strip_ax.set_title("Top marker and cell-type-specific footprint signatures", fontsize=9, fontweight="bold", pad=4)
    fig.text(0.035, 0.965, "A", fontsize=9, fontweight="bold", fontfamily="Arial", va="top")

    bottom_grid = outer[1].subgridspec(3, 4, hspace=0.38, wspace=0.28, width_ratios=[1.0, 1.0, 1.0, 0.07])
    axes = [fig.add_subplot(bottom_grid[row, col]) for row in range(3) for col in range(3)]
    cax = fig.add_subplot(bottom_grid[:, 3])
    draw_knn_umap_panel(fig, axes[: 1 + len(markers)], annotations, knn_scores, markers, cax=cax)
    for ax in axes[1 + len(markers) :]:
        ax.set_visible(False)
    fig.text(0.035, 0.61, "B", fontsize=9, fontweight="bold", fontfamily="Arial", va="top")
    fig.suptitle("Single-cell footprinting", fontsize=9, fontweight="bold", y=0.992)
    save_illustrator_svg(fig, output_prefix, tight=False)
    plt.close(fig)


def score_all_motif_per_cell_heatmap(
    annotations: pd.DataFrame,
    fragments: Path,
    h5ad: Path,
    motif_table: pd.DataFrame,
    output_prefix: Path,
    *,
    batch_size: int,
    max_sites_per_motif: int | None,
    flank: int,
    center_half_width: int,
    flank_inner: int,
    flank_outer: int,
    knn: int,
    create_index: bool,
    top_per_cell_type: int,
    top_min_specificity: float,
    marker_scores: pd.DataFrame | None,
    fig4_output_prefix: Path | None,
    fig4_markers: list[str],
    all_tf_review_prefix: Path | None = None,
    all_tf_review_panels_per_page: int = 12,
) -> None:
    ordered = ordered_cell_annotations(annotations)
    cell_order = ordered["barcode"].astype(str).tolist()
    annotation_columns = annotations["barcode"].astype(str).tolist()
    neighbors = knn_indices(annotations, h5ad, knn)

    matrices = []
    metadata_rows = []
    motif_ids = motif_table["motif_id"].astype(str).tolist()
    for batch_index, batch_ids in enumerate(batched(motif_ids, batch_size), start=1):
        batch_meta = motif_table.set_index("motif_id").loc[batch_ids].reset_index()
        sites = {}
        for _, row in batch_meta.iterrows():
            motif_id = str(row["motif_id"])
            motif_sites = read_motif_bed_sites(Path(row["bed_path"]), max_sites_per_motif)
            if motif_sites:
                sites[motif_id] = motif_sites
        if not sites:
            continue
        print(f"Scoring all-motif batch {batch_index}: {len(sites)} motifs", flush=True)
        profiles = count_fragment_profiles(fragments, annotations, sites, flank, create_index=create_index)
        batch_scores = score_knn_profiles(profiles, neighbors, list(sites), center_half_width, flank_inner, flank_outer)
        batch_matrix = (
            batch_scores.pivot_table(index="tf", columns="cell_index", values="knn_footprint_score", aggfunc="mean")
            .reindex(index=list(sites), columns=range(annotations.shape[0]))
        )
        batch_matrix.columns = annotation_columns
        batch_matrix = batch_matrix.reindex(columns=cell_order)
        z_values = row_zscore(batch_matrix.to_numpy(dtype=float))
        z_matrix = pd.DataFrame(z_values, index=batch_matrix.index, columns=batch_matrix.columns)
        matrices.append(z_matrix)

        for row_index, motif_id in enumerate(z_matrix.index):
            motif_values = z_values[row_index]
            group_means = {
                cell_type: float(np.nanmean(motif_values[ordered["cell_type"].to_numpy() == cell_type]))
                for cell_type in CELL_TYPES
            }
            dominant = max(group_means, key=group_means.get)
            meta_row = batch_meta[batch_meta["motif_id"] == motif_id].iloc[0].to_dict()
            meta_row.update(
                {
                    "dominant_cell_type": dominant,
                    "dynamic_range": float(np.nanmax(motif_values) - np.nanmin(motif_values)),
                    **{f"{cell_type}_mean_z": group_means[cell_type] for cell_type in CELL_TYPES},
                }
            )
            metadata_rows.append(meta_row)

    if not matrices:
        raise SystemExit("No motif batches produced per-cell scores.")

    matrix = pd.concat(matrices, axis=0)
    metadata = pd.DataFrame(metadata_rows).set_index("motif_id").loc[matrix.index]
    metadata.insert(0, "motif_id", metadata.index)
    metadata = metadata.reset_index(drop=True)
    group_order = {cell_type: index for index, cell_type in enumerate(CELL_TYPES)}
    metadata["_group_order"] = metadata["dominant_cell_type"].map(group_order)
    metadata = metadata.sort_values(["_group_order", "dynamic_range", "motif_id"], ascending=[True, False, True]).drop(columns="_group_order")
    matrix = matrix.loc[metadata["motif_id"]]
    metadata.insert(0, "rank", np.arange(1, len(metadata) + 1))
    output = pd.concat([metadata.reset_index(drop=True), matrix.reset_index(drop=True)], axis=1)
    output.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)
    plot_all_motif_per_cell_heatmap(matrix, metadata, ordered, output_prefix)
    top_matrix, top_metadata = write_top_motif_heatmap(
        matrix,
        metadata,
        ordered,
        output_prefix.with_name("pbmc5k_top_motif_per_cell_footprint_signature_heatmap"),
        per_cell_type=top_per_cell_type,
        min_specificity=top_min_specificity,
        marker_scores=marker_scores,
    )
    if fig4_output_prefix is not None and marker_scores is not None:
        plot_fig4_single_cell_footprinting(annotations, marker_scores, top_matrix, top_metadata, fig4_output_prefix, fig4_markers)
    if all_tf_review_prefix is not None:
        plot_all_tf_signature_review_pdfs(
            annotations,
            matrix,
            metadata,
            all_tf_review_prefix,
            panels_per_page=all_tf_review_panels_per_page,
        )


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
    parser.add_argument("--all-motif-bindetect-dir", help="Optional BINDetect output directory containing */beds/*_all.bed files for all-motif per-cell heatmap scoring.")
    parser.add_argument("--all-motif-results", help="BINDetect results table used to order and annotate all-motif heatmap rows.")
    parser.add_argument("--all-motif-score-table", help="Existing all-motif per-cell heatmap TSV to redraw all/top heatmaps without rescoring fragments.")
    parser.add_argument("--marker-score-table", help="Existing KNN marker score table used to orient selected marker rows in top heatmaps and Figure 4.")
    parser.add_argument("--all-motif-batch-size", type=int, default=50, help="Number of motif signatures to score per batch for the all-motif heatmap.")
    parser.add_argument("--max-sites-per-motif", type=int, default=200, help="Maximum motif instances per motif for all-motif heatmap scoring; use 0 for all sites.")
    parser.add_argument("--max-motifs", type=int, help="Optional all-motif smoke-test limit.")
    parser.add_argument("--top-motif-signatures-per-cell-type", type=int, default=40, help="Top cell-type-specific all-motif signatures to keep per broad cell type (default: 40).")
    parser.add_argument("--top-motif-min-specificity", type=float, default=0.5, help="Minimum dominant-vs-next cell-type mean z-score difference for top all-motif heatmap rows (default: 0.5).")
    parser.add_argument("--fig4-output-prefix", default="Fig4", help="Output prefix for the combined letter-size Figure 4 SVG when all-motif heatmap data are available.")
    parser.add_argument("--all-tf-review-prefix", default="pbmc5k_all_tf_footprint_signature_umaps", help="Output prefix for three multi-page all-TF signature review PDFs grouped by dominant broad cell type.")
    parser.add_argument("--all-tf-review-panels-per-page", type=int, default=12, help="Number of TF signature UMAP panels per all-TF review PDF page (default: 12).")
    parser.add_argument("--skip-all-tf-review-pdfs", action="store_true", help="Do not write the three all-TF signature review PDFs.")
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
    marker_scores_for_figures = pd.read_csv(args.marker_score_table, sep="\t") if args.marker_score_table else knn_scores

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
    plot_per_cell_signature_heatmap(
        annotations,
        knn_scores,
        summary_markers,
        outdir / "pbmc5k_per_cell_footprint_signature_heatmap",
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
    if args.all_motif_score_table:
        ordered = ordered_cell_annotations(annotations)
        matrix, metadata = read_all_motif_score_table(Path(args.all_motif_score_table))
        matrix = matrix.reindex(columns=ordered["barcode"].astype(str).tolist())
        metadata = add_motif_specificity(metadata)
        plot_all_motif_per_cell_heatmap(
            matrix,
            metadata,
            ordered,
            outdir / "pbmc5k_all_motif_per_cell_footprint_signature_heatmap",
        )
        top_matrix, top_metadata = write_top_motif_heatmap(
            matrix,
            metadata,
            ordered,
            outdir / "pbmc5k_top_motif_per_cell_footprint_signature_heatmap",
            per_cell_type=args.top_motif_signatures_per_cell_type,
            min_specificity=args.top_motif_min_specificity,
            marker_scores=marker_scores_for_figures,
        )
        plot_fig4_single_cell_footprinting(
            annotations,
            marker_scores_for_figures,
            top_matrix,
            top_metadata,
            outdir / args.fig4_output_prefix,
            markers,
        )
        if not args.skip_all_tf_review_pdfs:
            plot_all_tf_signature_review_pdfs(
                annotations,
                matrix,
                metadata,
                outdir / args.all_tf_review_prefix,
                panels_per_page=args.all_tf_review_panels_per_page,
            )
    if args.all_motif_bindetect_dir or args.all_motif_results:
        if not args.all_motif_bindetect_dir or not args.all_motif_results:
            raise SystemExit("--all-motif-bindetect-dir and --all-motif-results must be provided together.")
        max_sites_per_motif = None if args.max_sites_per_motif == 0 else args.max_sites_per_motif
        motif_table = read_bindetect_motif_table(
            Path(args.all_motif_results),
            Path(args.all_motif_bindetect_dir),
            args.max_motifs,
        )
        score_all_motif_per_cell_heatmap(
            annotations,
            Path(args.fragments),
            Path(args.h5ad),
            motif_table,
            outdir / "pbmc5k_all_motif_per_cell_footprint_signature_heatmap",
            batch_size=args.all_motif_batch_size,
            max_sites_per_motif=max_sites_per_motif,
            flank=args.flank,
            center_half_width=args.center_half_width,
            flank_inner=args.flank_inner,
            flank_outer=args.flank_outer,
            knn=args.knn,
            create_index=not args.no_create_fragment_index,
            top_per_cell_type=args.top_motif_signatures_per_cell_type,
            top_min_specificity=args.top_motif_min_specificity,
            marker_scores=marker_scores_for_figures,
            fig4_output_prefix=outdir / args.fig4_output_prefix,
            fig4_markers=markers,
            all_tf_review_prefix=None if args.skip_all_tf_review_pdfs else outdir / args.all_tf_review_prefix,
            all_tf_review_panels_per_page=args.all_tf_review_panels_per_page,
        )
    print(f"Wrote per-cell PBMC5k signature plots and score tables to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
