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


MARKERS = ("PAX5", "CEBPB", "TCF7")
CELL_TYPES = ("B_cell", "Monocyte", "T_NK_cell")
MARKER_GROUPS = {
    "PAX5": "B_cell",
    "CEBPB": "Monocyte",
    "TCF7": "T_NK_cell",
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


def label_groups(ax: plt.Axes, annotations: pd.DataFrame, fontsize: int = 7) -> None:
    for label, group in annotations.groupby("cell_type", sort=True):
        text = ax.text(
            group["umap_1"].median(),
            group["umap_2"].median(),
            str(label),
            fontsize=fontsize,
            weight="bold",
            ha="center",
            va="center",
            color="#111827",
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
        default="PAX5:B_cell,CEBPB:Monocyte,TCF7:T_NK_cell",
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

    plot_score_grid(
        annotations,
        [("KNN footprint", knn_scores, "knn_footprint_oriented_z", "oriented z")],
        markers,
        outdir / "pbmc5k_knn_footprint_signature_umap",
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
