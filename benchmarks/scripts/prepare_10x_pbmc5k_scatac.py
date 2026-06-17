#!/usr/bin/env python
"""Prepare the 10x PBMC5k scATAC dataset used by the scPrinter tutorial."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import snapatac2 as snap


HG38_CHROM_SIZES = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
    "chrY": 57_227_415,
    "chrM": 16_569,
}

BROAD_CELL_TYPE = {
    "Memory B": "B_cell",
    "Naive B": "B_cell",
    "CD14 Mono": "Monocyte",
    "CD16 Mono": "Monocyte",
    "cDC": "Monocyte",
    "pDC": "Monocyte",
    "CD4 Memory": "T_NK_cell",
    "CD4 Naive": "T_NK_cell",
    "CD8 Memory": "T_NK_cell",
    "CD8 Naive": "T_NK_cell",
    "MAIT": "T_NK_cell",
    "NK": "T_NK_cell",
}

REGION_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


def link_or_copy(source: Path, dest: Path, copy: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        return
    if copy:
        shutil.copy2(source, dest)
    else:
        os.symlink(source, dest)


def write_chrom_sizes(path: Path, chroms: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for chrom in chroms:
            handle.write(f"{chrom}\t{HG38_CHROM_SIZES[chrom]}\n")


def write_selected_regions(adata: ad.AnnData, output: Path, chroms: set[str] | None, selected_only: bool = True) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    selected = adata.var["selected"].to_numpy() if selected_only else None
    with output.open("w", encoding="utf-8") as handle:
        for idx, name in enumerate(adata.var_names):
            if selected is not None and not bool(selected[idx]):
                continue
            match = REGION_RE.match(str(name))
            if match is None:
                continue
            chrom = match.group("chrom")
            if chroms is not None and chrom not in chroms:
                continue
            handle.write(f"{chrom}\t{match.group('start')}\t{match.group('end')}\t{name}\n")
            count += 1
    return count


def write_annotations_and_umap(adata: ad.AnnData, outdir: Path) -> pd.DataFrame:
    obs = adata.obs.copy()
    umap = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["umap_1", "umap_2"])
    annotations = obs.join(umap)
    annotations = annotations.reset_index(names="barcode")
    annotations["snap_cell_type"] = annotations["cell_type"].astype(str)
    annotations["cell_type"] = annotations["snap_cell_type"].map(BROAD_CELL_TYPE).fillna("Unassigned")
    columns = [
        "barcode",
        "cell_type",
        "snap_cell_type",
        "umap_1",
        "umap_2",
        "n_fragment",
        "tsse",
        "frac_dup",
        "frac_mito",
        "doublet_probability",
        "doublet_score",
    ]
    annotations[columns].to_csv(outdir / "pbmc5k_scprinter_broad_annotations.tsv", sep="\t", index=False)
    annotations[["barcode", "snap_cell_type"]].to_csv(outdir / "pbmc5k_snatac2_subtype_annotations.tsv", sep="\t", index=False)
    return annotations[columns]


def plot_umap(annotations: pd.DataFrame, column: str, output_prefix: Path, title: str) -> None:
    colors = {
        "B_cell": "#3B82F6",
        "Monocyte": "#D97706",
        "T_NK_cell": "#059669",
        "Unassigned": "#9CA3AF",
        "Memory B": "#2563EB",
        "Naive B": "#60A5FA",
        "CD14 Mono": "#B45309",
        "CD16 Mono": "#F59E0B",
        "cDC": "#92400E",
        "pDC": "#FBBF24",
        "CD4 Memory": "#047857",
        "CD4 Naive": "#34D399",
        "CD8 Memory": "#0F766E",
        "CD8 Naive": "#5EEAD4",
        "MAIT": "#65A30D",
        "NK": "#16A34A",
    }
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for label, group in annotations.groupby(column, sort=True):
        ax.scatter(group["umap_1"], group["umap_2"], s=7, alpha=0.75, linewidths=0, label=label, color=colors.get(label, None))
    for label, group in annotations.groupby(column, sort=True):
        if label == "Unassigned":
            continue
        ax.text(group["umap_1"].median(), group["umap_2"].median(), str(label), fontsize=8, weight="bold", ha="center", va="center")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), dpi=220)
    fig.savefig(output_prefix.with_suffix(".pdf"))
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default="data/public/raw/10x_pbmc5k_scatac")
    parser.add_argument("--outdir", default="data/public/processed/pseudobulk_pbmc5k_scatac")
    parser.add_argument("--chroms", default="chr1,chr2", help="Comma-separated chromosomes for compact demo outputs.")
    parser.add_argument("--copy-fragments", action="store_true", help="Copy the 1 GB fragment file instead of symlinking the SnapATAC2 cache file.")
    args = parser.parse_args(argv)

    raw_dir = Path(args.raw_dir)
    outdir = Path(args.outdir)
    plot_dir = outdir / "plots"
    raw_dir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    chroms = [chrom.strip() for chrom in args.chroms.split(",") if chrom.strip()]
    chrom_set = set(chroms) if chroms else None

    fragment_source = Path(snap.datasets.pbmc5k(type="fragment"))
    h5ad_source = Path(snap.datasets.pbmc5k(type="annotated_h5ad"))
    fragment_dest = raw_dir / "atac_pbmc_5k_nextgem_fragments.tsv.gz"
    h5ad_dest = raw_dir / "atac_pbmc_5k_annotated.h5ad"
    link_or_copy(fragment_source, fragment_dest, args.copy_fragments)
    link_or_copy(h5ad_source, h5ad_dest, args.copy_fragments)

    adata = ad.read_h5ad(h5ad_source, backed="r")
    try:
        annotations = write_annotations_and_umap(adata, outdir)
        n_regions = write_selected_regions(adata, raw_dir / "atac_pbmc_5k_snatac2_selected_bins.bed", chroms=None)
        n_demo_regions = write_selected_regions(adata, raw_dir / "atac_pbmc_5k_snatac2_selected_bins.demo.bed", chroms=chrom_set)
    finally:
        adata.file.close()

    write_chrom_sizes(outdir / "hg38.chrom.sizes", chroms)
    plot_umap(annotations, "cell_type", plot_dir / "pbmc5k_umap_broad_3celltypes", "PBMC5k scATAC broad annotations")
    plot_umap(annotations, "snap_cell_type", plot_dir / "pbmc5k_umap_snatac2_subtypes", "PBMC5k scATAC SnapATAC2 annotations")

    summary = {
        "dataset": "10x PBMC5k scATAC, same fragment dataset used by the scPrinter PBMC scATAC tutorial",
        "fragment_file": str(fragment_dest),
        "annotated_h5ad": str(h5ad_dest),
        "annotation_file": str(outdir / "pbmc5k_scprinter_broad_annotations.tsv"),
        "selected_bins_bed": str(raw_dir / "atac_pbmc_5k_snatac2_selected_bins.bed"),
        "demo_selected_bins_bed": str(raw_dir / "atac_pbmc_5k_snatac2_selected_bins.demo.bed"),
        "chromosomes": chroms,
        "n_cells": int(len(annotations)),
        "cell_type_counts": annotations["cell_type"].value_counts().to_dict(),
        "snap_cell_type_counts": annotations["snap_cell_type"].value_counts().to_dict(),
        "n_selected_bins": n_regions,
        "n_demo_selected_bins": n_demo_regions,
        "source_fragment_url": "https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fragments.tsv.gz",
        "source_h5ad_provider": "snapatac2.datasets.pbmc5k(type='annotated_h5ad')",
    }
    (outdir / "pbmc5k_scatac_preparation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.Series(summary["cell_type_counts"], name="n_cells").rename_axis("cell_type").reset_index().to_csv(outdir / "pbmc5k_broad_cell_type_counts.tsv", sep="\t", index=False)
    pd.Series(summary["snap_cell_type_counts"], name="n_cells").rename_axis("snap_cell_type").reset_index().to_csv(outdir / "pbmc5k_snatac2_cell_type_counts.tsv", sep="\t", index=False)

    print(f"Wrote PBMC5k scATAC annotations, UMAPs, and selected-bin BEDs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
