#!/usr/bin/env python
"""Prepare 10x PBMC pseudobulk inputs for the fp-tools scATAC example."""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from pathlib import Path

import pandas as pd

TENX_URLS = {
    "fragments": "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
    "fragment_index": "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi",
    "analysis": "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_analysis.tar.gz",
    "peaks": "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peaks.bed",
    "peak_annotation": "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv",
}

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

CLUSTER_TO_CELL_TYPE = {
    1: "CD4_T",
    2: "CD14_Monocyte",
    3: "CD14_Monocyte",
    4: "CD4_T",
    5: "CD4_T",
    6: "NK_T_cytotoxic",
    7: "CD4_T",
    8: "B_cell",
    9: "FCGR3A_Monocyte",
    10: "Mixed_myeloid",
    11: "NK_T_cytotoxic",
    12: "Dendritic_cell",
    13: "Platelet_lowRNA",
    14: "CD4_T",
}

TF_MOTIFS = {
    "PAX5": ["PAX5_"],
    "TCF7": ["TCF7_", "TCF7L2_"],
    "CEBPB": ["CEBPB_"],
    "CTCF": ["CTCF_"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(output: Path) -> None:
    rows = [{"dataset": "10x_pbmc_multiome", "asset": key, "url": url} for key, url in TENX_URLS.items()]
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


def write_annotations(clusters_csv: Path, output: Path) -> pd.DataFrame:
    clusters = pd.read_csv(clusters_csv)
    clusters["cell_type"] = clusters["Cluster"].map(CLUSTER_TO_CELL_TYPE).fillna("Unassigned")
    clusters = clusters.rename(columns={"Barcode": "barcode", "Cluster": "tenx_gex_graph_cluster"})
    clusters[["barcode", "tenx_gex_graph_cluster", "cell_type"]].to_csv(output, sep="\t", index=False)
    return clusters


def write_chrom_sizes(output: Path, chroms: list[str]) -> None:
    with output.open("w", encoding="utf-8") as handle:
        for chrom in chroms:
            handle.write(f"{chrom}\t{HG38_CHROM_SIZES[chrom]}\n")


def write_tf_beds(mapping_bed: Path, outdir: Path, chroms: set[str], max_sites: int) -> dict[str, int]:
    outdir.mkdir(parents=True, exist_ok=True)
    handles = {tf: (outdir / f"{tf}.motif_peaks.bed").open("w", encoding="utf-8") for tf in TF_MOTIFS}
    seen: dict[str, set[tuple[str, str, str]]] = {tf: set() for tf in TF_MOTIFS}
    try:
        with mapping_bed.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, start, end, motif, *_ = line.rstrip("\n").split("\t")
                if chrom not in chroms:
                    continue
                for tf, prefixes in TF_MOTIFS.items():
                    if not any(motif.startswith(prefix) for prefix in prefixes):
                        continue
                    key = (chrom, start, end)
                    if key in seen[tf] or len(seen[tf]) >= max_sites:
                        continue
                    seen[tf].add(key)
                    center = (int(start) + int(end)) // 2
                    handles[tf].write(f"{chrom}\t{max(0, center - 1)}\t{center + 1}\t{tf}\t0\t.\n")
    finally:
        for handle in handles.values():
            handle.close()
    return {tf: len(sites) for tf, sites in seen.items()}


def package_example_data(paths: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(output, "w:gz") as archive:
        for path in paths:
            if path.exists():
                archive.add(path, arcname=path.name)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", default="data/public/raw/10x_pbmc")
    parser.add_argument("--outdir", default="data/public/processed/pseudobulk_pbmc")
    parser.add_argument("--chroms", default="chr1,chr2", help="Comma-separated chromosomes for the compact example.")
    parser.add_argument("--max-sites-per-tf", type=int, default=1000)
    parser.add_argument("--write-example-archive", action="store_true")
    args = parser.parse_args(argv)

    raw_dir = Path(args.raw_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chroms = [chrom.strip() for chrom in args.chroms.split(",") if chrom.strip()]

    write_manifest(Path("benchmarks/manifests/10x_pbmc_pseudobulk.tsv"))
    annotations = write_annotations(raw_dir / "analysis/clustering/gex/graphclust/clusters.csv", outdir / "pbmc_10x_cell_annotations.tsv")
    write_chrom_sizes(outdir / "hg38.chrom.sizes", chroms)
    site_counts = write_tf_beds(raw_dir / "analysis/tf_analysis/peak_motif_mapping.bed", outdir / "tf_sites", set(chroms), args.max_sites_per_tf)

    summary = {
        "dataset": "10x PBMC granulocyte-sorted Multiome, Cell Ranger ARC 2.0.0",
        "chromosomes": chroms,
        "n_cells": int(len(annotations)),
        "cell_type_counts": annotations["cell_type"].value_counts().to_dict(),
        "tf_site_counts": site_counts,
        "source_urls": TENX_URLS,
    }
    (outdir / "pbmc_pseudobulk_preparation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.Series(summary["cell_type_counts"], name="n_cells").rename_axis("cell_type").reset_index().to_csv(outdir / "pbmc_cell_type_counts.tsv", sep="\t", index=False)

    if args.write_example_archive:
        archive = outdir / "fp-tools-pbmc-pseudobulk-example-inputs.tar.gz"
        package_example_data(
            [
                outdir / "pbmc_10x_cell_annotations.tsv",
                outdir / "hg38.chrom.sizes",
                outdir / "pbmc_cell_type_counts.tsv",
                outdir / "pbmc_pseudobulk_preparation_summary.json",
                *sorted((outdir / "tf_sites").glob("*.bed")),
            ],
            archive,
        )
        (outdir / "fp-tools-pbmc-pseudobulk-example-inputs.sha256").write_text(f"{sha256(archive)}  {archive.name}\n", encoding="utf-8")
        print(f"Wrote {archive}")

    print(f"Wrote annotations and TF site sets to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
