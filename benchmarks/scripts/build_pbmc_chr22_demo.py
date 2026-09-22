#!/usr/bin/env python3
"""Build the versioned, real PBMC5k chromosome-22 teaching dataset."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import pysam
import yaml


SOURCE_SHA256 = {
    "fragments": "5fe44c0f8f76ce1534c1ae418cf0707ca5ef712004eee77c3d98d2d4b35ceaec",
    "h5ad": "592f1551c27d0cfe4d81e7febad624d6b7d3ebf977b0c3ea64e06b3f3d76f078",
    "annotations": "54d9428733c49ed4978baf390dfa91851f7c6d56abf9309e4bb874f4e9efcf39",
    "genome": "5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51",
    "blacklist": "31c69342df43bbc19dd8ef2886611a8150cb53e90f341de14d30e742c9251737",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build(args):
    for key, expected in SOURCE_SHA256.items():
        if sha256(getattr(args, key)) != expected:
            raise ValueError(f"Source checksum mismatch: {key}; do not redefine demo v1")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=False)
    annotations = pd.read_csv(args.annotations, sep="\t")
    annotations = annotations.sort_values("barcode").groupby("cell_type", sort=True).head(100)
    annotations = annotations.sort_values("barcode").reset_index(drop=True)
    expected = {"B_cell": 100, "Monocyte": 100, "T_NK_cell": 100}
    if annotations.cell_type.value_counts().to_dict() != expected or not annotations.barcode.is_unique:
        raise ValueError("Expected 100 distinct cells in each of the three PBMC groups")
    barcodes = set(annotations.barcode)
    source = ad.read_h5ad(args.h5ad, backed="r")
    try:
        rows = source.obs_names.get_indexer(annotations.barcode)
        if (rows < 0).any():
            raise ValueError("Annotation barcodes absent from AnnData")
        columns = np.flatnonzero(source.var_names.str.startswith("chr22:"))
        # Slice the sparse count matrix directly; do not load fragment_paired.
        counts = source.X[rows, :][:, columns]
        subset = ad.AnnData(X=counts, obs=source.obs.iloc[rows].copy(), var=source.var.iloc[columns].copy())
        for key in ("X_spectral", "X_umap"):
            if key in source.obsm:
                subset.obsm[key] = np.asarray(source.obsm[key][rows])
        if not subset.var["selected"].any():
            raise ValueError("No selected chromosome-22 bins")
        subset.write_h5ad(out / "genomic_bin_counts.h5ad", compression="gzip")
    finally:
        source.file.close()
    annotations.to_csv(out / "cell_annotations.tsv", sep="\t", index=False)
    fragment_counts = Counter()
    with pysam.TabixFile(str(args.fragments)) as fragments, (out / "fragments.tsv").open("w") as handle:
        for record in fragments.fetch("chr22"):
            barcode = record.split("\t")[3]
            if barcode in barcodes:
                handle.write(record + "\n")
                fragment_counts[barcode] += 1
    if set(fragment_counts) != barcodes:
        raise ValueError("Some selected cells have no chromosome-22 fragments")
    pysam.tabix_compress(str(out / "fragments.tsv"), str(out / "fragments.tsv.gz"))
    pysam.tabix_index(str(out / "fragments.tsv.gz"), preset="bed")
    (out / "fragments.tsv").unlink()
    with pysam.FastaFile(str(args.genome)) as genome:
        sequence = genome.fetch("chr22")
    with (out / "hg38_chr22.fa").open("w") as handle:
        handle.write(">chr22\n")
        for start in range(0, len(sequence), 60):
            handle.write(sequence[start:start + 60] + "\n")
    pysam.faidx(str(out / "hg38_chr22.fa"))
    (out / "hg38_chr22.chrom.sizes").write_text(f"chr22\t{len(sequence)}\n")
    with (out / "selected_bins.bed").open("w") as handle:
        for name in subset.var_names[subset.var["selected"].astype(bool)]:
            chrom, interval = name.split(":")
            start, end = interval.split("-")
            handle.write(f"{chrom}\t{start}\t{end}\n")
    with Path(args.blacklist).open() as source_blacklist, (out / "hg38_chr22.blacklist.bed").open("w") as handle:
        for line in source_blacklist:
            if line.startswith("chr22\t"):
                handle.write(line)
    config = {"version": 1, "run_mode": "single", "samples": [{
        "sample_id": "pbmc_chr22", "tool": "sc-footprinting",
        "fragments": "fragments.tsv.gz", "annotations": "cell_annotations.tsv",
        "h5ad": "genomic_bin_counts.h5ad", "group_by": "cell_type",
        "genome_sizes": "hg38_chr22.chrom.sizes", "genome": "hg38_chr22.fa",
        "peaks": "selected_bins.bed", "blacklist": "hg38_chr22.blacklist.bed",
        "motif_db": "jaspar2026_vertebrates", "outdir": "results",
    }], "comparisons": []}
    (out / "workflow.yml").write_text(yaml.safe_dump(config, sort_keys=False))
    manifest = {
        "dataset": "fp-tools-pbmc-chr22-demo-v1", "assembly": "hg38", "chromosome": "chr22",
        "selection": "First 100 lexicographically sorted barcodes per broad cell type",
        "cells": len(annotations), "groups": expected, "bins": subset.n_vars,
        "selected_bins": int(subset.var["selected"].sum()),
        "fragment_records": sum(fragment_counts.values()),
        "sources": {
            "fragments": "https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fragments.tsv.gz",
            "annotated_h5ad": "https://osf.io/download/e9vc3/",
            "genome": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            "blacklist": "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz",
        },
        "source_sha256": SOURCE_SHA256,
        "limitations": "Teaching subset, not a binding-validation benchmark. Counts, selected-bin mask and embeddings are inherited from the full source dataset; embeddings were not recomputed. Selected bins are accessible-region proxies, not newly called peaks. No biological replicates.",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "README.txt").write_text(
        "PBMC chromosome-22 example v1\n\n"
        "Extract this entire folder, then open a terminal in it.\n"
        "Run: run-yaml-workflow --config workflow.yml\n"
        "For the desktop app, load workflow.yml in Config and replace input/output paths with full paths.\n"
        "First output: results/pseudobulk_footprint_manifest.tsv\n"
        "Per-cell figures: results/plots/single_cell_footprinting/\n"
        "See manifest.json for source hashes and interpretation limits.\n"
    )
    files = sorted(path for path in out.iterdir() if path.is_file())
    (out / "SHA256SUMS.txt").write_text("".join(f"{sha256(path)}  {path.name}\n" for path in files))
    archive = out.with_suffix(".zip")
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for path in sorted(out.iterdir()):
            info = zipfile.ZipInfo(f"{out.name}/{path.name}", date_time=(2026, 9, 21, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            bundle.writestr(info, path.read_bytes())
    archive.with_suffix(".zip.sha256").write_text(f"{sha256(archive)}  {archive.name}\n")
    print(json.dumps({"archive": str(archive), "bytes": archive.stat().st_size, **manifest}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fragments", "h5ad", "annotations", "genome", "blacklist", "outdir"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    build(parser.parse_args())
