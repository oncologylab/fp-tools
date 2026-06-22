# Data Availability Statement

`fp-tools` is openly available at https://github.com/oncologylab/fp-tools under the
MIT license. The repository contains all benchmark manifests, schemas, analysis
scripts, figure generators, and small test fixtures needed to reproduce the
workflows described in this manuscript.

Public input datasets are obtained programmatically with the included discovery and
download scripts:

- **Chromatin accessibility datasets:** ENCODE Project bulk ATAC-seq records
  used for K562/HepG2 replicate demonstrations and manifest-driven workflow
  checks.
- **Motif catalogs:** JASPAR 2026 CORE vertebrate non-redundant motifs, with
  optional HOCOMOCO-scale databases for larger motif stress tests.
- **PBMC5k pseudobulk and single-cell signature demonstration:** public 10x
  Genomics 5k PBMC single-cell ATAC fragments (`atac_pbmc_5k_nextgem`) are
  retrieved from the original 10x Genomics public dataset and prepared by
  `benchmarks/scripts/prepare_10x_pbmc5k_scatac.py`. The workflow groups cells
  into broad B-cell, monocyte, and T/NK labels, runs pseudobulk footprinting,
  and generates pairwise volcano plots, marker heatmaps, and KNN-smoothed
  footprint-signature UMAPs with the PBMC5k plotting scripts.
Large raw public inputs, full pseudobulk fragments, cut-site bigWigs, and full
benchmark outputs are intentionally **not** stored in the main code repository;
they are regenerated from committed, versioned manifests using `benchmarks/scripts/`.
Small result tables, figure source tables, manuscript figures, and reproducibility
scripts stay in the main repository. The reviewer-facing reproduction guide is
provided in `docs/reproduce-paper.md`, with Conda and Docker environments in
`environment.yml` and `Dockerfile`. Each benchmark result records random seeds,
tool versions, command lines, and the resolved manifest to support exact
reproduction.
