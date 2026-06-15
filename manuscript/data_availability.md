# Data Availability Statement

`fp-tools` is openly available at https://github.com/oncologylab/fp-tools under the
MIT license. The repository contains all benchmark manifests, schemas, analysis
scripts, figure generators, and small test fixtures needed to reproduce the
workflows described in this manuscript.

Public input datasets are obtained programmatically with the included discovery and
download scripts:

- **Chromatin accessibility / TF occupancy labels:** ENCODE Project bulk ATAC-seq
  and matched TF ChIP-seq / CUT&RUN peak calls (human, GRCh38), e.g. K562 CTCF
  ATAC-seq (ENCFF926KTI) and K562 CTCF ChIP-seq (ENCFF362OPG).
- **Motif catalogs:** JASPAR 2026 CORE vertebrate non-redundant motifs, with
  optional HOCOMOCO-scale databases for larger motif stress tests.
- **Pseudobulk demonstration:** public 10x PBMC Multiome ATAC fragments,
  official 10x clustering outputs, TF-analysis motif mappings, and peak files are
  retrieved through `benchmarks/manifests/compact/10x_pbmc_pseudobulk.tsv` and prepared
  by `benchmarks/scripts/prepare_10x_pbmc_pseudobulk.py`.
- **Variant benchmarks (later tiers):** public caQTL and allele-specific datasets
  will be retrieved by manifest-driven scripts as those tiers are finalized.

Large raw public inputs, full pseudobulk fragments, cut-site bigWigs, and full
benchmark outputs are intentionally **not** stored in the main code repository;
they are regenerated from committed, versioned manifests using `benchmarks/scripts/`.
Small result tables, figure source tables, manuscript figures, and reproducibility
scripts stay in the main repository. The reviewer-facing reproduction guide is
provided in `docs/reproduce-paper.md`, with Conda and Docker environments in
`environment.yml` and `Dockerfile`. Large reusable example-data archives are
stored separately as release assets in `oncologylab/fp-tools-data` (https://github.com/oncologylab/fp-tools-data/releases/tag/pbmc-pseudobulk-v1). Each benchmark
result records random seeds, tool versions, command lines, and the resolved
manifest to support exact reproduction.
