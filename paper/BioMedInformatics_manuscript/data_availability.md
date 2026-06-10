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
- **Variant / pseudobulk benchmarks (later tiers):** public caQTL and
  allele-specific datasets, and public 10x scATAC/multiome fragments with cell
  annotations.

Large raw and processed data and full benchmark outputs are intentionally **not**
stored in version control; they are regenerated from the committed, versioned
manifests using `benchmarks/scripts/` (discovery, resumable download with checksum
reports, label-overlap table building, metric/calibration/bootstrap computation,
and the benchmark pipeline runner). Each benchmark result records random seeds,
tool versions, command lines, and the resolved manifest to support exact
reproduction.
