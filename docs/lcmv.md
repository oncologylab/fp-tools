# LCMV CD8 public data

The expanded LCMV CD8 v2 collection strengthens the paired ATAC-seq and
RNA-seq analysis while preserving the original audited dataset as v1.

## Coverage

| Layer | Conditions | ATAC units | RNA samples | Comparisons |
| --- | ---: | ---: | ---: | ---: |
| Primary paired | 14 | 37 | 35 | 18 ATAC + 18 RNA |
| Supporting assay-only | 6 | 11 | 2 | 7 ATAC + 3 RNA |
| Complete v2 | 20 | 48 | 37 | — |

V2 resolves 87 GSM libraries to 95 runs and 159 FASTQ files. Every primary
condition has both assays from the same study. The assays are paired at the
condition level—not from the same mouse, cells, or aliquot.

## Added paired conditions

- Milner: day-7 terminal effector and memory-precursor cells from spleen.
- Scott-Browne: day-8 effector, day-35 memory, and day-20 exhausted cells.

These additions complement the existing Milner intestinal TRM, Guan
naïve/MPEC/SLEC, Scott-Browne naïve, and four Beltra exhaustion states.

Eight primary contrasts match study, infection, day, and tissue. Ten
within-study trajectory contrasts are retained with explicit annotations for
time, tissue, infection, or baseline confounding. Cross-study integration is
descriptive only; v2 does not report pooled differential tests.

## Reproducibility contracts

- Library selection: `benchmarks/manifests/compact/lcmv_cd8_libraries_v2.tsv`
- Comparison design: `benchmarks/manifests/compact/lcmv_cd8_comparisons_v2.tsv`
- Restartable runner: `benchmarks/scripts/run_lcmv_v2.py`
- Output validator: `benchmarks/scripts/validate_lcmv_outputs.py`

The original selection remains available as
`benchmarks/manifests/compact/lcmv_cd8_libraries.tsv`.

## Server locations

| Data | Workspace-relative location |
| --- | --- |
| Shared raw FASTQs, references, and ATAC inputs | `data/public/raw/lcmv_cd8_bulk/` |
| Versioned v2 raw metadata | `data/public/raw/lcmv_cd8_bulk/metadata/v2/` |
| Audited v1 results | `data/public/processed/lcmv_cd8_bulk_fp_rna/` |
| Expanded v2 results | `data/public/processed/lcmv_cd8_bulk_fp_rna_v2/` |
| V2 inventory and transfer lists | `data/public/processed/lcmv_cd8_bulk_fp_rna_v2/validation/` |

Generated transfer lists provide a downstream-only package and a full
reproducibility package containing FASTQs and raw ATAC inputs. Large public
data and generated results remain ignored and are never committed to GitHub.
