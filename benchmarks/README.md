# fp-tools Benchmarks

This directory contains lightweight benchmark scaffolding from `DEV_PLAN.md`.
Large public data and generated benchmark results are intentionally ignored by git.

## Layout

- `manifests/`: versioned public-data manifests and schema documentation.
- `scripts/`: data discovery, download, metrics, and figure-generation helpers.
- `results/`: ignored output directory for benchmark runs.
- `download_reports/`: ignored reports from resumable public-data downloads.

## First benchmark priority

Start with released human GRCh38 ENCODE bulk ATAC-seq experiments matched to TF ChIP-seq or CUT&RUN labels. Commit manifests and scripts, not downloaded BAMs, bigWigs, or full outputs.

## Deferred Benchmark Scaffolds

Motif-relaxed/motif-free recovery and supervised calibration scripts remain in
this directory as development scaffolds, but they are not primary supported
workflows. See `../DEV_PLAN.md` for the validation boundary.

## Public-Label Benchmark Tables

Use `scripts/build_label_overlap_benchmark.py` to convert scored prediction intervals and public TF-binding BED labels into the standard TSV consumed by the metric, calibration, and figure scripts.

## Benchmark Result Folders

Use `scripts/run_benchmark_pipeline.py` after label-overlap table creation to combine one or more labeled prediction TSVs into a reproducible result folder containing metrics, calibration reports, optional bootstrap confidence intervals, and PDF/SVG/PNG multi-panel figures.

## LCMV CD8 multimodal collection

`manifests/compact/lcmv_cd8_libraries.tsv` is the curated GSM-level selection
for the four LCMV CD8 studies. The associated helpers deliberately keep bulky
FASTQs and results under ignored `data/public/` directories:

```bash
python benchmarks/scripts/resolve_lcmv_cd8_collection.py
python benchmarks/scripts/build_lcmv_cd8_downstream.py
python benchmarks/scripts/summarize_lcmv_rna.py \
  --project data/public/processed/lcmv_cd8_bulk_fp_rna \
  --gtf data/public/raw/lcmv_cd8_bulk/reference/mm10/gencode_m25/gencode.vM25.annotation.gtf
Rscript benchmarks/scripts/analyze_lcmv_rna.R \
  data/public/processed/lcmv_cd8_bulk_fp_rna
python benchmarks/scripts/validate_lcmv_outputs.py --verify-checksums
```

ATAC and RNA entries are matched at the study/condition level, not as aliquots
from the same mouse. The primary RNA contrasts use the paper-specific Kallisto
layer within each study. Beltra contrasts also include TMM/voom/limma results.
The uniform 21-mer Kallisto layer accommodates Milner's 25-base reads and is
reserved for cross-study visualization and explicitly exploratory RUVr
comparisons. The validator checks accessions, ENA checksums, BAM/BED/bigWig
integrity, motif and comparison completeness, count matrices, and DE tables;
its QC flags require review rather than automatic sample removal.
