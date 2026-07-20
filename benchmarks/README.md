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

The audited v1 selection remains in `manifests/compact/lcmv_cd8_libraries.tsv`.
The expanded v2 selection and its explicit evidence tiers are in
`manifests/compact/lcmv_cd8_libraries_v2.tsv` and
`manifests/compact/lcmv_cd8_comparisons_v2.tsv`. Bulky FASTQs and results remain
under ignored `data/public/` directories.

```bash
python benchmarks/scripts/run_lcmv_v2.py --dry-run --stage resolve
python benchmarks/scripts/run_lcmv_v2.py
```

V2 contains 87 GSM libraries, including 14 primary condition-level ATAC/RNA
pairs and a separate assay-only supporting layer. Pairing never implies the
same mouse or aliquot. The 18 primary comparisons are strictly within-study;
eight are matched-context comparisons and ten carry explicit time, tissue, or
infection caveats. Guan and Beltra use paper-specific k=31 Kallisto count
layers, while Milner and Scott-Browne use TopHat2/HTSeq counts. Uniform k=21
Kallisto supports descriptive integration only. Pooled cross-study
differential testing is intentionally absent from v2.

The runner is restartable and fingerprints every stage. The validator checks
accessions, ENA byte counts and MD5s, BAM/BED/bigWig structure, motif and
comparison completeness, count matrices, and differential tables. QC warnings
are reported rather than silently excluding samples. See the public
[LCMV data page](../docs/lcmv.md) for the condition matrix and transfer layout.
