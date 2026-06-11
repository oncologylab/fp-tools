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
this directory as development scaffolds, but they are not first-version paper
requirements. See `../docs/FUTURE_PLANS.md` for the validation work needed before
those workflows are presented as supported capabilities.

## Public-Label Benchmark Tables

Use `scripts/build_label_overlap_benchmark.py` to convert scored prediction intervals and public TF-binding BED labels into the standard TSV consumed by the metric, calibration, and figure scripts.

## Benchmark Result Folders

Use `scripts/run_benchmark_pipeline.py` after label-overlap table creation to combine one or more labeled prediction TSVs into a reproducible result folder containing metrics, calibration reports, optional bootstrap confidence intervals, and PDF/SVG/PNG multi-panel figures.
