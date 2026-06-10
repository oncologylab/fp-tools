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

## Motif-Removal Benchmarks

For motif-relaxed and motif-free recovery experiments, use `scripts/build_motif_removal_benchmark.py` to turn scored site tables into long-form benchmark predictions. The resulting TSV can be evaluated with `scripts/compute_binary_metrics.py` and rendered with `../paper/scripts/plot_benchmark_panels.py`.

## Calibration Benchmarks

Use `scripts/compute_calibration.py` with supervised prediction probabilities to produce reliability bins, ECE/MCE, and Brier summaries. Render those outputs with `../paper/scripts/plot_calibration_panels.py` for manuscript-ready calibration figures.

## Public-Label Benchmark Tables

Use `scripts/build_label_overlap_benchmark.py` to convert scored prediction intervals and public TF-binding BED labels into the standard TSV consumed by the metric, calibration, and figure scripts.

## Benchmark Result Folders

Use `scripts/run_benchmark_pipeline.py` after label-overlap table creation to combine one or more labeled prediction TSVs into a reproducible result folder containing metrics, calibration reports, optional bootstrap confidence intervals, and PDF/SVG/PNG multi-panel figures.
