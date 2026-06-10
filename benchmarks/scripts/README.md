# Benchmark Scripts

Planned helpers from `DEV_PLAN.md`:

- `build_encode_manifest.py`: query ENCODE and write a public-data manifest without downloading files.
- `download_manifest.py`: resumable downloads plus checksum and path reports.
- `compute_binary_metrics.py`: AUROC, AUPRC, recall@FDR, Brier score summaries, and optional bootstrap confidence intervals from scored labels.
- `compute_calibration.py`: reliability-bin, expected calibration error, maximum calibration error, and Brier summaries from probability-like predictions.
- `build_label_overlap_benchmark.py`: convert scored BED-like prediction intervals plus ChIP/CUT&RUN label BEDs into metrics-ready binary label/score tables.
- `build_motif_removal_benchmark.py`: create long-form motif-removal recovery benchmark tables from baseline, motif-free, supervised, or reranked site scores.
- `run_benchmark_pipeline.py`: combine labeled prediction TSVs, compute metrics/calibration/bootstrap summaries, and write PDF/SVG/PNG benchmark figures.
- `paper/scripts/plot_benchmark_panels.py`: PDF/SVG/PNG multi-panel benchmark figures for the BioMedInformatics manuscript.
- `paper/scripts/plot_calibration_panels.py`: PDF/SVG/PNG reliability curves and ECE panels.
- `paper/scripts/plot_multiscale_npz.py`: PDF/SVG/PNG multiscale tensor summary figures from `score-footprints --output-multiscale-npz`.

## Matched Public-Label Benchmark Tables

After downloading public ATAC and TF-binding label data, turn scored prediction intervals into the standard `label`, `score`, `method`, `tf`, and `cell` table consumed by the metric and calibration scripts:

```bash
python benchmarks/scripts/build_label_overlap_benchmark.py \
  --predictions benchmarks/results/ctcf_reranked_sites.bed \
  --labels-bed data/public/labels/ctcf_chip_peaks.bed \
  --score-col rank_score \
  --min-overlap-bp 1 \
  --method fp-tools-reranked \
  --tf CTCF \
  --cell K562 \
  --metadata-cols name motif_family \
  --out benchmarks/results/ctcf_labeled_predictions.tsv
```

The output can be passed directly to `compute_binary_metrics.py`, `compute_calibration.py`, and the paper figure scripts.

## End-to-End Benchmark Result Folder

After creating one or more labeled prediction TSVs, run the summary pipeline to create a reproducible result folder with combined predictions, metrics, calibration summaries, optional bootstrap CIs, and manuscript-ready figure panels:

```bash
python benchmarks/scripts/run_benchmark_pipeline.py \
  --predictions benchmarks/results/ctcf_labeled_predictions.tsv benchmarks/results/irf1_labeled_predictions.tsv \
  --outdir benchmarks/results/public_tfbs_benchmark \
  --bootstrap 1000 \
  --bins 10 \
  --title "fp-tools public TFBS benchmark"
```

The figure outputs are written under `<outdir>/figures/` as PDF, SVG, and PNG files, and `<outdir>/benchmark_run_summary.md` lists every generated artifact.

## Motif-Removal Recovery Benchmark

Use this scaffold after generating candidate, model, or reranked predictions. It simulates removing a motif ID or motif family from the known-motif catalog, zeroes the strict motif baseline for those sites by default, and compares recovery scores in the same metric/figure path used by the main benchmark.

```bash
python benchmarks/scripts/build_motif_removal_benchmark.py \
  --predictions benchmarks/results/ctcf_site_predictions.tsv \
  --remove-col motif_family \
  --remove-values CTCF \
  --baseline-score-col motif_score \
  --recovery-score-cols rank_score binding_probability candidate_score \
  --out-long benchmarks/results/ctcf_motif_removed_predictions.tsv \
  --out-summary benchmarks/results/ctcf_motif_removed_summary.tsv

python benchmarks/scripts/compute_binary_metrics.py \
  --predictions benchmarks/results/ctcf_motif_removed_predictions.tsv \
  --score-col score \
  --group-cols removal_target method tf cell \
  --out benchmarks/results/ctcf_motif_removed_metrics.tsv

python paper/scripts/plot_benchmark_panels.py \
  --metrics benchmarks/results/ctcf_motif_removed_metrics.tsv \
  --out-prefix paper/figures/figure_motif_removal_ctcf
```


Add bootstrap confidence intervals to a binary metric run when preparing paper tables:

```bash
python benchmarks/scripts/compute_binary_metrics.py \
  --predictions benchmarks/results/ctcf_test_predictions.tsv \
  --score-col binding_probability \
  --group-cols tf cell method \
  --bootstrap 1000 \
  --seed 2026 \
  --out benchmarks/results/ctcf_metrics.tsv \
  --out-bootstrap benchmarks/results/ctcf_metric_ci.tsv
```

## Calibration Reports

For supervised TFBS prediction tables with probability-like scores, compute reliability bins and render paper-ready calibration panels:

```bash
python benchmarks/scripts/compute_calibration.py \
  --predictions benchmarks/results/ctcf_test_predictions.tsv \
  --score-col binding_probability \
  --group-cols tf cell method \
  --bins 10 \
  --out-bins benchmarks/results/ctcf_calibration_bins.tsv \
  --out-summary benchmarks/results/ctcf_calibration_summary.tsv

python paper/scripts/plot_calibration_panels.py \
  --bins benchmarks/results/ctcf_calibration_bins.tsv \
  --summary benchmarks/results/ctcf_calibration_summary.tsv \
  --out-prefix paper/figures/figure_ctcf_calibration
```
