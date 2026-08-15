# Benchmark Scripts

Benchmark and validation helpers:

- `build_encode_manifest.py`: query ENCODE and write a public-data manifest without downloading files.
- `download_manifest.py`: resumable downloads plus checksum and path reports.
- `compute_binary_metrics.py`: AUROC, AUPRC, recall@FDR, Brier score summaries, and optional bootstrap confidence intervals from scored labels.
- `compute_calibration.py`: reliability-bin, expected calibration error, maximum calibration error, and Brier summaries from probability-like predictions.
- `build_label_overlap_benchmark.py`: convert scored BED-like prediction intervals plus ChIP/CUT&RUN label BEDs into metrics-ready binary label/score tables.
- `build_motif_removal_benchmark.py`: create long-form motif-removal recovery benchmark tables from baseline, motif-free, supervised, or reranked site scores.
- `run_benchmark_pipeline.py`: combine labeled prediction TSVs, compute metrics/calibration/bootstrap summaries, and write PDF/SVG/PNG benchmark figures.
- `benchmark_footprint_kernel.py`: run `call-footprints` with the legacy and fast footprint kernels, measure wall time, and compare output bigWigs and candidate BEDs.
- `manuscript/scripts/plot_benchmark_panels.py`: PDF/SVG/PNG multi-panel benchmark figures for the BioMedInformatics manuscript.
- `manuscript/scripts/plot_calibration_panels.py`: PDF/SVG/PNG reliability curves and ECE panels.
- `manuscript/scripts/plot_multiscale_npz.py`: PDF/SVG/PNG multiscale tensor summary figures from `call-footprints --output-multiscale-npz`.

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

## Footprint Kernel Speed and Consistency

Use this helper after changing the Cython footprint-scoring path. It runs the same corrected bigWig and BED regions through `call-footprints --footprint-kernel legacy` and `--footprint-kernel fast`, then writes timing and output-difference summaries.

```bash
python benchmarks/scripts/benchmark_footprint_kernel.py \
  --signal test_data/Bcell_corrected.bw \
  --regions test_data/merged_peaks.bed \
  --outdir /tmp/fptools_footprint_kernel_benchmark \
  --cores 1
```

The output directory contains the two bigWigs, two candidate BEDs, command logs, `kernel_benchmark_summary.tsv`, and `kernel_benchmark_summary.json`. For the full B-cell fixture used during development, the fast kernel ran about 3.6x faster than the legacy kernel while keeping bigWig differences below `2e-5` absolute error and candidate BED overlap above 99.99%.

To carry the comparison through motif matching and differential footprinting, provide a second corrected signal plus genome and motif inputs:

```bash
python benchmarks/scripts/benchmark_footprint_kernel.py \
  --signal test_data/Bcell_corrected.bw \
  --workflow-second-signal test_data/Tcell_corrected.bw \
  --workflow-first-name Bcell \
  --workflow-second-name Tcell \
  --workflow-first-condition Bcell \
  --workflow-second-condition Tcell \
  --regions test_data/merged_peaks.bed \
  --genome test_data/genome.fa.gz \
  --motifs test_data/individual_motifs/MA0050.2.jaspar \
  --outdir /tmp/fptools_footprint_kernel_workflow_check \
  --cores 4
```

This optional mode writes separate fast and legacy `match-motifs` and `diff-footprints` folders and reports the final motif-table row and numeric differences in the summary files.

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

python manuscript/scripts/plot_benchmark_panels.py \
  --metrics benchmarks/results/ctcf_motif_removed_metrics.tsv \
  --out-prefix manuscript/figures/figure_motif_removal_ctcf
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

python manuscript/scripts/plot_calibration_panels.py \
  --bins benchmarks/results/ctcf_calibration_bins.tsv \
  --summary benchmarks/results/ctcf_calibration_summary.tsv \
  --out-prefix manuscript/figures/figure_ctcf_calibration
```
### Eight-cell ENCODE ATAC project

`run_encode_atac_project.py` prepares the complete two-replicate ENCODE ATAC
design for GM12878, HCT116, HepG2, IMR-90, K562, MCF-7, PC-3, and Panc1. It
streams released GRCh38 alignment BAMs one at a time, retains compact
differential-ready outputs, and independently verifies the exact 16-sample set.

```bash
.venv/bin/python benchmarks/scripts/run_encode_atac_project.py preflight
.venv/bin/python benchmarks/scripts/run_encode_atac_project.py run
.venv/bin/python benchmarks/scripts/run_encode_atac_project.py verify
```

### Seven-line ENCODE cancer-cell project

The pairwise Q95 runner builds the cancer-cell resource from 17 biological
replicates in A549, HCT116, HepG2, K562, MCF-7, PC-3, and Panc1. K562 and
HepG2 use the same three-replicate ENCODE experiments as the preserved
standalone demonstration. Each comparison uses its own union of every released
GRCh38 IDR-thresholded peak file from the two selected experiments.

```bash
.venv/bin/python benchmarks/scripts/run_encode_cancer_pairwise_q95.py preflight
.venv/bin/python benchmarks/scripts/run_encode_cancer_pairwise_q95.py download-inputs --workers 3
.venv/bin/python benchmarks/scripts/run_encode_cancer_pairwise_q95.py run --download --cores 16
.venv/bin/python benchmarks/scripts/run_encode_cancer_pairwise_q95.py verify
```

The workflow scales corrected cut-site tracks to the across-sample median q95,
scores footprints from those scaled tracks, and runs differential analysis
without an additional normalization step. The preserved K562-HepG2 payload is
copied exactly; the remaining 20 unordered comparisons use the same method.
Pair-local generated bigWigs are removed only after the compact report payload
and result table validate. Downloaded BAMs are retained for the separate Box
archive and are never deleted by the runner.

```bash
.venv/bin/python benchmarks/scripts/build_encode_cancer_q95_site.py build
.venv/bin/python benchmarks/scripts/build_encode_cancer_q95_site.py verify
```

The dependency-free browser reads compact gzip payloads using the same report
schema as `diff-footprints`. It contains every motif-level result and supports
both directions of all 21 unordered comparisons through two selectors.
