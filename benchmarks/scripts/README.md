# Benchmark Scripts

Benchmark and validation helpers:

- `build_encode_manifest.py`: query ENCODE and write a public-data manifest without downloading files.
- `download_manifest.py`: resumable downloads plus checksum and path reports.
- `compute_binary_metrics.py`: AUROC, AUPRC, recall@FDR, Brier score summaries, and optional row or genomic-block bootstrap confidence intervals from scored labels.
- `compute_calibration.py`: reliability-bin, expected calibration error, maximum calibration error, and Brier summaries from probability-like predictions.
- `build_label_overlap_benchmark.py`: convert scored BED-like prediction intervals plus ChIP/CUT&RUN label BEDs into metrics-ready binary label/score tables.
- `build_motif_removal_benchmark.py`: create long-form motif-removal recovery benchmark tables from baseline, motif-free, supervised, or reranked site scores.
- `run_benchmark_pipeline.py`: combine labeled prediction TSVs, compute metrics/calibration/bootstrap summaries, and write PDF/SVG/PNG benchmark figures.
- `benchmark_footprint_kernel.py`: run `call-footprints` with the legacy and fast footprint kernels, measure wall time, and compare output bigWigs and candidate BEDs.
- `build_footprint_detectability_atlas.py`: collapse repeated ENCODE and nutrient comparisons to independent biological contexts and rank expression-supported weak aggregate-shape hypotheses.
- `build_footprint_site_labels.py`: create summit-supported positive, distant negative, and explicitly indeterminate motif-site labels with optional matched controls.
- `downsample_bam_fragments.py`: create deterministic pair-preserving BAM depth subsets that remain nested for a fixed seed.
- `build_footprint_ablation_plan.py`: write the depth, correction, and method task matrices from the locked study specification.
- `run_footprint_ablation_plan.py`: execute the signal plan with dependency, resume, and expected-output checks.
- `summarize_footprint_ablation.py`: collapse depth randomizations and report depth plateaus and correction gains.
- `classify_footprint_failure_modes.py`: apply prespecified diagnostic rules to matched-label correction/scoring ablations without interpreting low scores as TF absence.
- `evaluate_footprint_promotion.py`: compare a frozen candidate with the current method under the prespecified development or locked-holdout gates.
- `evaluate_nutrient_footprint_replication.py`: apply local cross-cell-line, RNA, external recovery, and occupancy replication tiers.
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

## Footprint Detectability Atlas

Build the first-stage audit from the seven-line ENCODE comparison resource,
the three nutrient-stress projects, and the matched nutrient RNA table:

```bash
.venv/bin/python benchmarks/scripts/build_footprint_detectability_atlas.py \
  --outdir benchmarks/results/footprint_detectability_atlas
```

The workflow ranks motifs within each source report, then collapses the six
pairwise appearances of each ENCODE cell line and repeated nutrient control
rows to one independent biological context. This is required because the raw
score scales differ between the Q95-normalized ENCODE resource and the
nutrient projects. The output folder contains:

- `detectability_context_scores.tsv.gz`: one motif row per independent context;
- `detectability_motif_summary.tsv`: all motifs with cohort-specific rank and expression summaries;
- `detectability_candidates.tsv`: expression-supported motifs consistently in the configured low percentile;
- `detectability_input_manifest.tsv`: portable paths, sizes, and SHA-256 hashes for every input;
- `detectability_atlas.html`: searchable standalone report;
- `detectability_metadata.json`: schema version, thresholds, and context counts.

These candidates are weak aggregate-shape hypotheses. They are not evidence
that the TF is absent, that bias correction failed, or that ATAC-seq contains
no occupancy information. Those diagnoses require matched orthogonal
occupancy labels plus correction and depth ablations.

After producing one long-form metrics table with one row per
cell/TF/motif/method, classify the matched-label tasks:

```bash
.venv/bin/python benchmarks/scripts/classify_footprint_failure_modes.py \
  --metrics benchmarks/results/footprint_method_ablation_metrics.tsv \
  --current-method "fp-tools footprint" \
  --raw-method "raw footprint" \
  --out benchmarks/results/footprint_failure_modes.tsv
```

Required columns are `cell`, `tf`, `motif_id`, `method`, and `auroc`;
`auprc`, `positive_sites`, `coverage_pass`, `depth_plateau`,
`protein_supported`, `motif_ambiguous`, and `bias_residual` add stronger
diagnostic evidence. Existing tables using `motif` and `chip_positive_sites`
are accepted as aliases. The `atac_information_limited` status is emitted only
when adequate orthogonal labels, protein support, and a depth plateau are all
recorded.

## Locked Footprint Improvement Study

Validate the preregistered cells, chromosome splits, depth series, method arms,
diagnostic thresholds, promotion gates, and external nutrient datasets:

```bash
.venv/bin/python benchmarks/scripts/validate_footprint_study.py
```

Create motif-site labels from a motif BED and an IDR peak file. Only sites
inside a peak and within the configured distance of its summit are positive.
Sites close to a peak without summit support remain indeterminate.

```bash
.venv/bin/python benchmarks/scripts/build_footprint_site_labels.py \
  --sites data/public/processed/K562/CTCF_sites.bed \
  --chip-peaks data/public/raw/encode/K562.CTCF.narrowPeak.gz \
  --features data/public/processed/K562/CTCF_site_features.tsv \
  --match-columns motif_score accessibility gc mappability tss_distance regulatory_class \
  --out benchmarks/results/footprint_detectability_v1/K562_CTCF_labels.tsv \
  --matched-out benchmarks/results/footprint_detectability_v1/K562_CTCF_matched.tsv \
  --indeterminate-out benchmarks/results/footprint_detectability_v1/K562_CTCF_indeterminate.tsv
```

The ablation sample TSV has five required columns: `sample`, `cell`, `bam`,
`peaks`, and the number of usable `fragments`. Build and inspect the executable
plan before starting large jobs:

```bash
.venv/bin/python benchmarks/scripts/build_footprint_ablation_plan.py \
  --samples benchmarks/results/footprint_detectability_v1/ablation_samples.tsv \
  --genome data/public/reference/hg38.fa.gz \
  --blacklist data/public/reference/hg38-blacklist.v2.bed \
  --outdir benchmarks/results/footprint_detectability_v1/ablation \
  --cores 16 \
  --check-paths

.venv/bin/python benchmarks/scripts/run_footprint_ablation_plan.py \
  --plan benchmarks/results/footprint_detectability_v1/ablation/ablation_commands.tsv \
  --dry-run
```

The same query-name hash and seed are reused across correction arms; increasing
depths for one seed are nested. The validated `fragments` count in the sample
table is passed to the downsampler so each large BAM needs only one streaming
pass; omit `--available-fragments` when using the downsampler directly if the
count is not already known. The plan includes raw signal, PWM and DWM bias
models, and a full-depth bias model reused at lower depths. It writes a separate
evaluation matrix for fp-tools and the locked comparator methods.

Summarize site-label metrics only after all methods have been evaluated on the
same sites:

```bash
.venv/bin/python benchmarks/scripts/summarize_footprint_ablation.py \
  --metrics benchmarks/results/footprint_detectability_v1/ablation_metrics.tsv \
  --outdir benchmarks/results/footprint_detectability_v1/diagnostics

.venv/bin/python benchmarks/scripts/compute_binary_metrics.py \
  --predictions benchmarks/results/footprint_detectability_v1/site_predictions.tsv \
  --group-cols cell tf method \
  --block-cols chrom \
  --bootstrap 1000 \
  --out benchmarks/results/footprint_detectability_v1/site_metrics.tsv \
  --out-bootstrap benchmarks/results/footprint_detectability_v1/site_metric_ci.tsv
```

Candidate development uses K562 and HepG2 only. After its code and parameters
are frozen, unlock the MCF-7, A549, HCT116, and Panc1 holdout exactly once:

```bash
.venv/bin/python benchmarks/scripts/evaluate_footprint_promotion.py \
  --metrics benchmarks/results/footprint_detectability_v1/frozen_method_metrics.tsv \
  --candidate fp-tools-candidate \
  --baseline fp-tools \
  --negative-controls benchmarks/results/footprint_detectability_v1/naked_dna_false_positives.tsv \
  --split locked_holdout \
  --unlock-holdout \
  --outdir benchmarks/results/footprint_detectability_v1/promotion
```

Evaluate the label-free footprint/PWM percentile-fusion hypothesis in the same
order. The first command exposes development validation chromosomes only. Add
`--unlock-development-test` only after freezing the fusion rule, and add
`--unlock-holdout` once for the final cell-line holdout. The evaluator records
paired candidate-minus-current metrics and optional chromosome-block bootstrap
support; it never changes the package default.

```bash
.venv/bin/python benchmarks/scripts/evaluate_site_evidence_fusion.py \
  --site-scores benchmarks/results/public_chip_site_scores.tsv.gz \
  --outdir benchmarks/results/footprint_detectability_v1/evidence_fusion/development \
  --bootstrap 1000
```

The first frozen candidate tested on 2026-08-30 used
`1 - (1 - footprint percentile) * (1 - PWM percentile)`. It passed the
available K562/HepG2 development slice but failed the locked holdout gates,
principally because MYC regressed. It therefore remains an opt-in reranking
experiment and must not replace the default footprint score.

The nutrient application stays outside model training. Its locked external
resources are GSE144833 (SUIT-2 non-adapted, adapted, and reverse-adapted ATAC
and RNA) and GSE137034/GSE137031/GSE137032 (full, low, and no-arginine ATAC
with ATF4 and CEBPB occupancy). The strongest evidence tier from
`evaluate_nutrient_footprint_replication.py` requires local three-cell-line
directionality, RNA concordance, external stress and recovery, and concordant
orthogonal occupancy.

## End-to-End Benchmark Result Folder

After creating one or more labeled prediction TSVs, run the summary pipeline to create a reproducible result folder with combined predictions, metrics, calibration summaries, optional bootstrap CIs, and manuscript-ready figure panels:

```bash
python benchmarks/scripts/run_benchmark_pipeline.py \
  --predictions benchmarks/results/ctcf_labeled_predictions.tsv benchmarks/results/irf1_labeled_predictions.tsv \
  --outdir benchmarks/results/public_tfbs_benchmark \
  --bootstrap 1000 \
  --block-cols chrom \
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
