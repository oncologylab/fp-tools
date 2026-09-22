# Benchmark Scripts

Benchmark and validation helpers:

The locked naked-DNA negative-control runs are listed in
`benchmarks/manifests/naked_dna_gse164997.tsv`; keep the three runs separate so
false-footprint rates can be checked for replicate consistency before pooling.

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
- `downsample_bam_depth_matrix.py`: create all missing depth/seed BAM subsets in one source-BAM scan while preserving the same deterministic nested-fragment rule.
- `evaluate_functional_template_transfer.py`: measure shape-only same-cell, cross-cell TF, leave-TF-out family, and global functional-template ceilings without opening locked holdouts.
- `evaluate_functional_depth_matrix.py`: apply 10M-frozen TF shape settings across deterministic 10M/25M/50M/full matrices, classify depth limits, and render per-TF depth/aggregate panels.
- `combine_frozen_functional_depth_matrices.py`: integrity-check and combine separately completed frozen depth runs, allowing the full endpoint to be added without rebuilding completed 10M/25M/50M profiles.
- `evaluate_frozen_functional_policy.py`: apply immutable label-free TF models on test chromosomes without refitting and emit exact-site raw/DWM/parametric comparisons, bootstrap-banded profiles, block-bootstrap deltas, and auditable per-site scores.
- `evaluate_frozen_functional_information_ceiling.py`: fit evaluation-only, leave-chromosome-out elastic-net classifiers on raw/DWM/parametric signals, motif/accessibility covariates, and label-free functional PCs. The v2 report keeps raw cut geometry as an explicit guardrail and separates assay-limited, signal-combination-limited, and shape-model-limited TF failures.
- `evaluate_frozen_bias_shrinkage.py`: test whether partial subtraction of a frozen direct or control-calibrated parametric expectation preserves useful raw footprint geometry. It freezes global and TF-specific validation choices, applies them without refitting on test chromosomes, and reports raw/DWM comparisons, chromosome-block intervals, aggregate curves, and exact site scores.
- `summarize_frozen_parametric_tf_evidence.py`: integrity-check and collate chromosome-test, biological-replicate, depth/seed, naked-DNA, and supervised-ceiling evidence for every prespecified development TF, including underpowered tasks and an explicit raw-signal guardrail.
- `render_bias_shrinkage_before_after.py`: emit a concise one-page DWM-versus-frozen-parametric report only when the same TF has significant chromosome-block AUROC and AUPRC gains in at least two internal cell contexts and the exact frozen policy passes independent naked-DNA safety.
- `render_frozen_functional_before_after.py`: emit one-page DWM-versus-candidate PDFs only for eligible TFs with positive AUROC/AUPRC bootstrap intervals and a passing independent naked-DNA screen.
- `evaluate_shifted_atac_null_policy.py`: score frozen detectors after deliberately misaligning motif-centered cellular ATAC profiles, providing a label-free accessibility/nucleosome null that complements naked DNA.
- `evaluate_naked_dna_functional_policy.py`: score the frozen functional policy on naked-DNA profiles; use `--candidate-only` for an independent control replicate when no same-replicate DWM artifact exists.
- `freeze_functional_call_thresholds.py`: freeze per-TF candidate and DWM call thresholds from development matched negatives before any naked-DNA profile is read, retaining auditable site-level validation scores for calibrated ensembles.
- `evaluate_frozen_functional_consensus.py`: transform count and FDA scores against frozen matched-negative distributions, select a global operator on validation chromosomes with a one-standard-error/CTCF rule, and evaluate every predeclared blend on test chromosomes and independent naked DNA without refitting. Test evaluation can integrity-check a separate raw-score manifest and report all operators against raw and DWM on their exact common finite support.
- `summarize_frozen_detector_evidence.py`: integrity-check and collate the validation-frozen TF-detector result with raw/DWM test bootstraps, biological-replicate stability, depth/seed stability, and independent naked-DNA safety for every predeclared development TF.
- `evaluate_shape_guarded_functional_screen.py`: compare validation-frozen count, FDA, GP-hybrid, anchored, and residualized detector families on exact raw-matched test support while requiring central protection before a numerical gain is called a footprint gain. This is a secondary development diagnostic, never promotion evidence.
- `evaluate_frozen_functional_naked_dna.py`: apply serialized label-free models and their immutable thresholds to an independent naked-DNA library, retaining zero-cut sites as valid negative observations and reporting independent plus paired support.
- `evaluate_factorization_residual_safety.py`: freeze per-cell/TF residual cutoffs on chr16--18 matched negatives, then screen every residual on label-free naked-DNA replicate 2 before the parametric configuration is frozen.
- `freeze_safe_parametric_configuration.py`: combine development ranking, CTCF non-regression, and naked-DNA safety into the only v2 configuration accepted for test-chromosome scoring.
- `calibrate_dual_null_posteriors.py`: require posterior calls to exceed conservative per-TF thresholds from both naked DNA and motif-misaligned cellular ATAC.
- `apply_frozen_dual_null_thresholds.py`: apply already frozen per-TF dual-null thresholds to an independent score table without refitting or borrowing a comparator replicate.
- `ensemble_parametric_bias_models.py`: combine compatible seed fits into one checksummed geometric coefficient ensemble for prespecified crossed controls.
- `fit_combined_control_bias.py`: fit the locked naked-DNA-plus-mitochondrial grid, validate and reuse complete seed/model pairs after interruption, and fail closed on incomplete or checksum-invalid artifacts.
- `render_functional_aggregate_comparison.py`: render blinded DWM-versus-strand per-TF aggregate curves with bootstrap confidence bands and an explicit blinding key.
- `sample_label_free_motif_sites.py`: deterministically sample larger per-TF motif pools from raw motif scans on training chromosomes, record checksums, and refuse label-bearing inputs.
- `evaluate_strand_label_free_models.py`: fit positionwise-count, conditional-multinomial spline/GP, FDA, and hybrid mixtures on separate label-free strand artifacts and score only development validation labels.
- `derive_label_free_family_selection.py`: derive an integrity-checked count/FDA/hybrid family-specific winner policy from a completed validation grid without reading test labels.
- `analyze_strand_bias_factorial.py`: compare all four SELMA/log-linear and +4/−4/+4/−5 combinations at common detector settings, with per-TF model, shift, and interaction effects.
- `assemble_label_free_candidate_matrix.py`: place prior-free DWM shape metrics and strand-aware label-free metrics on one audited task/candidate table.
- `freeze_tf_dependent_label_free_policy.py`: freeze family routes with mean-gain, AUPRC, context non-regression, and unseen-family DWM fallback gates before holdout labels are opened.
- `select_count_models_by_unlabeled_likelihood.py`: choose a TF's bias/shift/window count mixture by marginal likelihood on held-out unlabeled training chromosomes, then open development labels only for post-selection evaluation.
- `build_footprint_ablation_plan.py`: write the depth, correction, and method task matrices from the locked study specification.
- `run_footprint_ablation_plan.py`: execute the signal plan with dependency, resume, and expected-output checks.
- `summarize_footprint_ablation.py`: collapse depth randomizations and report depth plateaus and correction gains.
- `classify_footprint_failure_modes.py`: apply prespecified diagnostic rules to matched-label correction/scoring ablations without interpreting low scores as TF absence.
- `evaluate_footprint_promotion.py`: compare a frozen candidate with the current method under the prespecified development or locked-holdout gates.
- `evaluate_site_evidence_fusion.py`: test a frozen label-free footprint/PWM site-ranking candidate in the locked chromosome/cell sequence.
- `evaluate_bigwig_site_scores.py`: extract base-resolution scores from correction/scoring ablation bigWigs at fixed ChIP-labeled motif centers.
- `discover_encode_chip_peaks.py`: select unperturbed replicate-aware GRCh38 IDR ChIP peaks for the locked TF tasks and extract the required motif subset.
- `build_encode_tf_site_matrix.py`: combine cell-specific motif scans with conservative ENCODE ChIP summit labels.
- `search_tf_footprint_models.py`: run staged per-TF correction, geometry, normalization, and symmetry searches using train/validation chromosome separation.
- `match_tf_sites_on_accessibility.py`: optimally match positive and negative motif sites on motif score and local raw ATAC coverage.
- `compare_frozen_tf_candidates.py`: compare frozen TF candidates with legacy scores on identical finite sites.
- `plot_frozen_tf_profiles.py`: plot matched ChIP-positive and ChIP-negative aggregate profiles for frozen candidates.
- `render_tf_before_after_report.py`: render concise one-page, paired legacy-versus-frozen TF reports with held-out ROC/PR curves, aggregate profiles, bootstrap intervals, replicate and optional naked-DNA evidence, and explicit research-only scope.
- `evaluate_tf_geometry_naked_dna.py`: freeze per-method score thresholds on development ChIP-negative sites and test an exact TF-specific geometry plus the conventional footprint score on naked-DNA motif sites, reporting each method's independent finite support and their paired overlap.
- `evaluate_tf_geometry_external_transfer.py`: apply one frozen development-cell TF geometry without retuning to external labeled cells, compare it with the conventional footprint score on identical sites, and report chromosome-block, covariate-residual, and replicate evidence.
- `render_tf_external_transfer_report.py`: render a concise paired aggregate and metric PDF for a frozen TF geometry transferred without retuning to two external cells.
- `summarize_tf_footprint_search.py`: apply prespecified site-count, matching-balance, detectability, and point-gain statuses to frozen tests.
- `evaluate_tf_correction_transfer.py`: hold frozen TF geometry fixed while transferring it across raw, PWM, and DWM signals to isolate correction sensitivity.
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

`downsample_bam_depth_matrix.py` now records source, BAM and index checksums
with an identity for each target's sampling metadata. It resumes only complete
matching targets and preserves their modification times. An interrupted target
or historical metadata without the required checksums fails explicitly;
previously frozen data can still be reused through their original verified
registries. The helper does not rewrite old metadata or accept a BAM merely
because the file exists.

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

For signal arms that already produced footprint-score bigWigs, create a small
manifest with `cell`, `method`, and `signal`, then score the exact same motif
centers without rerunning motif discovery:

```bash
.venv/bin/python benchmarks/scripts/evaluate_bigwig_site_scores.py \
  --sites benchmarks/results/footprint_detectability_v1/development_sites.tsv.gz \
  --signals benchmarks/results/footprint_detectability_v1/ablation_signals.tsv \
  --chromosomes chr17 chr18 \
  --baseline-method raw \
  --bootstrap 1000 \
  --outdir benchmarks/results/footprint_detectability_v1/correction_metrics
```

Metrics use the common finite site set across methods. With `--bootstrap`, the
script also reports paired method-minus-baseline confidence intervals using
chromosomes as resampling blocks.

For TF-specific geometry experiments, first discover the locked ENCODE ChIP
resources and build the cell-specific label matrix, then run the staged search.
The search never reads test chromosomes unless `--profiles-only
--chromosome-splits test` is explicitly requested after candidates are frozen.

```bash
.venv/bin/python benchmarks/scripts/discover_encode_chip_peaks.py \
  --study benchmarks/manifests/footprint_detectability_v1.spec.json \
  --outdir benchmarks/results/footprint_detectability_v1/encode_labels \
  --motif-database data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
  --download

.venv/bin/python benchmarks/scripts/search_tf_footprint_models.py \
  --sites benchmarks/results/footprint_detectability_v1/encode_tf_site_labels.tsv.gz \
  --signals benchmarks/results/footprint_detectability_v1/correction_signals.tsv \
  --study benchmarks/manifests/footprint_detectability_v1.spec.json \
  --outdir benchmarks/results/footprint_detectability_v1/per_tf_search
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

### Functional footprint research

`evaluate_parametric_bias.py` runs the enzyme-bias stage before functional
model selection. It compares +4/-5 with +4/-4, fits the SELMA-style 10-mer and
81-bp log-linear models, and evaluates sample-specific, cross-cell pooled, and
pooled-plus-adapted fits using held-out control likelihood. A configuration is
not retained unless it beats the uniform conditional control on validation
chromosomes. Control-window caches include every extraction parameter and
input file identity.

```bash
.venv/bin/python benchmarks/scripts/evaluate_parametric_bias.py \
  --study benchmarks/manifests/footprint_functional_v1.spec.json \
  --sample K562=<K562-coordinate-sorted.bam> \
  --sample HepG2=<HepG2-coordinate-sorted.bam> \
  --genome <hg38.fa> \
  --peaks <merged-peaks.bed> \
  --blacklist <hg38-blacklist.bed> \
  --source gc_matched_low_signal_nonpeak \
  --outdir benchmarks/results/footprint_parametric_v1/nonpeak
```

Use matched unfiltered BAMs with `--source mitochondrial`; the fixed chrM
windows are deterministically divided into development and validation
partitions. Training-depth arms use binomial thinning with recorded seeds. The
benchmark writes strand-aligned cut motifs, likelihood/calibration metrics,
safe model artifacts, and the two configurations eligible for the downstream
functional screen. It never reads the locked test chromosomes.
After the initial pooling comparison, use `--pooled-only` for the larger
regularization, depth, and seed grid so rejected sample-specific fits are not
recomputed unnecessarily.

Before a retained model enters the functional screen,
`evaluate_bias_motif_leakage.py` measures its sequence-only response at every
eligible development TF motif against a nearby local-sequence control. It
reports response curves, bootstrap uncertainty, broad center/flank effects,
and cross-cell/model concordance without reading ChIP labels.

```bash
.venv/bin/python benchmarks/scripts/evaluate_bias_motif_leakage.py \
  --study benchmarks/manifests/footprint_functional_v1.spec.json \
  --motif-sites K562=<K562-unlabeled-motif-sites.tsv.gz> \
  --motif-sites HepG2=<HepG2-unlabeled-motif-sites.tsv.gz> \
  --model K562:loglinear81=<K562-model.npz> \
  --model HepG2:loglinear81=<HepG2-model.npz> \
  --genome <hg38.fa> \
  --outdir benchmarks/results/footprint_parametric_v1/motif_response
```

A flag means the bias model responds reproducibly to motif sequence and needs
review; it is not automatically called occupancy leakage, because genuine Tn5
sequence preference can overlap a TF motif. Promotion additionally requires
control-likelihood gain, cross-cell transfer, naked-DNA behavior, and no loss
of held-out TF signal.

`evaluate_functional_footprints.py` compares label-free spline, FDA,
Gaussian-process-equivalent, and hybrid footprint models on identical matched
motif sites. Its deployable candidates train on motif pools that contain no
ChIP/label columns; labels are used only for evaluation and the separately
reported supervised information ceiling. Test scoring requires
`--unlock-test` after the frozen-model table has been reviewed.

```bash
.venv/bin/python benchmarks/scripts/evaluate_functional_footprints.py \
  --study benchmarks/manifests/footprint_functional_v1.spec.json \
  --development-sites <development-sites.tsv.gz> \
  --unlabeled-sites K562=<K562-motif-sites.tsv.gz> \
  --unlabeled-sites HepG2=<HepG2-motif-sites.tsv.gz> \
  --tracks <raw-expected-track-manifest.tsv> \
  --outdir benchmarks/results/footprint_functional_v1/development
```

For a retained parametric model, add a `parametric_model` row beside the raw
track in the track manifest and supply `--genome`. Expected site profiles are
then derived from the frozen model and retain each site's observed total:

```text
cell   model            track              signal
K562   loglinear81_v1   raw                <K562-uncorrected.bw>
K562   loglinear81_v1   parametric_model   <K562-loglinear81.npz>
```

`build_strand_functional_profiles.py` is the strand-aware research path. It
extracts forward and reverse cuts directly from the BAM, predicts each strand
with the frozen bias model, and correctly reverses coordinates *and swaps cut
strands* for reverse-oriented motifs. Its NPZ contains plus/minus observed and
expected profiles, combined signed-deviance residuals, and shared and
antisymmetric strand residuals.

```bash
.venv/bin/python benchmarks/scripts/build_strand_functional_profiles.py \
  --sites <motif-sites.tsv.gz> \
  --cell K562 \
  --bam <coordinate-sorted.bam> \
  --genome <hg38.fa> \
  --bias-model <retained-bias-model.npz> \
  --read-shift 4,-4 \
  --out-prefix benchmarks/results/footprint_functional_v1/K562.strand_profiles
```

The artifact uses only numeric/string arrays with a checksummed JSON sidecar;
it does not use arbitrary pickle. ChIP labels are not used to construct any
profile channel.

`evaluate_strand_functional_templates.py` compares combined, shared-strand,
antisymmetric-strand, and multichannel shape detectors for every eligible TF.
It reports same-cell information ceilings separately from cross-cell TF and
leave-TF-out family/global transfer, and refuses artifacts that do not certify
label-free profile construction.

`evaluate_differential_functional_profiles.py` tests the entire condition
difference curve from replicate-aware strand-profile artifacts. A manifest
contains `sample`, `condition`, `replicate`, `profiles_npz`, and `sites_tsv`.
The result includes pointwise and simultaneous 95% bands, a global functional
permutation test, replicate-profile consistency, multiple-testing adjustment,
and changes in depth, width, shoulders, asymmetry, and periodicity.

```bash
.venv/bin/python benchmarks/scripts/evaluate_differential_functional_profiles.py \
  --manifest <nutrient-strand-profile-manifest.tsv> \
  --contrast stress,control \
  --channel combined \
  --outdir benchmarks/results/footprint_functional_v1/nutrient_frozen
```

At least two biological replicates per condition are required. Nutrient-stress
profiles are not evaluated until model choices and thresholds are frozen.

The output records input hashes, frozen model choices, per-site probabilities,
per-TF metrics, and explicit assay-, bias-, correction-, shape-, and
power-limitation diagnoses. Large profile caches and fitted benchmark outputs
remain ignored under `benchmarks/results/`.

It also writes `functional_aggregate_profiles.tsv.gz` with deterministic 95%
bootstrap bands, `functional_profile_descriptors.tsv` with center depletion,
width, shoulder distance, asymmetry, and periodicity, and
`functional_phenotype_clusters.tsv`. The clustering is based on the observed
functional shape rather than motif-family labels. The multi-page
`functional_aggregate_panels_blinded.pdf` compares matched groups across every
correction; its group key is stored separately so visual review can be done
before labels are revealed.

`evaluate_functional_promotion.py` is the fail-closed final audit. It requires
paired locked-holdout metrics plus functional-separation, naked-DNA,
bias-motif-response, replicate/seed/depth stability, runtime/memory/model-size,
uncertainty-coverage, and GP-versus-spline evidence. Missing evidence is a
failed gate. Its original specification required GM12878/IMR-90; those cells
have now been examined and provide diagnostic transfer evidence only. The
current frozen parametric study reserves SK-N-SH/GM23338 as new promotion
holdouts. Do not use this historical evaluator's original holdout designation
to claim new untouched validation.

```bash
.venv/bin/python benchmarks/scripts/evaluate_functional_promotion.py \
  --study benchmarks/manifests/footprint_functional_v1.spec.json \
  --metrics <locked-functional-metrics.tsv> \
  --candidate <correction:detector> \
  --baseline DWM:spline \
  --task-split locked_holdout \
  --unlock-holdout \
  --descriptors <profile-descriptors.tsv> \
  --negative-controls <naked-dna-fpr.tsv> \
  --resources <resource-metrics.tsv> \
  --uncertainty <uncertainty-coverage.tsv> \
  --stability <replicate-seed-depth-stability.tsv> \
  --leakage <bias-motif-review.tsv> \
  --complexity <gp-versus-spline.tsv> \
  --outdir benchmarks/results/footprint_functional_v1/promotion
```

### Completed ChromBPNet reference audit

`audit_chrombpnet_comparison.py` reads completed checksum-valid comparisons;
it never retrains a model or reruns prediction/evaluation. Run it from the
research worktree, pointing `--data-root` to the checkout containing the ignored
data. Relative study paths can resolve from the research worktree. Use a new
output directory for each audit.

```bash
python benchmarks/scripts/audit_chrombpnet_comparison.py \
  --data-root /path/to/data-checkout \
  --comparison-manifest /path/to/data-checkout/benchmarks/results/frozen_parametric_factorization_v1_run2/chrombpnet_comparison/HepG2_complete/chrombpnet_comparison_manifest.json \
  --comparison-manifest /path/to/data-checkout/benchmarks/results/frozen_parametric_factorization_v1_run2/chrombpnet_comparison/K562_complete/chrombpnet_comparison_manifest.json \
  --outdir benchmarks/results/chrombpnet_development_audit
```

The output includes per-TF metrics, aggregate PDF panels for every task,
signed shape descriptors, positive/negative retention and covariate balance,
and a checksummed manifest. Genomic coordinates, strand and TF identify sites;
the historical `site_id` column contains motif names and is not a unique locus.
The audit accepts only previously opened K562/HepG2 development comparisons.

The original HepG2 and K562 comparisons use different dispersions. DWM in these
comparisons means total-matched DWM expectations followed by signed deviance
and fixed geometry; it is not the production corrected-track footprint command.
ChromBPNet regulatory/no-bias arms score predicted profiles, whereas the three
bias arms score residual observed cuts. Their roles and aggregate units are
reported separately. The earlier factorization-residual safety table does not
certify direct parametric residuals. See `benchmarks/CHROMBPNET_AUDIT_20260908.md`
for the current scientific interpretation.

### Exact direct-bias candidate safety

`audit_direct_bias_safety.py` freezes code, recorded inputs and per-cell
dispersion, calibrates thresholds on chr16–18 negatives, then scores only
the existing naked-DNA replicate-2 library. The supplied plan fixes both
comparison manifests and input registries by checksum. It rejects final
replicate 3 and new holdout cells. Completed checksum-valid runs are verified
and skipped. Use the same code revision when verifying an existing freeze.

```bash
python benchmarks/scripts/audit_direct_bias_safety.py \
  --plan benchmarks/manifests/chrombpnet_direct_bias_safety_v1.json \
  --data-root /path/to/data-checkout \
  --outdir benchmarks/results/direct_bias_safety

python benchmarks/scripts/audit_direct_bias_safety_coverage.py \
  --safety-manifest benchmarks/results/direct_bias_safety/safety_manifest.json \
  --data-root /path/to/data-checkout \
  --outdir benchmarks/results/direct_bias_safety_coverage
```

Both stages are required. The archived DWM control track covers CTCF regions;
zero-filled profiles elsewhere cannot establish relative safety. The coverage
stage reuses existing scores and thresholds and compares observed count shapes
after removing a constant normalization scale. Read its `qualified_safety.tsv`
and `coverage_manifest.json`, retaining absolute and DWM-comparable denominators.
The first stage's provisional screen field does not certify DWM coverage.
These are development candidate checks; neither establishes final untouched
safety, no occupancy leakage, or calibrated functional uncertainty.

### One-page direct-bias comparison report

The report renderer requires the completed comparison, its shape/balance audit,
and the coverage-qualified direct-operator safety manifest. Legacy factorization
safety TSVs are rejected. It checks exact comparison identity, significant
AUROC/AUPRC gains over DWM and raw, canonical central protection, and a passing
candidate FPR screen. Safety values appear only on the candidate's metric row.

```bash
python benchmarks/scripts/render_chrombpnet_parametric_report.py \
  --comparison-manifest /path/to/data-checkout/benchmarks/results/frozen_parametric_factorization_v1_run2/chrombpnet_comparison/HepG2_complete/chrombpnet_comparison_manifest.json \
  --safety-manifest benchmarks/results/direct_bias_safety_coverage/coverage_manifest.json \
  --audit-manifest benchmarks/results/chrombpnet_development_audit/audit_manifest.json \
  --data-root /path/to/data-checkout \
  --cell HepG2 --tf CTCF \
  --outdir benchmarks/results/reports/HepG2_CTCF_direct_bias
```

The new output folder contains a one-page A4 PDF, preview, metrics TSV,
limitations README, portable provenance manifest and `SHA256SUMS`. Existing
reports are never overwritten. The report states the actual cell, dispersion,
covariate imbalance and safety denominators; it makes no equivalence claim
about ChromBPNet and does not imply final holdout safety or general promotion.

### Direct-bias development stability

`audit_direct_bias_stability.py` uses the exact recorded direct-bias operator
on previously examined chr19–22/X. The locked input registries contain five
10M/25M/50M seeds of replicate 1 and all three full-depth biological replicates
per development cell. The old matrices use tuning chromosomes; the script
reuses the two existing 10M test profiles and extracts only the missing test
profiles from checksum-verified BAMs and cached strand-specific propensities.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/scripts/audit_direct_bias_stability.py \
  --plan benchmarks/manifests/direct_bias_stability_v1.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/direct_bias_stability
```

The configuration, code and inputs freeze before new scoring. A per-cell
stage freezes sites with positive DWM expectations in every run and retains
complete existing covariate-matched pairs without rematching. Zero observed
counts remain in the fixed cohorts. Reports distinguish the original fixed
full-depth DWM expectation from each sample/depth-matched DWM expectation.
All runs have per-TF metrics, coverage and signed aggregate diagnostics;
chromosome-block metric intervals cover every full replicate and the 10M
seed-2026 reference. Five randomizations of one library are not five biological
replicates. Replicates 2 and 3 do not yet have a depth/seed factorial panel.

Completed stages are checksum-verified and reused. `--verify-only` checks the
complete result without scoring; after code changes, add `--code-ref` with the
immutable calculation commit. No new holdout labels, final naked-DNA control,
model refitting or production integration is allowed by this diagnostic.

`render_direct_bias_stability_report.py` creates a one-page CTCF summary only
after all 36 point-estimate directions, all six full-replicate intervals,
canonical protection and exact-operator candidate safety pass its report
checks. It preserves the observed-count balance drift, separate dispersions,
sparse-control caveat and missing general-promotion evidence. Reports use a
new immutable directory, portable source identities and SHA-256 checksums.

```bash
python benchmarks/scripts/render_direct_bias_stability_report.py \
  --stability-manifest benchmarks/results/direct_bias_stability/stability_manifest.json \
  --safety-manifest benchmarks/results/direct_bias_safety_coverage/coverage_manifest.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/reports/CTCF_direct_bias_stability
```

### Completing biological-replicate depth inputs

`prepare_direct_bias_replicate_depths.py` prepares the 60 missing combinations
for HepG2/K562 replicates 2 and 3: 10M/25M/50M with seeds 2026–2030. The
existing 36-run experiment and its manifests remain unchanged. Preparation
does not evaluate ChIP labels, refit the candidate or open a final control.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/scripts/prepare_direct_bias_replicate_depths.py \
  --plan benchmarks/manifests/direct_bias_replicate_depths_v1.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/direct_bias_replicate_depth_inputs \
  --phase freeze

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python benchmarks/scripts/prepare_direct_bias_replicate_depths.py \
  --plan benchmarks/manifests/direct_bias_replicate_depths_v1.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/direct_bias_replicate_depth_inputs \
  --phase all
```

The freeze verifies the four full-depth source BAMs, their indexes, the genome,
peaks, blacklist and the archived replicate-1 correction command. It records
code, the production source-tree identity, native extensions and versions.
Counting uses the original fragment eligibility rule through the pinned
`pysam` SAM filter, tested against all 4,096 flag combinations; it does not
introduce a MAPQ filter. The bit operations follow the
[SAMtools filter-expression specification](https://www.htslib.org/doc/samtools.html#FILTER_EXPRESSIONS).
The process count is at most two and a kernel lock prevents duplicate launches.

Stages `count`, `sample` and `correct` can be run individually. Sampling
preserves query-name pair decisions and nested depths, records actual selected
fragment counts, and rejects a target above the available depth. Each subset
gets a freshly estimated conventional DWM using the recorded replicate-1
command; the full-depth reused bias model is not substituted. Completion
requires all four tracks, the DWM artifact and diagnostic PDF with checksums.
Interrupted output directories fail explicitly and remain available for review.

`--phase verify --code-ref COMMIT` checks the final 60-row input registry,
all source files and nested stage artifacts without preparing data again.
Completed checksum-valid stages are reused. The generated registry is input
preparation evidence; a subsequent 96-run evaluation must freeze positive
expectation coverage across every run, preserve complete existing matched
pairs and retain zero observed counts before scoring any new runs.

### Complete 96-run development stability

`extend_direct_bias_stability.py` combines the archived 36-run experiment
with the 60 new inputs from `prepare_direct_bias_replicate_depths.py`. Its
plan is `manifests/direct_bias_stability_96_v1.json`. It requires the completed
preparation manifest and verifies source, sampling, correction and output
checksums, including the exact frozen DWM command. A temporary BAM or an
incomplete correction does not satisfy the prerequisite.

Both cells' positive expected coverage and complete existing matched pairs
are frozen before any new extraction or scoring. All 36 archived observed
and expected matrices are reused after checking coordinates, labels and the
matching DWM profile. Their metrics are reused only if the sites and both
cohort memberships are identical. If coverage changes, metrics and aggregate
bands are recomputed on the new support. Zero observed counts remain ties.
The existing candidate, geometry, dispersions and metric settings remain fixed.

The completed matrix has three biological replicates at each depth, five
nested seeds per replicate at 10M/25M/50M and one full-depth endpoint per
replicate. Principal chromosome-block intervals cover all six full endpoints
and all six 10M seed-2026 runs. The outputs include per-TF metrics, aggregate
profiles and pointwise bands, canonical-shape audits, coverage/balance tables,
a multipage PDF and nested completion manifests.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python benchmarks/scripts/extend_direct_bias_stability.py \
  --plan benchmarks/manifests/direct_bias_stability_96_v1.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/direct_bias_stability_96_20260908_v1
```

The plan pins the evaluator code; `--plan-sha256` additionally binds a queued
launch to the reviewed plan bytes. `--verify-only --code-ref COMMIT` verifies
historical code, all recorded input files, both supports and all 96 nested
stages without rerunning them. A kernel lock prevents simultaneous evaluators
from writing the same output. Keep the code and plan unchanged while queued
or running. These are post-hoc development results on examined chr19–22/X;
the new holdout labels and final naked-DNA control are outside this workflow.

The checked-in `direct_bias_stability_96_queue_v1.json` waits for the current
dated preparation output, then launches the pinned evaluator. It expects
the original data checkout at the sibling `../fp-tools` path. After verifying
the existing preparation PID, start one controller with that PID:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python benchmarks/scripts/run_background_experiment_queue.py \
  --plan benchmarks/manifests/direct_bias_stability_96_queue_v1.json \
  --watch-pid PREPARATION_PID --poll-seconds 30 --timeout-seconds 0
```

Do not start a second preparation or evaluator while an existing controller
is active. The queue fails if its watched preparation exits before producing
the prerequisite. Run the evaluator's `--verify-only` check before interpreting
the completed results; a queue completion record alone does not audit every
nested scientific input.

### Durable execution of long research jobs

`run_logged_experiment.py` supervises one checksum-pinned Linux command. Keep
each attempt's JSON specification, logs and receipts under ignored benchmark
results. The specification has schema `fp-tools-logged-experiment-job-v1`,
`name`, a `command` argument array starting with an absolute executable,
absolute `cwd` and `outdir`, and an `environment` mapping of explicit numerical
thread settings or `PYTHONPATH`. It never records the full inherited environment.

```bash
python benchmarks/scripts/run_logged_experiment.py \
  --spec /absolute/path/to/attempt/job.json \
  --spec-sha256 <sha256-of-job-json>
```

For detached execution, start that argument array with `subprocess.Popen`,
`start_new_session=True`, `stdin=subprocess.DEVNULL`, and stdout/stderr directed
to a file. The supervisor retains a byte-identical `specification.json`,
`launch.json`, `stdout.log` and `exit.json`. Launch evidence records both
supervisor and child PIDs, start ticks and boot identity. Monitor the actual
child identity when connecting a dependent queue. Reusing an attempt directory
is rejected before another command starts; preserve its evidence and choose
a new attempt directory after confirming that its processes are terminal.

SIGTERM/SIGINT sent to the supervisor are forwarded to its child group. If
writing a receipt fails, it stops the unrecorded child. Exit status is process
evidence only: verify the experiment's complete nested scientific manifests
before interpreting results. An interrupted sampling run still requires
review and preservation of its partial files, followed by the unchanged
frozen preparation command; a supervisor cannot make partial BAMs complete.

### Sensitivity to current-depth cut balance

`audit_direct_bias_depth_balance.py` uses only the completed 36-run archive.
The original matched cohorts balance the reference library; observed-cut
balance can drift at other depths or in other biological replicates. This
secondary diagnostic keeps an original pair only when its current log1p cut
totals differ by at most 0.2 standard deviations, using the already frozen
per-cell/TF scales from chr1–15. It does not find new partners or choose a
caliper from the resulting metrics. All 36 selections are frozen before new
scoring; every source TF appears in the retention table, including zero-pair
and underpowered outcomes.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python benchmarks/scripts/audit_direct_bias_depth_balance.py \
  --plan benchmarks/manifests/direct_bias_depth_balance_v1.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/direct_bias_depth_balance_20260908_v1
```

Use `--freeze-only` to stop after configuration and cohort selection, or
`--verify-only --code-ref COMMIT` to verify completed inputs and nested outputs.
The candidate, geometry, dispersions and statistical settings are inherited
unchanged from the archived study. Counts, expectations and primary metrics
are reused; only the secondary comparisons are newly evaluated. Outputs keep
the primary and secondary evidence roles explicit and report retained-pair
counts and current/reference covariate balance alongside performance and shape.
This is post-hoc conditioning on observed ATAC counts, which can introduce
selection effects. It does not establish a causal correction benefit,
calibrated uncertainty, new biological replication or promotion evidence.

`render_direct_bias_depth_balance.py` verifies the complete source graph and
all nested secondary outputs before drawing all TFs on shared vertical axes:

```bash
python benchmarks/scripts/render_direct_bias_depth_balance.py \
  --result benchmarks/results/direct_bias_depth_balance_20260908_v1/balance_manifest.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/direct_bias_depth_balance_20260908_v1/plots \
  --calculation-commit e8986ed
```

The figure set includes underpowered and non-estimable tasks, a multipage
PDF, two CTCF PNG previews and a checksum manifest. See
`benchmarks/DIRECT_BIAS_DEPTH_BALANCE_20260908.md` for the completed result.

### Direct-bias sequence inference resources

`benchmark_direct_bias_inference.py` measures the archived production Cython
DWM scorer and the frozen parametric scorer at identical genomic positions.
The plan selects 512 loci per development cell from existing qualified
coverage using coordinate hashes. ChIP labels and cut counts are not inputs
to timing. Six repetitions alternate method order in fresh processes pinned
to one CPU, with one thread and an identical sequence-input cache policy.

```bash
python benchmarks/scripts/benchmark_direct_bias_inference.py \
  --plan benchmarks/manifests/direct_bias_inference_benchmark_v1.json \
  --data-root /path/to/original/data-checkout \
  --outdir benchmarks/results/direct_bias_inference
```

Use `--freeze-only` to freeze inputs and code before timing. Existing complete
workers are verified and reused; partial workers fail explicitly. Use
`--verify-only`, optionally with the immutable `--code-ref`, to verify all
inputs, code and worker outputs without repeating measurements. Plan entries
with `root: research` resolve against this worktree; other inputs resolve
against `--data-root`.

The timing measures sequence inference after model loading and warmup. Peak
memory is the absolute Linux `VmHWM` in each fresh worker, including the shared
import baseline, input bundle, model and output. It is not a difference of
high-water marks or a parent's cumulative child peak. Candidate output must
reproduce the archived strand propensities, and each method must reproduce
its numerical output across repetitions.

These are component diagnostics. They exclude model fitting, BAM extraction,
FASTA I/O, DWM local regression, expected-count normalization, residual and
footprint scoring, and track output. They cannot certify the full workflow
runtime or memory promotion gates. The DWM pickles are checksum-pinned local
production artifacts; this is not an entry point for untrusted model files.

### Covariate-matched direct-bias diagnostic

`evaluate_chrombpnet_covariate_diagnostic.py` tests whether existing K562/HepG2
comparisons retain their discrimination and protection after tighter covariate
matching. It never refits the candidate or changes its recorded dispersion.
All three policies are declared in the plan: motif only; motif plus local
accessibility (primary); and those two plus observed cut total.

```bash
python benchmarks/scripts/evaluate_chrombpnet_covariate_diagnostic.py \
  --plan benchmarks/manifests/chrombpnet_covariate_diagnostic_v1.json \
  --data-root /path/to/data-checkout \
  --outdir benchmarks/results/chrombpnet_covariate_diagnostic
```

Covariate standard deviations come only from chr1–15. Accessibility is the
201-position raw-track sum used by the archived matching tables; the pinned
10M uncorrected tracks supply the missing training column and must reproduce
the archived test values. Observed cut total comes from the direct BAM-derived
profiles, preserving the distinction between that count and normalized track
coverage. Accessibility and cut total use `log1p` for matching.

Pairs share a chromosome and lie within 0.2 training standard deviations on
every matching variable. The assignment first maximizes the number of valid
1:1 pairs, then minimizes squared distance, with deterministic coordinate order.
All selected cohorts and balance diagnostics are frozen before model score
columns are read. No policy is chosen after seeing its AUROC or AUPRC.

The diagnostic reuses recorded scores, reconstructs profiles for the new
subsets and checks every reconstructed score against its original value.
It preserves checkpoints for matched cohorts, profile reconstruction and each
cell/policy evaluation. Completed checksum-valid stages are skipped on resume.
Outputs include all policy metrics, chromosome-block intervals for powered
tasks, signed shape summaries, balance and power failures, and primary-policy
aggregate PDFs. The `test` split field names the historical chromosome fold;
all these results are post-hoc development sensitivity, not untouched evidence.
Conditioning on accessibility can also remove biological occupancy signal;
this is not a causal adjustment or a new promotion result.

For an existing run after a verification-only code update, use `--verify-only`
with the immutable calculation-code commit. This checks historical inputs,
every recorded code blob, nested cohort files, profile caches and evaluation
outputs without recreating a configuration or rerunning calculations. Mutable
branch names are rejected for `--code-ref`.

```bash
python benchmarks/scripts/evaluate_chrombpnet_covariate_diagnostic.py \
  --plan benchmarks/manifests/chrombpnet_covariate_diagnostic_v1.json \
  --data-root /path/to/data-checkout \
  --outdir benchmarks/results/chrombpnet_covariate_diagnostic_20260908_v2 \
  --verify-only --code-ref c39dcd9
```

### Per-factor ledger for the completed 96-run experiment

`build_direct_bias_gate_ledger.py` summarizes the unchanged development result.
It preserves all 21 cell/TF tasks and all 48 runs per task, including unscored
and underpowered tasks. It checks candidate/comparator class counts, recomputes
paired metric gains, retains shape failures and checks the prespecified
principal interval matrix. It does not calculate a new locked mean-gain gate,
retune a model, read holdout labels or certify safety or promotion.

```bash
python benchmarks/scripts/build_direct_bias_gate_ledger.py \
  --manifest benchmarks/results/direct_bias_stability_96_20260908_v1/stability_manifest.json \
  --manifest-sha256 52fc25ca38d8a83b301ace2980a4947583a5605743cf729a66c6d9c8051da798 \
  --outdir benchmarks/results/direct_bias_gate_ledger_20260922_v1
```

The output `gate_ledger.json` records the source identities, script checksum,
per-task ranges, measured intervals, current-count imbalance and individual
shape exceptions. The reader verifies all four source summary-table hashes.
Exact repeated output is accepted; a different existing ledger is never
overwritten. Full nested-input verification remains the separate
`extend_direct_bias_stability.py --verify-only --code-ref c9ecd82` operation.
The ledger's missing safety/resource/uncertainty statuses cannot be interpreted
as passes. Consult the exact-operator safety report separately.

### Complete-workflow resource measurements

`../manifests/end_to_end_resource_protocol_v1.json` locks the premeasurement
design. Execution must additionally freeze the exact commands, inputs, code,
runtime, expected case matrix and validated measurement implementation. No
complete-workflow measurement exists yet. Candidate and production DWM must
write equivalent products; the same-operator DWM residual comparison is a
separate diagnostic, not a substitute for the production baseline.

The promotion evaluator's resource TSV now requires `cpu_seconds` alongside
`runtime_seconds`, `peak_memory_mb` and `model_size_mb`. It also requires
`benchmark_id` and `workload_id` (SHA-256 identities), `repetition` (0–5),
integer `cores`, and these explicit scope fields:

| Field | Required value |
| --- | --- |
| `resource_scope` | `complete_pipeline` |
| `pipeline_mode` | `fit_and_apply` |
| `cache_policy` | `fresh_process_warm_os_cache` |
| `cpu_scope` | `process_tree` |
| `memory_scope` | `absolute_process_tree_peak` |

Every candidate run must pair with a baseline run having the same workload,
repetition, cores, cache policy and pipeline mode. All six repetitions are
required. CPU ratios use within-workload means; memory ratios use absolute
within-workload maxima. Every workload must pass. The producer must verify
that the entire frozen case matrix is present: declarations in a TSV cannot
prove that its inputs, commands or measurement instrumentation were valid.
Older sequence-component tables remain diagnostic evidence and fail this gate.

On a Linux host with a systemd user manager and delegated CPU/memory cgroup
controllers, qualify process-tree accounting before scientific measurement:

```bash
python benchmarks/scripts/probe_process_tree_resources.py \
  --outdir benchmarks/results/process_tree_probe_local
```

This command runs synthetic sequential, concurrent, failed and orphaned-child
cases only. It requires a new output directory and retains logs and receipts.
Each workload receives a fresh kernel cgroup, outside the supervisor's memory
and CPU accounting. It reads the absolute simultaneous `memory.peak` and the
entire subtree's `cpu.stat`, including children that have already exited.
Memory is cgroup-accounted bytes (reported as decimal MB), including charged
file cache and kernel memory; it is not RSS. The frozen benchmark must use the
same accounting and cache policy for both methods. A failed qualification is
not permission to substitute sampled peaks or parent-only memory.

The probe stops and clears only its uniquely named temporary user services.
It neither reads scientific data nor supplies a computational promotion pass.

### Research extraction correctness

Read `../STRAND_EXTRACTION_AUDIT_20260922.md` before reusing historical strand
profiles. The corrected reader visits each requested chromosome once to avoid
double-counting long alignments across disjoint fetch windows and to retain
soft-clipped cuts outside aligned query spans. Old files remain preserved.

`audit_strand_extraction.py` freezes code, data and label-free training-chr1
coordinates before comparing the old and corrected reader against independent
sorted-cut counting. It never fits or scores a candidate. From the research
root, use a new ignored output directory:

```bash
python benchmarks/scripts/audit_strand_extraction.py \
  --data-root /path/to/original-data-checkout \
  --outdir benchmarks/results/strand_extraction_audit_local --freeze-only
python benchmarks/scripts/audit_strand_extraction.py \
  --data-root /path/to/original-data-checkout \
  --outdir benchmarks/results/strand_extraction_audit_local
```

The source registry is checksum-pinned to the completed 96-run configuration.
The script rejects changed inputs/code/runtime and an existing completed audit.
It is a correctness diagnostic, not a new evaluation or promotion result.

The fixed-candidate reassessment of all 96 archived development runs uses
`reassess_corrected_strand_counts.py`. Use a new research output directory;
never point it at the historical experiment. Freeze before running:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python benchmarks/scripts/reassess_corrected_strand_counts.py \
  --data-root /path/to/original-data-checkout \
  --outdir benchmarks/results/corrected_count_reassessment_local --freeze-only
python benchmarks/scripts/reassess_corrected_strand_counts.py \
  --data-root /path/to/original-data-checkout \
  --outdir benchmarks/results/corrected_count_reassessment_local
```

Append `--verify-only` to verify all corrected inputs and completed outputs
without extraction or scoring. The original candidate settings and cohorts
cannot change. Observed counts are freshly extracted; existing DWM expectations
and label-free sequence propensities retain their checked identities. Completed
stages are verified and reused; incomplete files are preserved and rejected.
The new five summary tables record metrics, shape, coverage, profiles and
bootstrap intervals. These are post-hoc implementation-reassessment results,
not untouched validation or authority to tune on examined chromosomes.
