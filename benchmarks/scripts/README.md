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
- `evaluate_shifted_atac_null_policy.py`: score frozen detectors after deliberately misaligning motif-centered cellular ATAC profiles, providing a label-free accessibility/nucleosome null that complements naked DNA.
- `evaluate_naked_dna_functional_policy.py`: score the frozen functional policy on naked-DNA profiles; use `--candidate-only` for an independent control replicate when no same-replicate DWM artifact exists.
- `calibrate_dual_null_posteriors.py`: require posterior calls to exceed conservative per-TF thresholds from both naked DNA and motif-misaligned cellular ATAC.
- `apply_frozen_dual_null_thresholds.py`: apply already frozen per-TF dual-null thresholds to an independent score table without refitting or borrowing a comparator replicate.
- `ensemble_parametric_bias_models.py`: combine compatible seed fits into one checksummed geometric coefficient ensemble for prespecified crossed controls.
- `render_functional_aggregate_comparison.py`: render blinded DWM-versus-strand per-TF aggregate curves with bootstrap confidence bands and an explicit blinding key.
- `evaluate_strand_label_free_models.py`: fit count spline/GP, FDA, and hybrid mixtures on separate label-free strand artifacts and score only development validation labels.
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
- `evaluate_tf_geometry_naked_dna.py`: freeze per-method score thresholds on development ChIP-negative sites and test an exact TF-specific geometry plus the conventional footprint score on common naked-DNA motif sites.
- `evaluate_tf_geometry_external_transfer.py`: apply one frozen development-cell TF geometry without retuning to external labeled cells, compare it with the conventional footprint score on identical sites, and report chromosome-block, covariate-residual, and replicate evidence.
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
failed gate. A development-only audit can explain what remains, but cannot
return a promotion pass because locked GM12878/IMR-90 validation is mandatory.

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
