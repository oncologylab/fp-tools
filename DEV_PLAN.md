# fp-tools Development Plan

Last updated: 2026-09-08

## Current Baseline

fp-tools is a command-first Python 3.11–3.13 package for bulk and pseudobulk
ATAC-seq footprinting, motif analysis, differential reports, and browser/YAML
wrappers. Scientific workflow logic belongs in
`src/fp_tools/tools/` and shared helpers in `src/fp_tools/utils/`; CLI and GUI
layers remain thin.

The supported console commands are declared in `pyproject.toml`. Public docs,
examples, package metadata, and smoke tests must use that declaration as the
source of truth. Removed TOBIAS-style console aliases must not be restored.

The current project layout, YAML compatibility, output filenames, and command
names are stable public contracts. Internal cleanup should preserve them unless
a separately planned breaking release explicitly changes the contract.

## Completed Capabilities

- Linux CLI/container raw-read ATAC preparation with modern and legacy
  processing profiles.
- A managed, versioned external-tool runtime downloaded on first use, with
  Linux raw-read components, cross-platform optional MEME components, and the
  complete Linux container retained as an explicit backend.
- Cross-platform one-command bulk analysis from coordinate-sorted BAM/BAI and
  matching peak BED files through portable multi-comparison HTML reports.
- Separate Linux CLI/container FASTQ preparation through `prepare-atac`, which
  writes the BAM/BAI and peak BED sample table accepted by the downstream bulk
  workflow.
- Tn5 bias correction, background scaling, footprint scoring, candidate calls,
  known-motif matching, de novo motif preparation, and aggregate plotting.
- Replicate-aware differential footprint reports, including per-sample motif
  matrices and empirical-Bayes residual-variance moderation with biological
  samples as the inferential units. Bundle and self-contained multi-comparison
  reviews share one browser implementation; standalone reports preserve input
  order with an exact-record selector and remove aggregate-only controls when
  profiles are unavailable. The shared browser also supports aggregate-neutral
  volcano highlighting, lightweight user-selected TF labels, and a two-mode
  ranked-motif waterfall with reciprocal score/significance labels, directional
  blue/red color strength, a stable square volcano viewport, and comparison
  labels embedded in exported waterfall, volcano, and combined-panel SVGs.
- Region-set differential footprinting for one sample or paired biological
  replicates, with equal region weighting, optional matching strata,
  resampling-based uncertainty, paired empirical-Bayes inference, and the same
  portable interactive report format.
- A compact HepG2 HNF4A/FOXA2 region-set example with four mutually exclusive,
  baseline-matched region classes, six pairwise comparisons, three biological
  replicates, complete motif statistics, and curated all-site aggregate views.
- Pseudobulk fragment/BAM generation and single-cell signature reporting.
- Streamlit GUI whose saved YAML runs through `run-yaml-workflow`.
- Local tests, package builds, GitHub CI, GitHub Pages deployment, and manual
  PyPI publication.
- A complete amd64/arm64 container with the external genomics toolchain and a
  browser GUI as its default entry point.
- BAM-first self-contained GUI executables for Windows x64 and Apple Silicon,
  with branded operating-system icons, a bundled native application window,
  and frozen-safe command and workflow dispatch. Every GUI starts from BAM/BAI
  plus peak BED files. Linux and Intel macOS use the browser-based Python
  package; the complete container remains available on amd64 and arm64.
- The Apple Silicon download remains an explicitly unsigned preview. Its DMG
  includes the same one-line quarantine-removal instruction shown on the
  installation page, and CI verifies the ad-hoc bundle signature, ARM64
  executable, DMG contents, quarantine removal, and native-window launch.
- Branded GitHub and MkDocs presentation with responsive embedded report and
  GUI demos plus automated desktop/mobile browser audits.
- Responsive GUI validation output with readable wrapped paths and locally
  scrolling YAML previews, verified in the native Windows and macOS audits.
  Validation uses concise field names, preserves whole words where possible,
  and failed desktop audits retain screenshots as CI artifacts for diagnosis.
- Responsive GUI navigation and forms across every destination: the mobile
  sidebar has named open/close controls, the persistent sidebar header stays
  compact, workspace details remain available on demand, loaded YAML refreshes
  all form controls, and main-panel labels use explicit high-contrast colors.
  Frozen Windows and Apple Silicon audits exercise every page at desktop and
  narrow widths and retain screenshots for visual review.
- The documentation, browser GUI, static demos, and native Windows/macOS
  windows use one explicit light theme regardless of the operating-system
  preference. Browser audits emulate a dark system preference to guard this
  presentation contract.
- A task-oriented GUI home page, neutral first-open forms, centralized field
  labels, compact advanced options, clearly disabled launch controls, and
  responsive form/run panels. Every route now reuses the Home page's light
  hero-card, typography, borders, controls, and action colors; form controls
  remain visibly bounded on white cards, and the form/run layout is compact at
  desktop widths before stacking cleanly at narrower breakpoints. Fresh desktop
  installations no longer open with repository-relative example paths
  presented as user inputs.
- Browser-verified plot controls and SVG exports for aggregate and
  aggregate-free reports, including explicit subplot-bound checks.
- Command-aware logging across shared analysis engines, so `match-motifs` and
  `diff-footprints` retain their public names in direct and wrapper-run logs.
- A manifest-pinned, storage-conscious ENCODE workflow for 17 biological
  replicates from seven cancer cell lines. Its resumable runner uses a
  pair-specific union of released IDR-thresholded peaks, peak-q95 scaling, and
  the same differential workflow as the preserved three-replicate
  K562-HepG2 report. The dependency-free static browser has two directional
  selectors and reads compact canonical report payloads.
- Benchmark-only footprint detectability auditing that ranks motifs within
  source analyses, collapses repeated pairwise/control rows to independent
  biological contexts, joins exact-name RNA evidence, records hashed inputs,
  and emits both machine-readable tables and a searchable HTML report. A
  separate matched-label classifier applies prespecified correction, scoring,
  coverage, ambiguity, and information-limit diagnoses without treating low
  aggregate scores as evidence of TF absence.
- A locked footprint-improvement study with 35 ENCODE cell/TF tasks,
  development/holdout and chromosome partitions, summit-supported motif-site
  labels, matched negative controls, deterministic nested depth subsets,
  raw/PWM/DWM/reused-bias correction arms, genomic-block uncertainty,
  depth/correction summaries, one-time holdout promotion gates, and tiered
  nutrient-stress replication. The external nutrient resources are pinned to
  GSE144833 and GSE137034 (including GSE137031/GSE137032).
- A label-free percentile-rank evidence-fusion primitive and locked evaluator.
  Its first footprint/PWM soft-OR candidate passed the available K562/HepG2
  development slice but failed the one-time cell-line holdout because it did
  not meet mean-gain or strong-positive non-regression gates. It remains
  experimental; the production footprint score is unchanged.
- An opt-in aggregate-shape detectability mode with per-site outer-flank RMS
  normalization, site-level 95% confidence bands, quantitative depletion
  diagnostics, explicit underpowered/not-detected states, and functional
  per-panel or grouped y-axis scaling. This mode improves visualization and
  failure identification without changing the production footprint score or
  presenting aggregate shape as proof of occupancy.
- An opt-in dual-geometry `call-footprints --score hybrid` arm that adds a
  low-weight, locally standardized wide symmetric depletion channel. It
  reproducibly improved CTCF and REST on locked K562/HepG2 chromosomes but
  regressed some JUND/MAX tasks, so it remains experimental and cannot replace
  the production `footprint` default.
- A staged per-TF footprint research engine that searches correction arm,
  center/shoulder geometry, local normalization, center statistic, and
  asymmetry penalty on train chromosomes, freezes candidates on validation
  chromosomes, and evaluates identical finite sites only after test unlock.
  Reproducible ENCODE peak discovery, cell-specific motif scans, summit-based
  labels, optimal motif/accessibility matching, matched legacy comparisons,
  and bound/unbound aggregate figures support the experiment. The first
  10-million-fragment test produced large CTCF gains and selective MEF2A,
  MEF2D, ARID3A, and MYC gains, but also clear TF/context regressions and
  underpowered tasks. It remains research-only and does not change the default
  scorer.
- A second locked research specification for a fast ChromBPNet-inspired
  alternative. The implementation now includes conditional parametric
  sequence-bias models, pooled sample adaptation, calibrated Poisson/NB
  residuals, weighted functional PCA, spline and 25-knot Matérn
  Gaussian-process-equivalent footprint mixtures, an exact-GP reference,
  hybrid FDA-GP scoring, and replicate-level functional differential tests.
  Model files use checksummed NPZ plus JSON rather than executable pickle.
  These APIs remain research-only until the GM12878/IMR-90, naked-DNA,
  uncertainty, runtime, and non-regression gates pass.
- A label-free parametric-bias benchmark now extracts GC-matched low-signal
  nonpeak or mitochondrial control windows, tests +4/-5 and +4/-4, fits
  SELMA-equivalent adjacent-simplex and 81-bp interaction models, and compares
  sample, pooled, and pooled-plus-adapted estimates. Validation likelihood must
  improve over the uniform conditional control before any model advances. The
  functional benchmark can consume a frozen NPZ model directly and derive
  total-preserving expected motif-site profiles without producing genome-wide
  intermediate tracks.
- The functional research path now has a strand-aware profile artifact. It
  preserves plus/minus observed and expected cuts, swaps strands when reverse
  motifs are orientation-aligned, and exposes combined, shared-strand, and
  antisymmetric signed-deviance channels. This enables genuine TF-specific
  asymmetry experiments rather than inferring asymmetry from a combined track.
- Functional evaluation now produces bootstrap-banded aggregate profiles,
  explicit depth/width/shoulder/asymmetry/periodicity descriptors, blinded
  multi-correction PDF panels, and unsupervised phenotype clusters learned
  from normalized residual curves rather than motif-family names. Numeric
  discrimination is therefore reviewed alongside visible footprint shape.
- Replicate-aware differential functional orchestration now consumes frozen
  strand-profile artifacts and reports the complete condition-difference
  curve, pointwise and simultaneous uncertainty, global permutation testing,
  BH-adjusted significance, replicate consistency, and changes in footprint
  depth, width, shoulders, asymmetry, and periodicity. Nutrient-stress data
  remain locked until model selection is complete.
- Bias-model leakage diagnostics now score every eligible development JASPAR
  motif against nearby local-sequence controls without reading ChIP labels.
  They report bootstrap-banded response functions, broad center/flank effects,
  and cross-cell/model concordance. A sequence response is a review flag—not
  automatic evidence of occupancy leakage—until naked-DNA and held-out TF
  behavior distinguish intrinsic Tn5 preference from overlearned chromatin.
- A fail-closed functional promotion auditor now implements the complete
  preregistered decision: locked external metrics, functional separation,
  positive-control non-regression, clustered bootstrap evidence, naked-DNA
  false positives, motif-response review, stability, CPU/memory/model size,
  uncertainty coverage, and GP-versus-spline justification. Development-only
  evidence or any missing table cannot produce a promotion pass.
- The frozen parametric factorization experiment now has checksum-locked
  SK-N-SH and GM23338 holdout manifests, LOG21/LOG41/LOG81 and SELMA10
  coefficient grids, bounded sample bias-strength calibration, flank-only
  accessibility splines, hierarchical latent footprint profiles, five
  residual ablations, and immutable resumable stage manifests. Independent
  enzyme controls are compared with paired block bootstraps against a safe
  TOBIAS-style all-pairs DWM reference. Candidate-screening and final
  naked-DNA libraries are extracted without fitting, and motif-centered
  observed-minus-predicted residuals have an explicit leakage gate. All new
  artifacts use checksummed NPZ plus JSON and remain research-only. Continuous
  residual safety is evaluated with cutoffs frozen from chr16--18 matched
  ChIP-negative sites before they are applied to label-free naked-DNA
  replicate 2; both finite-site and nonzero-cut false-positive rates must pass
  before a residual can enter the immutable configuration freeze.
- The first frozen conditional-multinomial checkpoint conditions explicitly on
  each motif site's cut total and learns only profile shape. Across the nine
  internally eligible K562/HepG2 tasks it improved mean AUROC by 0.0139,
  relative AUPRC by 4.38%, and functional separation by 28.7% over DWM, while
  all 14 evaluated task thresholds passed naked-DNA replicate-2 safety. The
  gain is not general: the same candidate was 0.0184 AUROC and 0.73% relative
  AUPRC below the exact raw-signal guardrail on average. CTCF is the clear
  internal exception, improving AUROC/AUPRC over raw by 0.0647/15.95% in HepG2
  and 0.0421/9.43% in K562 with visibly stronger canonical aggregate
  protection. FOXA1's smaller numerical gain does not yet have a convincing
  canonical aggregate shape. A count/FDA rank-consensus screen also failed the
  10% AUPRC gate, and FDA-heavy mixtures produced naked-DNA false positives.
  These results remain internal research evidence: full-depth confirmation,
  control-source qualification, the official ChromBPNet comparison, and
  locked external validation are still required before any promotion.
- Independent control-source qualification now favors mitochondrial-only
  training over combined mitochondrial plus naked-DNA training for both cut
  conventions. The retained internal +4/-5 reference remains the five-seed
  mitochondrial LOG21 ensemble; adding the naked-DNA library did not improve
  its paired held-out control likelihood. A validation-frozen partial-bias
  experiment identifies one reproducible TF-specific exception to the broader
  failure: CTCF shrinkage improves both AUROC and AUPRC over raw signal in
  HepG2 and K562 at 10M, 25M, and 50M fragments for every seed, and across all
  three biological replicates in each cell. Mean full-replicate AUROC gains
  over raw are 0.0410 in HepG2 and 0.0378 in K562. The same policy makes zero
  calls in the independent naked-DNA replicate-2 CTCF panels, with a 1.88%
  Wilson upper bound. For other eligible TFs, large gains over DWM generally
  return to raw-signal performance rather than add information; depth reports
  therefore distinguish genuine gains over raw from recovery after DWM
  overcorrection. This is a CTCF-specific internal research route, not evidence
  for a general correction replacement, and the package default is unchanged.
- A separate integrity-checked detector ledger now prevents TF-specific shape
  gains from being conflated with bias-correction gains. The globally frozen
  `count_only` operator passes independent naked-DNA replicate-2 safety for all
  14 evaluable task panels. FOXA1 has a reproducible occupancy-classification
  gain: on the exact common chromosome-test support it improves AUROC by
  0.0624 and relative AUPRC by 5.92% over raw signal, both metrics improve in
  all three HepG2 biological replicates, and 25M/50M seeds are stable. Visual
  and numeric shape auditing rejects this as a footprint improvement, however:
  the candidate has central enrichment rather than protection (depletion
  -0.0423). It is therefore classified as an occupancy-signal diagnostic, not
  a FOX footprint model. The detector also improves CTCF over DWM, but its
  AUROC comparison with raw signal changes with support/depth and therefore
  does not supersede the simpler CTCF shrinkage result. The other seven
  internally eligible tasks fail the raw guardrail, and 12 tasks remain
  underpowered. No general detector or correction promotion follows.
- A five-family secondary screen now compares validation-frozen count, FDA,
  anchored FDA, residualized FDA, and hybrid FDA-GP policies against raw and
  DWM scores on exact common test support. It joins independent naked-DNA
  safety and requires central protection, a negative ChIP-positive-versus-
  matched-negative center contrast, and pointwise aggregate evidence before a
  discrimination gain is called a footprint gain. Only HepG2 CTCF passes the
  strong shape guard: count improves AUROC by 0.0725 and relative AUPRC by
  15.90% over raw (bootstrap lower bounds 0.0406 and 8.96%), while anchored
  FDA improves AUROC by 0.0545 and relative AUPRC by 3.76% (lower bounds
  0.0408 and 0.80%). K562 CTCF remains support-sensitive, and no difficult TF
  family produces a safety-passing significant shape-guarded gain over raw.
  This confirms that the current FDA/GP variants do not yet solve hard-TF
  footprint detection and prevents accessibility-driven occupancy gains from
  being reported as footprints.
- A canonical-protection conditional model now constrains the learned center
  to depletion and the shoulders to enrichment. Validation-frozen models were
  evaluated for all 14 K562/HepG2 task contexts and passed independent
  naked-DNA replicate-2 safety. On exact chromosome-test support, CTCF gains
  over raw were 0.0589 AUROC/14.01% relative AUPRC in HepG2 and 0.0508/8.94%
  in K562 with visibly canonical protection. The model beats DWM on both
  metrics in all six biological replicates, but its full-replicate AUROC is
  below raw in all six; no other TF has a defensible shape-guarded gain. It is
  therefore a useful detector ablation, not a promotion candidate.
- Both complete 120-fit mitochondrial and naked-DNA sequence-bias grids are
  now collated into five-seed geometric ensembles. The naked-DNA grid ranks
  LOG21, +4/-5, L2 0.01 first on held-out conditional likelihood; cross-source
  qualification still retains the mitochondrial LOG21 ensemble for the
  internal footprint experiment. Resumption is fail-closed: NPZ/JSON
  checksums and the full source, shift, feature, penalty, depth, sample, seed,
  window, and cut identity must match before an existing fit is rescored.
- The official ChromBPNet reference is pinned to commit
  `09938fdb4397ec0006510e5251e48920a505d4de` and a digest-locked TensorFlow
  2.8 GPU container. On identical K562 test sites and conventional geometry,
  frozen CPU parametric bias scores CTCF at AUROC 0.7421/AUPRC 0.7180, close
  to ChromBPNet bias at 0.7509/0.7229 and above DWM at 0.6810/0.6455. Its
  paired gains over DWM have positive chromosome-block intervals; its AUPRC
  difference from ChromBPNet includes zero. Neither bias method is general:
  both regress eligible MYC, while the parametric model gives only a small
  ZNF384 gain over DWM and remains below raw. The checksum-locked HepG2 bias
  model, regulatory training, predictions, and same-site comparison are now
  checksum-verified complete. The HepG2 comparison used dispersion zero,
  whereas the recorded K562 comparison used 1.6144318781609228; their results
  must remain separate experiments rather than a pooled common-operator claim.
- A completed-reference audit now verifies comparison inputs, outputs and
  canonical identities and reports all 21 development TF tasks with signed
  protection, shoulders, asymmetry, label-dependent site retention, motif-score
  balance and accessibility balance. Only CTCF has significant parametric gains
  over both DWM and raw with canonical protection in both recorded comparisons.
  HepG2 CTCF scores 0.7149 AUROC/0.7388 AUPRC versus DWM 0.6477/0.6250 and
  ChromBPNet bias 0.7253/0.7570. The parametric mean gains across eligible tasks
  are only 0.0112 AUROC/1.99% relative AUPRC in HepG2 and 0.0053/1.25% in
  K562, below the general promotion gates. CTCF motif scores remain imbalanced
  (standardized differences 0.482 HepG2 and 0.463 K562). The older deep-reference
  report cited factorization-residual naked-DNA safety, which does not certify
  the direct parametric scoring operator. A new direct-operator audit freezes
  code, inputs and validation-only thresholds before scoring existing replicate-2
  profiles. Its coverage qualification reuses those scores and thresholds:
  only CTCF has sufficient archived DWM control coverage. On DWM-comparable
  informative sites, HepG2 has 0/121 false calls (Wilson upper 3.08%); K562
  has 1/114 (upper 4.80%, increase over DWM 0.88 percentage points). Each
  excludes one informative site lacking a DWM expectation; the absolute
  candidate FPR also passes with those sites retained. This is retrospective
  development candidate screening, not final untouched safety. Sparse counts
  and the absence of a large residual flag do not establish no occupancy
  leakage. Other TFs cannot pass relative safety from this CTCF-region DWM
  archive. See
  `benchmarks/CHROMBPNET_AUDIT_20260908.md` for the evidence and remaining work.
- The one-page direct-bias renderer now rejects legacy factorization safety,
  requires exact comparison identity plus coverage-qualified direct safety,
  and gates on significant paired gains over DWM/raw and canonical protection.
  It attaches safety only to the candidate, reports the actual cell/dispersion
  and covariate imbalance, and preserves historical reports. The new HepG2
  CTCF report is development evidence with candidate-screen safety; no general
  correction, final untouched safety or package promotion is claimed. Its six
  files were uploaded to the new Box folder
  `Yaoxiang/fp-tools/2026-09-08_HepG2_CTCF_direct_bias_candidate_screen`
  on 2026-09-08; `rclone check --download --one-way` reported six matches and
  zero differences. Historical report folders were preserved.
- A frozen post-hoc covariate sensitivity now covers all 21 development tasks
  with three declared matching policies and unchanged scores/dispersions. The
  primary motif/accessibility match retains 634 HepG2 and 609 K562 CTCF pairs
  with covariate SMDs below 0.01. Parametric gains over DWM remain significant
  in both cells (AUROC +0.0477/+0.0513, relative AUPRC +12.70%/+10.73%) with
  canonical protection, but K562's intervals against raw include zero. A third
  fixed policy also balances direct observed cut total and keeps that pattern.
  Fourteen tasks become underpowered; balanced FOXA1/MYC remain near chance
  and powered ZNF384 tasks remain weak without canonical protection. General
  mean-gain and difficult-family gates still fail. This is sensitivity on
  examined chromosomes, not new holdout evidence or a new candidate. See
  `benchmarks/CHROMBPNET_COVARIATE_DIAGNOSTIC_20260908.md` for complete results.
  Historical verification now checks nested cohort artifacts and can validate
  recorded calculation-code blobs against their immutable Git commit without
  rerunning completed experiments.
- Direct-bias stability now has a fixed research plan for the 36 available
  development BAMs: five seeds at 10M/25M/50M for replicate 1 plus all three
  full-depth biological replicates per cell. Existing partial-bias matrices
  use tuning chromosomes and cannot supply the missing test-site evidence.
  Preflight also found zero DWM expectations at 65 HepG2 and 66 K562 CTCF
  sites in the historical finite-support comparison. The new diagnostic
  freezes positive expected coverage across all runs, preserves complete
  existing matched pairs, and keeps zero observed-count sites. It separates
  the original full-depth DWM expectation from sample/depth-matched DWM.
  All 36 runs are complete and verified from calculation commit `c453a37`.
  CTCF retains canonical protection and positive point gains over raw and both
  DWM comparators throughout. All six full biological replicates have positive
  AUROC/AUPRC bootstrap intervals against raw; matched AUROC gains are
  0.0592–0.0627 in HepG2 and 0.0575–0.0666 in K562. No other powered TF passes
  the canonical-shape screen. Current cut-total SMD drifts above 0.1 in 11/36
  matched CTCF runs (maximum 0.134), so fixed matching does not eliminate
  accessibility confounding. General relative-AUPRC and difficult-family gates
  still fail; full runtime/memory, uncertainty and untouched safety remain
  incomplete. See `benchmarks/DIRECT_BIAS_STABILITY_20260908.md` for the full
  evidence, expected-coverage exclusions and limitations.
- A concise CTCF stability report now requires all declared run directions,
  full-replicate paired intervals, canonical protection and matching direct
  safety identities before rendering. Its 36-row metrics TSV retains balance
  drift and clearly distinguishes candidate-screen safety from untouched final
  safety. The six report files were placed in the new Box folder
  `Yaoxiang/fp-tools/2026-09-08_CTCF_direct_bias_stability_coverage`; download
  verification reported six matches and zero differences on 2026-09-08.
- A separate frozen sequence-inference benchmark compares the archived
  production Cython DWM scorer with the unchanged parametric model on 512
  coordinate-selected loci per development cell. Six paired repetitions use
  alternating method order and fresh single-CPU processes. It validates the
  candidate against archived strand propensities and measures absolute Linux
  worker peak memory. Its component ratios cannot certify the full workflow
  runtime/memory gates; fitting, BAM/FASTA I/O, local correction, detector
  scoring and track output remain outside this measurement.
  All 24 workers completed: candidate/DWM mean CPU ratios are 0.47385 in
  HepG2 and 0.47654 in K562, with maximum worker-memory ratios 0.99247 and
  0.99202. Candidate propensities reproduce the archived values; every method
  is numerically stable across timing repetitions. See
  `benchmarks/DIRECT_BIAS_INFERENCE_20260908.md` for the frozen protocol,
  checksums, validation and explicit disposition of the remaining gates.
- Nested depth-matrix resume now verifies source identity, sampling metadata,
  BAM and index checksums. Previously it silently skipped any existing BAM,
  including interrupted outputs. Regression fixtures reproduced six failures
  before the fix; the corrected implementation also rejects a source that
  differs from an explicitly frozen checksum. Historical data are preserved;
  this stricter path prepares missing inputs without rewriting old manifests.
- The missing depth/seed inputs for HepG2/K562 replicates 2 and 3 now have a
  separate frozen preparation plan. It declares all 60 missing combinations,
  the original pair-preserving hash selection and the same refitted DWM
  command used for replicate 1. A native fragment counter is verified against
  the original predicate for all 4,096 SAM flag combinations, including MAPQ
  zero and unpaired reads. Preparation checks source and nested output hashes,
  pins the production source tree/runtime, limits concurrency to two jobs and
  uses a kernel lock. This stage prepares data only; subsequent evaluation
  must freeze support across all 96 runs and preserve complete matched pairs.
  The four fragment-count stages are complete: HepG2 replicates 2/3 contain
  96,480,203/105,920,129 usable fragments, and K562 replicates 2/3 contain
  76,393,585/114,471,886. All support the requested 50M depth. Preparation
  is in progress; see `benchmarks/DIRECT_BIAS_REPLICATE_DEPTHS_20260908.md`
  for the frozen identities, tested counting rule and remaining requirements.
- The significant, safety-passing internal CTCF partial-bias result has a
  concise one-page report, metrics, checksums, and limitations in Box under
  `Yaoxiang/fp-tools/2026-09-02_CTCF_frozen_parametric_bias_shrinkage`.
  This is explicitly a research artifact; it does not change the package
  default and is not evidence for a general TF correction.
- The K562 deep-reference comparison has a second checksum-verified Box report
  under
  `Yaoxiang/fp-tools/2026-09-02_K562_CTCF_parametric_vs_ChromBPNet`.
  It contains a one-page same-site aggregate and metric comparison, its source
  metrics, and explicit limitations; the report does not claim general
  promotion from the CTCF result.

## Current ENCODE Resource

- All 21 pairwise comparisons are available in the static browser, and every
  aggregate profile uses the complete motif-site set rather than bound-only
  sites.
- The preserved K562-HepG2 comparison retains its reference scientific
  payload, and the browser loads compact report data, profile shards, and motif
  logos on demand.
- The public report page embeds this browser with K562 versus HepG2 as the
  default comparison.

## Near-Term Priorities

1. Keep the seven-line ENCODE cancer resource reproducible and
   storage-conscious. Preserve all 1,019 motifs, 17 biological replicates, and
   21 prespecified contrasts; resource membership remains independent of
   observed differential results, and the K562-HepG2 comparison must retain its
   exact preserved scientific payload.
2. Keep README, MkDocs, examples, CLI help, package metadata, and the static
   cancer-cell-line browser synchronized.
3. Add focused regression tests for every user-visible bug or command-contract
   change.
4. Reduce avoidable warnings and multiprocessing fragility without changing
   scientific output.
5. Improve GUI cancellation and output previews only through existing command
   and YAML interfaces.
6. Keep benchmark claims proportional to validated public data and recorded
   metrics.
7. Keep public pages free of broken assets, browser errors, responsive
   overflow, inaccessible navigation, and broken live embeds.
8. Keep region-set examples outcome-independent: define groups from external
   annotations, match baseline signal before testing, report all motifs, and
   use explicit display motifs only to configure the initial browser view.
9. Maintain BAM/peak scientific I/O parity on Windows, macOS, and Linux, plus
   desktop-app parity on Windows x64 and Apple Silicon. Test Linux raw-read and
   cross-platform MEME runtime artifacts separately, and keep the complete
   multi-architecture Linux container as a reproducible alternative backend.
10. Continue the locked footprint detectability matrix beyond the completed
    K562/HepG2 10-million-fragment, seed-2026 raw/PWM/DWM correction slice.
    Add the remaining depths, seeds, replicates, and naked-DNA controls before
    changing the default correction. The first label-free
    footprint/PWM fusion candidate has been rejected after its one-time
    MCF-7/A549/HCT116/Panc1 evaluation and must not be retuned on that holdout.
    Keep
    nutrient data as a prospective application rather than model-training data,
    and do not promote a new scorer until held-out performance, calibration,
    naked-DNA false-discovery control, and strong-positive non-regression gates
    pass.
11. Run the `footprint_functional_v1` stages in order: establish the cut-site
    convention and control-trained parametric bias models; compare residual
    formulations; freeze spline/FDA/GP/hybrid models on K562/HepG2; use the
    already examined GM12878/IMR-90 only as transfer diagnostics. Keep the
    current DWM correction and footprint score as defaults unless every
    scientific and CPU performance gate passes.
12. Complete `frozen_parametric_factorization_v1` without opening its new
    holdout labels out of order: finish the control-source grid, retain at most
    two frozen ensembles, select one residual on chr16--18, open chr19--22/X
    only after the configuration hash is written, then run depth/replicate,
    independent naked-DNA, diagnostic transfer, ChromBPNet reference, and the
    new SK-N-SH/GM23338 holdout stages. Nutrient-stress analysis and Box
    before/after reports follow only after a candidate passes its safety and
    significance gates. Do not merge this research branch into `main`.

## Deferred or Experimental Work

- Broad supervised TFBS prediction and calibration.
- Motif-removal, motif-relaxed, or motif-free recovery claims.
- Variant-scoring case studies beyond the existing optional utility.
- Footprint competition or nucleosome-decomposition claims.
- Cross-study RNA/ATAC causal interpretation where study and cell state are
  confounded.
- Detectability-aware occupancy probabilities, abstention statuses, and
  TF-family hierarchical models. The benchmark schemas are available, but
  these remain experimental until the matched-label and perturbation gates in
  the footprint detectability study are satisfied.

These areas may have scaffolding under `benchmarks/` or internal tool modules,
but they are not promoted as primary package workflows until their datasets,
metrics, validation, tests, and documentation are complete.

## Maintenance Gates

Before pushing broad changes, run the full unittest suite, all console-script
help checks, YAML dry runs, `pip check`, strict MkDocs, release artifact checks,
and relevant focused regressions. Release, GitHub Actions, Pages, and PyPI
instructions belong only in `RELEASE_CHECKLIST.md`.

The GUI remains a thin wrapper included in the standard installation. Future
GUI changes must keep direct CLI use primary, retain reusable YAML, avoid
hosted-service assumptions, enforce the BAM/BAI plus peak BED starting point,
and use the current static demo and layout as the visual baseline.
