# fp-tools Development Plan

Last updated: 2026-09-01

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
    formulations; freeze spline/FDA/GP/hybrid models on K562/HepG2; then score
    GM12878/IMR-90 exactly once. Keep the current DWM correction and footprint
    score as defaults unless every scientific and CPU performance gate passes.
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
