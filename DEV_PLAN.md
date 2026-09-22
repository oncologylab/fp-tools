# fp-tools Development Plan

Last updated: 2026-09-22

## Next-major-version research goal

**Corrected safety protocol:**
`benchmarks/manifests/corrected_count_safety_protocol_v1.json` fixes the next
negative-control reassessment before reading corrected control scores. Both
calibration and naked-DNA replicate-2 profiles used the affected extractor:
the archived inputs contain 36,931/26,647 development sites and 2,400/1,800
control sites for HepG2/K562, respectively. Sequence log-bias arrays are
available, so fresh observed counts can be combined with the same model
without fitting or new labels. Fixed historical thresholds isolate extraction
effects; independently recalibrated thresholds use only original tuning
negatives and the original rule. Both analyses are mandatory, with the latter
designated the primary development screen before outcomes are known.
Zero baseline expectation, profile disagreements and changed denominators
must remain explicit. This protocol is not an execution freeze or safety
result. `benchmarks/scripts/rebuild_corrected_safety_profiles.py` now provides
separate freeze, extraction and verification stages for the four archived
development/control inputs. It preserves site order, strand orientation,
original validity and fresh validity separately, and refuses changed inputs
or incomplete output reuse. Five focused tests passed, including invalid
sequence rows and corruption detection. The real input freeze/extraction and
corrected safety scorer remain outstanding.

The subsequent real four-input freeze completed successfully. Configuration
SHA-256: `dc053ee9d0e0d6154cabe0e3915a6a78ba39887ffe28e2850a5059db4aa3a9e9`.
Fresh extraction is now supervised with receipts under
`benchmarks/results/corrected_safety_jobs_20260922/execute_attempt1/`, with
that configuration hash checked before execution. Outputs are under
`benchmarks/results/corrected_safety_profiles_20260922_v1/`. Eighteen combined
safety-profile, corrected-count and extraction regression tests passed in
3.56 seconds. Completion and independent verification remain pending; these
profiles are not themselves safety scores or a passed safety gate.

**Current execution, September 22:** the corrected-count freeze completed
successfully (exit 0; 761.65 seconds), with configuration SHA-256
`59dd9a562d31e49ebd7bd2303d1fc458412f8170dad174d7292696ccc13d7b32`.
The unchanged-candidate 96-run reassessment was launched at 08:06:38 UTC
under the durable experiment supervisor. Its command checks that exact
configuration hash before executing, and the runner checks frozen calculation
code and runtime. Job specification SHA-256:
`bcf41ee29fbdc7635648eab950e59ef2a7cf57f1d0ee6ac017a0a6a83bac8f74`.
Receipts are in
`benchmarks/results/corrected_count_jobs_20260922/execute_attempt1/`;
scientific outputs are in
`benchmarks/results/corrected_count_reassessment_20260922_v1/`.
Launch is not completion: require the exit receipt, independent
`--verify-only` verification, paired comparison, and aggregate-plot review
before interpreting results. This remains post-hoc implementation reassessment,
not candidate tuning or untouched validation. Main remains clean at `636db5f`.

**September 22 extraction qualification:** the research strand-profile reader
had duplicate-window and soft-clipping count defects. They are now fixed and
confirmed against an independent training-only oracle. Historical stability,
reference and safety computations remain preserved, but their scientific
interpretation requires corrected-count reassessment; checksums alone do not
establish extraction correctness. See
`benchmarks/STRAND_EXTRACTION_AUDIT_20260922.md` before reusing their conclusions.

Work exclusively in the research worktree on
`research/footprint-improvement-20260830`. Improve footprint detection and
interpretable protection across transcription factor (TF) families while
measuring complete-workflow CPU time and peak memory. Main, released packages,
the published manuscript and the production dinucleotide weight matrix (DWM)
default remain unchanged. No integration or publication is authorized here.

September 22 status audit: research HEAD is `8a4a958`, with 18 local commits
after remote research HEAD `888c6e2`. Main is clean at `636db5f`. No relevant
experiment process was found. The completed HepG2 ChromBPNet reference and
96-run study supersede the older queued-stage handoff. The three primary
96-run manifests and 198 referenced stage/summary files match their recorded
checksums. This is a manifest/stage-summary audit, not a fresh verification of
every nested BAM/profile input or a scientific rerun. The original holdout
freeze matches its recorded hash and all 13 `labels_read` flags remain false;
no holdout label contents were opened.

The completed 96-run matrix is post-hoc development evidence on already
examined chromosomes. CTCF point gains over raw and sample-specific DWM are
positive throughout, but canonical protection passes 94/96 runs and some
low-depth bootstrap intervals cross zero. No other powered matched TF has
canonical protection. Independent exact-operator negative-control support is
limited to CTCF. Sequence inference uses about 47% of DWM CPU time, but this
excludes extraction, fitting, normalization, scoring and output. End-to-end
resource gates, calibrated functional uncertainty, occupancy-leakage exclusion,
two difficult families and untouched validation remain unestablished.

Proceed in this order:

1. Build a per-TF gate ledger from the existing reports, preserving original
   task denominators, coverage exclusions and evidence splits. Reverify nested
   inputs and the recorded runtime before reuse; do not rerun valid stages.
2. Freeze a fair complete-workflow resource protocol with identical inputs,
   outputs, core budgets and cache conditions. Measure fitting separately from
   reuse, CPU and elapsed time, and absolute peak memory in isolated processes.
   Require numerical/shape equivalence for implementation-only optimizations.
3. Diagnose the low-depth CTCF exceptions and difficult-family failures.
   Existing test results are descriptive only. Fit new candidates on chr1–15
   and select on chr16–18; never retune using examined test or transfer labels.
   Separate correction from scoring geometry in a frozen factorial comparison.
   Prefer simple parametric or hierarchical-spline models before more complex
   Gaussian-process models. Treat insufficient labels as a power limitation.
4. Freeze each candidate's code, data, model, thresholds and evaluation before
   scoring. Audit signed central protection, shoulders, asymmetry, coverage,
   motif/accessibility balance and exact-operator negative-control behavior.
   Require the existing scientific, uncertainty and computational gates in full;
   narrow CTCF gains and faster sequence inference cannot waive them.
5. Keep SK-N-SH/GM23338 labels and final naked-DNA replicate 3 locked until the
   prespecified prerequisites pass. Record failed candidates without promotion.
   Any later production integration requires separate authorization.

The exact current scientific results and limitations are in
`benchmarks/DIRECT_BIAS_STABILITY_96_20260912.md`; component resource results
are in `benchmarks/DIRECT_BIAS_INFERENCE_20260908.md`. No new scientific
experiment was launched during this status-and-goal setup.

### September 22 gate-ledger implementation

The new benchmark-only ledger validates the checksummed summary tables and
preserves all 21 tasks and 48 runs per task. HepG2 HMGA1 has one matched site
per class in the coverage table and no scored rows; it is retained explicitly
instead of disappearing from the denominator. Principal interval availability,
power, canonical protection and unresolved resource/uncertainty/safety gates
are separate fields. No new scientific score or candidate selection is made.

The ledger reproduces 46/48 canonical CTCF runs in HepG2 and 48/48 in K562.
The two failures are HepG2 replicate 2, 10M, seeds 2026 and 2028: both have one
significant central point, with central contrasts -0.054764 and +0.003386.
Both retain positive shoulders. K562 has 5/6 principal sample-DWM comparisons
with both lower bounds positive; HepG2 has 6/6. These observations describe
examined test data and cannot be used to tune a new candidate.

Forty-two focused ledger, stability-extension and report tests passed in
26.28 seconds. The real ledger and its exact immutable rerun passed, with
SHA-256 `33d65e9025ccb23ef09d640bf2e8757107547ba62656297f9c5c360693ed4d88`.
Output: `benchmarks/results/direct_bias_gate_ledger_20260922_v1/gate_ledger.json`.
The full historical nested-input verifier subsequently exited 0 and verified
all 96 runs, frozen inputs, archived reuse and outputs without rerunning a
stage (`--verify-only --code-ref c9ecd82`). The browser audit also passed all
32 pages at three viewport sizes with a dark-system preference. Historical
verification checks recorded code; it does not certify an identical current
scientific runtime for a new experiment.
No scientific computation, holdout access or production change occurred.

### Complete-workflow resource evidence contract

The research promotion checker previously accepted a generic runtime/memory
table without end-to-end scope or matched workload identities. Twelve new
regressions failed before correction. Resource evidence now requires total
process-tree CPU time, absolute process-tree peak memory, fitting plus
application, frozen benchmark/workload identities, matching core/cache settings
and all six repetitions per workload. It rejects nonfinite/negative values,
duplicate or unmatched runs and missing repetition pairs. Worst workload
ratios prevent a fast case from hiding a slow one. Wall time remains recorded
but cannot substitute for CPU time. Legacy component tables fail this gate.
These checks validate table consistency, not the authenticity of declarations;
the benchmark producer must still verify its frozen matrix and artifacts.

The premeasurement protocol is
`benchmarks/manifests/end_to_end_resource_protocol_v1.json`, SHA-256
`eeff0c5ecf2bee1741c8719205a8072b9f736e93d08fb0ef230c9df5b38adc2b`.
It specifies both development cells, all three biological replicates at 10M
and full depth, six timing repetitions and a shared four-core budget. It
includes fitting, input reading, correction, scoring and equivalent output
products. Frozen-model reuse and same-operator correction comparisons remain
separate diagnostics. This locks the protocol, not an executable benchmark:
exact command/input/code/runtime hashes, output-contract parity and validated
process-tree instrumentation are required before measurement. No runtime
ratio or computational promotion pass is claimed. The host's cgroup root is
not writable by the current user; a suitable isolated memory-accounting path
must be verified rather than substituting sampled or parent-only peaks.
After adding the six-repetition guard, all 31 focused resource, promotion and
ledger tests passed in 1.00 second. All 20 documentation-contract tests and the
strict MkDocs build passed. No timing or scientific model run was started.

### Process-tree measurement qualification

The cgroup root is not writable, but the existing user systemd manager can
create a dedicated delegated service. A benchmark-only synthetic probe now
moves its supervisor outside a fresh measured subtree, places the command in
that subtree before exec and reads the kernel's `cpu.stat` and `memory.peak`
after the command and descendants exit. It rejects leftover descendants and
failed commands, retains their measurements and cleans up only its own units.
The measured memory is kernel-accounted cgroup memory, including charged file
cache and kernel memory, not resident-set size. Both benchmark arms must use
this same metric and the frozen cache policy.

The qualified four-case probe used no biological data or model. Two sequential
32-MiB child allocations used 0.498529 CPU seconds and peaked at 42.209280 MB;
two simultaneous allocations used 0.497108 CPU seconds and peaked at 79.167488
MB. An exit-7 command and a parent that left a child running were both rejected.
All probe services were removed. These results distinguish descendant CPU and
simultaneous memory from parent-only or cumulative measurements; they are not
scientific speed or memory improvements. All 23 focused instrumentation and
resource-gate tests passed in 0.47 seconds.

Evidence: `benchmarks/results/process_tree_probe_20260922_v2/probe.json`,
SHA-256 `74fa8f93de827bffc2a948a0da66c7db80c3fd504d771d134c29912e46c2fce0`.
The earlier successful three-case probe remains preserved separately. The
complete-workflow execution freeze still requires scientific command adapters,
equivalent output products, numerical parity, input identities and the exact
runtime. No candidate timing, fitting, holdout reading or promotion occurred.

### Corrected research cut extraction

The BAM-to-score parity work exposed a long alignment being counted twice
when fetched by disjoint windows, and missed soft-clipped cuts outside queried
alignment spans. Three boundary regressions failed before the fix. The reader
now traverses each requested chromosome once and assigns each cut once per
containing site, preserving legitimate overlap and positional row order.
Forty focused regressions passed. No production code was changed.

A frozen, label-free chromosome-1 audit compared old and corrected extraction
on both development cells against independent sorted-cut counting. Exact
oracle parity passed at 2,381 HepG2 and 1,943 K562 sites. Counts changed at 35
and 15 sites, removing 173 and 45 duplicate assignments respectively. This
does not quantify the impact on test metrics or safety. Archived artifacts
remain unchanged, and no model was retuned. The new audit report records input,
code and result hashes and the requirements for a separately frozen corrected
reassessment. Runtime measurement must use the corrected reader.

### Fixed-candidate corrected-count reassessment

`reassess_corrected_strand_counts.py` freezes a new 96-run computation using
the original matrix, exact site/label order, complete matched pairs, model
propensities, dispersions, scoring geometry and bootstrap settings. Both cells'
contexts are frozen before corrected extraction. Every BAM is re-extracted;
no archived observed counts are reused. Existing valid DWM expectations are
reused with verified coordinate order. New outputs live separately from the
historical experiment and retain explicit post-hoc/non-promotion status.

Fifty-one focused reassessment, extraction and stability tests passed in
28.70 seconds. Fixtures cover all 96 runs, verification without recomputation,
corrupted output rejection, preserved incomplete artifacts, changed settings
and matrix rejection, and prohibitions on retuning, holdout reads and old-count
reuse. The new runner requires a configuration freeze before execution and
binds resumed stages to their configuration, run identity and checksums.
It does not establish a new method or close the independent-control gates.
The supervised freeze was launched with durable receipts under
`benchmarks/results/corrected_count_jobs_20260922/freeze_attempt1/` and a new
scientific output root `benchmarks/results/corrected_count_reassessment_20260922_v1/`.
Its job specification SHA-256 is
`3cdaec83fa33c2a613ba5af8e084933caacfcd58a1243c099c853b07842e35e4`.
Check the exit receipt before launching execution; launch evidence alone does
not certify the freeze. At the launch audit, library versions matched the
historical run and only the Python 3.12.3 build timestamp differed; the new
configuration records the current runtime. All 20 documentation-contract
tests and the strict MkDocs build passed after adding the runner.

The before/after reporter `compare_corrected_count_reassessment.py` is separate
from the pinned calculation. It requires the full corrected run matrix and
paired metric/shape/coverage/bootstrap identities, rejects changed class counts
or eligibility, retains unscored tasks via the coverage denominator, and emits
paired tables plus corrected aggregate plots. Six focused tests passed in
1.85 seconds. Reported reader-effect differences are descriptive; they are not
confidence intervals for the implementation change or a new method result.

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
  76,393,585/114,471,886. All support the requested 50M depth. All four sampling
  stages completed on 2026-09-09 UTC: independent checks passed for the exact
  60-target matrix and all 180 BAM/index/metadata files, including checksums,
  parent identities, metadata and indexed alignment counts. Conventional DWM
  correction and the 96-run scientific evaluation completed on September 10.
  Independent verification passed for all 60 corrections and all 96 runs on
  September 12; no completed scientific stage was rerun.
  See `benchmarks/DIRECT_BIAS_REPLICATE_DEPTHS_20260908.md` for the frozen
  identities, actual selected-depth ranges and remaining requirements.
- The first sampling process terminated with exit 143 before completing any
  sampling stage; the dependent queue correctly failed without starting an
  evaluation. The termination cause is unknown. All four source BAM/BAI pairs
  and count stages reverified, and 30 interrupted BAMs (54,098,858,454 bytes)
  were preserved with checksums and unchanged file identities. No partial BAM
  was accepted as complete. A separate logged-job supervisor now pins the
  command specification, records supervisor/child identities and terminal
  status, and forwards stop signals to its own child group. Detached launch
  survives the short launcher process; receipt failure stops an unrecorded
  child. Resume retains the original preparation/evaluator code, source tree,
  runtime, inputs, scientific settings and completed counts.
- The complete 96-run evaluator now extends the archived 36-run calculation
  with the 60 prepared replicate/depth inputs. It verifies nested sampling and
  correction provenance, pins the existing scientific code and runtime, and
  freezes both cells' all-run expectation coverage before new extraction or
  scoring. All 36 archived profiles are reused after coordinate/label and
  DWM-track checks. Metrics are reused only on identical original and matched
  cohorts; changed coverage requires a fresh comparison on the new support.
  The plan pins the evaluator checksum, and queued execution also pins the
  plan checksum. A small 96-run fixture exercises unchanged and changed
  coverage, complete pairs, zero-cut sites, reuse, verification and nested
  output corruption. The real calculation is complete: 36 archived
  extractions were reused and all 96 evaluations used the newly frozen common
  support. Matched CTCF AUROC/AUPRC point gains over raw and sample DWM are
  positive in all runs, but protection passes only 94/96 and some low-depth
  bootstrap intervals cross zero. No other powered matched TF has canonical
  protection. All 21 aggregate pages were reviewed. The result does not pass
  promotion gates; holdouts and final naked-DNA controls remain unopened.
  See `benchmarks/DIRECT_BIAS_STABILITY_96_20260912.md` for verified hashes,
  per-TF results, balance limitations and the historical-runtime audit.
- A separate current-depth balance diagnostic addresses the observed cut-total
  drift in the 36-run archive. Its fixed secondary rule retains original
  partners only within the existing 0.2-SD caliper on current log1p cut totals,
  using the original chr1–15 per-cell/TF scales. All selections precede new
  scores; underpowered and zero-pair outcomes remain explicit. Primary results
  and the queued 96-run evaluation are unchanged. This is a post-hoc
  sensitivity check, not a new model or causal validation. The completed
  36-run calculation balances all measured CTCF covariates below 0.04 SMD.
  All six full-depth CTCF runs retain significant gains over both DWM
  comparators, but HepG2 replicate 1's AUROC interval against raw includes
  zero, and several lower-depth raw/protection checks fail. FOXA1/MYC lose
  power after filtering; excluding them must not change the task denominator
  used to claim a promotion gate. CTCF remains the only powered canonical
  task. See `benchmarks/DIRECT_BIAS_DEPTH_BALANCE_20260908.md` for the complete
  result, 21-page figure audit, frozen identities and limitations.
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
