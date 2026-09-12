# fp-tools Development Plan

Last updated: 2026-09-12

## Current Baseline

### v0.2.3 release preparation

The user authorized a new release of main for retesting on 2026-09-12.
Version 0.2.3 includes the fixes for #65–#68, frozen motif-discovery environment
repair, and the audience-oriented documentation update. Package, desktop,
citation, example, and managed-runtime versions are synchronized. Research
methods and the published manuscript are unchanged. The two documented
single-cell source issues remain known limitations of this maintenance release.

Release artifacts are built from the versioned commit through the existing
desktop, runtime, container, and wheel workflows. Public downloads remain
marked as a prerelease until the required build and consumer checks pass.
The website's desktop links stay on the last verified release until the new
assets are available and verified.

Local release preflight passed: 529 tests, one skip, 26 warnings in 206.21
seconds; all console-script smoke checks; 18 YAML dry runs; `pip check`;
strict MkDocs; and source-distribution build/`twine check`. The source archive
contains 192 entries with no manuscript, environment, or agent-state
directories. A scan of 121 tracked package/packaging files found no private
keys, credential-bearing URLs, token patterns, or personal workspace paths.

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
- Checksum-verified managed hg38 and mm10 analysis references, including
  assembly-matched blacklists and FASTA indexes, shared by BAM-first bulk
  analysis and Linux preprocessing. Custom FASTA inputs never infer a
  blacklist.
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
- Reliable GUI background-job reconciliation based on retained child-process
  handles and batch results, including recovery after an application restart.
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
- Validated conversion of existing differential HTML reports and cache-safe
  aggregate-profile generation, with actionable failures instead of silently
  successful reports when requested aggregate data are unavailable.
- Frozen-safe motif-discovery dispatch, format-aware JASPAR/PFM-to-MEME
  conversion for Tomtom, verified TLS roots, and the exact Matplotlib raster,
  PDF, and SVG backends required by packaged applications.
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

### Unreleased maintenance: GitHub #68

YAML job expansion, core validation and GUI validation now share defaults
merging, with explicit per-job overrides. Required inputs, paths, choices and
GUI restrictions therefore validate the values that will actually run.
Five new regression cases pass; all 24 existing CLI/configuration tests pass.
No public command or scientific default changes.

### Unreleased maintenance: GitHub #66

Config text, uploads and path loading now report expected input errors without
replacing the previous configuration or truncating Save/Run controls. Failed
loads disable execution until a configuration is applied successfully.
Structural YAML errors use ValueError; unexpected implementation exceptions
remain visible. All 12 new GUI loading tests and five defaults tests pass.
The existing run-control unit test now uses the dictionary-style session-state
fixture, so absent errors behave like actual Streamlit state. All 39 GUI
launcher/loading tests pass, including the launch guard.

### Unreleased maintenance: GitHub #67

The bulk form labels the required input "Comparisons TSV", consistent with
existing CLI and config validation. The focused label/validation regression
reproduced the old mismatch and passes with the corrected label.

### Unreleased maintenance: GitHub #65

Legacy differential-report conversion and batch merging preserve selected
A/C/G/T motif matrices. Aggregate cards and logo-panel exports render their
information-content SVG logos, with existing SVG/PNG fallbacks. Invalid or
missing matrices do not invalidate aggregate profiles. Fourteen focused
tests pass, including five Chromium rendering/export cases; all ten existing
aggregate batch tests also pass. The bundled legacy report (1,019 motifs,
200 positions, pre-rendered PNG logos) retains its six displayed logos and
curves with no browser errors. Matrix-only behavior is verified separately
with synthetic fixtures because that bundled report has no motif matrices.
SVG logos share the image height constraint; browser regressions and visual
inspection confirm they fit inside the card instead of being cropped.

### Maintenance validation, 2026-09-12

Linux CI installs pytest explicitly and runs the complete pytest suite,
including the new function-based regressions and existing unittest cases.
This avoids importing pytest-dependent tests into an environment without the
test runner and ensures the new cases execute on Linux as well as Windows.

The Windows CI follow-up exposed AppTest leaking its temporary `__main__`
module into subsequent multiprocessing tests. GUI test teardown now restores
the original process entrypoint. The failure was reproduced locally by running
a GUI test followed by parallel normalization with the spawn start method;
the same sequence passes after isolation. Production multiprocessing code is
unchanged.

The full pytest suite passed: 506 passed, one skipped, 26 warnings in 210.95
seconds. All 24 aggregate tests passed after the final SVG sizing correction.
Console-script smoke checks, the call-footprints YAML dry run and pip check
passed. Strict MkDocs passed in 3.98 seconds; the browser audit passed all
32 pages at three viewport sizes. These are Linux/source checks; newly built
native macOS/Windows frozen applications have not been validated.

With the user's subsequent authorization, the maintenance fixes were pushed
through `7d54f47` and #65–#68 were closed after verification. All ten CI jobs
passed at that commit, as did the documentation build and Pages deployment.
GitHub was rechecked on 2026-09-12 and had no open issues. These remain
unreleased maintenance changes; the package version is still 0.2.2. Research
remains in its separate clean worktree for the next major version; the
published manuscript is unchanged.

### Audience-oriented documentation review, 2026-09-12

Reviewed all 29 website Markdown pages, including every public command guide,
plus README, example instructions, contribution/security guidance, and the
interactive demos. The guides now explain when to use a command, required
inputs, concrete output paths, and how to review results. Minimal bulk sample
and comparison tables remain source-neutral; the single-cell workflow now
describes required annotation columns and AnnData counts/features. Corrected
discovery execution instructions, signature flags and SVG output names,
normalization QC paths, and project-mode table/output descriptions. The API
reference is regenerated from the revised guides and installed command help.

The GUI demo now identifies example-only behavior, provides valid copyable
YAML, and uses source-neutral comparison names. Eight failing example checks
were reproduced before correction. Twenty-two new regression cases validate
all eleven displayed configurations, parse the resulting command arguments,
and dry-run without reading input files or creating run state. A small CSS
adjustment prevents long paths and YAML from widening mobile pages. Report
demos use clearer controls and errors; the region example's initial title
now matches its actual dataset. Scientific payloads, production algorithms,
command interfaces, and defaults are unchanged.

Validation: 45 focused documentation/demo tests passed in 2.17 seconds. The
full pytest suite passed with 529 passed, one skipped, and 26 warnings in
200.76 seconds. All 18 primary command examples parsed across the 17 guides;
all 18 repository YAML examples passed dry runs. Console-script smoke checks,
`pip check`, JavaScript syntax checks, and `git diff --check` passed. Strict
MkDocs passed in 1.03 seconds. The GUI demo's twelve routes passed desktop
and mobile browser checks with no JavaScript errors. The site-wide browser
audit passed all 33 pages at three viewport sizes, including all 17 command
guides; screenshots of the homepage, single-cell workflow, output report, and
mobile GUI demo were visually reviewed.
The 0.2.2 source distribution built successfully and passed `twine check`;
its 192 entries contain no manuscript, environment, or agent-state directories.
No package release or version change was made.

CI caught one stale release-checklist assertion after the final documentation
edit changed the test command from unittest discovery to pytest. The failure
was reproduced locally; the contract now requires `pytest -q` and retains all
other release-gate checks. All 62 documentation, demo, and release-metadata
tests passed in 30.27 seconds after this correction.

Two pre-existing single-cell behaviors need a separate source fix:

- `pseudobulk_footprints.py` sets a motif-database default before resolving
  custom motifs, so custom-only `sc-footprinting` inputs also include the
  default database. The revised guide makes no custom-only claim for this
  command.
- `find_signature_fp.py` reads annotations without checking `snap_cell_type`,
  then requires that column when attaching annotations. The guides now state
  the actual required columns; early validation still needs correction.

These are Linux/source checks. No native frozen desktop application was built
or validated during this documentation review.

### Ongoing priorities

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
10. Continue footprint-improvement experiments on the research branch for the
    next major version, following the current locked manifests on that branch.
    Keep nutrient data as a prospective application rather than training data.
    Retain production defaults until all prespecified scientific and
    computational promotion gates pass.

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

Before pushing broad changes, run the full pytest suite, all console-script
help checks, YAML dry runs, `pip check`, strict MkDocs, release artifact checks,
and relevant focused regressions. Release, GitHub Actions, Pages, and PyPI
instructions belong only in `RELEASE_CHECKLIST.md`.

The GUI remains a thin wrapper included in the standard installation. Future
GUI changes must keep direct CLI use primary, retain reusable YAML, avoid
hosted-service assumptions, enforce the BAM/BAI plus peak BED starting point,
and use the current static demo and layout as the visual baseline.
