# fp-tools Development Plan

Last updated: 2026-09-18

## Current Baseline

### v0.2.7 release preparation

Version metadata, desktop metadata, managed-runtime archive names, and example
version pins are synchronized to 0.2.7. The bulk manual explains the unique
condition-pair requirement and the separate-report option for repeated labels.
Download links remain on the verified prior release until new assets pass.

The seven GitHub workflows were reviewed. Removed the redundant `build`
reinstallation and the separate Windows I/O/launcher test invocation (those
tests remain in the full Windows suite). Retained all regression tests and
platform, frozen-app, managed-runtime, container, wheel, documentation-browser,
and public-consumer checks: each still covers a supported interface or artifact.
The opt-in slow correction regression remains available for scientific changes.
Release verification results will be recorded after the corresponding jobs finish.

Initial release CI found one host-dependent pseudobulk command assertion: it
expected eight cores on a four-core runner despite the new shared cap. That
explicit-limit test now controls the available budget at sixteen cores; the
separate automatic-selection and excessive-budget regressions remain enabled.

Windows CI also exposed CRLF expansion across nested workflow wrappers. Console
decoding now normalizes line endings incrementally before writing through the
platform's text stream; raw log bytes remain unchanged. A Linux reproduction
with a split CRLF failed before the fix and was added alongside the native
Windows nested-wrapper regression.

The tagged macOS desktop passed command dispatch and automatic-core checks but
exposed an audit helper that clicked an already-open example dropdown closed.
The helper now selects a visible option directly and only clicks to open a
closed dropdown. Browser regressions cover both Streamlit interaction styles.
This change affects validation infrastructure only; packaged source is unchanged.

### Unreleased issue #76: comparison preflight and report preservation

Bulk workflows with a combined review now reject repeated unordered condition
pairs (including reversed pairs) and self-comparisons before reference resolution
or output creation, including dry runs. Duplicate errors identify both comparison
IDs and their physical TSV line numbers. Separate reports selected with
`--review-format none` still permit repeated pairs and explicit sample subsets.

Static review bundles are built in temporary storage before replacing report
output. Invalid payloads, duplicate pairs, filename collisions, and invalid
default comparison/motif/plot options therefore preserve existing output and
leave no new partial report. This protects against input-validation failures;
it does not claim atomic replacement during disk or filesystem failures.

All 26 initial regressions failed before the fix. Final related regressions:
**57 passed in 2.35 seconds**. Full Linux/source suite: **621 passed, 1 skipped,
26 warnings in 251.98 seconds**. Both affected console commands passed help
checks, and `git diff --check` passed. Native desktop builds were not rerun.
No scientific calculations, research files, or manuscript files changed.

### Unreleased workflow usability improvements

Managed reference resolution now offers opt-in progress messages for downloads,
checksum verification, cache reuse, and index preparation. The existing hg38
and mm10 manifests and assembly-specific blacklists are unchanged. Custom
references still never infer a blacklist; explicit overrides and disabling
filtering retain their existing behavior. Reference regressions: **11 passed
in 0.38 seconds**, including corrupted downloads, cache reuse, overrides,
disabled filtering, and download-free dry runs. This is fixture-based source
verification; it does not claim a new full-genome download or native desktop run.

Workflow execution now resolves omitted core limits from process-visible logical
processors, including affinity restrictions. Bulk, single-cell, preprocessing,
YAML, and graphical forms share the automatic setting; explicit limits remain
available. Sample concurrency stays within the resolved total budget. Clearing
a loaded graphical core limit saves a portable null value and restores automatic
selection. Single-cell annotation validation still precedes output creation.

Bulk, single-cell, preprocessing, and YAML subprocess diagnostics stream live
while retaining their existing log files. Stages report commands, completion,
failures, and resume skips. Binary pipeline output and machine-readable counts
remain separate from console diagnostics. Frozen command dispatch flushes its
Python streams, and child processes inherit the launcher's process group.

Validation on Linux/source:

- Full regression suite: **593 passed, 1 skipped, 26 warnings in 252.42 seconds**.
  After the final process-group compatibility adjustment, execution, graphical
  job, runtime, and platform regressions passed **47 tests in 3.07 seconds**.
- A synthetic two-sample bulk workflow completed all five stages with automatic
  32 cores and with one core. Four signal tracks and differential numerical
  values agreed at relative tolerance 1e-5 and absolute tolerance 1e-6. Resume
  skipped all five stages without rewriting tracks.
- Read-only UCSC checks confirmed both pinned reference FASTA checksums against
  upstream checksum lists. Both downloaded blacklist checksums matched their
  manifests; all 636 hg38 and 3,435 mm10 intervals fit their assembly's chromosome
  sizes. No full reference genome was downloaded for this check.
- All console-script checks, **18 example YAML dry runs**, and `pip check` passed.
- The source graphical browser audit passed, including actual signature runs;
  an additional browser check passed clearing and saving a loaded core limit.

These are unreleased maintenance changes. Native Windows/macOS frozen execution,
full-scale single-cell analysis, and full raw-read preprocessing were not rerun.
Production scientific defaults, research work, and manuscript files are unchanged.

Documentation and both distributed bulk YAML examples now show managed hg38
without a fixed core count or a mismatched custom blacklist. Normal examples
omit fixed core limits; the complete command reference retains the optional
limit. Workflow guides explain automatic resources and live diagnostics.
Documentation/contracts/demo regressions passed **44 tests in 2.52 seconds**;
strict MkDocs and the browser audit of **33 pages at three viewports**, including
the dark color preference, passed. The command reference was regenerated from
the current parsers.

### Post-v0.2.5 display maintenance

The user's official Windows v0.2.5 retest passed #70–#73, actual example-sized
bulk/differential/signature/discovery runs, 543 tests with 22 skips and 69
subtests, and the GUI audit on rerun. This is user-reported Windows evidence,
not macOS or full-scale single-cell certification.

Issue #74 reproduced in the source GUI: the branding Markdown container's
negative 16-pixel bottom margin left the title bottom at 37.58px while Workspace
started at 30.84px. Removing that margin for the branding wrapper restores its
natural height. A browser regression failed before the fix and passed afterward,
checking containment and a minimum 4px gap on Home/Config, collapsed/expanded
Workspace, 1280/1920 widths, and 1/1.25/2 device scale and CSS zoom factors.

Issue #75 reproduced in the shared report formatter: exponent zero became an
incomplete `e+` suffix. JavaScript's native exponential output now retains the
zero. Formatting tests cover 0, 1, 0.05, 1e-10, NaN, infinity, and invalid text.
The focused regression failed before the fix; all 47 related report/GUI tests
passed afterward (6.75 seconds). The browser audit passed FDR=1 cards, volcano
tooltips, and downloaded motif-card SVG text, plus both existing aggregate and
aggregate-free plot-control fixtures. Numerical results are unchanged.

Both website demo formatter copies are synchronized. Native desktop smoke now
generates an FDR=1 report through the frozen command and checks cards, tooltips,
and the downloaded SVG in a real browser. The same command path passed with
the source adapter before native builds. Documentation/demo contracts passed
37 tests in 24.94 seconds. These fixes are being prepared for v0.2.6; the
published v0.2.5 artifacts are unchanged.

Release preflight passed 567 tests, one skip, and 26 warnings in 261.46 seconds.
The complete source GUI audit passed, including branding geometry and real
signature runs. Strict MkDocs and the 33-page, three-viewport docs audit passed.
Version metadata and runtime artifact names are synchronized to v0.2.6 while
external-tool pins and scientific implementations remain unchanged.

v0.2.6 is tagged at 7b37433. Its isolated editable build and dependency check
passed, followed by 22 installed release/report tests (28.37 seconds), all
console-script checks, and 18 YAML dry runs. The sdist passed twine and hygiene
checks (197 entries); 122 tracked package/packaging files passed the secret and
private-path scan.

CI run 34903980310 passed all ten jobs. The full Linux job reported 549 passed,
16 skipped, 26 warnings, and 87 subtests; each Windows Python 3.11–3.13 job
reported 537 passed, 28 skipped, one warning, and 69 subtests. Counts are kept
separate for each environment. Desktop run 34904025645 passed both platforms,
including frozen report generation, FDR=1 browser/SVG checks, branding geometry,
actual GUI signature runs, and macOS signing/DMG validation.

Managed runtimes (34904025641), containers (34904025640), and Publish
(34904026255) passed. PyPI has all 15 wheels and one sdist. All 13 public
binary/archive downloads match their checksums. Desktop SHA-256 values:

- Windows EXE: `21e1f07199f79477bcb5a90f7e13249b1ee736b49d9a9b9e2fa02041c65a5ba4`
- Apple Silicon DMG: `6479e871e1d4be068b1cf1c692cf5b5f6381b0987b655bcb13b45d7e9bd7117b`

Fresh PyPI checks passed console scripts, FDR browser/export regressions,
identical direct/list-YAML signature scores with SVG/PDF output, annotation
preflight, and first-use/cached Linux managed discovery. Release smoke run
34905697537 passed all three public-download jobs: macOS and Windows default
and extra-argument discovery, including fresh/cached runtimes and all five
required output types. No scientific implementation or shared scientific
helper changed relative to v0.2.5.

Final download-page contracts passed 37 tests in 29.99 seconds, followed by
strict MkDocs and the 33-page, three-viewport browser audit. Both live website
demo JavaScript files were independently checked for the corrected formatter.
Existing report files retain embedded JavaScript; regenerate their review HTML
to update formatting without repeating statistical analysis. The existing
optional STREME HTML-helper warning remains documented. Native checks use small
fixtures and do not certify full-scale single-cell analysis. Research work and
the published manuscript were not modified.

### v0.2.5 maintenance work

Issue #70: wrapper runtime replacement and path translation now stop at
`--extra-args`, preserving external arguments verbatim. Both focused
regressions failed before the fix; all 15 runtime tests passed afterward.
Native public Windows discovery remains a release gate.

Issue #71: managed WSL commands now initialize the guest executable PATH in a
noninteractive shell, with user arguments passed positionally. The delegation
regression failed before the fix; all 16 runtime tests passed afterward.

Issue #72: shared TSV/CSV schema validation requires `snap_cell_type` alongside
barcode, cell type, and UMAP coordinates. GUI preflight checks headers; built-in
signature and sc-footprinting commands reject before output creation. Custom
signature scripts keep their own schema. Three initial regressions reproduced
the defect; the expanded annotation/pseudobulk suite passed 17 tests. The GUI
browser audit now includes rejection and recovery before its real signature run.

Issue #73: sc-footprinting now selects JASPAR only when neither a database nor
custom motifs are supplied. GUI defaults expose custom motifs and leave the
database empty; explicit saved database choices remain explicit. Four motif
selection regressions failed before the fix and passed afterward, checking
both resolver inputs and generated differential commands.

Additional verification: 40 annotation/configuration tests passed (6 warnings,
65.35 seconds), and eight motif command/YAML regressions passed (2.11 seconds).
Release smoke now requires public Windows EXE discovery in two isolated runners
(default and extra-argument cases), with fresh WSL imports and cached repeats.

The first full run found two maintenance inconsistencies (559 passed, one skip):
the packaged single-cell example still carried the removed default, and the
command guide needed to show its newly documented custom-motif option. Both
were corrected; all 64 focused runtime/GUI/docs tests passed in 2.49 seconds.
The source GUI audit passed including annotation rejection/recovery and
identical real signature scores. The docs audit passed 33 pages at three
viewports; console smoke, dependency checks, and strict MkDocs also passed.

v0.2.5 release preparation synchronizes package, desktop, citation, example,
and runtime versions while preserving all external-tool pins. The isolated
editable build succeeded, all 18 YAML dry runs passed, and the source archive
passed twine validation and hygiene checks (197 entries). A scan of 122 tracked
package/packaging files found no secret or private-path patterns. Native and
public-consumer verification remain required before stable release.

Final local preflight passed 566 tests, one skip, and 26 warnings in 252.65
seconds. Both the complete source GUI audit (including dedicated-form and
Config annotation rejection/recovery) and the 33-page, three-viewport docs
audit passed. Version v0.2.5 is tagged immutably at e509a3a.

CI run 34887790047 passed all ten jobs: Linux full-suite results were 548
passed, 16 skipped, 26 warnings, and 87 subtests; each Windows Python 3.11–3.13
full-suite result was 536 passed, 28 skipped, one warning, and 69 subtests.
Counts are reported separately for each environment, including skips and subtests.
Desktop run 34888066236 passed both native platforms, including actual GUI
signature runs, parent/child completion, annotation preflight, SVG/PDF output,
and macOS signing/DMG verification. Public desktop checksums were verified:

- Windows EXE: `ae58dd25f9502760178378ac724868f2e40b2642a89544d1a15e780437a6a9ff`
- Apple Silicon DMG: `37f738d2d40a53e837b593991593957cc3bc27cdfa5ba6703a9913f9f157c8cd`

Managed runtimes run 34888066224 passed, including Windows WSL import.
Container run 34888066250 and Publish run 34888067597 passed. The complete
PyPI inventory contains 15 wheels and one sdist. All 13 public binary/archive
downloads (nine runtimes, two containers, two desktops) match their checksums.
Fresh PyPI first-use and cached Linux managed discovery passed all five
outputs without reinstalling the runtime. Fresh PyPI direct and list-YAML
signature runs produced identical score tables and SVG/PDF outputs; invalid
annotations were rejected before output creation.

Release smoke run 34889899040 passed all three public-download jobs: macOS
first-use/cached discovery, and two isolated Windows EXE jobs using default
parameters and STREME extra arguments. Each Windows job verified the EXE
checksum, fresh WSL import, cached reuse, paths with spaces, both runtime flag
spellings, and all five output types. This is end-to-end discovery evidence,
separate from the runtime's absolute-path STREME version check.
Download-page documentation/release contracts passed 37 tests in 29.44 seconds.
These tests use small deterministic fixtures; no full large-data single-cell
certification or new scientific performance claim is made. Production scoring
defaults, research work, and the published manuscript remain unchanged.

Nonblocking diagnostic: in the public Windows runs with paths containing
spaces, STREME's optional upstream HTML helper warned that its XML input was
not found (with additional Perl locale warnings). STREME text output, converted
motifs, Tomtom TSV, and fp-tools' own summary TSV/HTML were present and nonempty.
The optional STREME HTML file is not included in the five-output acceptance
claim; no claim of warning-free vendor-tool execution is made.

The v0.2.5 download links passed strict MkDocs and the final 33-page,
three-viewport browser audit before publication. Website deployment and issue
closure are tracked in GitHub rather than changing the immutable release tag.

### v0.2.4 marker-list maintenance release

Issue #69 was reproduced on clean main at 9a072cb: the GUI's YAML marker list
passed validation but expanded into multiple values for the single-value
`find-signature-fp --markers` argument. Six focused regressions failed before
the fix, with eleven controls passing. The shared serializer now joins only
this tool/flag's list into one comma-separated argument. Scalar strings,
defaults, explicit overrides, empty-list omission, other list arguments, and
the supplied configuration are preserved. No scientific implementation changed.

All 65 focused configuration, GUI, example, and wrapper tests passed in 64.18
seconds, including real list-YAML execution with byte-identical primary score
TSVs compared with direct scalar CLI execution. Both paths produced SVG/PDF
plots. The native smoke and browser audits now require this marker-list path,
including successful parent/child statuses from dedicated-form and Config runs.
The user selected a full v0.2.4 release; native/public-consumer verification
remains required before stable designation and closing #69.
The complete source GUI audit also passed against Streamlit 1.63, including
both real marker-list launches and identical scores. All console-script checks,
18 YAML dry runs, and `pip check` passed.
Release preflight passed 550 tests, one skip, and 26 warnings in 267.24 seconds;
strict MkDocs and the 33-page/three-viewport documentation audit also passed.
The source archive passed `twine check` with 194 entries and no manuscript,
environment, or agent-state directories. A scan of 119 tracked package and
packaging files found no secret/private-path patterns. Version metadata and
all nine runtime artifact names are synchronized to 0.2.4; external-tool pins
are unchanged, and website downloads remain on verified v0.2.3 until new assets
pass their checks.
The normal isolated editable build succeeded; the initial non-isolated attempt
used local Cython 3.2.5 outside the declared build range and was discarded.
The installed 0.2.4 package passed dependency and compiled-extension import
checks, 35 release/marker tests in 90.89 seconds, and 14 rebuilt-kernel regression
tests with one skip in 89.24 seconds. No build-requirement changes were needed.

All release artifacts use immutable tag v0.2.4 at eee84dd. CI run 34873541255
passed all ten jobs. Desktop run 34873565672 passed both native platforms,
including frozen list-YAML versus scalar-CLI scores, both real GUI launch
paths, parent/child completion, SVG/PDF output, and Mac DMG checks. Public
desktop SHA-256 values were independently verified:

- Apple Silicon DMG: `76ad18f42d51ba9446b0631592d9a7dcb7fba9ac06707ea139cb2e7e828d9cd1`
- Windows x64 EXE: `f4efbb56ae7d2b30e86f92b555f6ecc46e4d33339c804c9f311f45c15447b012`

Managed runtimes run 34873565698 passed, including a real Windows WSL2 import;
all nine public archives match their sidecars. Container run 34873565689 passed
both architectures, and both public archives match their checksum manifest.
Release smoke run 34874636136 passed real cold-cache and cached motif discovery
from the public Apple Silicon DMG.

Publish run 34873604321 delivered all 15 wheels and one source distribution.
A fresh Linux wheel-only PyPI installation passed dependency/console/YAML
checks, real marker-list YAML with byte-identical scalar CLI score tables and
SVG/PDF output, and cold-cache/cached managed motif discovery with an unchanged
cache marker. Research stayed at 8a4a958; the published manuscript was untouched.
The final v0.2.4 download links passed 37 documentation/release contracts in
29.43 seconds, strict MkDocs, and the 33-page audit at three viewport sizes.

### v0.2.3 maintenance release

The user authorized a new release of main for retesting on 2026-09-12.
Version 0.2.3 includes the fixes for #65–#68, frozen motif-discovery environment
repair, and the audience-oriented documentation update. Package, desktop,
citation, example, and managed-runtime versions are synchronized. Research
methods and the published manuscript are unchanged. The two documented
single-cell source issues remain known limitations of this maintenance release.

Release artifacts were built through the existing desktop, runtime, container,
and wheel workflows. The release remained a prerelease while native checks
were completed. Both public desktop downloads are now checksum-verified, and
the website links point to v0.2.3.

Local release preflight passed: 529 tests, one skip, 26 warnings in 206.21
seconds; all console-script smoke checks; 18 YAML dry runs; `pip check`;
strict MkDocs; and source-distribution build/`twine check`. The source archive
contains 192 entries with no manuscript, environment, or agent-state
directories. A scan of 121 tracked package/packaging files found no private
keys, credential-bearing URLs, token patterns, or personal workspace paths.

The first native desktop audit found three stale selectors for the pre-#67
"Comparisons TSV (optional)" label. The app correctly shows the required
"Comparisons TSV" input. A new consumer-contract regression reproduced the
mismatch; all 19 focused GUI/default/loading tests pass with the corrected
audit. This validation-only follow-up changes no packaged source. Desktop
bundles are rebuilt with the existing manual `publish_version=0.2.3` workflow;
the release tag and package/runtime sources remain immutable.

A second macOS audit reached configuration editing and exposed a timing-sensitive
audit read: Enter implicitly submitted the form before the scripted button click,
and the audit sampled values again after its retrying assertions had passed.
The audit now commits text on blur, submits once, waits for the new widget
revision, and uses retrying assertions without a second instantaneous read.
A real-browser regression reproduced the unintended submission. All 20 focused
checks and the complete source-GUI browser audit pass after correction; no
application code changed.

The Windows audit also exposed an asynchronous example-selector race. The
browser audit now waits for and selects the exact example option, verifies the
selection, and uses normal actionable button clicks for all three loaders.
The complete source-GUI audit and three repeated load/edit/example/upload
cycles pass; the 20 focused checks pass in 2.20 seconds. Native desktop checks
must still pass before these downloads are marked stable.

Final main regression verification passed: 531 tests, one skip, 26 warnings in
200.19 seconds. PyPI now contains all 15 version-0.2.3 platform wheels and the
source archive. A fresh Linux wheel-only installation passed `pip check`, all
console-script checks, a YAML dry run, and real managed STREME discovery with
JASPAR conversion, Tomtom, and TSV/HTML summaries on both first-use and cached
runs; the runtime cache marker was unchanged on reuse. All nine public runtime
archives and both container archives match their published SHA-256 checksums.

The next native audit isolated a dependency-specific test interaction:
Streamlit 1.63 uses a click-triggered React Aria dropdown, while the original
local environment used Streamlit 1.58. Filling text alone leaves the new menu
closed. This was reproduced against the fresh public PyPI install and in a
browser regression (one failure before correction). Explicitly opening the menu
fixes the test interaction; 21 focused checks and three repeated loader cycles
against Streamlit 1.63 pass, as does the complete GUI browser audit against the
fresh PyPI installation. Packaged application source remains identical to
the immutable v0.2.3 tag.

The Apple Silicon job in Desktop bundles run 34699598672 passed frozen command
dispatch, native startup, the full GUI audit, signing/quarantine checks, and DMG
verification/mounting. Its public DMG SHA-256 is
`1b40068d76aab3727bfe8dc5e1a5217f7e92f85e1289f234d1adb0f54409d3bd`.
Release smoke run 34700091438 passed real first-use and cached managed motif
discovery from that public DMG. All ten CI jobs in run 34699559584 passed.
The Windows job in the same desktop run also passed frozen command dispatch,
native startup, and the complete GUI audit. Both desktop bundles use commit
fa6e8b4 for validation, with no package-source differences from the release tag.
Desktop run 34699598672 completed successfully and attached both executables
and their checksum manifest. The public Windows executable SHA-256 is
`86c84896c6c3aea5c4e7dc58aaa34c4d37c77bb126fc01a6a41e4f5bc6f50ed2`.
The release contains 24 assets, including nine runtime/checksum pairs, two
container archives plus checksums, and the two desktop downloads plus checksums.
The updated desktop links passed 37 documentation/release contracts in 28.10
seconds, strict MkDocs, and the 33-page browser audit at three viewport sizes.
Final release rechecks found no open GitHub issues; the research worktree stayed
at 8a4a958 and the published manuscript was not modified.

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
