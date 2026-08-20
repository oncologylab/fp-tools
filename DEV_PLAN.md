# fp-tools Development Plan

Last updated: 2026-08-20

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

- Raw-read ATAC preparation with modern and legacy processing profiles.
- A managed, versioned external-tool runtime downloaded on first use, with
  native Linux/macOS archives, a private WSL2 runtime for Windows, and the
  complete container retained as an explicit backend.
- One-command bulk analysis from FASTQ files or public run accessions through
  portable multi-comparison HTML reports.
- Tn5 bias correction, background scaling, footprint scoring, candidate calls,
  known-motif matching, de novo motif preparation, and aggregate plotting.
- Replicate-aware differential footprint reports, including per-sample motif
  matrices and empirical-Bayes residual-variance moderation with biological
  samples as the inferential units, plus bundle and aggregate-free standalone
  multi-comparison reviews.
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
- Self-contained GUI executables for Windows x64, macOS Intel and Apple
  Silicon, and Linux x64 and ARM64, with frozen-safe command and workflow
  dispatch.
- Branded GitHub and MkDocs presentation with responsive embedded report and
  GUI demos plus automated desktop/mobile browser audits.
- A manifest-pinned, storage-conscious ENCODE workflow for 17 biological
  replicates from seven cancer cell lines. Its resumable runner uses a
  pair-specific union of released IDR-thresholded peaks, peak-q95 scaling, and
  the same differential workflow as the preserved three-replicate
  K562-HepG2 report. The dependency-free static browser has two directional
  selectors and reads compact canonical report payloads.

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
9. Maintain binary-wheel and desktop-app scientific I/O parity on Windows,
   macOS, and Linux. Test managed-runtime artifacts and keep the complete
   multi-architecture container as a reproducible alternative backend.

## Deferred or Experimental Work

- Broad supervised TFBS prediction and calibration.
- Motif-removal, motif-relaxed, or motif-free recovery claims.
- Variant-scoring case studies beyond the existing optional utility.
- Footprint competition or nucleosome-decomposition claims.
- Cross-study RNA/ATAC causal interpretation where study and cell state are
  confounded.

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
hosted-service assumptions, and use the current static demo and layout as the
visual baseline.
