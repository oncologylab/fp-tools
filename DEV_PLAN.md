# fp-tools Development Plan

Last updated: 2026-07-21

## Current Baseline

fp-tools is a command-first Python 3.12 package for bulk and pseudobulk
ATAC-seq footprinting, motif analysis, differential reports, and optional
browser/YAML wrappers. Scientific workflow logic belongs in
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
- Tn5 bias correction, background scaling, footprint scoring, candidate calls,
  known-motif matching, de novo motif preparation, and aggregate plotting.
- Replicate-aware differential footprint reports and multi-comparison reviews.
- Pseudobulk fragment/BAM generation and single-cell signature reporting.
- Optional Streamlit GUI whose saved YAML runs through `run-workflow`.
- Local tests, package builds, GitHub CI, GitHub Pages deployment, and manual
  PyPI publication.
- Branded GitHub and MkDocs presentation with responsive standalone report and
  GUI demos plus automated desktop/mobile browser audits.

## Near-Term Priorities

1. Keep README, MkDocs, examples, CLI help, and package metadata synchronized.
2. Add focused regression tests for every user-visible bug or command-contract
   change.
3. Reduce avoidable warnings and multiprocessing fragility without changing
   scientific output.
4. Improve GUI cancellation and output previews only through existing command
   and YAML interfaces.
5. Keep benchmark claims proportional to validated public data and recorded
   metrics.
6. Keep public pages free of broken assets, browser errors, responsive
   overflow, inaccessible navigation, and stale screenshots.

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

The GUI remains an optional wrapper. Future GUI changes must keep direct CLI
use primary, retain reusable YAML, avoid hosted-service assumptions, and use
the current static demo and layout as the visual baseline.
