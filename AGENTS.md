# Repository Guidelines

## Project Structure

Core package code lives in `src/fp_tools/`. Command parsers are in
`src/fp_tools/parsers.py`, command entry points are thin wrappers, scientific
implementations live under `src/fp_tools/tools/`, and shared helpers live under
`src/fp_tools/utils/`. Small fixtures live in `test_data/`; example YAML configs
and lightweight examples live under `examples/`; public-data and benchmark
helpers live under `benchmarks/`; website documentation lives under `docs/`.

The `manuscript/` directory is intentionally local and ignored. Do not commit it
unless the user explicitly asks to publish the manuscript source.

## Architecture Rules

- Keep the package command-first. The GUI must call existing commands or stable
  Python entry points rather than moving workflow logic into Streamlit code.
- Direct CLI remains primary. YAML configs and the GUI are wrapper paths, not
  replacements for direct command use.
- GUI-saved YAML must remain runnable with `run-yaml-workflow`.
- Use the current public command names in docs and examples:
  `prepare-atac`, `bulk-footprinting`, `atac-correct`, `call-footprints`,
  `match-motifs`, `diff-footprints`, `normalize-bigwig`, `plot-aggregate`,
  `review-multi-comparisons`, `run-yaml-workflow`, `fp-tools-gui`,
  `discover-motifs`, `summarize-motifs`, `pseudobulk-fragments`,
  `find-signature-fp`, and `sc-footprinting`.
- Do not reintroduce the removed TOBIAS-style console aliases. Use the current
  command names above in code, tests, docs, and examples.

## Common Development Commands

Run commands from the repository root. Prefer ignored output directories under
`examples/`, `benchmarks/results/`, or `/tmp` for exploratory work.

### Environment And Status

```bash
pwd
git status --short --branch
git diff --stat
git diff --check
.venv/bin/python -m pip show fp-tools-bio
.venv/bin/python -m pip check
.venv/bin/python scripts/smoke_console_scripts.py
```

### Focused Smoke Examples

```bash
.venv/bin/atac-correct \
  --bams test_data/Bcell.bam \
  --genome test_data/genome.fa.gz \
  --peaks test_data/merged_peaks.bed \
  --blacklist test_data/blacklist.bed \
  --outdir examples/atacorrect/atac_correct_test_local \
  --cores 1

.venv/bin/call-footprints \
  --signals examples/atacorrect/atac_correct_test_local/Bcell_corrected.bw \
  --regions test_data/merged_peaks.bed \
  --outdir examples/scorebigwig \
  --score footprint \
  --cores 1

.venv/bin/match-motifs \
  --signals examples/scorebigwig/Bcell_footprints.bw \
  --genome test_data/genome.fa.gz \
  --peaks test_data/merged_peaks_annotated.bed \
  --outdir examples/motif_matches/Bcell \
  --motif-db jaspar2026_vertebrates \
  --cores 1
```

### YAML And GUI

```bash
.venv/bin/run-yaml-workflow --config examples/gui_configs/call_footprints_single.yml --dry-run
.venv/bin/fp-tools-gui --host 0.0.0.0 --port 8891 --run-dir examples/gui_runs
```

Keep YAML configs portable: paths and sample lists should be explicit, and the
same config should work through `run-yaml-workflow` without GUI-only state.

For release, package-build, full-test, GitHub Actions, website-deployment, and
PyPI procedures, use `RELEASE_CHECKLIST.md` as the source of truth. Do not
duplicate those steps here.

For documentation changes, preserve the existing four-page MkDocs navigation
and the two standalone demos. Run the strict MkDocs build and, when Playwright
is installed, `python scripts/audit_docs.py --site-dir site`.

## Coding Style

Follow the existing codebase style before introducing new patterns. Python
modules use `snake_case`; CLI functions should stay thin and delegate to tool
modules. Use 4 spaces for new Python code. Some legacy files use tabs; do not
reformat unrelated blocks. Put shared logic in `src/fp_tools/utils/` or
`src/fp_tools/tools/` rather than in wrappers.

Use `apply_patch` for manual file edits. `AGENTS.md` and `DEV_PLAN.md` are
tracked project guidance and must stay current. Do not commit transient agent
state such as `.claude/`, `.codex`, or local scratch plans.

## Testing Guidance

The test suite covers CLI smoke behavior, config expansion, public aliases,
docs contracts, regression fixtures, benchmark helpers, multiscale scoring,
motif discovery, pseudobulk generation, replicate reports, variant scoring, and
GUI/config helpers. Add tests in proportion to risk:

- narrow helpers: unit tests
- command contracts: parser or CLI smoke tests
- output behavior: fixture or regression tests
- docs/API changes: docs contract tests and `mkdocs build --clean --strict`
- package/release changes: build artifacts, `twine check`, and a fresh wheel
  console-script smoke test

## Data And Output Hygiene

Keep the project root limited to source, packaging metadata, small fixtures,
examples, benchmark scripts, and website docs. Store manual-run outputs in
ignored example directories, `benchmarks/results/`, `data/public/`, or `/tmp`.
Do not commit bulky public data, manuscript drafts, generated paper files, or
local GUI run outputs unless explicitly requested.
