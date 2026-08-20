# fp-tools Release Checklist

Use this checklist before publishing `fp-tools-bio` or preparing paper benchmark artifacts.

## 1. Environment

- Use a clean Python 3.12 environment for release orchestration. Wheels target
  Python 3.11–3.13 on Windows, macOS, and Linux.
- Confirm editable install:

```bash
.venv/bin/python -m pip show fp-tools-bio
.venv/bin/python -m pip check
```

## 2. Test Suite

Run the full local test suite:

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/smoke_console_scripts.py
```

Required coverage before release:

- public console scripts
- YAML config expansion and dry-run behavior
- command `--help` smoke checks
- core-count handling
- progress logging behavior
- stable fixture summaries for existing bigWig/BED test data

## 3. CLI Smoke Checks

Primary current API checks:

```bash
.venv/bin/prepare-atac --help
.venv/bin/bulk-footprinting --help
.venv/bin/atac-correct --help
.venv/bin/call-footprints --help
.venv/bin/match-motifs --help
.venv/bin/diff-footprints --help
.venv/bin/normalize-bigwig --help
.venv/bin/plot-aggregate --help
.venv/bin/review-multi-comparisons --help
.venv/bin/run-yaml-workflow --help
.venv/bin/fp-tools-gui --help
.venv/bin/discover-motifs --help
.venv/bin/summarize-motifs --help
.venv/bin/pseudobulk-fragments --help
.venv/bin/find-signature-fp --help
.venv/bin/sc-footprinting --help
.venv/bin/run-yaml-workflow --config examples/gui_configs/call_footprints_single.yml --dry-run
```

## 4. Build Artifacts

Release wheels are built by `cibuildwheel` in the manual `Publish` workflow.
The required artifact set is:

- CPython 3.11, 3.12, and 3.13
- Windows AMD64
- macOS x86_64 and arm64
- manylinux x86_64 and aarch64

Build an sdist locally only for preflight inspection:

```bash
.venv/bin/python -m pip install build twine
.venv/bin/python -m build --sdist
.venv/bin/python -m twine check dist/*
```

After uploading, verify a fresh install from PyPI:

```bash
python -m venv /tmp/fp-tools-pypi-smoke
/tmp/fp-tools-pypi-smoke/bin/python -m pip install --upgrade pip
/tmp/fp-tools-pypi-smoke/bin/python -m pip install --only-binary=:all: "fp-tools-bio==<version>"
/tmp/fp-tools-pypi-smoke/bin/atac-correct --help >/dev/null
/tmp/fp-tools-pypi-smoke/bin/plot-aggregate --help >/dev/null
/tmp/fp-tools-pypi-smoke/bin/fp-tools-gui --help >/dev/null
```

The `Desktop bundles` workflow must produce and smoke-test Windows x64, macOS
Intel, macOS Apple Silicon, Linux x64, and Linux ARM64 executables. Its live
health check verifies that the bundled Streamlit application starts, and its
command check verifies frozen child-process dispatch. Release assets include a
SHA-256 checksum manifest.

The `Container` workflow builds and tests the complete environment, then
publishes `linux/amd64` and `linux/arm64` images to
`ghcr.io/oncologylab/fp-tools` for tagged releases.

The manual GitHub Actions `Publish` workflow uses the repository
`PYPI_API_TOKEN` secret. Do not paste PyPI tokens into chat, shell history, or
committed files. Rotate any token that was exposed outside a secret manager.

The manual GitHub Actions `Publish` workflow builds and tests every wheel,
checks all artifacts, then uploads the complete set with the configured token.
Do not publish when any platform wheel is absent.

## 5. Metadata And Docs

- Confirm `pyproject.toml` version is correct.
- Confirm `project.urls` point to `https://github.com/oncologylab/fp-tools`.
- Confirm README renders correctly on GitHub.
- Validate the MkDocs site locally with `.venv/bin/mkdocs build --clean --strict`, then push documentation changes to `main`. The GitHub Actions `Docs` workflow deploys GitHub Pages from `main`; do not use `mkdocs gh-deploy` or create a `gh-pages` branch for this repository.
- When Playwright is installed, run
  `.venv/bin/python scripts/audit_docs.py --site-dir site` to check every
  documentation and demo page at desktop and mobile widths.
- Confirm `LICENSE`, `CITATION.cff`, `environment.yml`, and `Dockerfile` are present and current.
- Confirm `pyproject.toml`, `src/fp_tools/__init__.py`, `CITATION.cff`, release
  tag, and example version pins agree.

## 6. Data Hygiene

Do not commit:

- downloaded public data under `data/public/raw/` or `data/public/processed/`
- benchmark result directories under `benchmarks/results/`
- generated paper figures/tables except intentional manuscript previews and small examples
- BAM/BAI fixtures beyond existing local-only test data

## 7. Paper/Benchmark Gate

Before using outputs in a manuscript:

- freeze the exact public data manifest
- save command logs and environment versions
- validate benchmark manifests with `python benchmarks/scripts/validate_manifests.py --manifest-dir benchmarks/manifests`
- save metrics tables used by each figure
- label chromosome-4 benchmark results as pilot evidence unless whole-genome or chromosome-held-out validation has been completed
- generate both vector and PNG figure outputs
- write Data Availability and Code Availability notes
