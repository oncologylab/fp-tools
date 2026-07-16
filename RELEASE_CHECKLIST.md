# fp-tools Release Checklist

Use this checklist before publishing `fp-tools-bio` or preparing paper benchmark artifacts.

## 1. Environment

- Use Python 3.12 in the project virtualenv.
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
.venv/bin/atac-correct --help
.venv/bin/call-footprints --help
.venv/bin/match-motifs --help
.venv/bin/diff-footprints --help
.venv/bin/normalize-bigwig --help
.venv/bin/plot-aggregate --help
.venv/bin/review-multi-comparisons --help
.venv/bin/plot-motif-aggregate-grid --help
.venv/bin/run-workflow --help
.venv/bin/fp-tools-gui --help
.venv/bin/motif-discovery --help
.venv/bin/motif-summary --help
.venv/bin/fp-tools-score-variants --help
.venv/bin/pseudobulk-fragments --help
.venv/bin/find-signature-fp --help
.venv/bin/pseudobulk-footprints --help
.venv/bin/run-workflow --config examples/gui_configs/call_footprints_single.yml --dry-run
```

## 4. Build Artifacts

Build source and wheel artifacts. On Linux releases, install `auditwheel` and
`patchelf` first so `scripts/build_release.sh` can repair platform wheels to
manylinux tags accepted by PyPI:

```bash
.venv/bin/python -m pip install build twine auditwheel patchelf
./scripts/build_release.sh
```

The build script uses isolated `python -m build`, removes unrepaired
`linux_x86_64` wheels when a repaired manylinux wheel is produced, and leaves
the upload-ready files in `dist/`. Validate metadata before upload:

```bash
.venv/bin/python -m twine check dist/*
```

After uploading, verify a fresh install from PyPI:

```bash
python -m venv /tmp/fp-tools-pypi-smoke
/tmp/fp-tools-pypi-smoke/bin/python -m pip install --upgrade pip
/tmp/fp-tools-pypi-smoke/bin/python -m pip install "fp-tools-bio[gui]==<version>"
/tmp/fp-tools-pypi-smoke/bin/atac-correct --help >/dev/null
/tmp/fp-tools-pypi-smoke/bin/plot-aggregate --help >/dev/null
/tmp/fp-tools-pypi-smoke/bin/fp-tools-gui --help >/dev/null
```

The manual GitHub Actions `Publish` workflow uses the repository
`PYPI_API_TOKEN` secret. Do not paste PyPI tokens into chat, shell history, or
committed files. Rotate any token that was exposed outside a secret manager.

The preferred publish path is the manual GitHub Actions `Publish` workflow. It
repairs the Linux wheel, checks metadata, smoke-tests every declared console
script, and uploads with the configured token.

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
