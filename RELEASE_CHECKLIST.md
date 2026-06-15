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
```

Required coverage before release:

- public console script aliases and legacy aliases
- YAML config expansion and dry-run behavior
- command `--help` smoke checks
- core-count handling
- progress logging behavior
- stable fixture summaries for existing bigWig/BED test data

## 3. CLI Smoke Checks

Primary current API checks:

```bash
.venv/bin/atac-correct --help
.venv/bin/call-footprints --help
.venv/bin/match-motifs --help
.venv/bin/diff-footprints --help
.venv/bin/plot-aggregate --help
.venv/bin/plot-aggregate-batch --help
.venv/bin/run-workflow --help
.venv/bin/motif-discovery --help
.venv/bin/motif-summary --help
.venv/bin/pseudobulk-fragments --help
.venv/bin/run-workflow --config examples/gui_configs/plotaggregate_single.yml --dry-run
```

Compatibility alias checks:

```bash
.venv/bin/ATACorrect --help
.venv/bin/FootprintScores --help
.venv/bin/ScoreBigwig --help
.venv/bin/BINDetect --help
.venv/bin/PlotAggregate --help
```

## 4. Build Artifacts

Build source and wheel artifacts:

```bash
./scripts/build_release.sh
```

The build script uses isolated `python -m build`. Validate metadata when `twine` is available:

```bash
.venv/bin/python -m twine check dist/*
```

## 5. Metadata And Docs

- Confirm `pyproject.toml` version is correct.
- Confirm `project.urls` point to `https://github.com/oncologylab/fp-tools`. The currently published PyPI `0.1.7` metadata still points to the old repository until a new release is published.
- Confirm README renders on GitHub, especially the feature comparison table.
- Confirm `LICENSE`, `CITATION.cff`, `.zenodo.json`, `environment.yml`, and `Dockerfile` are present and current.

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
