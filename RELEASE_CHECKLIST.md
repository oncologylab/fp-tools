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

```bash
.venv/bin/ATACorrect --help
.venv/bin/FootprintScores --help
.venv/bin/BINDetect --help
.venv/bin/PlotAggregate --help
.venv/bin/fp-tools-run --help
.venv/bin/fp-tools-gui --help
.venv/bin/fp-tools-run --config examples/gui_configs/plotaggregate_single.yml --dry-run
```

## 4. Build Artifacts

Build source and wheel artifacts:

```bash
./scripts/build_release.sh
```

Validate metadata when `twine` is available:

```bash
.venv/bin/python -m twine check dist/*
```

## 5. Metadata And Docs

- Confirm `pyproject.toml` version is correct.
- Confirm `project.urls` point to `https://github.com/oncologylab/fp-tools`.
- Confirm README renders on GitHub, especially the feature comparison table.
- Confirm `DEV_PLAN.md` remains local-only and that the former research report stays consolidated in `DEV_PLAN.md`.

## 6. Data Hygiene

Do not commit:

- downloaded public data under `data/public/raw/` or `data/public/processed/`
- benchmark result directories under `benchmarks/results/`
- generated paper figures/tables except intentional small examples
- BAM/BAI fixtures beyond existing local-only test data

## 7. Paper/Benchmark Gate

Before using outputs in a manuscript:

- freeze the exact public data manifest
- save command logs and environment versions
- save metrics tables used by each figure
- generate both vector and PNG figure outputs
- write Data Availability and Code Availability notes
