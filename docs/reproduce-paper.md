# Reproduce the Manuscript Analyses

This guide records the manuscript-level reproduction path for reviewers. It is
intended for the tagged source state `manuscript-revision-2026-06-22`.

## 1. Prepare the Environment

From the repository root:

```bash
python3.12 -m venv .venv
.venv/bin/pip install -e .
.venv/bin/python -m pip check
```

Optional container and conda-style environments are documented in `Dockerfile`
and `environment.yml`.

## 2. Run Fast Checks

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/atac-correct --help
.venv/bin/call-footprints --help
.venv/bin/match-motifs --help
.venv/bin/diff-footprints --help
.venv/bin/fp-tools-pseudobulk --help
.venv/bin/fp-tools-gui --help
```

These checks verify the installed command surface and compact regression suite.

## 3. Recreate Public Inputs

Large public inputs are not stored in Git. The repository tracks manifests,
schemas, and preparation scripts used to recreate them. The main public data
sources are ENCODE bulk ATAC-seq, 10x Genomics PBMC Multiome, 10x Genomics
PBMC5k scATAC, hg38 resources, and JASPAR 2026 CORE vertebrate motifs.

Representative entry points:

```bash
.venv/bin/python benchmarks/scripts/prepare_10x_pbmc_pseudobulk.py --write-example-archive
.venv/bin/python benchmarks/scripts/build_encode_manifest.py --help
.venv/bin/python benchmarks/scripts/download_manifest.py --help
```

Downloaded files and large generated outputs should remain under ignored
`data/public/`, `benchmarks/results/`, or example-output directories.

## 4. Rebuild Manuscript Figures

Figure builders and source tables are stored under `manuscript/scripts/`,
`manuscript/tables/`, and `manuscript/figures/`. The current main figures are
edited SVG/PNG manuscript assets. Figure 5 can be regenerated from its source
table with:

```bash
.venv/bin/python manuscript/scripts/plot_fig5_engineering_gui.py
```

Other figure-specific builders are retained for auditability in the manuscript
and benchmark script directories.

## 5. Compile the Manuscript

```bash
cd manuscript
latexmk -pdf -shell-escape -interaction=nonstopmode main.tex
rg -n "Undefined|Citation.*undefined|Reference.*undefined|LaTeX Warning" main.log
```

The expected manuscript output is `manuscript/main.pdf`.

## 6. Scope Notes

The manuscript is a software and reproducibility paper. The public-data examples
and local engineering measurements are intended to demonstrate workflow
operation, output provenance, and biological plausibility. They are not an
exhaustive cross-tool benchmark.
