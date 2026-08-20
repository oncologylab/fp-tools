<p align="center">
  <img src="docs/assets/fp_tools_logo_horizontal.svg" alt="fp-tools — regulatory footprinting" width="560">
</p>

<p align="center">
  <a href="https://pypi.org/project/fp-tools-bio/"><img alt="PyPI" src="https://img.shields.io/pypi/v/fp-tools-bio?color=1f9d55"></a>
  <a href="https://github.com/oncologylab/fp-tools/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/oncologylab/fp-tools/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://github.com/oncologylab/fp-tools/blob/main/LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-1967b3"></a>
</p>

<p align="center">
  <a href="https://oncologylab.github.io/fp-tools/"><strong>Documentation</strong></a>
  ·
  <a href="https://oncologylab.github.io/fp-tools/ENCODE-Cancer-Cell-lines-Footprinting/"><strong>Output demo with ENCODE cancer cell lines</strong></a>
  ·
  <a href="https://oncologylab.github.io/fp-tools/demos/gui/fp-tools-gui-static-demo.html"><strong>GUI demo</strong></a>
  ·
  <a href="https://pypi.org/project/fp-tools-bio/"><strong>PyPI</strong></a>
</p>

`fp-tools` is a command-first toolkit for footprinting Tn5-based chromatin
profiling data, including ATAC-seq, CUT&Tag, and related assays. It provides bias
correction, motif analysis, replicate-aware comparisons, and single-cell
footprint signatures. The GUI and YAML runner call the same commands.

## Install

Choose one route:

| Route | Best for | Start |
| --- | --- | --- |
| Desktop app | Windows or Apple Silicon macOS | [Download](https://github.com/oncologylab/fp-tools/releases) |
| Python package | Windows, macOS, or Linux with Python 3.11–3.13 | `python -m pip install fp-tools-bio` |
| Container | Complete reproducible environment | `docker build -t fp-tools:latest https://github.com/oncologylab/fp-tools.git#main` |

Python package example:

```bash
python -m pip install --pre fp-tools-bio
fp-tools-gui
```

External genomics programs are downloaded into a private, versioned cache on
first use. They do not need to be installed separately. Docker remains an
optional reproducible backend.

## Bulk ATAC-seq

`bulk-footprinting` can run from FASTQ files or public run accessions through
the final interactive comparison report.

```bash
bulk-footprinting \
  --reads-table reads.tsv \
  --comparison-table comparisons.tsv \
  --genome hg38 \
  --outdir project \
  --cores 8
```

The wrapper runs `atac-correct`, `call-footprints`, `match-motifs`,
`diff-footprints`, and `review-multi-comparisons`. Each command can also be run
directly. `diff-footprints --comparison-axis regions` compares matched genomic
region sets within one sample or across biological replicates.

## Single-cell ATAC-seq

`sc-footprinting` groups fragments, runs pseudobulk footprinting, and produces
per-cell KNN footprint-signature heatmaps and UMAPs.

```bash
sc-footprinting \
  --fragments fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --h5ad embedding.h5ad \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --outdir project/single_cell
```

## Main commands

| Area | Commands |
| --- | --- |
| Core analysis | `prepare-atac`, `atac-correct`, `call-footprints`, `match-motifs`, `diff-footprints`, `normalize-bigwig` |
| Workflows | `bulk-footprinting`, `sc-footprinting`, `run-yaml-workflow`, `fp-tools-gui`, `fp-tools-runtime` |
| Reports | `plot-aggregate`, `review-multi-comparisons` |
| De novo motifs | `discover-motifs`, `summarize-motifs` |
| Single-cell utilities | `pseudobulk-fragments`, `find-signature-fp` |

Use `<command> --help` for complete options. Practical examples and the API
reference are available in the [documentation](https://oncologylab.github.io/fp-tools/).
