# Tool overview

Choose a command by the task you want to perform. Each guide explains its
inputs, gives an example command, and lists the files it writes. See the
[Command reference](../api.md) for all command options.

If you are new to fp-tools, start with the [bulk ATAC-seq workflow](workflows/bulk-atac-seq.md):
it includes minimal sample and comparison tables and explains where to find
your results. For fragments and cell annotations, use the
[single-cell workflow](workflows/single-cell.md). Linux users starting from
FASTQ files can prepare BAM and peak files with `prepare-atac` first.

## Core analysis

Examples in command guides are templates: replace input paths with your files
and run from the folder containing those paths. Output paths are relative to
that same working folder. For a complete analysis, follow a workflow guide,
then use each command's “Main outputs” section to locate the result.

- [`atac-correct`](commands/atac-correct.md) — correct ATAC-seq cut-site signal for Tn5 sequence bias.
- [`call-footprints`](commands/call-footprints.md) — calculate footprint-score tracks from corrected signal.
- [`match-motifs`](commands/match-motifs.md) — scan motifs and summarize motif-associated footprint scores.
- [`diff-footprints`](commands/diff-footprints.md) — compare conditions or user-defined region sets, with replicate-aware statistics.
- [`normalize-bigwig`](commands/normalize-bigwig.md) — normalize corrected cut-site signals over shared background regions.

## Visualization and review

- [`plot-aggregate`](commands/plot-aggregate.md) — plot signal around motif sites or user-defined BED regions; also export motif-by-comparison PDF grids.
- [`review-multi-comparisons`](commands/review-multi-comparisons.md) — combine differential reports into one static comparison browser.

## Workflow and interface

- [`bulk-footprinting`](commands/bulk-footprinting.md) — run the complete bulk workflow from BAM/BAI and peak BED inputs.
- [`sc-footprinting`](commands/sc-footprinting.md) — run pseudobulk and per-cell single-cell ATAC-seq footprinting.
- [`run-yaml-workflow`](commands/run-yaml-workflow.md) — run saved command settings from a YAML file.
- [`fp-tools-gui`](commands/fp-tools-gui.md) — launch the browser interface.
- [`fp-tools-runtime`](commands/fp-tools-runtime.md) — check or install the external tools used for read preparation and motif discovery.

## Linux preprocessing

- [`prepare-atac`](commands/prepare-atac.md) — prepare FASTQ inputs as filtered BAM, peak, alignment coverage, and QC outputs from the Linux CLI or Linux container.

## De Novo Motif Discovery

- [`discover-motifs`](commands/discover-motifs.md) — discover motifs from footprint candidates.
- [`summarize-motifs`](commands/summarize-motifs.md) — summarize discovered motifs and known-motif matches.

## Single-cell ATAC-seq utilities

- [`pseudobulk-fragments`](commands/pseudobulk-fragments.md) — group fragments by cell annotation.
- [`find-signature-fp`](commands/find-signature-fp.md) — plot per-cell footprint-signature heatmaps and UMAPs.

## Practical glossary

| Term | Meaning in this guide |
| --- | --- |
| GUI | Graphical user interface: forms and buttons for setting up analyses. |
| CLI | Command-line interface: commands entered in a terminal. |
| YAML | YAML Ain't Markup Language: a text format used to save analysis settings and paths. |
| BAM / BAI | Binary alignment/map file containing aligned reads, and its matching index for efficient access. |
| BED | Browser extensible data: genomic intervals, with zero-based starts and end positions excluded. |
| bigWig | A genome-wide numerical signal track. |
| FASTA | A sequence-file format; reference DNA must match your inputs' genome assembly. |
| AnnData / h5ad | Annotated data stored in an HDF5 file; single-cell analysis needs matching barcodes and genomic-bin counts. |
| Pseudobulk | Fragments combined from cells in the same annotated group. |
| KNN smoothing | K-nearest-neighbor smoothing: combining information from nearby cells to reduce sparse measurements. |
| UMAP | Uniform manifold approximation and projection: a low-dimensional view of similarities between cells. |
| TF | Transcription factor; a motif match alone does not establish which protein is bound. |
| FDR | False discovery rate: the expected proportion of false discoveries among selected statistical results. |
| DWM | Dinucleotide weight matrix: the sequence-bias model used by the default correction. |
| Managed runtime | Supporting analysis programs that fp-tools downloads and prepares when needed. |

Logging verbosity is a number from `0` (silent) to `5` (all diagnostic messages).
Help calls level `5` “spam”; enter the number, not that word.
