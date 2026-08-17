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

`fp-tools` is a command-first toolkit for ATAC-seq bias correction, footprint
scoring, motif analysis, replicate-aware comparisons, and single-cell
footprint signatures. The optional GUI saves YAML that remains runnable with
`run-workflow`.

## Install

```bash
pip install fp-tools-bio
```

For the browser GUI:

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

Raw-read processing also requires the genomics tools included in
`environment.yml` or the project Docker image.

## Standard workflow

```text
prepare-atac → atac-correct → call-footprints → match-motifs → diff-footprints
```

For processed BAM and peak files, use a sample table:

```text
sample	condition	bam	peaks
A1	conditionA	A1.bam	A1_peaks.bed
B1	conditionB	B1.bam	B1_peaks.bed
```

```bash
atac-correct \
  --sample-table samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project

call-footprints \
  --sample-table samples.tsv \
  --regions project/peaks/merged_peaks_filtered.bed \
  --outdir project

match-motifs \
  --sample-table samples.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project

diff-footprints \
  --sample-table samples.tsv \
  --comparison-table comparisons.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

Repeated condition labels define biological replicates. Output includes motif
tables, replicate statistics, aggregate profiles, and portable HTML/SVG
reports.

## Single-cell workflow

Use `pseudobulk-fragments` for grouping only, or `pseudobulk-footprints` for
the complete grouped workflow. `find-signature-fp` produces per-cell footprint
signature heatmaps and UMAPs.

```bash
pseudobulk-footprints \
  --fragments fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project/pseudobulk
```

## Optional de novo motif discovery

Candidate footprints from `call-footprints` can be passed to
`motif-discovery`, then summarized with `motif-summary` or included in
`match-motifs` and `diff-footprints`.

```bash
motif-discovery \
  --candidates candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo
```

## Commands

| Command | Purpose |
| --- | --- |
| `prepare-atac` | Prepare raw ATAC-seq reads. |
| `atac-correct` | Correct Tn5 sequence bias. |
| `call-footprints` | Create footprint-score tracks. |
| `match-motifs` | Scan motifs and call bound sites. |
| `diff-footprints` | Compare conditions or replicates. |
| `normalize-bigwig` | Scale signal tracks. |
| `plot-aggregate` | Plot motif-centered profiles. |
| `review-multi-comparisons` | Review multiple comparison reports. |
| `plot-motif-aggregate-grid` | Export aggregate comparison grids. |
| `motif-discovery` | Run optional de novo motif discovery. |
| `motif-summary` | Summarize discovered motifs. |
| `pseudobulk-fragments` | Group single-cell fragments. |
| `pseudobulk-footprints` | Run the complete pseudobulk workflow. |
| `find-signature-fp` | Plot single-cell footprint signatures. |
| `run-workflow` | Run a saved YAML configuration. |
| `fp-tools-gui` | Open the optional browser GUI. |
| `fp-tools-score-variants` | Score sequence variants in motifs. |

Check command syntax directly:

```bash
prepare-atac --help
atac-correct --help
call-footprints --help
match-motifs --help
diff-footprints --help
normalize-bigwig --help
plot-aggregate --help
review-multi-comparisons --help
plot-motif-aggregate-grid --help
run-workflow --help
fp-tools-gui --help
motif-discovery --help
motif-summary --help
fp-tools-score-variants --help
pseudobulk-fragments --help
find-signature-fp --help
pseudobulk-footprints --help
```

Example YAML files are in `examples/gui_configs/` and run identically through
the GUI or `run-workflow`.
