# Get Started

`fp-tools` provides command-line and optional browser workflows for ATAC-seq
footprinting, motif analysis, replicate comparisons, and single-cell
footprint signatures.

## Install

```bash
pip install fp-tools-bio
```

For the GUI:

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

## Bulk ATAC-seq workflow

Starting from FASTQ or archive accessions:

```bash
prepare-atac \
  --samples metadata.tsv \
  --genome hg38 \
  --outdir project/raw
```

For processed data, create a tab-separated sample table:

```text
sample	condition	bam	peaks
A1	conditionA	A1.bam	A1_peaks.bed
B1	conditionB	B1.bam	B1_peaks.bed
```

Then run the standard commands:

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

Repeated condition labels define biological replicates. The comparison output
includes motif statistics, aggregate profiles, and a portable interactive HTML
report. See the [complete command reference](api.md) for all options.

## Single-cell workflow

`pseudobulk-fragments` groups fragments by cell annotation.
`pseudobulk-footprints` runs the complete grouped workflow, and
`find-signature-fp` creates per-cell signature heatmaps and UMAP reports.

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

```bash
motif-discovery \
  --candidates candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo
```
