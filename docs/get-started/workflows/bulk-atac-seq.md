# Bulk ATAC-seq workflow

Use this workflow for individual samples, biological replicates, or
condition-level comparisons.

## Starting from FASTQ files

[`prepare-atac`](../commands/prepare-atac.md) can download public runs or use
local FASTQ files and produces the BAM, peak, coverage, and QC files required
for footprint analysis.

```bash
prepare-atac \
  --samples metadata.tsv \
  --genome hg38 \
  --outdir project/raw
```

## Starting from aligned data

Create a tab-separated sample table. Repeated condition names identify
biological replicates.

```text
sample	condition	bam	peaks
A1	conditionA	A1.bam	A1_peaks.bed
B1	conditionB	B1.bam	B1_peaks.bed
```

## Main analysis

Run the complete workflow with [`bulk-footprinting`](../commands/bulk-footprinting.md):

```bash
bulk-footprinting \
  --sample-table samples.tsv \
  --comparison-table comparisons.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project
```

## Run step by step

<div class="fp-command-chain" markdown="1">

[`atac-correct`](../commands/atac-correct.md)
<span>→</span>
[`call-footprints`](../commands/call-footprints.md)
<span>→</span>
[`match-motifs`](../commands/match-motifs.md)
<span>→</span>
[`diff-footprints`](../commands/diff-footprints.md)

</div>

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

The workflow produces motif statistics, aggregate profiles, and a static
comparison browser. See the [Bulk output example](../output-examples/bulk-atac-seq.md).
