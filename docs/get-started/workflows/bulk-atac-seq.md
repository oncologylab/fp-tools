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

The comparison produces motif statistics, aggregate profiles, and a portable
interactive HTML report. See the [Bulk output example](../output-examples.md#bulk-atac-seq).
