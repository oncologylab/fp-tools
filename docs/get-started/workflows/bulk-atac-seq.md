# Bulk ATAC-seq workflow

This guide starts with either paired FASTQ files or aligned BAM and peak files,
then compares motif-associated footprint scores between biological conditions.
The worked example uses three HepG2 and three K562 biological replicates from
ENCODE experiments [ENCSR291GJU](https://www.encodeproject.org/experiments/ENCSR291GJU/)
and [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/).

## Before you start

Install fp-tools and verify the external programs used by raw-read processing:

```bash
prepare-atac --doctor --profile modern
```

For a custom genome, provide a matching FASTA, Bowtie2 index, chromosome names,
blacklist, and peak files. All BAMs, peaks, and the FASTA must use the same
assembly. The examples below use hg38.

## Starting from FASTQ files

[`prepare-atac`](../commands/prepare-atac.md) accepts a tab-separated or
comma-separated sample sheet.

For paired-end data:

```text
sample	condition	fastq_1	fastq_2
HepG2_rep1	HepG2	fastq/HepG2_rep1_R1.fastq.gz	fastq/HepG2_rep1_R2.fastq.gz
K562_rep1	K562	fastq/K562_rep1_R1.fastq.gz	fastq/K562_rep1_R2.fastq.gz
```

For single-end data, omit `fastq_2`:

```text
sample	condition	fastq_1
HepG2_rep1	HepG2	fastq/HepG2_rep1.fastq.gz
K562_rep1	K562	fastq/K562_rep1.fastq.gz
```

[Download the ENCODE FASTQ sheet](../../demos/data/encode/encode_hepg2_k562_fastq_urls.tsv)
or use the [local FASTQ template](../../demos/data/encode/local_fastq_template.tsv).
Use the columns shown for your library layout. The ENCODE download includes
optional provenance columns that may be left out.

```bash
prepare-atac --samples encode_hepg2_k562_fastq_urls.tsv --genome hg38 --outdir project --cores 8
```

The downstream aligned-data table is written to
`project/metadata/samples.tsv`; use that generated table for the remaining
commands.

## Starting from aligned ENCODE data

The sample table requires `sample`, `condition`, `bam`, and `peaks`:

```text
sample	condition	bam	peaks
HepG2_rep1	HepG2	data/bams/HepG2_rep1.bam	data/peaks/HepG2_peaks.bed
K562_rep1	K562	data/bams/K562_rep1.bam	data/peaks/K562_peaks.bed
```

Use the complete [HepG2/K562 ENCODE sample table](../../demos/data/encode/encode_hepg2_k562_bams.tsv)
with its [download helper](../../demos/data/encode/download_encode_hepg2_k562.sh),
or copy the
[local BAM/peak template](../../demos/data/encode/local_bam_peak_template.tsv)
for your own data.

Repeated `condition` values define biological replicates. Every BAM needs a
`.bai` index next to it.

## Define the comparison

The comparison direction is `cond1 - cond2`. For HepG2 versus K562:

```text
comparison	cond1	cond2
HepG2_vs_K562	HepG2	K562
```

[Download the HepG2/K562 comparison table](../../demos/data/encode/encode_hepg2_k562_comparisons.tsv).
Condition names must exactly match the sample table; the `comparison` value
becomes the output-directory name.

## Run the complete aligned-data workflow

```bash
bulk-footprinting --sample-table encode_hepg2_k562_bams.tsv --comparison-table encode_hepg2_k562_comparisons.tsv \
  --genome hg38.fa.gz --blacklist hg38.blacklist.bed --motif-db jaspar2026_vertebrates --outdir project --cores 8
```

Preview the expanded commands without running them:

```bash
bulk-footprinting --sample-table encode_hepg2_k562_bams.tsv --comparison-table encode_hepg2_k562_comparisons.tsv \
  --genome hg38.fa.gz --blacklist hg38.blacklist.bed --outdir project --dry-run
```

The wrapper runs `atac-correct`, `call-footprints`, `match-motifs`,
`diff-footprints`, and `review-multi-comparisons`. See the
[`bulk-footprinting` guide](../commands/bulk-footprinting.md) for its exact
output tree.

## Full seven-cell-line design

The full example contains 17 samples from seven cancer cell lines:

- [17-sample ENCODE BAM/peak table](../../demos/data/encode/encode_cancer_7line_bams.tsv)
- [21-comparison table](../../demos/data/encode/encode_cancer_7line_comparisons.tsv)
- [Interactive ENCODE cancer-cell-line output](../../reports.md)

## Run the core commands separately

The same analysis can be run step by step:

<div class="fp-command-chain" markdown="1">

[`atac-correct`](../commands/atac-correct.md)
<span>→</span>
[`call-footprints`](../commands/call-footprints.md)
<span>→</span>
[`match-motifs`](../commands/match-motifs.md)
<span>→</span>
[`diff-footprints`](../commands/diff-footprints.md)

</div>

Use the exact project-layout paths documented on each command page. The
[Bulk output example](../output-examples/bulk-atac-seq.md) shows how to inspect
the final interactive report.

## Example QC files

- [ATACorrect diagnostic PDF](../../demos/qc/encode/A549_rep1_atacorrect.pdf)
- [Normalization QC table](../../demos/qc/encode/A549_normalize_bigwig_qc.tsv)
- [Replicate motif-score matrix](../../demos/qc/encode/A549_motif_score_matrix.tsv)
