# Bulk ATAC-seq workflow

`bulk-footprinting` runs a complete bulk ATAC-seq analysis from FASTQ files or
prepared BAM and peak files. It produces footprint scores, differential motif
results, and an interactive report.

This example compares three HepG2 and three K562 biological replicates from
ENCODE experiments [ENCSR291GJU](https://www.encodeproject.org/experiments/ENCSR291GJU/)
and [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/).

## Start from FASTQ files

Download the [ENCODE FASTQ sample sheet](../../demos/data/encode/encode_hepg2_k562_fastq_urls.tsv)
and the [HepG2/K562 comparison file](../../demos/data/encode/encode_hepg2_k562_comparisons.tsv).
For your own data, begin with the [local FASTQ template](../../demos/data/encode/local_fastq_template.tsv).

```bash
bulk-footprinting --reads-table encode_hepg2_k562_fastq_urls.tsv \
  --comparison-table encode_hepg2_k562_comparisons.tsv --genome hg38 \
  --outdir project --cores 8
```

## Start from aligned files

Download the [ENCODE BAM and peak sample sheet](../../demos/data/encode/encode_hepg2_k562_bams.tsv)
and its [download helper](../../demos/data/encode/download_encode_hepg2_k562.sh).
For your own data, use the [local BAM and peak template](../../demos/data/encode/local_bam_peak_template.tsv).
The same [comparison file](../../demos/data/encode/encode_hepg2_k562_comparisons.tsv)
works for these aligned inputs.

```bash
bulk-footprinting --sample-table encode_hepg2_k562_bams.tsv \
  --comparison-table encode_hepg2_k562_comparisons.tsv \
  --genome hg38.fa.gz --blacklist hg38.blacklist.bed \
  --motif-db jaspar2026_vertebrates --outdir project --cores 8
```

Repeated condition names in a sample sheet define biological replicates. BAM,
peak, blacklist, and genome files must use the same genome assembly.

## Review the results

Open the generated interactive report to review differential motifs and
aggregate footprint profiles. See the [bulk output example](../output-examples/bulk-atac-seq.md)
or explore the [ENCODE cancer-cell-line output demo](../../reports.md).

For additional options and output paths, see the
[`bulk-footprinting` guide](../commands/bulk-footprinting.md) or the
[complete API reference](../../api.md#bulk-footprinting).
