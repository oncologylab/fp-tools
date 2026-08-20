# Bulk ATAC-seq workflow

`bulk-footprinting` runs a complete bulk ATAC-seq analysis from prepared BAM/BAI
and peak BED files. It produces footprint scores, differential motif results,
and an interactive report.

This example compares three HepG2 and three K562 biological replicates from
ENCODE experiments [ENCSR291GJU](https://www.encodeproject.org/experiments/ENCSR291GJU/)
and [ENCSR868FGK](https://www.encodeproject.org/experiments/ENCSR868FGK/).

## Run from aligned files

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

??? note "Optional FASTQ-to-BAM preparation"

    On Linux, [`prepare-atac`](../commands/prepare-atac.md) can prepare FASTQ
    files as BAM/BAI and peak BED inputs. Run it separately before
    `bulk-footprinting`:

    ```bash
    prepare-atac --samples reads.tsv --genome hg38 --outdir prepared_project
    ```

    Use the generated `metadata/samples.tsv` with `bulk-footprinting`. This
    preprocessing step is not part of the bulk footprinting wrapper.

## Review the results

Open the generated interactive report to review differential motifs and
aggregate footprint profiles. See the [bulk output example](../output-examples/bulk-atac-seq.md)
or explore the [ENCODE cancer-cell-line output demo](../../reports.md).

For additional options and output paths, see the
[`bulk-footprinting` guide](../commands/bulk-footprinting.md) or the
[complete API reference](../../api.md#bulk-footprinting).
