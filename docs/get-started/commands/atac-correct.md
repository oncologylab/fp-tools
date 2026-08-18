# `atac-correct`

Estimate and correct Tn5 sequence bias in ATAC-seq cut-site signal before
footprint scoring.

## Example command

```bash
atac-correct \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project
```

## Primary inputs

- One or more ATAC-seq BAM files or a project sample table.
- Reference genome FASTA.
- Shared or per-sample peak BED files and an optional blacklist.

## Main outputs

- Bias-corrected cut-site bigWig tracks.
- Optional uncorrected, expected, and bias tracks.
- Bias-correction QC figure, merged peaks, and run logs.

Continue with [`call-footprints`](call-footprints.md), or see the
[complete `atac-correct` reference](../../api.md#atac-correct).
