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

- `--sample-table` — sample names, BAM files, and peak BED files for one or more libraries.
- `--genome` — reference genome FASTA.
- `--blacklist` — genomic regions excluded from signal correction.
- `--outdir` — project directory for corrected tracks and QC outputs.

## Main outputs

- Bias-corrected cut-site bigWig tracks.
- Optional uncorrected, expected, and bias tracks.
- Bias-correction QC figure, merged peaks, and run logs.

Continue with [`call-footprints`](call-footprints.md), or see the
[complete `atac-correct` reference](../../api.md#atac-correct).
