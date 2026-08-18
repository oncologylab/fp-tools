# `normalize-bigwig`

Normalize multiple bigWig tracks with scale estimates calculated from the same
background regions.

## Example command

```bash
normalize-bigwig \
  --sample-table project/metadata/samples.tsv \
  --background project/peaks/merged_peaks_filtered.bed \
  --outdir project \
  --method background-scale \
  --stat q95 \
  --target median
```

## Primary inputs

- A project sample table or explicit bigWig paths.
- Shared background BED regions.
- Normalization method, summary statistic, and across-sample target.

## Main outputs

- One normalized bigWig per sample.
- Background statistics and scaling-factor table.
- Manifest of normalized tracks for downstream commands.

See the [Bulk ATAC-seq workflow](../workflows/bulk-atac-seq.md) and the
[complete `normalize-bigwig` reference](../../api.md#normalize-bigwig).
