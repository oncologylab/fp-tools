# [`normalize-bigwig`](../../api.md#normalize-bigwig)

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

- `--sample-table` — samples and input bigWig tracks.
- `--background` — shared BED regions used to estimate scaling factors.
- `--outdir` — project directory for normalized tracks and statistics.
- `--method` — normalization method; the example uses `background-scale`.
- `--stat` — background summary statistic; the example uses `q95`.
- `--target` — across-sample target; the example uses the median.

## Main outputs

- One normalized bigWig per sample.
- Background statistics and scaling-factor table.
- Manifest of normalized tracks for downstream commands.

See the [Bulk ATAC-seq workflow](../workflows/bulk-atac-seq.md) and the
[complete `normalize-bigwig` reference](../../api.md#normalize-bigwig).
