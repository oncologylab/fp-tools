---
core_nav:
  previous:
    title: diff-footprints
    url: get-started/commands/diff-footprints/
---

# [`normalize-bigwig`](../../api.md#normalize-bigwig)

Put corrected cut-site tracks on a comparable signal scale using the same
background regions for every sample. This is an optional step after
`atac-correct`; use the resulting tracks for downstream scoring or plotting
when you need this normalization.

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

- `--sample-table` — TSV with `sample` and `condition` columns; reuse the table from `atac-correct`. Tracks are read from `{project}/samples/{sample}/atac_correct/`.
- `--background` — shared BED intervals used to calculate comparable background statistics.
- `--outdir` — project directory represented by `{project}` below.
- `--method` — transformation; `background-scale` multiplies each signal by a shared-target scale factor.
- `--stat` — within-sample background statistic; the example uses the 95th percentile.
- `--target` — across-sample target for the selected statistic; the example uses the median.

## Main outputs

| Path | Meaning |
| --- | --- |
| `{project}/samples/{sample}/normalize/{sample}_corrected_q95_scaled.bw` | Q95-scaled bias-corrected cut-site signal bigWig for one sample. |
| `{project}/normalize_bigwig_qc.tsv` | Background statistics, selected statistic, target, and scale factor for every sample. |
| `{project}/normalize_bigwig_manifest.tsv` | Sample-to-input/output signal mapping for downstream use. |

Check the QC table to see how much each track was scaled. The manifest lists
the output paths to pass to the next command. Scaling changes signal
magnitude; it does not by itself provide evidence of TF binding.

In custom layout, default outputs use
`{outdir}/{input_stem}.background_scale_{stat}.bw`, plus the two QC tables in
`{outdir}`. `background-zscore` instead writes standardized signal and uses a
method-specific filename suffix.

[Download a representative ENCODE normalization QC table](../../demos/qc/encode/A549_normalize_bigwig_qc.tsv).
See the [Bulk ATAC-seq workflow](../workflows/bulk-atac-seq.md) and the
[complete `normalize-bigwig` reference](../../api.md#normalize-bigwig).
