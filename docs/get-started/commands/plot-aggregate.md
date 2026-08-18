# `plot-aggregate`

Plot average signal around motif sites or other genomic regions as a static
figure or interactive HTML report.

## Example command

```bash
plot-aggregate \
  --sample-table project/metadata/samples.tsv \
  --motifs SPIB CEBPB \
  --site-set bound \
  --outdir project
```

## Primary inputs

- `--sample-table` — samples, conditions, and signal bigWig tracks.
- `--motifs` — motif names or identifiers to plot.
- `--site-set` — motif-site set; the example uses bound sites.
- `--outdir` — project directory containing motif results and receiving plots.

## Main outputs

- Static PDF or interactive HTML aggregate plots.
- Optional aggregated signal and summary tables.

See the [Bulk output example](../output-examples.md#bulk-atac-seq) and the
[complete `plot-aggregate` reference](../../api.md#plot-aggregate).
