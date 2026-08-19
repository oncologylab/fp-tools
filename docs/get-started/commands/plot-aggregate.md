# [`plot-aggregate`](../../api.md#plot-aggregate)

Plot average signal around motif sites or other genomic regions as a static
figure or interactive HTML report.

Multiple user-defined BED files are supported through `--TFBS`. Multiple
`--regions` BED files can restrict or compare distinct regions of interest.

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
- A multipage motif-by-comparison PDF when `--motif-grid` is used with a
  `review-multi-comparisons` bundle.

```bash
plot-aggregate \
  --input-html project/reports/review_multi_comparisons/index.html \
  --motif-grid \
  --output project/reports/motif_aggregate_grid.pdf
```

See the [Bulk output example](../output-examples/bulk-atac-seq.md) and the
[complete `plot-aggregate` reference](../../api.md#plot-aggregate).
