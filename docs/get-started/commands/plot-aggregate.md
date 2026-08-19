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

- `--sample-table` — samples, conditions, and bias-corrected cut-site signal bigWigs used for the aggregate profiles.
- `--motifs` — motif names or identifiers to plot.
- `--site-set` — motif-site set; the example uses bound sites.
- `--outdir` — project directory containing motif results and receiving plots.

## Main outputs

- `{project}/reports/plot_aggregate.html` — default project-layout interactive aggregate report with motif-centered signal profiles.
- the exact `--output` path — static PDF/PNG/SVG or interactive HTML in custom layout.
- the exact `--output-txt` path — optional per-position aggregate values.
- the exact `--output-aggregated-signals`, `--output-aggregated-scores`, and `--output-aggregated-stats` paths — optional source tables when requested.
- the exact `--output` path in `--motif-grid` mode — multipage motif-by-comparison PDF built from a review bundle.

When both signal types are available, use footprint score bigWigs for motif
statistics and bias-corrected cut-site signal bigWigs for observed aggregate
profiles; label the chosen signal explicitly in figure captions.

```bash
plot-aggregate \
  --input-html project/reports/review_multi_comparisons/index.html \
  --motif-grid \
  --output project/reports/motif_aggregate_grid.pdf
```

See the [Bulk output example](../output-examples/bulk-atac-seq.md) and the
[complete `plot-aggregate` reference](../../api.md#plot-aggregate).
