# [`plot-aggregate`](../../api.md#plot-aggregate)

Plot the average cut-site signal around motif sites to inspect the shape of a
footprint. Start with corrected bigWigs from `atac-correct` and motif results
from `match-motifs`. Outputs can be static figures or an interactive HTML
report.

## Example command

```bash
plot-aggregate \
  --sample-table project/metadata/samples.tsv \
  --motifs SPIB CEBPB \
  --site-set bound \
  --outdir project
```

## Primary inputs

- `--sample-table` — TSV with `sample` and `condition` columns; the command reads corrected tracks and motif results from `{project}/samples/{sample}/`.
- `--motifs` — motif names or identifiers to plot.
- `--site-set` — sites to average; the example uses sites classified as bound in each sample.
- `--outdir` — project directory containing motif results and receiving plots.

Replace `SPIB CEBPB` with motifs present in your motif results.

## Main outputs

- `{project}/reports/plot_aggregate.html` — default project-layout interactive aggregate report with motif-centered signal profiles.
- the exact `--output` path — static PDF/PNG/SVG or interactive HTML in custom layout.
- the exact `--output-txt` path — optional per-position aggregate values.
- the exact `--output-aggregated-signals`, `--output-aggregated-scores`, and `--output-aggregated-stats` paths — optional source tables when requested.
- the exact `--output` path in `--motif-grid` mode — multipage motif-by-comparison PDF built from a review bundle.

When both signal types are available, use footprint score bigWigs for motif
statistics and bias-corrected cut-site signal bigWigs for observed aggregate
profiles; label the chosen signal explicitly in figure captions.

Look for central cut-site depletion relative to the flanking signal. Also
compare the number of sites and the sample-to-sample consistency. Bound-site
plots can use different sites in each sample; choose `--site-set all` to
inspect all matched motif sites instead.

For your own site sets, supply BED files with `--TFBS`. Use `--regions` to
restrict or compare regions of interest.

## Export a grid from a review report

```bash
plot-aggregate \
  --input-html project/reports/review_multi_comparisons/index.html \
  --motif-grid \
  --output project/reports/motif_aggregate_grid.pdf
```

See the [Bulk output example](../output-examples/bulk-atac-seq.md) and the
[complete `plot-aggregate` reference](../../api.md#plot-aggregate).
