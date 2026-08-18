# `plot-motif-aggregate-grid`

Export motif-by-comparison aggregate profiles from a multi-comparison review
report.

## Example command

```bash
plot-motif-aggregate-grid \
  --outdir project \
  --output project/reports/motif_aggregate_grid.pdf \
  --source-tsv project/reports/motif_aggregate_grid_source.tsv \
  --rows-per-page 16 \
  --fill-missing-profiles
```

## Primary inputs

- A `review-multi-comparisons` HTML report or project directory.
- Optional motif ordering, page size, flank width, and missing-profile settings.
- Optional RNA expression tables and motif-to-gene mapping.

## Main outputs

- A multipage PDF with motifs as rows and comparisons as columns.
- A source TSV containing plotted statistics, profile provenance, and site counts.

See the [complete `plot-motif-aggregate-grid` reference](../../api.md#plot-motif-aggregate-grid).
