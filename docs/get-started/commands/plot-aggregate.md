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

- Signal bigWigs.
- Motif-site BED files, other region BED files, or a `match-motifs` directory.
- Optional labels, site-set selection, flank size, and normalization mode.

## Main outputs

- Static PDF or interactive HTML aggregate plots.
- Optional aggregated signal and summary tables.

See the [Bulk output example](../output-examples.md#bulk-atac-seq) and the
[complete `plot-aggregate` reference](../../api.md#plot-aggregate).
