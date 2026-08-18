# `diff-footprints`

Compare motif-associated footprint scores across conditions, including
replicate-supported contrasts.

## Example command

```bash
diff-footprints \
  --sample-table project/metadata/samples.tsv \
  --comparison-table project/metadata/comparisons.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

## Primary inputs

- `--sample-table` — samples, conditions, footprint tracks, and reusable motif-result folders.
- `--comparison-table` — condition pairs to compare.
- `--genome` — reference genome FASTA.
- `--peaks` — accessible-region BED file.
- `--motif-db` — built-in motif database name.
- `--outdir` — project directory for statistics, figures, and HTML reports.

## Main outputs

- Motif-level differential footprint statistics.
- Replicate score matrices and diagnostics when replicates are present.
- A portable interactive HTML report with volcano and aggregate views.

See the [Bulk output example](../output-examples.md#bulk-atac-seq) and the
[complete `diff-footprints` reference](../../api.md#diff-footprints).
