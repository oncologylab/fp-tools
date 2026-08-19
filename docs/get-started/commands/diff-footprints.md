# [`diff-footprints`](../../api.md#diff-footprints)

Compare motif-associated footprint scores across conditions or between
user-defined region sets measured in the same sample(s).

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
- Region-set analyses also report confidence intervals, motif prevalence,
  region counts, per-replicate effects, and matching balance.

For a region-set comparison, use `--comparison-axis regions`, provide two or
more BED files with `--regions`, and name them with `--region-labels`. An
optional `--region-strata-column` preserves accessibility or other matching
strata during resampling. One sample uses a stratified label-permutation test;
two or more biological replicates use a paired empirical-Bayes model.

See the [Bulk output example](../output-examples/bulk-atac-seq.md), the
[region-set comparison example](../output-examples/region-set-comparison.md), and the
[complete `diff-footprints` reference](../../api.md#diff-footprints).
