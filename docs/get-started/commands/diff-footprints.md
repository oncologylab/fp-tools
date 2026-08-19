---
core_nav:
  previous:
    title: match-motifs
    url: get-started/commands/match-motifs/
  next:
    title: normalize-bigwig
    url: get-started/commands/normalize-bigwig/
---

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

Each comparison is written below
`{project}/comparisons/{comparison}/`, where `{comparison}` is taken from the
comparison table and `{prefix}` defaults to `diff_footprints`:

| Path | Meaning |
| --- | --- |
| `{prefix}_results.txt` | Tab-separated motif-level differential footprint statistics; change direction is `cond1 - cond2`. |
| `{prefix}_results.xlsx` | Excel copy of the result table unless `--skip-excel` is used. |
| `{prefix}_distances.txt` | Motif distances used for clustering related motifs. |
| `{prefix}_{cond1}_{cond2}.html` | Portable interactive report with volcano, motif, and embedded aggregate-profile views. |
| `{prefix}_replicate_report.tsv` | Long-form per-replicate diagnostic data when replicate reporting is active. |
| `{prefix}_replicate_summary.tsv` | Motif-level replicate agreement summary. |
| `{prefix}_replicate_report.png` | Replicate diagnostic figure. |
| `{prefix}_figures.pdf` and `{prefix}_clusters.pdf` | Optional static summaries written with `--static-plots`. |
| `{motif}/beds/{motif}_{condition}_bound.bed` | Motif instances classified as bound for a condition when full motif outputs are required. |

Region-set analyses use the same result/report patterns and add confidence
intervals, motif prevalence, region counts, per-replicate effects, and matching
balance to the result tables.

For a region-set comparison, use `--comparison-axis regions`, provide two or
more BED files with `--regions`, and name them with `--region-labels`. An
optional `--region-strata-column` preserves accessibility or other matching
strata during resampling. One sample uses a stratified label-permutation test;
two or more biological replicates use a paired empirical-Bayes model.
Use `--plot-aggregate-motifs` to choose an ordered aggregate panel without
limiting the motifs tested, and `--default-aggregate-plots` to set its initial
size.

See the [Bulk output example](../output-examples/bulk-atac-seq.md), the
[region-set comparison example](../output-examples/region-set-comparison.md), and the
[complete `diff-footprints` reference](../../api.md#diff-footprints).
