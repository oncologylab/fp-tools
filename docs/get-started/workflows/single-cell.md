# Single-cell workflow

The single-cell workflow groups fragments into pseudobulk samples and can then
calculate per-cell footprint signatures.

## Main commands

<div class="fp-command-chain" markdown="1">

[`pseudobulk-fragments`](../commands/pseudobulk-fragments.md)
<span>→</span>
[`pseudobulk-footprints`](../commands/pseudobulk-footprints.md)
<span>→</span>
[`find-signature-fp`](../commands/find-signature-fp.md)

</div>

- [`pseudobulk-fragments`](../commands/pseudobulk-fragments.md) performs grouping only.
- [`pseudobulk-footprints`](../commands/pseudobulk-footprints.md) runs grouping, correction, footprint scoring, motif analysis, and reports.
- [`find-signature-fp`](../commands/find-signature-fp.md) creates per-cell signature heatmaps and UMAP figures.

## Example

```bash
pseudobulk-footprints \
  --fragments fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project/pseudobulk
```

See the [Single-cell output example](../output-examples.md#single-cell-atac-seq).
