# Single-cell ATAC-seq workflow

`sc-footprinting` groups fragments into pseudobulk samples, runs the footprint
analysis, and calculates per-cell footprint signatures.

## Main commands

<div class="fp-command-chain" markdown="1">

[`sc-footprinting`](../commands/sc-footprinting.md)

</div>

- [`sc-footprinting`](../commands/sc-footprinting.md) runs the complete workflow.
- [`pseudobulk-fragments`](../commands/pseudobulk-fragments.md) and [`find-signature-fp`](../commands/find-signature-fp.md) remain available as focused utilities.

## Example

```bash
sc-footprinting \
  --fragments fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --h5ad cell_embedding.h5ad \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project/pseudobulk
```

See the [Single-cell output example](../output-examples/single-cell-atac-seq.md).
