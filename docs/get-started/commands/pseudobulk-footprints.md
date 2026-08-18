# `pseudobulk-footprints`

Run grouping, bias correction, footprint scoring, motif analysis, and optional
per-cell signature reporting for single-cell ATAC-seq data.

## Example command

```bash
pseudobulk-footprints \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project/pseudobulk
```

## Primary inputs

- Single-cell fragments or a barcode-tagged BAM plus cell annotations.
- Grouping column, genome, genome sizes, and accessible peaks.
- Optional motif database, signature-site directory, and H5AD embedding.

## Main outputs

- Pseudobulk fragments, pseudo-BAMs, and corrected/footprint bigWigs.
- Motif-aware differential reports and aggregate plots.
- Optional per-cell footprint-signature heatmaps and UMAP figures.

See the [Single-cell workflow](../workflows/single-cell.md) and the
[complete `pseudobulk-footprints` reference](../../api.md#pseudobulk-footprints).
