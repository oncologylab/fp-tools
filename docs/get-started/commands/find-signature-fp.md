# `find-signature-fp`

Calculate and plot per-cell footprint signatures from completed pseudobulk or
motif analyses.

## Example command

```bash
find-signature-fp \
  --annotations cell_annotations.tsv \
  --fragments pbmc_fragments.tsv.gz \
  --h5ad pbmc_embedding.h5ad \
  --tf-site-dir marker_motif_sites \
  --all-motif-results project/pseudobulk/pseudobulk_diff_footprints_results.txt \
  --outdir project/pseudobulk/signature_fp
```

## Primary inputs

- Cell annotations, fragments, and a single-cell H5AD containing an embedding.
- Motif-site directories and result tables from prior footprint analysis.
- Optional marker definitions and KNN settings.

## Main outputs

- Per-cell footprint-signature heatmaps.
- Cell-type and footprint-signature UMAP figures.
- Source tables for the plotted signatures.

See the [Single-cell output example](../output-examples.md#single-cell-atac-seq)
and the [complete `find-signature-fp` reference](../../api.md#find-signature-fp).
