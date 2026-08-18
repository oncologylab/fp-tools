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

- `--annotations` — barcode-level cell annotation table.
- `--fragments` — indexed single-cell fragment file.
- `--h5ad` — single-cell object containing the spectral or UMAP embedding.
- `--tf-site-dir` — motif-site directories from the footprint analysis.
- `--all-motif-results` — completed motif-level differential result table.
- `--outdir` — directory for per-cell scores, heatmaps, and UMAP figures.

## Main outputs

- Per-cell footprint-signature heatmaps.
- Cell-type and footprint-signature UMAP figures.
- Source tables for the plotted signatures.

See the [Single-cell output example](../output-examples.md#single-cell-atac-seq)
and the [complete `find-signature-fp` reference](../../api.md#find-signature-fp).
