# [`find-signature-fp`](../../api.md#find-signature-fp)

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

Under `{outdir}` the default names include:

| Path | Meaning |
| --- | --- |
| `knn_footprint_signature_scores.tsv` | Per-cell KNN-smoothed footprint protection scores for selected TFs. |
| `knn_footprint_orientation_summary.tsv` | Direction/orientation checks used to make marker scores comparable. |
| `chromvar_like_motif_activity_scores.tsv` | Companion accessibility-derived motif activity scores. |
| `knn_footprint_signature_umap.svg` and `.pdf` | Per-marker footprint-signature UMAP panels. |
| `per_cell_footprint_signature_heatmap.svg` and `.pdf` | Selected-marker per-cell heatmap. |
| `single_cell_footprinting_summary.svg` and `.pdf` | Combined heatmap and representative UMAP summary. |
| `all_motif_per_cell_footprint_signature_heatmap.tsv` | Optional all-motif score matrix and metadata when all-motif inputs are supplied. |

Additional top-motif and all-TF review files use their requested output prefix.

See the [Single-cell output example](../output-examples/single-cell-atac-seq.md)
and the [complete `find-signature-fp` reference](../../api.md#find-signature-fp).
