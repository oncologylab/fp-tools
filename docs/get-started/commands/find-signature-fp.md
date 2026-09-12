# [`find-signature-fp`](../../api.md#find-signature-fp)

Score selected motif sites in individual cells and plot the results on a UMAP
and in heatmaps. The command pools signal from nearby cells to reduce sparsity,
so the scores describe relative footprint signatures rather than independent
binding calls for every cell.

## Example command

```bash
find-signature-fp --annotations cell_annotations.tsv --fragments pbmc_fragments.tsv.gz --h5ad pbmc_embedding.h5ad \
  --all-motif-diff-dir project/pseudobulk/diff_footprints \
  --all-motif-results project/pseudobulk/diff_footprints/pseudobulk_diff_footprints_results.txt \
  --outdir project/pseudobulk/signature_fp
```

## Primary inputs

- `--annotations` — TSV with `barcode`, `cell_type`, `snap_cell_type`, `umap_1`, and `umap_2` columns.
- `--fragments` — single-cell fragment file with chromosome, start, end, and barcode in its first four columns; the command creates a missing Tabix index by default.
- `--h5ad` — AnnData file with matching cell names, genomic-bin counts, and a boolean `selected` column in `var`. Bin names must use `chromosome:start-end`; the default bin size is 500 bases.
- `--all-motif-diff-dir` — completed `diff-footprints` directory containing motif-site BED files.
- `--all-motif-results` — motif-level results from that directory; supply it together with `--all-motif-diff-dir`.
- `--outdir` — directory for per-cell scores, heatmaps, and UMAP figures.

Use `cell_type` for broad labels and `snap_cell_type` for detailed labels; they
can be the same if you have one annotation level. AnnData cell names must match
the annotation barcodes. Smoothing uses `obsm['X_spectral']` if available,
otherwise the annotation UMAP coordinates. The companion activity scores need
the genomic-bin counts, so an embedding-only AnnData file is not sufficient.

## Main outputs

Under `{outdir}` the default names include:

| Path | Meaning |
| --- | --- |
| `knn_footprint_signature_scores.tsv` | Per-cell KNN-smoothed footprint protection scores for selected TFs. |
| `knn_footprint_orientation_summary.tsv` | Direction/orientation checks used to make marker scores comparable. |
| `chromvar_like_motif_activity_scores.tsv` | Companion accessibility-derived motif activity scores. |
| `knn_footprint_signature_umap.svg` | Per-marker footprint-signature UMAP panels. |
| `per_cell_footprint_signature_heatmap.svg` | Selected-marker per-cell heatmap. |
| `single_cell_footprinting_summary.svg` | Combined heatmap and representative UMAP summary. |
| `all_motif_per_cell_footprint_signature_heatmap.tsv` | Optional all-motif score matrix and metadata when all-motif inputs are supplied. |

Open the summary SVG first, then use the score tables to inspect individual
cells. Additional top-motif plots and all-TF review PDFs are produced when
all-motif inputs are supplied.

## Choose marker TFs

Choose TFs with `--markers TF1,TF2`. The default markers are
`STAT6,FOSB,CEBPA,IRF8,RELA,ZNF683,NR4A1,SMAD3`; each selected TF must have motif
sites in your inputs. For selected-marker reports only, use `--tf-site-dir`
and omit `--all-motif-diff-dir` and `--all-motif-results` from the example.
This alternative directory must contain files named `{TF}.motif_hits.bed` or
`{TF}.motif_peaks.bed`.

See the [Single-cell output example](../output-examples/single-cell-atac-seq.md)
and the [complete `find-signature-fp` reference](../../api.md#find-signature-fp).
