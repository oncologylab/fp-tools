# [`sc-footprinting`](../../api.md#sc-footprinting)

Analyze single-cell ATAC-seq by combining cells into pseudobulk groups, calling
footprints and motif differences between groups, then mapping selected
footprint signatures back to individual cells.

## Example command

```bash
sc-footprinting --fragments pbmc_fragments.tsv.gz --annotations cell_annotations.tsv --h5ad cell_embedding.h5ad \
  --group-by cell_type --genome-sizes hg38.chrom.sizes --genome hg38.fa.gz --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates --outdir project/pseudobulk
```

## Primary inputs

- `--fragments` — TSV or TSV.GZ with chromosome, start, end, and cell barcode in its first four columns.
- `--annotations` — cell annotation TSV with `barcode`, `cell_type`, `snap_cell_type`, `umap_1`, and `umap_2`, plus any additional grouping columns.
- `--h5ad` — AnnData file containing the same cells and genomic-bin counts used for the companion motif-activity scores. See the requirements below.
- `--group-by` — annotation column used to define pseudobulk groups.
- `--genome-sizes` — two-column chromosome-name and length file used to write grouped signal tracks.
- `--genome` — reference genome FASTA matching the fragments and peak coordinates.
- `--peaks` — accessible-region BED file.
- `--motif-db` — built-in motif database name.
- `--outdir` — directory for pseudobulk tracks, motif results, and reports.

The per-cell reports need annotation barcodes that match the AnnData cell names.
Use `cell_type` for the broad cell labels and `snap_cell_type` for detailed
labels; these can be the same if you have only one annotation level. UMAP
coordinates go in `umap_1` and `umap_2`.

The AnnData file needs genomic-bin names such as `chr1:0-500`, a count matrix,
and a boolean `selected` column in its feature annotations (`var`). The
companion activity calculation uses 500-base bins. Nearest-neighbor smoothing
uses `obsm['X_spectral']` if present, otherwise the annotation UMAP coordinates.
An embedding-only AnnData file is not sufficient.

## Main outputs

Paths below are relative to `--outdir`; `{group}` is a retained pseudobulk group:

| Path | Meaning |
| --- | --- |
| `pseudobulk/{group}.fragments.tsv.gz` and `.tbi` | Indexed fragments for each retained cell group. |
| `pseudobulk/{group}.cutsites.cpm.bw` | Group cut-site signal bigWig. |
| `pseudobulk/{group}.pseudo_pairs.sorted.bam` and `.bai` | Pseudo-paired alignment used for bias correction. |
| `atacorrect/{group}/{group}_corrected.bw` | Bias-corrected cut-site signal per group. |
| `footprints/{group}_footprints.bw` | Footprint score signal per group. |
| `diff_footprints/pseudobulk_diff_footprints_results.txt` | Optional motif-level group comparison results. |
| `plots/single_cell_footprinting/` | Per-cell score tables, heatmaps, and UMAP figures from `find-signature-fp`. |
| `pseudobulk_footprint_manifest.tsv` | Group paths and workflow completion state. |
| `pseudobulk_footprint_commands.sh` | Exact generated commands for reproducibility. |
| `logs/{stage}.stdout.log` and `{stage}.stderr.log` | Captured output for each stage. |

Check the manifest for completed groups, then review the motif comparison and
per-cell heatmaps. Per-cell scores combine signal from neighboring cells to
reduce sparsity; interpret them as relative signatures, not independent
binding calls for every cell. Use `--single-cell-signature-markers` to choose
the TFs shown in the marker reports.

See the [Single-cell workflow](../workflows/single-cell.md) and the
[complete `sc-footprinting` reference](../../api.md#sc-footprinting).
