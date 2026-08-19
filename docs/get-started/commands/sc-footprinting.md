# [`sc-footprinting`](../../api.md#sc-footprinting)

Run grouping, bias correction, footprint scoring, motif analysis, and per-cell
signature reporting for single-cell ATAC-seq data.

## Example command

```bash
sc-footprinting \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --h5ad cell_embedding.h5ad \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project/pseudobulk
```

## Primary inputs

- `--fragments` — single-cell fragment file.
- `--annotations` — barcode-level cell annotation table.
- `--h5ad` — AnnData file containing the cell embedding used for KNN smoothing.
- `--group-by` — annotation column used to define pseudobulk groups.
- `--genome-sizes` — chromosome sizes used to write grouped signal tracks.
- `--genome` — reference genome FASTA.
- `--peaks` — accessible-region BED file.
- `--motif-db` — built-in motif database name.
- `--outdir` — directory for pseudobulk tracks, motif results, and reports.

## Main outputs

- Pseudobulk fragments, pseudo-BAMs, and corrected/footprint bigWigs.
- Motif-aware differential reports and aggregate plots.
- Per-cell footprint-signature heatmaps and UMAP figures.

See the [Single-cell workflow](../workflows/single-cell.md) and the
[complete `sc-footprinting` reference](../../api.md#sc-footprinting).
