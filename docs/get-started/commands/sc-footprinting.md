# [`sc-footprinting`](../../api.md#sc-footprinting)

Run grouping, bias correction, footprint scoring, motif analysis, and per-cell
signature reporting for single-cell ATAC-seq data.

## Example command

```bash
sc-footprinting --fragments pbmc_fragments.tsv.gz --annotations cell_annotations.tsv --h5ad cell_embedding.h5ad \
  --group-by cell_type --genome-sizes hg38.chrom.sizes --genome hg38.fa.gz --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates --outdir project/pseudobulk
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

`{outdir}` contains a complete staged workflow:

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

See the [Single-cell workflow](../workflows/single-cell.md) and the
[complete `sc-footprinting` reference](../../api.md#sc-footprinting).
