# [`pseudobulk-fragments`](../../api.md#pseudobulk-fragments)

Group single-cell ATAC fragments by a cell-annotation column to create
pseudobulk inputs.

## Example command

```bash
pseudobulk-fragments \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --write-cutsite-bigwigs \
  --outdir project/pseudobulk/fragments
```

## Primary inputs

- `--fragments` — single-cell fragment TSV or TSV.GZ file.
- `--annotations` — barcode-level cell annotation table.
- `--group-by` — annotation column used to define pseudobulk groups.
- `--genome-sizes` — chromosome sizes used to write signal tracks.
- `--write-cutsite-bigwigs` — write cut-site bigWigs for retained groups.
- `--outdir` — directory for grouped fragments, tracks, and QC outputs.

## Main outputs

- One fragment file per retained cell group.
- Group manifest and QC summary.
- Optional indexed fragments, cut-site bigWigs, and pseudo-BAMs.

Continue with [`sc-footprinting`](sc-footprinting.md), or see the
[complete `pseudobulk-fragments` reference](../../api.md#pseudobulk-fragments).
