# `pseudobulk-fragments`

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

- Single-cell fragment TSV or TSV.GZ file.
- Cell annotation table and grouping column.
- Genome sizes when cut-site bigWigs or pseudo-BAMs are requested.

## Main outputs

- One fragment file per retained cell group.
- Group manifest and QC summary.
- Optional indexed fragments, cut-site bigWigs, and pseudo-BAMs.

Continue with [`pseudobulk-footprints`](pseudobulk-footprints.md), or see the
[complete `pseudobulk-fragments` reference](../../api.md#pseudobulk-fragments).
