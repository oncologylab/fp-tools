# [`pseudobulk-fragments`](../../api.md#pseudobulk-fragments)

Combine single-cell ATAC-seq fragments into groups such as cell types or
donor–cell-type pairs. Each group becomes a pseudobulk sample for downstream
analysis.

## Example command

```bash
pseudobulk-fragments --fragments pbmc_fragments.tsv.gz --annotations cell_annotations.tsv --group-by cell_type \
  --genome-sizes hg38.chrom.sizes --write-cutsite-bigwigs --outdir project/pseudobulk/fragments
```

## Primary inputs

- `--fragments` — TSV or TSV.GZ with chromosome, start, end, and barcode in its first four columns; a fifth column can contain fragment counts.
- `--annotations` — TSV or CSV with `barcode` and the column named by `--group-by` (for example, `cell_type`).
- `--group-by` — annotation column used to define pseudobulk groups; use `donor,cell_type` to keep donors separate within each cell type.
- `--genome-sizes` — two-column chromosome-name and length file matching the fragments; required for the example's signal tracks.
- `--write-cutsite-bigwigs` — write cut-site bigWigs for retained groups.
- `--outdir` — directory for grouped fragments, tracks, and QC outputs.

## Main outputs

`{group}` is the annotation value converted to a filename-safe name. Paths below
are relative to `{outdir}`:

| Path | Meaning |
| --- | --- |
| `{group}.fragments.tsv` or `{group}.fragments.tsv.gz` | Fragments assigned to the group; compressed/indexed form is controlled by the command options. |
| `{group}.fragments.tsv.gz.tbi` | Optional Tabix index for random genomic access. |
| `{group}.cutsites.cpm.bw` | Optional CPM-normalized cut-site signal bigWig written by `--write-cutsite-bigwigs`. |
| `{group}.pseudo_pairs.sorted.bam` and `.bai` | Optional pseudo-paired alignment written with `--write-pseudo-bams`; use `atac-correct --read_shift 0 0` on these files. |
| `pseudobulk_manifest.tsv` | Per-group paths, cell/fragment counts, and filter status. |
| `fp_tools_manifest.yml` | Machine-readable run settings and retained groups. |
| `pseudobulk_downstream_commands.sh` | Optional generated downstream command examples. |

Inspect `pseudobulk_manifest.tsv` for each group's cell and fragment counts.
Choose `--min-cells` and `--min-fragments` to exclude groups with too little data.
To run grouping and the full footprint analysis together, start with
[`sc-footprinting`](sc-footprinting.md). See the
[complete `pseudobulk-fragments` reference](../../api.md#pseudobulk-fragments).

## Match cell barcodes

By default, barcode matching ignores a trailing suffix such as `-1`. Add
`--no-strip-barcode-suffix` when suffixes distinguish cells in your dataset.
Use `--barcode-column` if your annotation barcode column has a different name.
