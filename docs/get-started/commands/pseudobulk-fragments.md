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

For each sanitized `{group}` under `{outdir}`:

| Path | Meaning |
| --- | --- |
| `{group}.fragments.tsv` or `{group}.fragments.tsv.gz` | Fragments assigned to the group; compressed/indexed form is controlled by the command options. |
| `{group}.fragments.tsv.gz.tbi` | Optional Tabix index for random genomic access. |
| `{group}.cutsites.cpm.bw` | Optional CPM-normalized cut-site signal bigWig written by `--write-cutsite-bigwigs`. |
| `{group}.pseudo_pairs.sorted.bam` and `.bai` | Optional pseudo-paired alignment used by `atac-correct` with read shift `0 0`. |
| `pseudobulk_manifest.tsv` | Per-group paths, cell/fragment counts, and filter status. |
| `fp_tools_manifest.yml` | Machine-readable run settings and retained groups. |
| `pseudobulk_downstream_commands.sh` | Optional generated downstream command examples. |

Continue with [`sc-footprinting`](sc-footprinting.md), or see the
[complete `pseudobulk-fragments` reference](../../api.md#pseudobulk-fragments).
