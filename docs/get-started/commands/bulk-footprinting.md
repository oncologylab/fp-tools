# [`bulk-footprinting`](../../api.md#bulk-footprinting)

Run the complete bulk ATAC-seq workflow from aligned BAM files, peak BED files,
and an explicit comparison table.

The [bulk workflow guide](../workflows/bulk-atac-seq.md) provides a runnable
six-sample ENCODE design and the complete seven-cell-line tables.

## Example command

```bash
bulk-footprinting --sample-table samples.tsv --comparison-table comparisons.tsv --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed --plot-aggregate all --review-format auto --outdir project --cores 8
```

## Primary inputs

- `--sample-table` — sample, condition, BAM, and peak BED columns.
- `--comparison-table` — comparison, condition 1, and condition 2 columns.
- `--genome` — reference FASTA.
- `--blacklist` — optional blacklist BED.
- `--outdir` — project output directory.
- `--cores` — total worker cores.
- `--plot-aggregate` — `sig`, `all` (default), `top`, or `off`.
- `--review-format` — `auto` (default), `bundle`, `standalone`, or `none`.

With `--review-format auto`, aggregate-free runs produce one standalone HTML;
other runs retain the static browser bundle.

```bash
bulk-footprinting \
  --sample-table samples.tsv \
  --comparison-table comparisons.tsv \
  --genome hg38.fa.gz \
  --plot-aggregate off \
  --review-format auto \
  --outdir project \
  --cores 8
```

## Main outputs

`{project}` is the `--outdir`, `{sample}` comes from the sample table, and
`{comparison}` comes from the comparison table:

| Path | Meaning |
| --- | --- |
| `{project}/samples/{sample}/atac_correct/{sample}_corrected.bw` | Bias-corrected cut-site signal. |
| `{project}/samples/{sample}/footprints/{sample}_footprints.bw` | Footprint score signal. |
| `{project}/samples/{sample}/match_motifs/motif_matches_results.txt` | Per-sample motif summary and binding calls. |
| `{project}/comparisons/{comparison}/diff_footprints_results.txt` | Motif-level differential statistics. |
| `{project}/comparisons/{comparison}/diff_footprints_{cond1}_{cond2}.html` | Portable interactive comparison report. |
| `{project}/reports/review_multi_comparisons/index.html` | Static browser combining every requested comparison. |
| `{project}/reports/review_multi_comparisons.html` | Aggregate-free portable review written when standalone HTML review mode is selected. |
| `{project}/logs/bulk_footprinting/bulk_footprinting_commands.sh` | Exact commands generated for all five stages. |
| `{project}/logs/bulk_footprinting/{stage}.stdout.log` and `{stage}.stderr.log` | Stage-specific logs for troubleshooting. |

The wrapper does not run `normalize-bigwig`; `--normalization` controls only
the differential stage. See the workflow guide for how this differs from the
pair-specific ENCODE demo method.

FASTQ preparation is a separate optional step with [`prepare-atac`](prepare-atac.md).
