# [`bulk-footprinting`](../../api.md#bulk-footprinting)

Run bulk ATAC-seq from raw reads or aligned inputs through interactive reports.

The [bulk workflow guide](../workflows/bulk-atac-seq.md) provides a runnable
HepG2-versus-K562 ENCODE example from FASTQ or aligned inputs.

## Example command

```bash
bulk-footprinting --reads-table reads.tsv --comparison-table comparisons.tsv --genome hg38 \
  --outdir project --cores 8
```

## Primary inputs

- `--reads-table` — sample, condition, and local FASTQ or public run-accession columns.
- `--comparison-table` — comparison, condition 1, and condition 2 columns.
- `--genome` — `hg38`, `mm10`, or a custom genome label.
- `--outdir` — project output directory.
- `--cores` — total worker cores.

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
| `{project}/logs/bulk_footprinting/bulk_footprinting_commands.sh` | Exact commands generated for the workflow stages. |
| `{project}/logs/bulk_footprinting/{stage}.stdout.log` and `{stage}.stderr.log` | Stage-specific logs for troubleshooting. |

Use `--sample-table` instead of `--reads-table` to start from prepared BAM and
peak files. `--runtime auto` downloads the pinned external tools only when raw
reads require them. `--runtime system` and `--runtime container` are optional
advanced backends.
