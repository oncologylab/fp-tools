# [`bulk-footprinting`](../../api.md#bulk-footprinting)

Run a complete bulk ATAC-seq analysis from aligned reads and peak regions to
footprint scores, motif comparisons, and interactive reports.

The [bulk workflow guide](../workflows/bulk-atac-seq.md) provides minimal sample
and comparison tables for a two-condition analysis.

## Example command

```bash
bulk-footprinting --sample-table samples.tsv --comparison-table comparisons.tsv --genome hg38 \
  --outdir project --cores 8
```

## Primary inputs

- `--sample-table` — TSV with `sample`, `condition`, `bam`, and `peaks` columns. Each BAM must be coordinate-sorted and have a matching BAI index.
- `--comparison-table` — TSV with `comparison`, `cond1`, and `cond2` columns. Use condition names from the sample table.
- `--genome` — managed `hg38` or `mm10` assembly, or a reference FASTA matching the BAM and peak coordinates.
- `--outdir` — project output directory.
- `--cores` — total worker cores.

Use the same genome assembly and chromosome names for every BAM and BED file.

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

Start by opening the comparison HTML report. Use the combined review to compare
motif results across all requested comparisons, and inspect aggregate profiles
alongside the statistics.

Add `--dry-run` to check the inputs and inspect the commands before starting.

## Reference and motif options

Choosing `hg38` or `mm10` downloads and verifies the matching reference and
blacklist as needed. Use `--reference-dir` to choose where they are cached.
`--blacklist` replaces a managed assembly's blacklist, while
`--no-blacklist` disables it. Custom FASTA inputs never infer a blacklist.

Choose a packaged motif database with `--motif-db`, provide custom files with
`--motifs`, or combine both options. With neither option, the workflow uses
`jaspar2026_vertebrates`. Run `bulk-footprinting --list-motif-dbs` to list the
packaged databases.

If you have FASTQ files on Linux, run
[`prepare-atac`](prepare-atac.md) first, then provide its generated
`metadata/samples.tsv` to this command.
