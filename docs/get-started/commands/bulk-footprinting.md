# [`bulk-footprinting`](../../api.md#bulk-footprinting)

Run the complete bulk ATAC-seq workflow from aligned BAM files, peak BED files,
and an explicit comparison table.

## Example command

```bash
bulk-footprinting \
  --sample-table samples.tsv \
  --comparison-table comparisons.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project \
  --cores 8
```

## Primary inputs

- `--sample-table` — sample, condition, BAM, and peak BED columns.
- `--comparison-table` — comparison, condition 1, and condition 2 columns.
- `--genome` — reference FASTA.
- `--blacklist` — optional blacklist BED.
- `--outdir` — project output directory.
- `--cores` — total worker cores.

## Main outputs

- Bias-corrected and footprint-score tracks for every sample.
- Motif and replicate-aware differential results for each requested comparison.
- One static multi-comparison browser under `reports/review_multi_comparisons/`.

FASTQ preparation is a separate optional step with [`prepare-atac`](prepare-atac.md).
