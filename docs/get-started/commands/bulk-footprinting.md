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
  --plot-aggregate all \
  --review-format auto \
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

- Bias-corrected and footprint-score tracks for every sample.
- Motif and replicate-aware differential results for each requested comparison.
- A static browser under `reports/review_multi_comparisons/`, or an
  aggregate-free `reports/review_multi_comparisons.html` file.

FASTQ preparation is a separate optional step with [`prepare-atac`](prepare-atac.md).
