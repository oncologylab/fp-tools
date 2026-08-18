# `prepare-atac`

Prepare public or local ATAC-seq reads as the filtered alignments, peaks,
coverage tracks, and QC files needed for footprint analysis.

## Example command

```bash
prepare-atac \
  --samples metadata.tsv \
  --genome hg38 \
  --outdir project/raw
```

## Primary inputs

- A TSV or CSV sample sheet containing public accessions, local FASTQ paths, or HTTPS FASTQ links.
- A named `hg38` or `mm10` reference, or a configured custom genome.
- Optional YAML settings for trimming, alignment, filtering, peak calling, QC, and compute resources.

## Main outputs

- Filtered BAM and BAI files, peak BED files, and RP10M coverage bigWigs.
- Per-sample QC files and command logs.
- Merged project peaks, resolved settings, and a downstream sample table.

Continue with [`atac-correct`](atac-correct.md), or see the
[complete `prepare-atac` reference](../../api.md#prepare-atac).
