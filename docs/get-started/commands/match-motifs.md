# [`match-motifs`](../../api.md#match-motifs)

Scan accessible regions for motif sites and summarize their footprint scores
for one or more samples.

## Example command

```bash
match-motifs \
  --signals A_footprints.bw B_footprints.bw \
  --sample-names A B \
  --genome hg38.fa.gz \
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --sample-output-root project/samples
```

## Primary inputs

- `--signals` — footprint score bigWig tracks.
- `--sample-names` — labels corresponding to the input tracks.
- `--genome` — reference genome FASTA.
- `--peaks` — accessible-region BED file searched for motif sites.
- `--motif-db` — built-in motif database name.
- `--sample-output-root` — root directory for per-sample motif results.

## Main outputs

- Per-sample motif binding summaries.
- Motif-site and score caches reusable by `diff-footprints`.
- Optional bound, unbound, and all-site BED files for each motif.

Continue with [`diff-footprints`](diff-footprints.md), or see the
[complete `match-motifs` reference](../../api.md#match-motifs).
