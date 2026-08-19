---
core_nav:
  previous:
    title: call-footprints
    url: get-started/commands/call-footprints/
  next:
    title: diff-footprints
    url: get-started/commands/diff-footprints/
---

# [`match-motifs`](../../api.md#match-motifs)

Scan accessible regions for motif instances, measure the footprint score at
each instance, and classify sample-specific bound and unbound sites. Run this
when motif locations and per-sample motif summaries are needed.

## Example command

```bash
match-motifs --signals A_footprints.bw B_footprints.bw --sample-names A B --genome hg38.fa.gz \
  --peaks merged_peaks.bed --motif-db jaspar2026_vertebrates --sample-output-root project/samples
```

## Primary inputs

- `--signals` — one footprint score bigWig per sample.
- `--sample-names` — sample labels in the same order as `--signals`.
- `--genome` — assembly-matched reference FASTA used to scan motif sequences.
- `--peaks` — accessible-region BED searched for motif instances.
- `--motif-db` — packaged motif collection; the example uses JASPAR 2026 vertebrates.
- `--sample-output-root` — root represented by `{sample_root}` below.

## Main outputs

For each `{sample}`, the default output directory is
`{sample_root}/{sample}/match_motifs/`:

| Path | Meaning |
| --- | --- |
| `motif_matches_results.txt` | Tab-separated motif summary with site counts and per-sample mean scores. |
| `motif_matches_results.xlsx` | Excel copy of the motif summary unless `--skip-excel` is used. |
| `motif_matches_distances.txt` | Motif-similarity distances used for motif clustering. |
| `motif_matches_replicate_motif_score_matrix.tsv` | Motif-by-sample footprint score matrix when multiple samples are analyzed together. |
| `cache/motif_sites.tsv.gz` | Compact scanned motif-site cache reusable by differential analysis. |
| `cache/background_scores.tsv.gz` | Compact background-score cache. |
| `{motif}/beds/{motif}_{sample}_all.bed` | All scanned instances for one motif. |
| `{motif}/beds/{motif}_{sample}_bound.bed` | Instances classified as bound in the sample. |
| `{motif}/beds/{motif}_{sample}_unbound.bed` | Instances classified as unbound in the sample. |

`{motif}` follows the selected `--naming` convention, such as
`CTCF_MA0139.2`. `--motif-outputs summary` omits the per-motif BED files but
keeps the summary and reusable caches.

[Download a representative ENCODE replicate score matrix](../../demos/qc/encode/A549_motif_score_matrix.tsv).
Continue with [`diff-footprints`](diff-footprints.md), or see the
[complete `match-motifs` reference](../../api.md#match-motifs).
