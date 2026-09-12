---
core_nav:
  previous:
    title: atac-correct
    url: get-started/commands/atac-correct/
  next:
    title: match-motifs
    url: get-started/commands/match-motifs/
---

# [`call-footprints`](../../api.md#call-footprints)

Calculate a footprint score at each position in accessible regions using the
corrected bigWigs from `atac-correct`. The score measures local cut-site
depletion relative to the surrounding signal; higher scores indicate stronger
footprint evidence.

## Example command

```bash
call-footprints \
  --signals A_corrected.bw B_corrected.bw \
  --sample-names A B \
  --regions merged_peaks.bed \
  --sample-output-root project/samples
```

## Primary inputs

- `--signals` — one bias-corrected cut-site signal bigWig per sample.
- `--sample-names` — labels in the same order as `--signals`.
- `--regions` — BED intervals in which scores are calculated; normally the project merged, filtered peaks.
- `--sample-output-root` — root represented by `{sample_root}` below.

Replace `A_corrected.bw` and `B_corrected.bw` with your corrected tracks. In a
project, these are under `project/samples/{sample}/atac_correct/`; use
`project/peaks/merged_peaks_filtered.bed` for the regions.

## Main outputs

| Path | Meaning |
| --- | --- |
| `{sample_root}/{sample}/footprints/{sample}_footprints.bw` | Base-resolution footprint score bigWig used by `match-motifs` and `diff-footprints`. |
| `{sample_root}/{sample}/footprints/{sample}_candidate_footprints.bed` | Optional local score maxima for de novo motif discovery; written with `--call-candidates`. |
| user-selected `*.npz` | Optional compressed scale-by-position score arrays when `--score multiscale` is used with an NPZ output option. |

In direct mode, `--output result.bw` writes exactly `result.bw`; multiple
signals written through `--outdir {outdir}` use
`{outdir}/{signal_stem}_footprints.bw`.

Pass the footprint score tracks to `match-motifs` to summarize evidence at
known motif sites. Keep the corrected cut-site tracks for aggregate plots,
which show the signal shape around those sites.

Continue with [`match-motifs`](match-motifs.md), or see the
[complete `call-footprints` reference](../../api.md#call-footprints).
