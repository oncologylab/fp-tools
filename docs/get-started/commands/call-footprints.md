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

Calculate a continuous footprint score from bias-corrected cut-site signal
within accessible regions. Higher local depletion relative to flanking signal
produces stronger footprint evidence.

## Example command

```bash
call-footprints --signals A_corrected.bw B_corrected.bw --sample-names A B \
  --regions merged_peaks.bed --sample-output-root project/samples
```

## Primary inputs

- `--signals` — one bias-corrected cut-site signal bigWig per sample.
- `--sample-names` — labels in the same order as `--signals`.
- `--regions` — BED intervals in which scores are calculated; normally the project merged, filtered peaks.
- `--sample-output-root` — root represented by `{sample_root}` below.

## Main outputs

| Path | Meaning |
| --- | --- |
| `{sample_root}/{sample}/footprints/{sample}_footprints.bw` | Base-resolution footprint score bigWig used by `match-motifs` and `diff-footprints`. |
| `{sample_root}/{sample}/footprints/{sample}_candidate_footprints.bed` | Optional local score maxima for de novo motif discovery; written with `--call-candidates`. |
| user-selected `*.npz` | Optional compressed scale-by-position score arrays when `--score multiscale` is used with an NPZ output option. |

In direct mode, `--output result.bw` writes exactly `result.bw`; multiple
signals written through `--outdir {outdir}` use
`{outdir}/{signal_stem}_footprints.bw`.

An experimental dual-geometry arm is available for method evaluation:

```bash
call-footprints \
  --signal sample_corrected.bw \
  --regions merged_peaks.bed \
  --score hybrid \
  --output sample_hybrid_footprints.bw
```

`hybrid` retains the standard footprint score and adds a low-weight,
locally standardized 33 bp central-depletion channel with symmetric 32 bp
shoulders. It improved wide CTCF/REST footprints in the locked K562/HepG2
experiment, but reduced some JUND/MAX metrics. It is therefore opt-in and is
not the production default.

Continue with [`match-motifs`](match-motifs.md), or see the
[complete `call-footprints` reference](../../api.md#call-footprints).
