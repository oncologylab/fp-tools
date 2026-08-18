# `call-footprints`

Calculate continuous footprint score tracks from corrected ATAC-seq cut-site
signal.

## Example command

```bash
call-footprints \
  --signals A_corrected.bw B_corrected.bw \
  --sample-names A B \
  --regions merged_peaks.bed \
  --sample-output-root project/samples
```

## Primary inputs

- `--signals` — one or more corrected cut-site bigWig tracks.
- `--sample-names` — labels corresponding to the input tracks.
- `--regions` — BED regions in which footprint scores are calculated.
- `--sample-output-root` — root directory for per-sample footprint outputs.

## Main outputs

- One footprint score bigWig per input signal.
- Optional candidate-footprint BED files for de novo motif discovery.
- Optional multiscale score arrays and candidate summaries.

Continue with [`match-motifs`](match-motifs.md), or see the
[complete `call-footprints` reference](../../api.md#call-footprints).
