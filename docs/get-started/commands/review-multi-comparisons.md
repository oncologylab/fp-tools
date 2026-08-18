# [`review-multi-comparisons`](../../api.md#review-multi-comparisons)

Combine multiple differential-footprint reports into one interactive review
page.

## Example command

```bash
review-multi-comparisons \
  --outdir project \
  --display-panels 8 \
  --aggregate-legends hide \
  --recompute-missing-aggregate-profiles \
  --cores 16
```

## Primary inputs

- `--outdir` — project directory containing differential-footprint reports.
- `--display-panels` — number of comparison panels shown initially.
- `--aggregate-legends` — initial aggregate-legend visibility.
- `--recompute-missing-aggregate-profiles` — calculate profiles absent from source reports.
- `--cores` — worker processes used for profile calculation.

## Main outputs

- A portable HTML report for selecting comparisons and motifs.
- Coordinated volcano, bar, logo, and aggregate-profile views with SVG export.

Continue with [`plot-motif-aggregate-grid`](plot-motif-aggregate-grid.md), or
see the [complete `review-multi-comparisons` reference](../../api.md#review-multi-comparisons).
