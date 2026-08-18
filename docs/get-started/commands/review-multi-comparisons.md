# `review-multi-comparisons`

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

- One or more `diff-footprints` HTML reports or a project directory.
- Optional comparison labels and initial panel count.
- Optional project data for completing aggregate profiles not embedded in the source reports.

## Main outputs

- A portable HTML report for selecting comparisons and motifs.
- Coordinated volcano, bar, logo, and aggregate-profile views with SVG export.

Continue with [`plot-motif-aggregate-grid`](plot-motif-aggregate-grid.md), or
see the [complete `review-multi-comparisons` reference](../../api.md#review-multi-comparisons).
