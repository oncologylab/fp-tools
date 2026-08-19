# [`review-multi-comparisons`](../../api.md#review-multi-comparisons)

Combine differential-footprint reports as a scalable browser bundle or one
self-contained HTML report.

## Example command

```bash
review-multi-comparisons \
  --inputs project/comparisons \
  --output-dir project/reports/review_multi_comparisons
```

## Primary inputs

- `--inputs` — report files or directories containing differential reports.
- `--output-dir` — destination for the complete static bundle.

## Main outputs

- `index.html` with local CSS, JavaScript, gzip data, profile shards, and logos.
- Or one portable HTML containing coordinated volcano, ranked-motif, logo,
  and SVG-export views. Aggregate controls appear only when profiles exist.

## Standalone output

Use `--output-html` instead of `--output-dir`. Aggregate profiles are optional,
and `--labels` keeps repeated condition pairs distinct:

```bash
review-multi-comparisons \
  --inputs baseline/report.html dose1/report.html dose2/report.html \
  --labels Baseline "Dose 1" "Dose 2" \
  --output-html review.html
```

Continue with the `--motif-grid` mode of [`plot-aggregate`](plot-aggregate.md), or
see the [complete `review-multi-comparisons` reference](../../api.md#review-multi-comparisons).
