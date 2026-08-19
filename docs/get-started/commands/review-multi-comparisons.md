# [`review-multi-comparisons`](../../api.md#review-multi-comparisons)

Combine multiple differential-footprint reports into one static browser with
two condition selectors.

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
- Coordinated volcano, bar, logo, and aggregate-profile views with SVG export.

Continue with the `--motif-grid` mode of [`plot-aggregate`](plot-aggregate.md), or
see the [complete `review-multi-comparisons` reference](../../api.md#review-multi-comparisons).
