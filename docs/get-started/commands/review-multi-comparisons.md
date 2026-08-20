# [`review-multi-comparisons`](../../api.md#review-multi-comparisons)

Combine differential-footprint reports as a scalable browser bundle or one
self-contained HTML report.

## Example command

```bash
review-multi-comparisons --inputs project/comparisons --output-dir project/reports/review_multi_comparisons \
  --default-comparison "HNF4A + FOXA2" "No HNF4A/FOXA2" \
  --default-aggregate-motifs MA1494.2 MA0484.3 MA0047.4 MA0148.5 MA0046.3 MA0153.2 MA0102.5 MA0466.4 \
  --default-aggregate-plots 8 --documentation-url https://oncologylab.github.io/fp-tools/
```

## Primary inputs

- `--inputs` — report files or directories containing differential reports.
- `--output-dir` — destination for the complete static bundle.
- `--default-comparison` — condition or region pair shown first.
- `--default-aggregate-motifs` — ordered motif panel shown first.
- `--default-aggregate-plots` — initial number of aggregate panels.
- `--documentation-url` — optional link back to the documentation site.

## Main outputs

`{bundle}` is the `--output-dir`:

| Path | Meaning |
| --- | --- |
| `{bundle}/index.html` | Browser entry point; open or publish this file together with the full bundle. |
| `{bundle}/app.js` and `{bundle}/styles.css` | Local application code and styling. |
| `{bundle}/data/metadata.json` | Comparison index and payload checksums. |
| `{bundle}/data/reports/{comparison}.json.gz` | Compact data for one comparison. |
| `{bundle}/data/profiles/` | Aggregate-profile shards loaded on demand. |
| `{bundle}/data/logos/` | Motif logo assets. |

Project mode defaults to
`{project}/reports/review_multi_comparisons/index.html`. The directory is a
portable unit; copying only `index.html` produces a broken report.

## Standalone output

Use `--output-html` instead of `--output-dir`. Aggregate profiles are optional,
and `--labels` keeps repeated condition pairs distinct. The exact output is the
path passed to `--output-html`; it is one portable HTML file with coordinated
volcano, ranked-motif, logo, and SVG-export views. Aggregate controls appear
only when profiles exist. One **Comparison** list selects the exact input record
in `--labels` order, so repeated condition pairs remain distinct.

```bash
review-multi-comparisons --inputs baseline/report.html dose1/report.html dose2/report.html \
  --labels Baseline "Dose 1" "Dose 2" --output-html review.html
```

Continue with the `--motif-grid` mode of [`plot-aggregate`](plot-aggregate.md), or
see the [complete `review-multi-comparisons` reference](../../api.md#review-multi-comparisons).
