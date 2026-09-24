# [`review-multi-comparisons`](../../api.md#review-multi-comparisons)

Review several `diff-footprints` comparisons in one place. Combine the
existing reports into a browser bundle or a single HTML file, then switch
between comparisons to explore motif statistics and aggregate profiles.

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

`{bundle}` is the `--output-dir`:

| Path | Meaning |
| --- | --- |
| `{bundle}/index.html` | Report entry page, served together with the other bundle files. |
| `{bundle}/app.js`, `{bundle}/plot_controls.js`, and `{bundle}/styles.css` | Local application code, shared plot behavior, and styling. |
| `{bundle}/data/metadata.json` | Comparison index and payload checksums. |
| `{bundle}/data/reports/{comparison}.json.gz` | Compact data for one comparison. |
| `{bundle}/data/profiles/` | Aggregate-profile shards loaded on demand. |
| `{bundle}/data/logos/` | Motif logo assets. |

Project mode defaults to
`{project}/reports/review_multi_comparisons/index.html`. Keep the entire
directory together when copying or publishing a bundle.

## Open the bundle

The bundle loads its data through a web server. To view it locally, run:

```bash
python -m http.server 8000 --bind 127.0.0.1 \
  --directory project/reports/review_multi_comparisons
```

Open [http://127.0.0.1:8000/](http://127.0.0.1:8000/) in your browser. Stop the
server with `Ctrl+C` when you finish. To open a report directly without a
server, use the standalone output below.

To choose the bundle's opening view, add `--default-comparison` with the two
condition or region labels, `--default-aggregate-motifs` with the motif names
or IDs in display order, and `--default-aggregate-plots` with the number of
panels. The selected labels and motifs must exist in the input reports.
Use `--documentation-url` to add a documentation link.

## Standalone output

Use `--output-html` instead of `--output-dir` to create one HTML file that you
can open directly or share. It includes volcano plots, ranked motifs, logos,
and SVG exports. Aggregate controls appear when the input reports contain
profiles. Use `--labels` to give each input report a distinct name in the
**Comparison** list, especially when reports compare the same condition pair.

In the ranked-motif plot, switch between differential footprint score and
`-log10(p-value)` to change the ranking view. Use the legend to interpret the
color scale and direction of change. In the volcano plot, use **Label TFs** to label
selected TF names, motif IDs, or output prefixes, separated by commas. Select
`(none)` in the highlight control to clear highlighting. SVG exports retain
the selected comparison label and plot settings.

```bash
review-multi-comparisons --inputs baseline/report.html dose1/report.html dose2/report.html \
  --labels Baseline "Dose 1" "Dose 2" --output-html review.html
```

## Adjust plots and compare results

Open **Plot ranges and display filters** to change the view:

- **Single comparison** shows one comparison with several selected motifs.
  **Side by side** shows two to eight comparisons. Click a selected motif to
  inspect that same motif across the comparison panels.
- **Shared volcano X/Y** and **Shared waterfall X** use the same scales across
  every comparison in the report, including comparisons not currently visible.
  Shared waterfall scaling follows the ranking, top-motif count and filter.
- Each plot has range sliders and numeric minimum/maximum inputs. A custom
  range overrides sharing for that plot; **Reset** returns it to automatic or
  shared scaling. A notice identifies ranges that clip data.
- **Shared aggregate Y** links the displayed aggregate footprint plots. You
  can adjust that shared range or override individual plots. These settings
  change the axes, not the signal values.

**Top motifs (total)** sets the number of waterfall bars, with positive and
negative changes represented when available. The initial display filter requires
**false discovery rate (FDR) ≤ 0.05** and **absolute footprint-score change
(|ΔFP|) ≥ 0.1**. Choose FDR or raw p-value, enter both cutoffs, and click
**Apply filters**. Both criteria must be met. Missing values do not pass.

Failing volcano points remain gray unless **Hide failing volcano points** is
checked. Only passing motifs enter the waterfall or receive differential
highlighting. You can still select a failing motif's aggregate for inspection;
the report labels it accordingly. The volcano Y-axis always shows
`−log10(raw p-value)`, including when filtering by FDR.

Axis settings are remembered while the report is open. Figure exports retain
the visible settings and filter criteria. Original statistics and full-result
table downloads are unchanged; these controls do not rerun the analysis.

Continue with the `--motif-grid` mode of [`plot-aggregate`](plot-aggregate.md), or
see the [complete `review-multi-comparisons` reference](../../api.md#review-multi-comparisons).
