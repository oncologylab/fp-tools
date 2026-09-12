# [`summarize-motifs`](../../api.md#summarize-motifs)

Turn completed motif-discovery results into a table of motif sequences,
significance values, and optional matches to known motifs. Add an HTML report
to review the results in a browser.

## Example command

```bash
summarize-motifs --meme-txt project/de_novo/sample/streme/streme.txt \
  --tomtom-tsv project/de_novo/sample/tomtom/tomtom.tsv \
  --out-tsv project/de_novo/sample/motif_summary.tsv --out-html project/de_novo/sample/motif_summary.html
```

## Primary inputs

- `--meme-txt` — discovery output such as `meme.txt`, `streme.txt`, or `dreme.txt`.
- `--tomtom-tsv` — optional Tomtom known-motif matches.
- `--out-tsv` — compact output table for discovered motifs and matches.
- `--out-html` — optional browser report.

Omit `--tomtom-tsv` if you did not run known-motif matching. This command reads
existing results; it does not rerun motif discovery.

## Main outputs

- The `--out-tsv` file — tab-separated discovered motif IDs, consensus sequences, significance values, and known-database matches when available.
- The `--out-html` file — optional portable HTML table containing the same summary and consensus-sequence logos when available.

Open `motif_summary.html` in a browser or import `motif_summary.tsv` into a
spreadsheet. Known-motif matches help identify candidate TF families; they do
not establish which TF is bound.

See the [de novo motif discovery workflow](../workflows/de-novo-motif-discovery.md)
and the [complete `summarize-motifs` reference](../../api.md#summarize-motifs).
