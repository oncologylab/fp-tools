# [`summarize-motifs`](../../api.md#summarize-motifs)

Summarize MEME, STREME, DREME, and Tomtom results in a compact report.

## Example command

```bash
summarize-motifs \
  --meme-txt project/de_novo/sample/streme/streme.txt \
  --tomtom-tsv project/de_novo/sample/tomtom/tomtom.tsv \
  --out-tsv project/de_novo/sample/motif_summary.tsv
```

## Primary inputs

- `--meme-txt` — MEME-compatible discovery output.
- `--tomtom-tsv` — optional Tomtom known-motif matches.
- `--out-tsv` — compact output table for discovered motifs and matches.

## Main outputs

- the exact `--out-tsv` path — tab-separated discovered motif IDs, consensus sequences, significance values, and known-database matches when available.
- the exact `--out-html` path — optional portable HTML table containing the same summary and motif logos when available.

The command does not rename the requested output prefix; in the example the
primary file is `project/de_novo/sample/motif_summary.tsv`.

See the [de novo motif discovery workflow](../workflows/de-novo-motif-discovery.md)
and the [complete `summarize-motifs` reference](../../api.md#summarize-motifs).
