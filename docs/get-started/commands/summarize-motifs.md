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

- Motif names and consensus sequences in TSV format.
- Known-database matches when Tomtom results are supplied.
- Optional portable HTML summary.

See the [de novo motif discovery workflow](../workflows/de-novo-motif-discovery.md)
and the [complete `summarize-motifs` reference](../../api.md#summarize-motifs).
