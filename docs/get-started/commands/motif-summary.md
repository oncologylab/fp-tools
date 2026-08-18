# `motif-summary`

Summarize MEME, STREME, DREME, and Tomtom results in a compact report.

## Example command

```bash
motif-summary \
  --meme-txt project/de_novo/sample/streme/streme.txt \
  --tomtom-tsv project/de_novo/sample/tomtom/tomtom.tsv \
  --out-tsv project/de_novo/sample/motif_summary.tsv
```

## Primary inputs

- Motif-discovery output in MEME-compatible text format.
- Optional Tomtom match table.
- Output table and optional HTML paths.

## Main outputs

- Motif names and consensus sequences in TSV format.
- Known-database matches when Tomtom results are supplied.
- Optional portable HTML summary.

See the [de novo motif discovery workflow](../workflows/de-novo-motif-discovery.md)
and the [complete `motif-summary` reference](../../api.md#motif-summary).
