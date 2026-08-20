# [`discover-motifs`](../../api.md#discover-motifs)

Prepare or run de novo motif discovery from candidate footprint intervals or
an existing FASTA file.

## Example command

```bash
discover-motifs --candidates project/samples/sample/footprints/sample_candidate_footprints.bed --genome hg38.fa.gz \
  --flank 75 --method streme --known-motif-db jaspar2026_vertebrates --outdir project/de_novo/sample --execute
```

## Primary inputs

- `--candidates` — candidate-footprint BED intervals.
- `--genome` — reference genome used to extract candidate sequences.
- `--flank` — bases included on each side of a candidate center.
- `--method` — discovery method; the example uses STREME.
- `--known-motif-db` — optional known-motif database for Tomtom matching.
- `--outdir` — directory for candidate FASTA files and discovery results.
- `--execute` — run discovery immediately using the managed MEME Suite runtime.

## Main outputs

`{outdir}` is the selected discovery directory:

| Path | Meaning |
| --- | --- |
| `{outdir}/candidate_sequences.fa` | Reference sequences extracted around candidate footprint intervals. |
| `{outdir}/run_motif_discovery.sh` | Reproducible MEME/DREME/STREME command plan. |
| `{outdir}/{method}/streme.txt` or the method-equivalent MEME output | De novo motif models when `--execute` is used. |
| `{outdir}/tomtom/tomtom.tsv` | Optional similarity matches to the selected known-motif database. |
| `{outdir}/motif_summary.tsv` and `motif_summary.html` | Summary targets written by the generated plan after discovery and matching complete. |

Continue with [`summarize-motifs`](summarize-motifs.md), or see the
[complete `discover-motifs` reference](../../api.md#discover-motifs).
