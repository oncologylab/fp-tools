# [`discover-motifs`](../../api.md#discover-motifs)

Find recurring DNA sequence patterns in candidate footprint regions, without
starting from a known motif list. Provide candidate intervals and a reference
genome, or use sequences you have already extracted into a FASTA file.

## Example command

```bash
discover-motifs --candidates project/samples/sample/footprints/sample_candidate_footprints.bed --genome hg38.fa.gz \
  --flank 75 --method streme --known-motif-db jaspar2026_vertebrates --outdir project/de_novo/sample --execute
```

## Primary inputs

- `--candidates` — candidate-footprint BED intervals.
- `--genome` — reference FASTA matching the candidate coordinates; required with `--candidates`.
- `--flank` — bases included on each side of a candidate center.
- `--method` — discovery method; the example uses STREME.
- `--known-motif-db` — optional known-motif database for Tomtom matching.
- `--outdir` — directory for candidate FASTA files and discovery results.
- `--execute` — run discovery immediately using the managed MEME Suite runtime.

In the example, `--flank 75` extracts up to 75 bases on each side of each
candidate center.
Without `--execute`, the command writes the FASTA and a command script for
inspection. It does not run discovery or create result summaries.

## Main outputs

`{outdir}` is the selected discovery directory:

| Path | Meaning |
| --- | --- |
| `{outdir}/candidate_sequences.fa` | Reference sequences extracted around candidate footprint intervals. |
| `{outdir}/run_motif_discovery.sh` | Commands for discovery, optional known-motif matching, and summary reports. |
| `{outdir}/streme/streme.txt`, `meme/meme.txt`, or `dreme/dreme.txt` | Discovered motif models from the selected method when `--execute` is used. |
| `{outdir}/tomtom/tomtom.tsv` | Optional similarity matches to the selected known-motif database. |
| `{outdir}/motif_summary.tsv` and `motif_summary.html` | Motif tables written after successful execution. |

Open `motif_summary.html` to review discovered motifs and any known-motif
matches. These matches suggest motif identities; similar TFs can share the
same motif. To regenerate a summary from existing results, use
[`summarize-motifs`](summarize-motifs.md). See the
[complete `discover-motifs` reference](../../api.md#discover-motifs).

## Choose your sequence input

Create candidate intervals with `call-footprints --output-bed`. If you already
have a sequence FASTA, pass it with `--fasta` instead of `--candidates` and
`--genome`. In that case the command uses your existing sequences directly.
