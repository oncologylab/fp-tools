# [`discover-motifs`](../../api.md#discover-motifs)

Prepare or run de novo motif discovery from candidate footprint intervals or
an existing FASTA file.

## Example command

```bash
discover-motifs \
  --candidates project/samples/sample/footprints/sample_candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo/sample
```

## Primary inputs

- `--candidates` — candidate-footprint BED intervals.
- `--genome` — reference genome used to extract candidate sequences.
- `--flank` — bases included on each side of a candidate center.
- `--method` — discovery method; the example uses STREME.
- `--known-motif-db` — optional known-motif database for Tomtom matching.
- `--outdir` — directory for candidate FASTA files and discovery results.

## Main outputs

- Candidate FASTA sequences.
- A runnable discovery plan or executed results.
- Optional Tomtom comparisons with known motifs.

Continue with [`summarize-motifs`](summarize-motifs.md), or see the
[complete `discover-motifs` reference](../../api.md#discover-motifs).
