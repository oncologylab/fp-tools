# `motif-discovery`

Prepare or run de novo motif discovery from candidate footprint intervals or
an existing FASTA file.

## Example command

```bash
motif-discovery \
  --candidates project/samples/sample/footprints/sample_candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo/sample
```

## Primary inputs

- Candidate BED intervals or candidate sequences in FASTA format.
- Reference genome when sequences must be extracted from BED intervals.
- MEME, STREME, or DREME and an optional known-motif database.

## Main outputs

- Candidate FASTA sequences.
- A runnable discovery plan or executed motif-discovery results.
- Optional Tomtom comparisons with known motifs.

Continue with [`motif-summary`](motif-summary.md), or see the
[complete `motif-discovery` reference](../../api.md#motif-discovery).
