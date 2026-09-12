# De novo motif discovery

This optional workflow discovers enriched sequence motifs from candidate
footprint intervals and compares them with known motif databases.

## Main commands

<div class="fp-command-chain" markdown="1">

[`call-footprints`](../commands/call-footprints.md)
<span>→</span>
[`discover-motifs`](../commands/discover-motifs.md)
<span>→</span>
[`summarize-motifs`](../commands/summarize-motifs.md)

</div>

First use [`call-footprints`](../commands/call-footprints.md) with `--output-bed`
to create candidate footprint intervals. Supply that BED file and a FASTA
reference with the same assembly and chromosome names:

```bash
discover-motifs \
  --candidates candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo \
  --execute
```

The command extracts sequences around the candidates, runs STREME, compares
the discovered motifs with JASPAR using Tomtom, and writes a summary. External
tools are installed automatically on first use. Without `--execute`, it only
prepares the sequence file and a script for later execution.

Open `project/de_novo/motif_summary.html` to review the motifs and their known
matches, or use `motif_summary.tsv` for further analysis. Similar TFs can share
a motif, so a database match suggests a motif identity rather than confirming
which protein binds. Use
[`summarize-motifs`](../commands/summarize-motifs.md) to rebuild a summary from
existing MEME, STREME, DREME, or Tomtom results.
