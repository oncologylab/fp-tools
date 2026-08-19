# Optional de novo motif discovery

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

First use [`call-footprints`](../commands/call-footprints.md) with candidate
calling enabled. Then run motif discovery:

```bash
discover-motifs \
  --candidates candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo
```

[`summarize-motifs`](../commands/summarize-motifs.md) converts MEME/STREME/DREME and
Tomtom results into a compact table and optional HTML report.
