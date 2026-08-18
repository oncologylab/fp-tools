# Tool overview

The command names below link to short practical guides. The
[API Reference](../api.md) contains the complete option listings.

## Core analysis

- [`prepare-atac`](commands/prepare-atac.md) — prepare public or local FASTQ files as filtered BAM, peak, coverage, and QC outputs.
- [`atac-correct`](commands/atac-correct.md) — correct ATAC-seq cut-site signal for Tn5 sequence bias.
- [`call-footprints`](commands/call-footprints.md) — calculate footprint score tracks from corrected signal.
- [`match-motifs`](commands/match-motifs.md) — scan motifs and summarize motif-associated footprint scores.
- [`diff-footprints`](commands/diff-footprints.md) — compare motif-associated footprint scores across conditions and replicates.
- [`normalize-bigwig`](commands/normalize-bigwig.md) — normalize bigWig tracks using shared background regions.

## Visualization and review

- [`plot-aggregate`](commands/plot-aggregate.md) — plot aggregate signal around motif sites or other region sets.
- [`review-multi-comparisons`](commands/review-multi-comparisons.md) — browse multiple differential-footprint reports together.
- [`plot-motif-aggregate-grid`](commands/plot-motif-aggregate-grid.md) — export motif-by-comparison aggregate profiles as a multipage PDF.

## Workflow and interface

- [`run-workflow`](commands/run-workflow.md) — run command-compatible jobs from a YAML configuration.
- [`fp-tools-gui`](commands/fp-tools-gui.md) — launch the optional browser interface.

## Motifs and variants

- [`motif-discovery`](commands/motif-discovery.md) — prepare or run de novo motif discovery from footprint candidates.
- [`motif-summary`](commands/motif-summary.md) — summarize discovered motifs and known-motif matches.
- [`fp-tools-score-variants`](commands/fp-tools-score-variants.md) — annotate variants with footprint, sequence, motif, or model score changes.

## Single-cell analysis

- [`pseudobulk-fragments`](commands/pseudobulk-fragments.md) — group single-cell ATAC fragments by annotation.
- [`find-signature-fp`](commands/find-signature-fp.md) — plot per-cell footprint-signature heatmaps and UMAPs.
- [`pseudobulk-footprints`](commands/pseudobulk-footprints.md) — run the complete pseudobulk footprinting workflow.
