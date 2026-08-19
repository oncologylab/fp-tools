# Tool overview

The guides below cover common inputs and outputs. See the
[API Reference](../api.md) for every option.

## Core analysis

- [`prepare-atac`](commands/prepare-atac.md) — prepare FASTQ inputs as filtered BAM, peak, coverage, and QC outputs.
- [`atac-correct`](commands/atac-correct.md) — correct ATAC-seq cut-site signal for Tn5 sequence bias.
- [`call-footprints`](commands/call-footprints.md) — calculate footprint-score tracks from corrected signal.
- [`match-motifs`](commands/match-motifs.md) — scan motifs and summarize motif-associated footprint scores.
- [`diff-footprints`](commands/diff-footprints.md) — compare conditions or user-defined region sets, with replicate-aware statistics.
- [`normalize-bigwig`](commands/normalize-bigwig.md) — normalize bigWig tracks over shared background regions.

## Visualization and review

- [`plot-aggregate`](commands/plot-aggregate.md) — plot signal around motif sites or user-defined BED regions; also export motif-by-comparison PDF grids.
- [`review-multi-comparisons`](commands/review-multi-comparisons.md) — combine differential reports into one static comparison browser.

## Workflow and interface

- [`bulk-footprinting`](commands/bulk-footprinting.md) — run the complete bulk workflow from aligned BAM and peak files.
- [`sc-footprinting`](commands/sc-footprinting.md) — run pseudobulk and per-cell single-cell ATAC-seq footprinting.
- [`run-yaml-workflow`](commands/run-yaml-workflow.md) — run command-compatible jobs from YAML.
- [`fp-tools-gui`](commands/fp-tools-gui.md) — launch the optional browser interface.

## De Novo Motif Discovery

- [`discover-motifs`](commands/discover-motifs.md) — discover motifs from footprint candidates.
- [`summarize-motifs`](commands/summarize-motifs.md) — summarize discovered motifs and known-motif matches.

## Single-cell ATAC-seq utilities

- [`pseudobulk-fragments`](commands/pseudobulk-fragments.md) — group fragments by cell annotation.
- [`find-signature-fp`](commands/find-signature-fp.md) — plot per-cell footprint-signature heatmaps and UMAPs.
