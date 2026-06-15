# Workflow Figure Revision Prompts

These prompts align the existing workflow schematics with the current supported fp-tools manuscript scope. They are intended for regenerating the bitmap figures, not for manuscript prose.

## Overall workflow figure (`manuscript/figures/fp-tools-workflow.png`)

Revise the figure as a supported-methods workflow, not a roadmap. Use the title `fp-tools workflow for ATAC-seq footprinting and motif-centered analysis`.

Show one central green core with five sequential steps:

1. `Tn5 bias correction` - corrected cut-site tracks and QC.
2. `Footprint calling` - center-versus-flank footprint scores and optional candidate BED output.
3. `Motif matching` - known-motif scanning, per-sample motif-site tables, and bound-site summaries.
4. `Differential footprint analysis` - two-condition, replicate-aware, or time-course comparisons.
5. `Aggregate reporting` - volcano plot, motif logo, replicate-aware aggregate profiles, and editable SVG export.

Keep optional modules as side branches, not as a large competing panel:

- `Pseudobulk fragments`: single-cell ATAC fragments plus annotations become group-level fragment files and cut-site tracks.
- `De novo motif discovery`: candidate BED/FASTA supports de novo-only discovery or known database plus de novo supplement.
- `Reproducible workflows`: CLI and YAML/batch execution record parameters, versions, and source tables.

Remove or avoid these items in the figure:

- `BINDetect` as a public-facing label.
- `Supervised TFBS prediction`, `variant scoring`, `competition decomposition`, and `multiscale / nucleosome-aware scoring` as current supported manuscript claims.
- `first-version`, `prototype`, `future`, `legacy`, or `scaffold` language.
- Any implication that samples are aligned at this stage. The comparison uses a shared peak/motif-site universe after peak merging and motif scanning.

Use consistent labels across the paper: `B cell`, `T cell`, `replicate`, `condition mean`, `motif sites`, `corrected cut-site signal`, and `differential footprint report`.

## Replicate-aware differential footprint figure

Use the title `Replicate-aware differential footprint analysis`.

Recommended flow:

1. `Shared analysis universe`: merged peaks and motif sites define common loci across all samples.
2. `Replicate inputs`: corrected footprint or cut-site tracks from B-cell Rep1/Rep2 and T-cell Rep1/Rep2.
3. `Optional normalization`: none, condition-quantile, or sample-quantile normalization; show that the same normalization is used for aggregate profiles.
4. `Differential summary`: condition means, replicate variation, log2 fold change or delta score, p-value, q-value, and replicate support.
5. `Interactive report`: volcano plot, selected motif summary, motif logo, and aggregate profiles with thin replicate lines plus thick condition means.

Do not use `BINDetect` in the title or panel labels. If backward compatibility must be shown, put it only in a small footnote outside the main graphic.

## De novo motif discovery figure

Use the title `De novo motif discovery as standalone analysis and database supplement`.

Show two parallel modes from the same candidate footprints:

- `De novo-only`: candidate BED -> candidate-centered FASTA -> STREME/MEME -> discovered motifs -> aggregate validation.
- `Database + de novo supplement`: known motif database, such as JASPAR2026, plus candidate-derived motifs -> differential footprint analysis -> added motif families in the same report.

Make clear that Tomtom provides motif-similarity labels, not definitive TF identity.

## Pseudobulk figure

Use the title `Pseudobulk ATAC fragments for motif-centered footprinting`.

Recommended flow:

1. 10x-style fragments and cell annotations.
2. Group by cell type or donor plus cell type.
3. Write group-level fragment files, manifest, QC table, and CPM-normalized cut-site bigWigs.
4. Scan JASPAR2026 motifs in peaks and aggregate around exact motif centers.
5. Show lineage examples with natural motif-hit counts, not fixed subsamples: SPIB, RUNX3, CEBPB, and CEBPA.

The workflow figure should show the pipeline; the aggregate plot should be presented as example output below or in a separate panel.
