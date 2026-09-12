---
hide:
  - toc
---

<div class="fp-home" markdown="1">
<div class="fp-home-logo">
  <img src="assets/fp_tools_logo_horizontal.svg" alt="fp-tools">
</div>
<div class="fp-home-copy" markdown="1">

# Tn5-based chromatin footprinting and regulatory motif analysis

Use `fp-tools` to analyze ATAC-seq and CUT&Tag data, score chromatin footprints,
and compare motif-associated signals between conditions. Run a complete
workflow or individual steps from the command line, or use the GUI to prepare
and run the same analyses.

<div class="fp-badges">
  <a href="https://pypi.org/project/fp-tools-bio/">PyPI</a>
  <a href="https://github.com/oncologylab/fp-tools">GitHub</a>
  <span>Python 3.11–3.13</span>
  <span>MIT license</span>
</div>

- Correct Tn5 sequence bias and score footprints from BAM and peak files.
- Compare motif-associated footprint scores across samples and replicates.
- Analyze pseudobulk and per-cell footprint signatures.
- Export static figures and portable interactive HTML reports.

Choose your starting point:

- **BAM/BAI and peak BED files:** follow the runnable [bulk ATAC-seq workflow](get-started/workflows/bulk-atac-seq.md).
- **Single-cell fragments:** follow the [single-cell workflow](get-started/workflows/single-cell.md).
- **Existing fp-tools outputs:** open the [output demo](reports.md) or choose a command in the [tool overview](get-started/tool-overview.md).

Linux command-line and container users can optionally prepare FASTQ files with
[`prepare-atac`](get-started/commands/prepare-atac.md).

[Install fp-tools](get-started/installation.md){ .fp-text-link }
&nbsp;·&nbsp;
[View the tool overview](get-started/tool-overview.md){ .fp-text-link }

</div>
</div>
