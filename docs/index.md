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

`fp-tools` provides command-line and browser workflows for bulk and single-cell
footprinting of Tn5-based chromatin profiling data, including ATAC-seq,
CUT&Tag, and related assays. It connects bias correction, footprint scoring,
motif analysis, replicate comparisons, and reusable reports.

<div class="fp-badges">
  <a href="https://pypi.org/project/fp-tools-bio/">PyPI</a>
  <a href="https://github.com/oncologylab/fp-tools">GitHub</a>
  <span>Python 3.11–3.13</span>
  <span>MIT license</span>
</div>

- Process raw reads or start from aligned Tn5-based chromatin data.
- Compare motif-associated footprint scores across samples and replicates.
- Analyze pseudobulk and per-cell footprint signatures.
- Export static figures and portable interactive HTML reports.

Choose your starting point:

- **FASTQ files:** prepare and QC reads with [`prepare-atac`](get-started/commands/prepare-atac.md).
- **BAM and peak files:** follow the runnable [bulk ATAC-seq workflow](get-started/workflows/bulk-atac-seq.md).
- **Single-cell fragments:** follow the [single-cell workflow](get-started/workflows/single-cell.md).
- **Existing fp-tools outputs:** open the [output demo](reports.md) or choose a command in the [tool overview](get-started/tool-overview.md).

[Install fp-tools](get-started/installation.md){ .fp-text-link }
&nbsp;·&nbsp;
[View the tool overview](get-started/tool-overview.md){ .fp-text-link }

</div>
</div>
