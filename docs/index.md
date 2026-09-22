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

`fp-tools` analyzes bulk and single-cell
ATAC-seq (assay for transposase-accessible chromatin using sequencing) data. It
corrects sequence bias, scores footprints around DNA motifs, and compares
samples or cell groups. Use the command-line interface (CLI) or graphical user
interface (GUI).

A footprint or motif match does not prove transcription factor (TF) binding.
CUT&Tag (cleavage under targets and tagmentation) needs assay-specific controls
and interpretation.

New to these file formats? See the [practical glossary](get-started/tool-overview.md#practical-glossary).

<div class="fp-badges">
  <a href="https://pypi.org/project/fp-tools-bio/">PyPI</a>
  <a href="https://github.com/oncologylab/fp-tools">GitHub</a>
  <span>Python 3.11–3.13</span>
  <span>MIT license</span>
</div>

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
