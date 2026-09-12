# Single-cell ATAC-seq output example

The PBMC5k example shows per-cell footprint signatures across B cells,
monocytes, and T/NK cells.

Use the heatmap to compare motif signatures across cell groups. In the UMAPs,
each panel places the same cells at the same coordinates; the colors show
either cell labels or the selected motif's signature score. Similar colors
within a group indicate a shared pattern of signal.

## Signature heatmap

<figure class="fp-output-figure">
  <a href="../../../assets/pbmc5k_signature_heatmap.png">
    <img
      src="../../../assets/pbmc5k_signature_heatmap.png"
      alt="Per-cell PBMC5k footprint-signature heatmap grouped into B-cell, monocyte, and T/NK-cell blocks">
  </a>
  <figcaption>Per-cell footprint-signature heatmap. Columns are cells grouped by broad cell type; rows are motif signatures.</figcaption>
</figure>

## Eight marker UMAPs

<figure class="fp-output-figure">
  <a href="../../../assets/pbmc5k_eight_marker_umaps.png">
    <img
      src="../../../assets/pbmc5k_eight_marker_umaps.png"
      alt="PBMC5k broad-cell-type UMAP and footprint-signature UMAPs for STAT6, FOSB, CEBPA, IRF8, RELA, ZNF683, NR4A1, and SMAD3">
  </a>
  <figcaption>Broad cell-type annotation and per-cell footprint-signature UMAPs for STAT6, FOSB, CEBPA, IRF8, RELA, ZNF683, NR4A1, and SMAD3.</figcaption>
</figure>

Select either figure to open the full-resolution image. The
[Single-cell workflow](../workflows/single-cell.md) describes the commands used
to produce these outputs.

The signatures are smoothed using neighboring cells to reduce sparse signal.
They describe relative motif-associated patterns, not direct measurements of
TF binding in individual cells.
