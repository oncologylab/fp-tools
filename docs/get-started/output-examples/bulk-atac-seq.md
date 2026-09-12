# Bulk ATAC-seq output example

The interactive ENCODE cancer-cell-line report supports condition selection,
motif search, volcano plots, motif logos, aggregate profiles, and SVG export.

Choose two conditions, then search for a motif or select it in the volcano
plot. Read its score difference and adjusted p-value (**FDR**) alongside the aggregate
cut-site profiles. Look for agreement between replicates and a central
reduction in cuts relative to the flanks when assessing footprint protection.
Use the export controls to save a plot for closer review.

Motif scores help compare signals; they do not by themselves identify the
bound protein, especially when several TFs share similar motifs.

<div class="fp-live-demo-wrap fp-output-example">
  <iframe
    class="fp-live-demo fp-report-demo"
    src="../../../ENCODE-Cancer-Cell-lines-Footprinting/"
    title="Interactive bulk ATAC-seq differential-footprint example"
    loading="eager">
  </iframe>
</div>

[Open the complete Bulk output demo](../../reports.md){ target="_blank" rel="noopener noreferrer" }
