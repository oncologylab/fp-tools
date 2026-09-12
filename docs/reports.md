---
hide:
  - navigation
  - toc
---

# Output Demo with ENCODE Cancer Cell Lines

Select any two of the seven cell lines to inspect differential footprints,
motif statistics, logos, and aggregate profiles.

1. Choose **Condition 1** and **Condition 2**. Positive changes indicate higher
   scores in Condition 1.
2. Search for a motif or select a point in the volcano plot to inspect it.
3. Compare the aggregate curves and replicate agreement alongside the motif's
   effect size and adjusted p-value (**FDR**). Use the plot controls to export figures.

In the plot labels, `_up` indicates higher scores in the named condition.
`n.s.` marks motifs that do not meet the report's significance and filtering
thresholds.

These reports summarize motif-associated footprint signals. A motif match or
score difference alone does not establish that a particular TF is bound.

<div class="fp-live-demo-wrap">
  <iframe
    class="fp-live-demo fp-report-demo"
    src="../ENCODE-Cancer-Cell-lines-Footprinting/"
    title="Interactive ENCODE cancer cell-line footprint report"
    loading="eager">
  </iframe>
</div>

<a href="../ENCODE-Cancer-Cell-lines-Footprinting/" target="_blank" rel="noopener noreferrer">Open the output demo in a full page</a>
