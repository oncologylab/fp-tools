---
hide:
  - toc
---

# Report Demos

`fp-tools` writes standalone HTML reports. They can be opened locally, shared with collaborators, or hosted on GitHub Pages.

## Interactive Differential Footprint Report

K562 versus HepG2 is shown initially. Use the condition selectors to browse
all 21 pairwise comparisons across seven ENCODE cancer cell lines.

<div class="fp-live-demo fp-live-demo-report">
  <iframe
    src="../ENCODE-Cancer-Cell-lines-Footprinting/"
    title="Interactive fp-tools differential footprint report"
    loading="eager"
    allowfullscreen></iframe>
</div>

<p class="fp-live-demo-link"><a href="../ENCODE-Cancer-Cell-lines-Footprinting/">Open the report in a full browser window</a></p>

## Typical Report Contents

<div class="fp-grid">
  <div class="fp-card">
    <h3>Volcano summary</h3>
    <p>Quickly identify motifs with condition-biased footprint scores.</p>
  </div>
  <div class="fp-card">
    <h3>Motif table</h3>
    <p>Search, sort, and review motif-level statistics.</p>
  </div>
  <div class="fp-card">
    <h3>Aggregate profiles</h3>
    <p>Inspect corrected cut-site signal around motif centers.</p>
  </div>
</div>

These reports do not need a Python server after they are created.

## Aggregate Browser

<figure class="fp-wide-image">
  <img src="../assets/interface_plot_aggregate_batch_html.png"
       alt="Multi-panel aggregate footprint browser comparing B-cell and T-cell profiles">
  <figcaption>Multi-sample aggregate reports support selectable panels,
  condition means, individual replicates, and editable SVG export.</figcaption>
</figure>
