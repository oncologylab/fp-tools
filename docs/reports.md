# Report Demos

`fp-tools` writes standalone HTML reports. They can be opened locally, shared with collaborators, or hosted on GitHub Pages.

## ENCODE Cancer Cell-line Footprinting

The [browser for seven ENCODE-profiled cancer cell lines included in this analysis](https://oncologylab.github.io/fp-tools/ENCODE-Cancer-Cell-lines-Footprinting/)
contains all 21 replicate-aware pairwise comparisons and all 1,019 JASPAR 2026
vertebrate motifs included in the harmonized, sample-quantile-normalized
analysis. Two condition selectors open each comparison in the same report layout
used by the standalone `diff-footprints` output.

## Differential Footprint Demo

<div class="fp-demo-callout">
  <h3>K562 vs HepG2 differential footprint report</h3>
  <p>Example output from a motif-aware condition comparison with volcano results, motif tables, logos, and aggregate footprint profiles.</p>
  <a class="fp-button primary" href="../demos/reports/diff_footprints_K562_HepG2.html">Open report demo</a>
</div>

<figure class="fp-wide-image">
  <a href="../demos/reports/diff_footprints_K562_HepG2.html">
    <img src="../assets/interface_diff_footprints_html.png"
         alt="Interactive differential footprint report with motif statistics, volcano plot, and aggregate profile">
  </a>
  <figcaption>Interactive motif selection, differential evidence, aggregate
  profiles, and editable SVG export in one portable report.</figcaption>
</figure>

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
