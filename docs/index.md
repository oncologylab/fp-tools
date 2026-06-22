<section class="fp-hero">
  <div>
    <img class="fp-lockup" src="assets/fp_tools_logo_horizontal.svg" alt="fp-tools regulatory footprinting">
    <p class="fp-eyebrow">ATAC-seq footprinting without custom scripting</p>
    <h1>Correct ATAC signal, call footprints, compare TF activity, and review interactive reports.</h1>
    <p class="fp-lede">
      fp-tools is a command-line and browser-guided toolkit for ATAC-seq footprint analysis.
      It supports bulk ATAC, motif-centered differential footprinting, aggregate plots, and
      pseudobulk single-cell ATAC workflows.
    </p>
    <div class="fp-actions">
      <a class="fp-button primary" href="workflow/">Start workflow</a>
      <a class="fp-button" href="reports/">View report demos</a>
      <a class="fp-button" href="gui/">Open GUI guide</a>
    </div>
    <div class="fp-install">pip install fp-tools-bio</div>
  </div>
</section>

## What fp-tools Does

<div class="fp-grid">
  <div class="fp-card">
    <h3>Correct ATAC signal</h3>
    <p>Convert BAM input into bias-corrected cut-site bigWigs.</p>
  </div>
  <div class="fp-card">
    <h3>Score footprints</h3>
    <p>Create footprint score tracks over peaks or candidate regions.</p>
  </div>
  <div class="fp-card">
    <h3>Compare conditions</h3>
    <p>Use motif-aware differential reports for replicates or time courses.</p>
  </div>
  <div class="fp-card">
    <h3>Review reports</h3>
    <p>Open standalone HTML reports in any browser.</p>
  </div>
  <div class="fp-card">
    <h3>Use the GUI</h3>
    <p>Load examples, edit paths, preview YAML, and launch runs.</p>
  </div>
  <div class="fp-card">
    <h3>Run pseudobulk ATAC</h3>
    <p>Group single-cell fragments into cell-type footprint profiles.</p>
  </div>
</div>

## Core Workflow

<div class="fp-command-chain">
  atac-correct -> call-footprints -> diff-footprints -> plot-aggregate-batch
</div>

For most condition comparisons, start with `diff-footprints`. It scans motifs, compares footprint scores, and writes a browser-ready HTML report.

## Quick Install

```bash
pip install fp-tools-bio
```

For the browser GUI:

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```
