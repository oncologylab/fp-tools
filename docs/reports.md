# Static Report Demos

`fp-tools` HTML reports are designed to be standalone static files. That makes them a good fit for GitHub Pages.

## Demo Gallery

<div class="fp-grid">
  <div class="fp-card">
    <h3>Differential footprint report</h3>
    <p>Volcano plot, motif logos, comparison table, and embedded aggregate profiles from a compact two-condition run.</p>
  </div>
  <div class="fp-card">
    <h3>Aggregate browser</h3>
    <p>Multi-sample and multi-TF aggregate profiles with searchable TF selection and configurable layouts.</p>
  </div>
  <div class="fp-card">
    <h3>Pseudobulk marker review</h3>
    <p>Marker-focused single-cell ATAC pseudobulk output built from toy or compact public example data.</p>
  </div>
</div>

Curated demo reports should be placed under `docs/demos/reports/`.

!!! note
    Avoid committing very large full benchmark review outputs to the docs site. For large reports, publish them as GitHub Release assets or in the companion data repository, then link to them from this page.

## Add A Report

1. Generate the report with a small fixture or curated public example.
2. Copy the standalone HTML file to `docs/demos/reports/`.
3. Add a link below.

Example link format:

```markdown
- [Differential footprint report](demos/reports/diff-footprints-demo.html)
```

## Current Links

- Demo report slots are ready; add curated standalone HTML files under `docs/demos/reports/`.

## Why Static Reports Work

The report HTML embeds the JavaScript and compressed JSON payload needed for review. GitHub Pages can serve it without a backend, database, or Python runtime.
