# Static Report Demos

`fp-tools` HTML reports are designed to be standalone static files. That makes them a good fit for GitHub Pages.

## Demo Gallery

Curated demo reports should be placed under:

```text
docs/demos/reports/
```

Recommended examples:

- a small `diff-footprints` report with volcano plot, motif logos, and aggregate profiles;
- a `plot-aggregate-batch` report showing multi-sample aggregate browsing;
- a pseudobulk marker report generated from toy or compact public example data.

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
