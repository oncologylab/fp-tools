---
hide:
  - navigation
  - toc
---

# GUI Demo

The GUI runs the same fp-tools commands and saves reusable YAML
configurations.

Bulk GUI workflows start from coordinate-sorted BAM/BAI files and matching
peak BED files. FASTQ-to-BAM preparation is available separately through the
Linux CLI or Linux container.

```bash
fp-tools-gui
```

<div class="fp-live-demo-wrap">
  <iframe
    class="fp-live-demo fp-gui-demo"
    src="../demos/gui/fp-tools-gui-static-demo.html"
    title="Interactive fp-tools GUI demonstration"
    loading="eager">
  </iframe>
</div>

<a href="../demos/gui/fp-tools-gui-static-demo.html" target="_blank" rel="noopener noreferrer">Open the GUI demonstration in a full page</a>

The GUI includes example configurations. Repository copies under
`examples/gui_configs/` can also be run directly:

```bash
run-yaml-workflow --config examples/gui_configs/call_footprints_single.yml
```
