---
hide:
  - toc
---

# GUI Demo

[fp-tools-gui](api.md#fp-tools-gui) is an optional browser interface for users who prefer forms over long command lines. It still runs the same fp-tools commands and saves reusable YAML configs.

## Open The GUI

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

For a workstation or cloud VM:

```bash
fp-tools-gui --host 0.0.0.0 --port 8891
```

Open the printed URL in a browser. Public access also requires the port to be allowed by your firewall or cloud security group.

## What The GUI Does

<div class="fp-grid">
  <div class="fp-card">
    <h3>Choose a command</h3>
    <p>Bulk ATAC, motif reports, variants, and pseudobulk tools are grouped in the sidebar.</p>
  </div>
  <div class="fp-card">
    <h3>Load examples</h3>
    <p>Use bundled YAML configs or paste your own paths.</p>
  </div>
  <div class="fp-card">
    <h3>Run and review</h3>
    <p>Launch jobs, inspect logs, and reuse the same YAML on the command line.</p>
  </div>
</div>

## Interactive Preview

This browser preview supports the navigation and form controls but does not
run analyses.

<div class="fp-live-demo fp-live-demo-gui">
  <iframe
    src="../demos/gui/fp-tools-gui-static-demo.html"
    title="Interactive fp-tools GUI preview"
    loading="eager"
    allowfullscreen></iframe>
</div>

<p class="fp-live-demo-link"><a href="../demos/gui/fp-tools-gui-static-demo.html">Open the GUI preview in a full browser window</a></p>

## Example Configs

Bundled configs are in `examples/gui_configs/`. The same file can be loaded in the GUI or run directly with [run-workflow](api.md#run-workflow):

```bash
run-workflow --config examples/gui_configs/call_footprints_single.yml
```
