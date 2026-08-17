---
hide:
  - toc
---

# GUI Demo

The optional GUI runs the same fp-tools commands and saves reusable YAML
configurations.

```bash
pip install "fp-tools-bio[gui]"
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

<a href="../demos/gui/fp-tools-gui-static-demo.html">Open the GUI demonstration in a full page</a>

Example configurations are in `examples/gui_configs/` and can also be run
directly:

```bash
run-workflow --config examples/gui_configs/call_footprints_single.yml
```
