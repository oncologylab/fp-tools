---
hide:
  - navigation
  - toc
---

# GUI Demo

Explore the interface below: choose a command, inspect its settings, and view
an example run. This demonstration does not run analyses or read your files.
Install fp-tools to use the working GUI, which runs the same commands as the
CLI and saves reusable YAML configurations.

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

In the installed GUI, load an example configuration and replace its paths with
your own input files. You can save the settings as YAML and run them from the
command line. From a repository checkout, try an example without launching
the analysis:

```bash
run-yaml-workflow --config examples/gui_configs/call_footprints_single.yml --dry-run
```

Remove `--dry-run` to run it. See the [GUI guide](get-started/commands/fp-tools-gui.md)
for setup, validation, and run-history instructions.
