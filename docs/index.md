# fp-tools

`fp-tools` is a command-first Python package for ATAC-seq footprinting, motif-aware differential footprint analysis, aggregate visualization, de novo motif-discovery preparation, and pseudobulk single-cell ATAC workflows.

The PyPI distribution is `fp-tools-bio`; the import package is `fp_tools`.

```bash
pip install fp-tools-bio
```

Optional GUI support:

```bash
pip install "fp-tools-bio[gui]"
```

Optional DESeq2-style fixed-site analysis:

```bash
pip install "fp-tools-bio[deseq2]"
```

![fp-tools regulatory footprinting framework](assets/fp-tools-workflow.png)

## What The Site Contains

- **Workflow**: end-to-end ATAC correction, footprint scoring, differential footprint analysis, aggregate plotting, and pseudobulk routes.
- **Commands**: concise CLI reference with the main entry points.
- **Reports**: static standalone HTML report demos that can be hosted directly on GitHub Pages.
- **GUI Demo**: guidance for the optional Streamlit GUI and external live-demo hosting.
- **API Reference**: generated Python reference for stable package modules.

## Quick Command Chain

```text
atac-correct -> call-footprints -> match-motifs / diff-footprints -> normalize-bigwig -> plot-aggregate
                                      |                 |         -> plot-aggregate-batch
                                      -> motif-discovery
```

For most two-condition analyses, start with `diff-footprints`; it scans motifs internally, handles repeated condition names as replicates, writes differential tables, and can embed aggregate profiles into standalone HTML reports.
