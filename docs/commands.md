# Commands

Install the package, then inspect any command with `--help`.

```bash
pip install fp-tools-bio
```

<div class="fp-grid">
  <div class="fp-card">
    <h3>Start with bulk ATAC</h3>
    <p>atac-correct, call-footprints, match-motifs, and diff-footprints cover the core footprint workflow.</p>
  </div>
  <div class="fp-card">
    <h3>Review reports</h3>
    <p>plot-aggregate generates static figures or standalone HTML aggregate browsers.</p>
  </div>
  <div class="fp-card">
    <h3>Add optional routes</h3>
    <p>Use pseudobulk and motif-discovery utilities when working with single-cell ATAC or candidate motifs.</p>
  </div>
</div>

## Core Commands

| Command | Purpose |
| --- | --- |
| `atac-correct` | Correct ATAC-seq cut-site signal for Tn5 sequence bias. |
| `call-footprints` | Calculate footprint scores from one or more bigWigs and optionally write candidate BED intervals. |
| `match-motifs` | Scan motifs for one or more footprint tracks and infer bound/unbound motif sites. |
| `diff-footprints` | Compare motif-associated footprint scores across conditions or replicates. |
| `normalize-bigwig` | Background-match corrected or footprint-score bigWigs before aggregate visualization. |
| `plot-aggregate` | Plot aggregate signal around motif sites as PDF/SVG-style output or HTML. |
| `plot-aggregate-batch` | Compatibility command for interactive aggregate HTML reports. |
| `run-workflow` | Run optional YAML batch configurations. |
| `fp-tools-gui` | Launch the optional Streamlit GUI. |

## Optional Utilities

| Command | Purpose |
| --- | --- |
| `motif-discovery` | Prepare candidate-centered de novo motif-discovery runs. |
| `motif-summary` | Summarize MEME/Tomtom outputs into TSV and HTML reports. |
| `fp-tools-score-variants` | Score variants with allele checks, candidate overlaps, sequence deltas, and optional motif/model deltas. |
| `pseudobulk-fragments` | Group single-cell ATAC fragments into pseudobulk fragment files. |
| `find-signature-fp` | Plot per-cell footprint-signature heatmaps and UMAP reports. |
| `pseudobulk-footprints` | Run the full pseudobulk fragment, correction, footprint, report, aggregate, and optional signature-reporting workflow. |

## Differential Footprint Defaults

For multi-sample analyses, score footprints from q95-scaled corrected bigWigs and keep `--normalization none` unless you are running an explicit sensitivity check.

## Standard Bulk Workflow

```text
atac-correct -> call-footprints -> match-motifs / motif-discovery -> diff-footprints -> plot-aggregate
```

`match-motifs` is the motif-site review step for one or more footprint tracks. `motif-discovery` is the optional de novo route from candidate footprint intervals. `diff-footprints` performs multi-condition comparisons and writes the differential report.

## Pseudobulk Workflow

```text
pseudobulk-fragments
  -> pseudobulk fragments and pseudo-BAMs
  -> atac-correct / call-footprints / match-motifs / diff-footprints
  -> find-signature-fp
  -> marker signature heatmaps and UMAP reports
```

Use `pseudobulk-footprints` when you want this route in one wrapper command. Provide `--motif-db` for motif-aware reports. Provide `--tf-site-dir` and `--single-cell-signature-h5ad` to add marker footprint-signature heatmaps and UMAP outputs.

## De Novo Motif Workflow

```text
call-footprints --output-bed
  -> motif-discovery
  -> motif-summary
  -> diff-footprints with --motifs alone, or --motif-db plus --motifs
```

Use de novo-only mode to test discovered motifs by themselves. Use database-plus-de-novo mode when discovered motifs should supplement a standard JASPAR or HOCOMOCO scan.
