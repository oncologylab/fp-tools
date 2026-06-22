# Commands

Install the package, then inspect any command with `--help`.

```bash
pip install fp-tools-bio
```

<div class="fp-grid">
  <div class="fp-card">
    <h3>Start with bulk ATAC</h3>
    <p>atac-correct, call-footprints, and diff-footprints cover the core footprint workflow.</p>
  </div>
  <div class="fp-card">
    <h3>Review reports</h3>
    <p>plot-aggregate and plot-aggregate-batch generate static figures and standalone HTML browsers.</p>
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
| `call-footprints` | Calculate footprint scores and optionally write ranked candidate BED intervals. |
| `match-motifs` | Scan motifs for one sample and infer bound/unbound motif sites. |
| `diff-footprints` | Compare motif-associated footprint scores across conditions or replicates. |
| `normalize-bigwig` | Background-match corrected or footprint-score bigWigs before aggregate visualization. |
| `plot-aggregate` | Plot static aggregate signal around TFBS or region sets. |
| `plot-aggregate-batch` | Create an interactive multi-sample, multi-TF aggregate HTML report. |
| `run-workflow` | Run optional YAML batch configurations. |
| `fp-tools-gui` | Launch the optional Streamlit GUI. |

## Optional Utilities

| Command | Purpose |
| --- | --- |
| `motif-discovery` | Prepare candidate-centered de novo motif-discovery runs. |
| `motif-summary` | Summarize MEME/Tomtom outputs into TSV and HTML reports. |
| `fp-tools-score-variants` | Score variants with allele checks, candidate overlaps, sequence deltas, and optional motif/model deltas. |
| `pseudobulk-fragments` | Group single-cell ATAC fragments into pseudobulk fragment files. |
| `pseudobulk-footprints` | Run the full pseudobulk fragment, correction, footprint, and report workflow. |

## Compatibility Aliases

Legacy aliases remain available:

```bash
ATACorrect --help
FootprintScores --help
ScoreBigwig --help
BINDetect --help
PlotAggregate --help
```

## Differential Footprint Defaults

`diff-footprints` uses the native BINDetect-style motif-aware comparison. For multi-sample analyses, score footprints from q95-scaled corrected bigWigs and keep `--normalization none` unless you are running an explicit sensitivity check.
