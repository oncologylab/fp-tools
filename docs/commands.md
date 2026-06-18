# Commands

Install the package, then inspect any command with `--help`.

```bash
pip install fp-tools-bio
```

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

## Fixed Motif-Site Statistical Backends

`diff-footprints` can reuse an existing motif-site reference:

- `--method deseq2-cutcount`: raw shifted Tn5 insertion counts over fixed motif-site windows, analyzed with PyDESeq2.
- `--method footprint-score`: continuous footprint-score signal over fixed motif-site windows, analyzed with an empirical-Bayes moderated test.

Install PyDESeq2 support with:

```bash
pip install "fp-tools-bio[deseq2]"
```
