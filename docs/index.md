# Get Started

`fp-tools` runs ATAC-seq footprint workflows from the command line or the optional browser GUI. The same YAML configuration can be saved from the GUI and rerun with `run-workflow`.

## Install

```bash
pip install fp-tools-bio
```

For the browser GUI:

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

## Standard Workflow

<div class="fp-command-chain">
  atac-correct -> call-footprints -> diff-footprints -> plot-aggregate-batch
</div>

### 1. Correct ATAC Signal

```bash
atac-correct \
  --bam sample.bam \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --blacklist hg38.blacklist.bed \
  --outdir results/atac_correct/sample
```

Output: corrected cut-site bigWigs and QC files.

### 2. Call Footprints

```bash
call-footprints \
  --signal results/atac_correct/sample/sample_corrected.bw \
  --regions peaks.bed \
  --output results/footprints/sample_footprints.bw
```

Output: a footprint score bigWig. Add `--output-bed` to write ranked candidate intervals.

### 3. Compare Conditions

```bash
diff-footprints \
  --signals \
    A_rep1_footprints.bw A_rep2_footprints.bw \
    B_rep1_footprints.bw B_rep2_footprints.bw \
  --aggregate-signals \
    A_rep1_corrected.bw A_rep2_corrected.bw \
    B_rep1_corrected.bw B_rep2_corrected.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --cond-names A A B B \
  --motif-db jaspar2026_vertebrates \
  --normalization none \
  --plot-aggregate sig \
  --outdir results/diff_footprints/A_vs_B
```

Repeated names in `--cond-names` define biological replicates. Output includes motif tables, BED files, volcano-style results, and a standalone HTML report.

### 4. Review Aggregate Plots

```bash
plot-aggregate-batch \
  --input-html results/diff_footprints/A_vs_B/diff_footprints_A_B.html \
  --output results/reports/aggregate_browser.html
```

Output: an interactive HTML browser for motif-centered aggregate profiles.

## Single-Sample Motif Matching

Use `match-motifs` when you have one footprint score track and want motif-site tables and bound/unbound motif calls.

```bash
match-motifs \
  --signals sample_footprints.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir results/motif_matches/sample
```

## Pseudobulk Single-Cell ATAC

Use these commands when starting from single-cell fragments and cell annotations:

```bash
pseudobulk-fragments --help
pseudobulk-footprints --help
```

`pseudobulk-footprints` groups fragments, runs ATAC correction, scores footprints, and can produce motif-aware reports.

## Where To Go Next

<div class="fp-grid">
  <div class="fp-card">
    <h3>Command List</h3>
    <p>See all primary commands and what each one does.</p>
    <p><a href="commands/">Open Commands</a></p>
  </div>
  <div class="fp-card">
    <h3>Command Manuals</h3>
    <p>Read full command help for options, inputs, and outputs.</p>
    <p><a href="api/">Open API Reference</a></p>
  </div>
  <div class="fp-card">
    <h3>Report Demo</h3>
    <p>Open a static differential-footprint report in the browser.</p>
    <p><a href="reports/">Open Reports</a></p>
  </div>
  <div class="fp-card">
    <h3>GUI Demo</h3>
    <p>Preview the browser interface and tutorial layout.</p>
    <p><a href="gui/">Open GUI Demo</a></p>
  </div>
</div>
