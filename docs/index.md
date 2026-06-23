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
  atac-correct -> call-footprints -> match-motifs -> diff-footprints -> plot-aggregate-batch
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

### 3. Match Motifs

```bash
match-motifs \
  --signals results/footprints/sample_footprints.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir results/motif_matches/sample
```

Output: motif-site tables, bound/unbound motif calls, and files that can be reviewed before multi-condition comparisons.

### 4. Compare Conditions

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

### 5. Review Aggregate Plots

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

Use `pseudobulk-fragments` to group single-cell fragments by annotation, then run the standard footprint commands on the grouped pseudobulk samples. Use `pseudobulk-footprints` when you want grouping, correction, scoring, motif reports, aggregate plots, and optional signature reporting in one wrapper command.

```bash
pseudobulk-footprints \
  --fragments pbmc_fragments.tsv.gz \
  --annotations cell_annotations.tsv \
  --group-by cell_type \
  --genome-sizes hg38.chrom.sizes \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --tf-site-dir marker_motif_sites \
  --single-cell-signature-h5ad pbmc_embedding.h5ad \
  --outdir results/pseudobulk
```

Standard outputs include pseudobulk fragments, pseudo-BAMs, corrected cut-site bigWigs, footprint-score bigWigs, motif-aware differential reports, aggregate plots, and optional single-cell footprint-signature heatmaps/UMAPs. The combined signature plot is written as `plots/single_cell_footprinting/single_cell_footprinting.svg` when `--single-cell-signature-h5ad` and `--tf-site-dir` are supplied.

For a lighter first step, use `pseudobulk-fragments` to only group fragments and write cut-site tracks. After motif-aware analysis, use `find-signature-fp` as a standalone reporting step for marker signature heatmaps and UMAPs.

```bash
find-signature-fp \
  --annotations cell_annotations.tsv \
  --fragments pbmc_fragments.tsv.gz \
  --h5ad pbmc_embedding.h5ad \
  --tf-site-dir marker_motif_sites \
  --all-motif-results results/pseudobulk/pseudobulk_diff_footprints_results.txt \
  --all-motif-diff-dir results/pseudobulk/diff_footprints \
  --outdir results/pseudobulk/signature_fp
```

## De Novo Motif Discovery

Use de novo motif discovery when you want to start from candidate footprint intervals rather than only a curated database.

```bash
call-footprints \
  --signal results/atac_correct/sample/sample_corrected.bw \
  --regions peaks.bed \
  --output results/footprints/sample_footprints.bw \
  --output-bed results/footprints/sample_candidates.bed

motif-discovery \
  --candidates results/footprints/sample_candidates.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir results/de_novo/sample
```

This writes candidate FASTA, a runnable MEME/STREME/DREME shell script, and optional Tomtom comparison against the selected known motif database.

Use the discovered motifs in either mode:

```bash
# De novo-only
diff-footprints ... --motifs results/de_novo/sample/streme/streme.txt

# Database plus de novo supplement
diff-footprints ... \
  --motif-db jaspar2026_vertebrates \
  --motifs results/de_novo/sample/streme/streme.txt
```

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
