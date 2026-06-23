# Get Started

`fp-tools` runs ATAC-seq footprint workflows from the command line or the optional browser GUI. The same YAML configuration can be saved from the GUI and rerun with [run-workflow](api.md#run-workflow).

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
  <a href="api/#atac-correct">atac-correct</a> -> <a href="api/#call-footprints">call-footprints</a> -> 3a. <a href="api/#match-motifs">match-motifs</a> / 3b. <a href="api/#motif-discovery">motif-discovery</a> -> <a href="api/#diff-footprints">diff-footprints</a> -> <a href="api/#plot-aggregate">plot-aggregate</a>
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
  --signals results/atac_correct/sample/sample_corrected.bw \
  --regions peaks.bed \
  --outdir results/footprints
```

Output: footprint score bigWig files. Use `--signals` with multiple corrected bigWigs to score several samples in one run. Add `--output-bed` or `--output-bed-dir` to write genomic coordinates for footprint peaks needed by de novo motif discovery.

### 3a. Match Motifs

```bash
match-motifs \
  --signals results/footprints/sample_footprints.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir results/motif_matches/sample
```

For multiple footprint tracks, pass all tracks to `--signals` and provide one sample name per track:

```bash
match-motifs \
  --signals \
    results/footprints/A_footprints.bw \
    results/footprints/B_footprints.bw \
  --sample-names A B \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir results/motif_matches/A_B
```

Output: for each sample, the summary table contains motif binding statistics. Each motif folder contains BED files such as `<motif>_all.bed`, `<motif>_<sample>_bound.bed`, and `<motif>_<sample>_unbound.bed`.

### 3b. Discover De Novo Motifs

Use de novo motif discovery when you want to start from candidate footprint intervals rather than only a curated database.

```bash
call-footprints \
  --signals results/atac_correct/sample/sample_corrected.bw \
  --regions peaks.bed \
  --outdir results/footprints \
  --output-bed-dir results/footprints/candidates

motif-discovery \
  --candidates results/footprints/candidates/sample_candidates.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir results/de_novo/sample
```

Output: candidate FASTA files, a runnable MEME/STREME/DREME command script, motif-discovery outputs, and optional Tomtom comparison against the selected known motif database. Discovered motifs can be used alone with `--motifs`, or added to a database run with `--motif-db` plus `--motifs`.

### 4. Compare Conditions

```bash
diff-footprints \
  --signals \
    A_rep1_footprints.bw A_rep2_footprints.bw \
    B_rep1_footprints.bw B_rep2_footprints.bw \
  --sample-names A_R1 A_R2 B_R1 B_R2 \
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

Repeated names in `--cond-names` define biological replicates. Output includes motif tables, BED files, volcano-style results, and [a standalone HTML report](demos/reports/diff_footprints_K562_HepG2.html).

### 5. Plot Aggregate

```bash
plot-aggregate \
  --match-dir results/motif_matches/sample \
  --signals results/atac_correct/sample/sample_corrected.bw \
  --motifs SPIB CEBPB \
  --site-set bound \
  --format html \
  --output results/reports/aggregate_browser.html
```

Output: static PDF/SVG-style aggregate plots or an interactive HTML subplot browser, depending on `--format` or the `--output` extension.

## Pseudobulk Single-Cell ATAC

Use [pseudobulk-fragments](api.md#pseudobulk-fragments) to group single-cell fragments by annotation, then run the standard footprint commands on the grouped pseudobulk samples. Use [pseudobulk-footprints](api.md#pseudobulk-footprints) when you want grouping, correction, scoring, motif reports, aggregate plots, and optional signature reporting in one wrapper command.

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

For a lighter first step, use [pseudobulk-fragments](api.md#pseudobulk-fragments) to only group fragments and write cut-site tracks. After motif-aware analysis, use [find-signature-fp](api.md#find-signature-fp) as a standalone reporting step for marker signature heatmaps and UMAPs.

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

## Where To Go Next

<div class="fp-grid">
  <div class="fp-card">
    <h3>Command Manuals</h3>
    <p>See the command overview and full help for every supported command.</p>
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
