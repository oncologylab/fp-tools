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
  <a href="api/#atac-correct">atac-correct</a> -> <a href="api/#call-footprints">call-footprints</a> -> <a href="api/#match-motifs">match-motifs</a> / <a href="api/#motif-discovery">motif-discovery</a> -> <a href="api/#diff-footprints">diff-footprints</a> -> <a href="api/#plot-aggregate">plot-aggregate</a>
</div>

### 1. Correct ATAC Signal

```bash
atac-correct \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project
```

Output: corrected cut-site bigWigs and QC files. With a sample table and `--outdir project`, fp-tools uses the project layout by default, reads a generic `sample	condition	bam	peaks` table, merges the sample peak BED files, writes `project/peaks/merged_peaks.bed`, and writes analysis-ready peaks to `project/peaks/merged_peaks.analysis.bed`.

Optional q95 scaling for multi-sample projects:

```bash
normalize-bigwig \
  --sample-table project/metadata/samples.tsv \
  --background project/peaks/merged_peaks.analysis.bed \
  --outdir project \
  --method background-scale \
  --stat q95 \
  --target median
```

### 2. Call Footprints

```bash
call-footprints \
  --sample-table project/metadata/samples.tsv \
  --regions project/peaks/merged_peaks.analysis.bed \
  --outdir project
```

Output: footprint score bigWig files and candidate footprint BEDs for de novo motif discovery. In project mode, fp-tools uses q95-scaled corrected tracks when present and otherwise falls back to corrected tracks.

### 3a. Match Motifs

```bash
match-motifs \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks.analysis.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

Output: for each sample, the summary table contains motif binding statistics. Each motif folder contains BED files such as `<motif>_all.bed`, `<motif>_<sample>_bound.bed`, and `<motif>_<sample>_unbound.bed`.

### 3b. Discover De Novo Motifs

Use de novo motif discovery when you want to start from candidate footprint intervals rather than only a curated database.

```bash
call-footprints \
  --signals project/samples/sample/atac_correct/sample_corrected.bw \
  --sample-names sample \
  --regions merged_peaks.bed \
  --sample-output-root project/samples

motif-discovery \
  --candidates project/samples/sample/footprints/sample_candidate_footprints.bed \
  --genome hg38.fa.gz \
  --flank 75 \
  --method streme \
  --known-motif-db jaspar2026_vertebrates \
  --outdir project/de_novo/sample
```

Output: candidate FASTA files, a runnable MEME/STREME/DREME command script, motif-discovery outputs, and optional Tomtom comparison against the selected known motif database. Discovered motifs can be used alone with `--motifs`, or added to a database run with `--motif-db` plus `--motifs`.

### 4. Compare Conditions

```bash
diff-footprints \
  --sample-table project/metadata/samples.tsv \
  --comparison-table project/metadata/comparisons.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks.analysis.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

`metadata/comparisons.tsv` uses generic columns: `comparison`, `cond1`, and `cond2`. Project-mode `diff-footprints` reuses cached motif-site tables and background scores from prior `match-motifs` runs, and samples with the same condition are treated as replicates. Output includes motif tables, BED files, volcano-style results, and [a standalone HTML report](demos/reports/diff_footprints_K562_HepG2.html).

### 5. Plot Aggregate

```bash
plot-aggregate \
  --sample-table project/metadata/samples.tsv \
  --motifs SPIB CEBPB \
  --site-set bound \
  --outdir project
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
  --peaks merged_peaks.bed \
  --motif-db jaspar2026_vertebrates \
  --tf-site-dir marker_motif_sites \
  --single-cell-signature-h5ad pbmc_embedding.h5ad \
  --outdir project/pseudobulk
```

Standard outputs include pseudobulk fragments, pseudo-BAMs, corrected cut-site bigWigs, footprint-score bigWigs, motif-aware differential reports, aggregate plots, and optional single-cell footprint-signature heatmaps/UMAPs. The combined signature plot is written as `plots/single_cell_footprinting/single_cell_footprinting.svg` when `--single-cell-signature-h5ad` and `--tf-site-dir` are supplied.

For a lighter first step, use [pseudobulk-fragments](api.md#pseudobulk-fragments) to only group fragments and write cut-site tracks. After motif-aware analysis, use [find-signature-fp](api.md#find-signature-fp) as a standalone reporting step for marker signature heatmaps and UMAPs.

```bash
find-signature-fp \
  --annotations cell_annotations.tsv \
  --fragments pbmc_fragments.tsv.gz \
  --h5ad pbmc_embedding.h5ad \
  --tf-site-dir marker_motif_sites \
  --all-motif-results project/pseudobulk/pseudobulk_diff_footprints_results.txt \
  --all-motif-diff-dir project/pseudobulk/diff_footprints \
  --outdir project/pseudobulk/signature_fp
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
