<section class="fp-hero">
  <div class="fp-hero-copy">
    <p class="fp-eyebrow">Reproducible regulatory genomics</p>
    <h1>ATAC-seq footprints, motifs, and reports from one command-first toolkit</h1>
    <p>Turn raw reads or processed tracks into bias-corrected signal,
    footprint scores, motif-aware comparisons, aggregate profiles, and
    portable HTML reports.</p>
    <div class="fp-actions">
      <a class="fp-button primary" href="#install">Install fp-tools</a>
      <a class="fp-button" href="reports/">Explore reports</a>
      <a class="fp-button" href="api/">Command reference</a>
    </div>
  </div>
  <img src="assets/fp_tools_logo_icon.svg" alt="" role="presentation">
</section>

`fp-tools` runs ATAC-seq footprint workflows from the command line or the
optional browser GUI. The same YAML configuration can be saved from the GUI
and rerun with [run-workflow](api.md#run-workflow).

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
  <a href="api/#prepare-atac">prepare-atac</a> -> <a href="api/#atac-correct">atac-correct</a> -> <a href="api/#call-footprints">call-footprints</a> -> <a href="api/#match-motifs">match-motifs</a> / <a href="api/#motif-discovery">motif-discovery</a> -> <a href="api/#diff-footprints">diff-footprints</a> -> <a href="api/#plot-aggregate">plot-aggregate</a>
</div>

### 0. Prepare Raw Reads

```bash
prepare-atac \
  --samples metadata.tsv \
  --genome hg38 \
  --outdir project/raw
```

The table needs an `ID` or `run_accession` column for SRR, ERR, or DRR data. It
can instead provide `fastq_1` and `fastq_2` paths or HTTPS URLs. fp-tools first
looks for checksum-verified compressed FASTQs at ENA and can fall back to the
NCBI SRA Toolkit.

The default pipeline trims adapters with fastp, aligns reads with Bowtie2,
keeps confidently mapped alignments, removes PCR duplicates, mitochondrial
reads, and blacklist regions, and calls peaks with MACS3. The `legacy-atac`
profile uses Trim Galore, Bowtie2 local alignment, Picard duplicate marking,
and HOMER coverage and peak calling. It accepts paired-end libraries and older
single-end ATAC libraries using the corresponding single-read command forms.
Both profiles write filtered
BAM/BAI files, peak BED files, RP10M bigWigs, QC summaries, command logs, merged
peaks, and a downstream `metadata/samples.tsv`.

Named `hg38` and `mm10` references are cached; custom FASTA, Bowtie2 index, and
blacklist files are also supported. Run `prepare-atac --doctor` before a
production job. Thread counts scale to `--cores`, and `--memory-gb` controls the
memory-aware sample scheduler and sort settings. Use `--write-default-config`
to create an editable YAML file with every processing option.

### 1. Correct ATAC Signal

```bash
atac-correct \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --blacklist hg38.blacklist.bed \
  --outdir project
```

Output: corrected cut-site bigWigs and QC files. With a sample table and `--outdir project`, fp-tools uses the project layout by default, reads a generic `sample	condition	bam	peaks` table, merges the sample peak BED files, writes `project/peaks/merged_peaks.bed`, and writes filtered peaks to `project/peaks/merged_peaks_filtered.bed` with mitochondrial chromosomes excluded. Downstream project-mode commands use this filtered BED when peak/background regions are omitted or when the project raw merged BED is passed.
For production runs where corrected bigWigs are the required output, add `--skip-qc` to skip the diagnostic PDF and pre/post correction bias-count summaries without changing the corrected signal tracks.

Optional q95 scaling for multi-sample projects:

```bash
normalize-bigwig \
  --sample-table project/metadata/samples.tsv \
  --background project/peaks/merged_peaks_filtered.bed \
  --outdir project \
  --method background-scale \
  --stat q95 \
  --target median
```

### 2. Call Footprints

```bash
call-footprints \
  --sample-table project/metadata/samples.tsv \
  --regions project/peaks/merged_peaks_filtered.bed \
  --outdir project
```

Output: footprint score bigWig files. Add `--call-candidates` in project mode, or `--output-bed`/`--output-beds` in custom mode, to also write genomic coordinates for footprint peaks used by de novo motif discovery. In project mode, fp-tools uses q95-scaled corrected tracks when present and otherwise falls back to corrected tracks. The default footprint kernel uses the faster prefix-sum implementation; use `--footprint-kernel legacy` only when exact historical floating-point behavior is required.

### 3a. Match Motifs

```bash
match-motifs \
  --sample-table project/metadata/samples.tsv \
  --genome hg38.fa.gz \
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

Output: for each sample, the summary table contains motif binding statistics, motif-site/background caches are written for fast reuse by `diff-footprints`, and per-motif BED folders are materialized by default with files such as `<motif>_all.bed`, `<motif>_<sample>_bound.bed`, and `<motif>_<sample>_unbound.bed`. In project mode, `match-motifs` uses one shared motif scan across samples by default and then writes normal per-sample folders; per-motif BED folders are written in the background after report-ready outputs. Cached `diff-footprints` comparisons may create internal per-motif shard caches on first reuse, then reuse those shards for later comparisons. Use `--match-scan-mode per-sample` only when independent sample scans are needed for debugging. Use `--motif-outputs summary` only when you want to skip permanent BED folders.

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
  --peaks project/peaks/merged_peaks_filtered.bed \
  --motif-db jaspar2026_vertebrates \
  --outdir project
```

`metadata/comparisons.tsv` uses generic columns: `comparison`, `cond1`, and `cond2`. Project-mode `diff-footprints` reuses compact motif-site tables and background scores from prior `match-motifs` runs, and samples with the same condition are treated as biological replicates. With at least two replicates per condition, the command writes a per-sample motif-score matrix and fits empirical-Bayes moderated condition contrasts using samples—not motif sites—as the inferential units. Output includes motif tables, volcano-style results, replicate diagnostics, and [a standalone HTML report](demos/reports/diff_footprints_K562_HepG2.html). Aggregate profiles in the report are capped by `--plot-aggregate-top-n` in `sig` and `top` modes; increase this value when you want more motif profiles in the HTML.

### Workflow Improvements

fp-tools preserves the interpretable TOBIAS-style center-versus-flank footprint score while improving the multi-sample workflow. q95 scaling can align corrected cut-site tracks over shared reference regions before scoring. The default footprint-scoring path uses optimized Cython-backed kernels, and project-mode motif analysis uses one shared motif scan plus compact caches that `diff-footprints` reuses for biological-replicate empirical-Bayes comparisons and HTML/SVG reports.

For a complete multi-condition shell-script template, see
`examples/nutrient_stress_project/run_ctrl_vs_10fbs.sh` in the GitHub
repository. It documents the expected raw-data layout, explains how to prepare
portable `samples.tsv` and `comparisons.tsv` files from an
`ATAC_Nutrients_hg38_*.txt` table, supports `CHECK_ONLY=1` input validation, and
then runs the full project workflow.

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
