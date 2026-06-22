# Workflow

`fp-tools` can be used from the command line or through the optional browser GUI. The same YAML configs work in both places.

## 1. Correct ATAC Signal

```bash
atac-correct \
  --bam sample.bam \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --blacklist hg38.blacklist.bed \
  --outdir results/atac_correct/sample
```

Output: corrected cut-site bigWigs and QC files.

## 2. Call Footprints

```bash
call-footprints \
  --signal results/atac_correct/sample/sample_corrected.bw \
  --regions peaks.bed \
  --output results/footprints/sample_footprints.bw
```

Output: a footprint score bigWig. Add `--output-bed` when you also want ranked candidate intervals.

## 3. Compare Conditions

```bash
diff-footprints \
  --signals A_rep1_footprints.bw A_rep2_footprints.bw B_rep1_footprints.bw B_rep2_footprints.bw \
  --aggregate-signals A_rep1_corrected.bw A_rep2_corrected.bw B_rep1_corrected.bw B_rep2_corrected.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --cond-names A A B B \
  --motif-db jaspar2026_vertebrates \
  --normalization none \
  --plot-aggregate sig \
  --outdir results/diff_footprints/A_vs_B
```

Repeated condition names define replicates. Output includes motif tables, BED files, volcano-style results, and a standalone HTML report.

## 4. Review Aggregate Plots

```bash
plot-aggregate-batch \
  --input-html results/diff_footprints/A_vs_B/diff_footprints_A_B.html \
  --output results/reports/aggregate_browser.html
```

Output: an interactive HTML browser for motif-centered aggregate profiles.

## Pseudobulk Single-Cell ATAC

Use these commands when starting from single-cell fragments and cell annotations:

```bash
pseudobulk-fragments --help
pseudobulk-footprints --help
```

`pseudobulk-footprints` runs grouping, ATAC correction, footprint scoring, and optional motif-aware reports.

## Browser GUI

```bash
pip install "fp-tools-bio[gui]"
fp-tools-gui
```

Use the sidebar to pick a command, load an example YAML config, preview the command-ready config, and launch the run.
