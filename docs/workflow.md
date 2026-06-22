# Workflow

`fp-tools` keeps the command line as the primary interface. The optional GUI and YAML runner orchestrate the same command surface.

## Core Bulk ATAC Workflow

```bash
atac-correct \
  --bam sample.bam \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --blacklist hg38-blacklist.bed \
  --outdir results/atacorrect/sample
```

```bash
call-footprints \
  --signal results/atacorrect/sample/sample_corrected.bw \
  --regions peaks.bed \
  --output results/footprints/sample_footprints.bw \
  --output-bed results/footprints/sample_candidate_footprints.bed
```

```bash
diff-footprints \
  --signals B_rep1_footprints.bw B_rep2_footprints.bw T_rep1_footprints.bw T_rep2_footprints.bw \
  --aggregate-signals B_rep1_corrected.bw B_rep2_corrected.bw T_rep1_corrected.bw T_rep2_corrected.bw \
  --genome hg38.fa.gz \
  --peaks peaks.bed \
  --cond-names Bcell Bcell Tcell Tcell \
  --normalization none \
  --plot-aggregate sig \
  --outdir results/diff_footprints/Bcell_vs_Tcell
```

Repeated condition names define biological replicates. By default, motif-aware commands use the bundled `jaspar2026_vertebrates` database; use `--motif-db hocomoco14_core` or `--motifs motifs.jaspar` to change inputs. When aggregate signals are supplied, `diff-footprints` writes a standalone interactive HTML report with volcano-style differential evidence, motif logos, and aggregate profiles.

## Aggregate Visualization

Use `normalize-bigwig` to make corrected cut-site tracks comparable before aggregate plotting:

```bash
normalize-bigwig \
  --bigwigs B_rep1_corrected.bw B_rep2_corrected.bw T_rep1_corrected.bw T_rep2_corrected.bw \
  --background merged_peaks.50bp_bins.bed \
  --method background-scale \
  --stat q95 \
  --target median \
  --outdir results/normalized_corrected_bigwigs
```

Then render static or interactive aggregate reports:

```bash
plot-aggregate-batch \
  --input-html results/diff_footprints/Bcell_vs_Tcell/diff_footprints_Bcell_Tcell.html \
  --output results/reports/aggregate_browser.html
```

## Pseudobulk Route

The pseudobulk commands group single-cell fragments, run footprint scoring, and generate marker-focused differential reports:

```bash
pseudobulk-fragments --help
pseudobulk-footprints --help
```

See the paper reproduction page for public-data workflows.

## Variant Scoring

Use `fp-tools-score-variants` when candidate footprints or motif databases should be summarized at variant alleles:

```bash
fp-tools-score-variants \
  --variants variants.bed \
  --genome hg38.fa.gz \
  --candidate-scores candidate_footprints.bed \
  --motif-db jaspar2026_vertebrates \
  --out variant_scores.tsv
```
