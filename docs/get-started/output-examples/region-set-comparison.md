# Region-set comparison

This example compares K562 accessible regions overlapping CTCF ChIP-seq peaks
with accessibility-matched accessible regions that do not overlap CTCF peaks.
The ATAC-seq data are ENCODE experiment
[ENCSR483RKN](https://www.encodeproject.org/experiments/ENCSR483RKN/)
(replicates ENCFF512VEZ and ENCFF987XOV; conservative peaks ENCFF695IGF), and
the CTCF ChIP-seq peaks are ENCFF362OPG from
[ENCSR000AKO](https://www.encodeproject.org/experiments/ENCSR000AKO/).

The two sets were restricted to their shared baseline-accessibility range,
divided into 50 baseline-signal strata, and sampled equally within every
stratum. Each enhancer contributes one mean motif score even when it contains
multiple instances of the same motif.

## Example command

```bash
diff-footprints \
  --comparison-axis regions \
  --signals K562_rep1_footprints.bw K562_rep2_footprints.bw \
  --sample-names K562_rep1 K562_rep2 \
  --cond-names K562 K562 \
  --regions CTCF_bound.bed matched_control.bed \
  --region-labels CTCF_bound matched_control \
  --region-strata-column 4 \
  --genome hg38.fa \
  --motif-db jaspar2026_vertebrates \
  --aggregate-signals K562_rep1_corrected.bw K562_rep2_corrected.bw \
  --outdir K562_CTCF_region_comparison
```

The primary effect is the matching-stratum-adjusted difference in enhancer-level
footprint score. With two K562 replicates, significance is calculated from the
paired replicate effects using the empirical-Bayes model.

[Motif-level results](../../demos/data/region_set_K562_CTCF_results.tsv) ·
[matching QC](../../demos/data/region_set_K562_CTCF_matching_qc.tsv) ·
[CTCF-bound BED](../../demos/data/region_set_K562_CTCF_bound.bed) ·
[matched-control BED](../../demos/data/region_set_K562_matched_control.bed)

<div class="fp-live-demo-wrap fp-output-example">
  <iframe
    class="fp-live-demo fp-report-demo"
    src="../../../demos/reports/region_set_K562_CTCF.html"
    title="K562 CTCF-bound versus accessibility-matched region report"
    loading="eager">
  </iframe>
</div>
