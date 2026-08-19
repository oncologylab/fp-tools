# Region-set comparison

This example compares four classes of accessible HepG2 regions defined by
HNF4A and FOXA2 ChIP-seq: HNF4A + FOXA2, HNF4A only, FOXA2 only, and no
HNF4A/FOXA2. The last group contains accessible regions overlapping neither
ChIP-seq peak set.

The classes were restricted to a shared baseline-accessibility range, divided
into 50 strata, and sampled equally within every stratum. Each class contains
4,786 regions. Three HepG2 ATAC-seq replicates are analyzed as paired
measurements of the same region classes.

## Example command

```bash
diff-footprints \
  --comparison-axis regions \
  --signals HepG2_rep1_footprints.bw HepG2_rep2_footprints.bw HepG2_rep3_footprints.bw \
  --sample-names "HepG2 rep 1" "HepG2 rep 2" "HepG2 rep 3" \
  --regions HNF4A_FOXA2.bed HNF4A_only.bed FOXA2_only.bed No_HNF4A_FOXA2.bed \
  --region-labels "HNF4A + FOXA2" "HNF4A only" "FOXA2 only" "No HNF4A/FOXA2" \
  --region-strata-column 4 \
  --genome hg38.fa \
  --motif-db jaspar2026_vertebrates \
  --aggregate-signals HepG2_rep1_corrected.bw HepG2_rep2_corrected.bw HepG2_rep3_corrected.bw \
  --plot-aggregate-motifs MA1494.2 MA0484.3 MA0047.4 MA0148.5 MA0046.3 MA0153.2 MA0102.5 MA0466.4 \
  --default-aggregate-plots 8 \
  --outdir HepG2_region_comparison
```

The primary effect is the matching-stratum-adjusted difference in region-level
footprint score. All 1,019 motifs are tested; the eight specified motifs only
set the initial aggregate display. Significance is calculated from paired
replicate effects using the empirical-Bayes model.

In the default view, HNF4A/HNF4G and FOXA1/FOXA2 footprints are stronger in
co-bound regions. The HNF4A-only versus FOXA2-only view separates the two
factor families.

[All motif results](../../demos/data/region_set_HepG2_HNF4A_FOXA2_results.tsv.gz) ·
[matching QC](../../demos/data/region_set_HepG2_matching_qc.tsv) ·
[matching summary](../../demos/data/region_set_HepG2_matching_summary.tsv) ·
[source manifest](../../demos/data/region_set_HepG2_source_manifest.tsv)

<div class="fp-live-demo-wrap fp-output-example">
  <iframe
    class="fp-live-demo fp-report-demo"
    src="../../../demos/reports/region_set_HepG2_HNF4A_FOXA2/"
    title="HepG2 HNF4A and FOXA2 region-set report"
    loading="eager">
  </iframe>
</div>
