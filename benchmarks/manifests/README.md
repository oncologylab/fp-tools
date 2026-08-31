# Public Data Manifests

Manifest files are tab-separated and versioned. They describe public datasets without storing bulky data in git.

Required columns:

| Column | Meaning |
|---|---|
| `source` | ENCODE, JASPAR, HOCOMOCO, caQTL, 10x, IGVF, etc. |
| `benchmark_tier` | smoke, bulk, depth, motif_removal, variant, pseudobulk |
| `cell_type` | Biosample or pseudobulk group label |
| `donor` | Donor identifier when available |
| `tf` | TF label for binding benchmarks |
| `assay` | ATAC-seq, TF ChIP-seq, CUT&RUN, scATAC, multiome, variant |
| `experiment_accession` | Public experiment accession |
| `file_accession` | Public file accession |
| `assembly` | Genome assembly, initially GRCh38 |
| `output_type` | BAM, bigWig, peaks BED, fragments, motif file, variant file |
| `file_format` | bam, bigWig, bed, tsv.gz, meme, jaspar, etc. |
| `url` | Download URL |
| `checksum` | MD5/SHA if available |
| `status` | released or equivalent public status |
| `local_path` | Planned path after download |
| `split` | train, validation, test, smoke, or not_applicable |
| `notes` | Filtering or curation notes |


## Manifest locations

Top-level `*.tsv` files in this directory use the full schema above and are
validated in CI. Compact source manifests used by helper scripts live in
`benchmarks/manifests/compact/` and use explicitly registered compact schemas.
Run validation with:

```bash
python benchmarks/scripts/validate_manifests.py --manifest-dir benchmarks/manifests
```

`footprint_functional_v1.spec.json` is the preregistered follow-up to the
original footprint-detectability study. It freezes K562/HepG2 development
chromosomes, the previously unscored GM12878/IMR-90 cell holdout, the +4/-5
versus +4/-4 shift comparison, parametric bias and functional-model grids, and
scientific/computational promotion gates. The old
`footprint_detectability_v1.spec.json` remains unchanged as the record of the
earlier study and its already-viewed holdout.

`compact/functional_holdout_atac_bams.tsv` pins the four coordinate-sorted
ENCODE ATAC-seq BAMs used to prepare the locked GM12878 and IMR-90 transfer
panel. File size and MD5 are checked before any holdout analysis. Downloading
these assay inputs does not open or score the preregistered TF ChIP labels.
The companion `compact/functional_holdout_atac_peaks.tsv` pins the pooled,
conservative IDR-thresholded ENCODE4 v1.9.1 peaks from the same two ATAC-seq
experiments. They define accessible motif candidates without using TF labels.

`functional_detector_policy_v1.freeze.json` records the policy freeze after the
complete DWM/SELMA/log-linear, shift, FDA, spline, GP, and hybrid development
factorial. The compact policy and holdout-route tables make the fail-closed
decision executable: only families passing every development gate use a new
route, and unseen or failed families retain DWM. The freeze fingerprints its
source metrics and explicitly precedes all locked GM12878/IMR-90 label scoring.

`functional_max_wtc11_validation_v1.freeze.json` preregisters a subsequent,
untouched WTC11/MAX transfer test after the GM12878 diagnostic exposed a
promising anchored LOG81 route. Its companion TSV pins two ENCODE ATAC-seq
BAMs, the label-free ATAC peak universe, and the independent MAX ChIP-seq peak
label. The candidate, DWM comparator, chromosome split, matching covariates,
and pass thresholds are fixed before the ChIP peak content is downloaded.
