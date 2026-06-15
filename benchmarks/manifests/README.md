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
