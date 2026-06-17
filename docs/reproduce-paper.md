# Reproducing the fp-tools manuscript

This repository supports two reproducibility paths:

- **Smoke path:** runs on committed fixtures and checks that the package, paper
  scripts, and LaTeX manuscript build correctly.
- **Full public-data path:** reruns the Buenrostro bulk ATAC, ENCODE K562/HepG2
  de novo motif validation, 10x PBMC pseudobulk, and benchmark scaffolds from downloaded
  public data under `data/public/`.

Large public downloads and generated benchmark outputs are not stored in git.
The committed source of truth is the code, manifests, command scripts, and small
source TSVs beside manuscript figures.

## Environment

Use either Conda/Mamba or Docker.

```bash
mamba env create -f environment.yml
mamba activate fp-tools
```

```bash
docker build -t fp-tools:paper .
docker run --rm -it -v "$PWD":/work -w /work fp-tools:paper
```

## Smoke checks

```bash
make test
make paper-smoke
```

The smoke path regenerates selected figure outputs and compiles
`manuscript/main.pdf`.

## Full public-data reruns

The full path expects public data under `data/public/` and may require substantial
CPU, memory, and disk space.

```bash
bash scripts/run_buenrostro_2x2_atac_replicate_demo.sh
bash scripts/run_encode_k562_hepg2_atac_demo.sh
STREME_NMOTIFS=250 \
THREADS=4 \
DIFF_PLOT_AGGREGATE=off \
OUT_DIR="$PWD/data/public/processed/encode_k562_hepg2_atac_replicates/fp_tools/denovo_motif_validation_maxcover_n250" \
  bash scripts/run_encode_k562_hepg2_denovo_motif_validation.sh
.venv/bin/python manuscript/scripts/plot_denovo_motif_validation.py \
  --validation-dir data/public/processed/encode_k562_hepg2_atac_replicates/fp_tools/denovo_motif_validation_maxcover_n250 \
  --out-prefix manuscript/figures/denovo_motif_validation
.venv/bin/python benchmarks/scripts/prepare_10x_pbmc_pseudobulk.py --write-example-archive
.venv/bin/python manuscript/scripts/prepare_pseudobulk_motif_sites.py \
  --peaks data/public/raw/10x_pbmc/pbmc_granulocyte_sorted_10k_atac_peaks.bed \
  --genome data/public/raw/genome/hg38.fa \
  --motifs data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
  --outdir data/public/processed/pseudobulk_pbmc/tf_sites_motif_centered \
  --summary data/public/processed/pseudobulk_pbmc/tf_sites_motif_centered/motif_centered_site_summary.tsv \
  --plot-sites-per-tf 0 \
  --site-selection all
.venv/bin/pseudobulk-footprints \
  --fragments data/public/raw/10x_pbmc/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz \
  --annotations data/public/processed/pseudobulk_pbmc/pbmc_10x_cell_annotations.tsv \
  --group-by cell_type \
  --min-cells 300 \
  --min-fragments 50000 \
  --genome-sizes data/public/processed/pseudobulk_pbmc/hg38.chrom.sizes \
  --genome data/public/raw/genome/hg38.fa \
  --peaks data/public/raw/10x_pbmc/pbmc_granulocyte_sorted_10k_atac_peaks.bed \
  --motifs data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
  --tf-site-dir data/public/processed/pseudobulk_pbmc/tf_sites_motif_centered \
  --site-summary data/public/processed/pseudobulk_pbmc/tf_sites_motif_centered/motif_centered_site_summary.tsv \
  --tfs auto \
  --outdir data/public/processed/pseudobulk_pbmc/footprints_full \
  --cores 8
.venv/bin/python manuscript/scripts/plot_pseudobulk_tf_aggregates.py \
  --manifest data/public/processed/pseudobulk_pbmc/footprints_full/pseudobulk_footprint_manifest.tsv \
  --tf-site-dir data/public/processed/pseudobulk_pbmc/tf_sites_motif_centered \
  --site-summary data/public/processed/pseudobulk_pbmc/tf_sites_motif_centered/motif_centered_site_summary.tsv \
  --out-prefix manuscript/figures/supp_pseudobulk_corrected_footprints \
  --signal-column footprint_bigwig \
  --value-column footprint_score \
  --ylabel "Footprint score" \
  --groups B_cell,CD4_T,CD14_Monocyte \
  --tfs SPIB,RUNX3,CEBPB \
  --flank 100
.venv/bin/python benchmarks/scripts/prepare_10x_pbmc5k_scatac.py --chroms chr1,chr2
.venv/bin/python manuscript/scripts/prepare_pseudobulk_motif_sites.py \
  --peaks data/public/raw/10x_pbmc5k_scatac/atac_pbmc_5k_snatac2_selected_bins.demo.bed \
  --genome data/public/raw/genome/hg38.fa \
  --motifs data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
  --outdir data/public/processed/pseudobulk_pbmc5k_scatac/tf_sites_motif_centered \
  --summary data/public/processed/pseudobulk_pbmc5k_scatac/tf_sites_motif_centered/motif_centered_site_summary.tsv \
  --candidates 'B_cell:PAX5,SPIB,POU2F2;T_NK:TCF7,ZBTB7B;Myeloid:CEBPB,CEBPA' \
  --chroms chr1,chr2 \
  --plot-sites-per-tf 1500 \
  --motif-pvalue 1e-4
.venv/bin/python benchmarks/scripts/plot_pbmc5k_per_cell_signatures.py \
  --annotations data/public/processed/pseudobulk_pbmc5k_scatac/pbmc5k_scprinter_broad_annotations.tsv \
  --fragments data/public/raw/10x_pbmc5k_scatac/atac_pbmc_5k_nextgem_fragments.tsv.gz \
  --h5ad data/public/raw/10x_pbmc5k_scatac/atac_pbmc_5k_annotated.h5ad \
  --tf-site-dir data/public/processed/pseudobulk_pbmc5k_scatac/tf_sites_motif_centered \
  --outdir data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/per_cell_signature_demo \
  --markers PAX5,CEBPB,TCF7,CEBPA,SPIB,ZBTB7B,POU2F2
.venv/bin/python benchmarks/scripts/plot_pbmc5k_pseudobulk_markers.py \
  --annotations data/public/processed/pseudobulk_pbmc5k_scatac/pbmc5k_scprinter_broad_annotations.tsv \
  --aggregate-screen data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/pseudobulk_footprint_aggregate_screen.tsv \
  --bindetect-results data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/bindetect/pseudobulk_bindetect_results.txt \
  --outdir data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/marker_demo \
  --markers PAX5,EBF1,POU2F2,POU2AF1,BCL6,SPIB,CEBPB,CEBPA,TCF7,LEF1,ZBTB7B,RUNX3,GATA3
cp data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/per_cell_signature_demo/pbmc5k_knn_footprint_signature_umap.pdf \
  manuscript/figures/pbmc5k_knn_footprint_signature_umap.pdf
cp data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/per_cell_signature_demo/pbmc5k_knn_footprint_signature_umap.png \
  manuscript/figures/pbmc5k_knn_footprint_signature_umap.png
cp data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/marker_demo/pbmc5k_volcano_pairwise_directional_markers.pdf \
  manuscript/figures/pbmc5k_volcano_pairwise_directional_markers.pdf
cp data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/marker_demo/pbmc5k_volcano_pairwise_directional_markers.png \
  manuscript/figures/pbmc5k_volcano_pairwise_directional_markers.png
cp data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/plots/marker_demo/pbmc5k_volcano_pairwise_directional_markers.tsv \
  manuscript/figures/pbmc5k_volcano_pairwise_directional_markers.tsv
```

After public-data outputs exist, regenerate the remaining manuscript figures with the scripts
under `manuscript/scripts/` and compile the paper:

```bash
make paper-pdf
```

## Benchmark manifests

Top-level TSV files in `benchmarks/manifests/` follow the full benchmark manifest
schema documented in `benchmarks/manifests/README.md`. Compact source manifests
used by helper scripts live under `benchmarks/manifests/compact/` and have their
own explicit schemas validated by `benchmarks/scripts/validate_manifests.py`.

## Engineering benchmarks

Use `benchmarks/scripts/run_engineering_benchmark.py` to record runtime and memory
metadata for future fp-tools, TOBIAS, HINT, or other external-tool comparisons:

```bash
python benchmarks/scripts/run_engineering_benchmark.py \
  --label fp-tools-demo \
  --cores 8 \
  --out benchmarks/results/engineering_runtime.tsv \
  -- diff-footprints --help
```
