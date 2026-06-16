#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${ROOT_DIR}"

FP_TOOLS_ENV=${FP_TOOLS_ENV:-${ROOT_DIR}/.venv}
MODE=${1:-full}

FRAGMENTS=${FRAGMENTS:-data/public/raw/10x_pbmc/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz}
ANNOTATIONS=${ANNOTATIONS:-data/public/processed/pseudobulk_pbmc/pbmc_10x_cell_annotations.tsv}
GENOME_SIZES=${GENOME_SIZES:-data/public/processed/pseudobulk_pbmc/hg38.chrom.sizes}
GENOME=${GENOME:-data/public/raw/genome/hg38.fa}
PEAKS=${PEAKS:-data/public/raw/10x_pbmc/pbmc_granulocyte_sorted_10k_atac_peaks.bed}
TF_SITE_DIR=${TF_SITE_DIR:-data/public/processed/pseudobulk_pbmc/tf_sites_motif_centered}
SITE_SUMMARY=${SITE_SUMMARY:-${TF_SITE_DIR}/motif_centered_site_summary.tsv}
CORES=${CORES:-4}
MOTIFS=${MOTIFS:-data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt}

if [[ "${MODE}" == "full" ]]; then
  OUTDIR=${OUTDIR:-data/public/processed/pseudobulk_pbmc/footprints_full}
  INCLUDE_CHROMS=${INCLUDE_CHROMS:-chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX}
  GROUP_NAMES=${GROUP_NAMES:-B_cell,CD4_T,NK_T_cytotoxic,CD14_Monocyte,FCGR3A_Monocyte,Dendritic_cell,Mixed_myeloid}
  TFS=${TFS:-SPIB,RUNX3,CEBPB,CEBPA}
else
  OUTDIR=${OUTDIR:-data/public/processed/pseudobulk_pbmc/footprint_demo}
  INCLUDE_CHROMS=${INCLUDE_CHROMS:-chr1,chr2}
  GROUP_NAMES=${GROUP_NAMES:-B_cell,CD4_T,CD14_Monocyte}
  TFS=${TFS:-SPIB,RUNX3,CEBPB}
fi

"${FP_TOOLS_ENV}/bin/pseudobulk-footprints" \
  --fragments "${FRAGMENTS}" \
  --annotations "${ANNOTATIONS}" \
  --group-by cell_type \
  --min-cells 300 \
  --min-fragments 50000 \
  --include-chroms "${INCLUDE_CHROMS}" \
  --groups "${GROUP_NAMES}" \
  --genome-sizes "${GENOME_SIZES}" \
  --genome "${GENOME}" \
  --peaks "${PEAKS}" \
  --motifs "${MOTIFS}" \
  --tf-site-dir "${TF_SITE_DIR}" \
  --site-summary "${SITE_SUMMARY}" \
  --tfs "${TFS}" \
  --outdir "${OUTDIR}" \
  --cores "${CORES}" \
  "${@:2}"
