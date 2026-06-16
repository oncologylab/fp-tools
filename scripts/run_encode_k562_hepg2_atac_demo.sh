#!/usr/bin/env bash
set -Eeuo pipefail

on_error() {
  local rc=$?
  printf '[%s] ERROR line %s rc=%s: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${BASH_LINENO[0]:-${LINENO}}" "$rc" "$BASH_COMMAND" >&2
}
trap on_error ERR

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="${ROOT_DIR}/data/public/raw/encode_k562_hepg2_atac"
OUT_DIR="${ROOT_DIR}/data/public/processed/encode_k562_hepg2_atac_replicates"
FP_DIR="${OUT_DIR}/fp_tools"
REF_DIR="${OUT_DIR}/reference"
LOG_DIR="${OUT_DIR}/logs"
STATUS_DIR="${OUT_DIR}/status"
PEAK_DIR="${OUT_DIR}/peaks"

GENOME="${ROOT_DIR}/data/public/raw/genome/hg38.fa"
JASPAR2026_MOTIFS="${ROOT_DIR}/data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
BLACKLIST="${REF_DIR}/hg38-blacklist.v2.bed"
CHROM_SIZES="${REF_DIR}/hg38.chrom.sizes"
MERGED_PEAKS="${PEAK_DIR}/merged_peaks.bed"
PEAK_BINS="${PEAK_DIR}/merged_peaks.50bp_bins.bed"
RUN_MANIFEST="${OUT_DIR}/run_manifest.tsv"

FP_TOOLS_ENV="${FP_TOOLS_ENV:-${ROOT_DIR}/.venv}"
ENV_DIR="${FP_TOOLS_ATAC_ENV:-/home/exouser/miniforge3/envs/fp-tools-atac}"
BIN_DIR="${ENV_DIR}/bin"
SAMTOOLS="${BIN_DIR}/samtools"
BEDTOOLS="${BIN_DIR}/bedtools"
ATAC_CORRECT="${FP_TOOLS_ENV}/bin/atac-correct"
CALL_FOOTPRINTS="${FP_TOOLS_ENV}/bin/call-footprints"
DIFF_FOOTPRINTS="${FP_TOOLS_ENV}/bin/diff-footprints"
NORMALIZE_BIGWIG="${FP_TOOLS_ENV}/bin/normalize-bigwig"
THREADS="${THREADS:-$(nproc)}"

mkdir -p "${RAW_DIR}/bam" "${RAW_DIR}/peaks" "${FP_DIR}" "${REF_DIR}" "${LOG_DIR}" "${STATUS_DIR}" "${PEAK_DIR}"

SAMPLES=(
  "K562_rep1 K562 ENCSR868FGK ENCFF077FBI https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam"
  "K562_rep2 K562 ENCSR868FGK ENCFF128WZG https://www.encodeproject.org/files/ENCFF128WZG/@@download/ENCFF128WZG.bam"
  "K562_rep3 K562 ENCSR868FGK ENCFF534DCE https://www.encodeproject.org/files/ENCFF534DCE/@@download/ENCFF534DCE.bam"
  "HepG2_rep1 HepG2 ENCSR291GJU ENCFF624SON https://www.encodeproject.org/files/ENCFF624SON/@@download/ENCFF624SON.bam"
  "HepG2_rep2 HepG2 ENCSR291GJU ENCFF926KFU https://www.encodeproject.org/files/ENCFF926KFU/@@download/ENCFF926KFU.bam"
  "HepG2_rep3 HepG2 ENCSR291GJU ENCFF990VCP https://www.encodeproject.org/files/ENCFF990VCP/@@download/ENCFF990VCP.bam"
)

PEAK_FILES=(
  "K562 ENCSR868FGK ENCFF135AEX https://www.encodeproject.org/files/ENCFF135AEX/@@download/ENCFF135AEX.bed.gz"
  "K562 ENCSR868FGK ENCFF223QDM https://www.encodeproject.org/files/ENCFF223QDM/@@download/ENCFF223QDM.bed.gz"
  "K562 ENCSR868FGK ENCFF433EPT https://www.encodeproject.org/files/ENCFF433EPT/@@download/ENCFF433EPT.bed.gz"
  "K562 ENCSR868FGK ENCFF771HDN https://www.encodeproject.org/files/ENCFF771HDN/@@download/ENCFF771HDN.bed.gz"
  "K562 ENCSR868FGK ENCFF948AFM https://www.encodeproject.org/files/ENCFF948AFM/@@download/ENCFF948AFM.bed.gz"
  "K562 ENCSR868FGK ENCFF993BAP https://www.encodeproject.org/files/ENCFF993BAP/@@download/ENCFF993BAP.bed.gz"
  "HepG2 ENCSR291GJU ENCFF161ZZX https://www.encodeproject.org/files/ENCFF161ZZX/@@download/ENCFF161ZZX.bed.gz"
  "HepG2 ENCSR291GJU ENCFF576UEM https://www.encodeproject.org/files/ENCFF576UEM/@@download/ENCFF576UEM.bed.gz"
  "HepG2 ENCSR291GJU ENCFF609BSU https://www.encodeproject.org/files/ENCFF609BSU/@@download/ENCFF609BSU.bed.gz"
  "HepG2 ENCSR291GJU ENCFF915FZC https://www.encodeproject.org/files/ENCFF915FZC/@@download/ENCFF915FZC.bed.gz"
  "HepG2 ENCSR291GJU ENCFF919RKQ https://www.encodeproject.org/files/ENCFF919RKQ/@@download/ENCFF919RKQ.bed.gz"
  "HepG2 ENCSR291GJU ENCFF935GLR https://www.encodeproject.org/files/ENCFF935GLR/@@download/ENCFF935GLR.bed.gz"
)

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_file() {
  if [[ ! -s "$1" ]]; then
    log "ERROR: missing required file: $1"
    exit 1
  fi
}

require_exec() {
  if [[ ! -x "$1" ]]; then
    log "ERROR: missing executable: $1"
    exit 1
  fi
}

done_marker() {
  [[ -s "${STATUS_DIR}/$1.done" ]]
}

mark_done() {
  date '+%Y-%m-%d %H:%M:%S' > "${STATUS_DIR}/$1.done"
}

run_step() {
  local name="$1"
  shift
  if done_marker "${name}"; then
    log "SKIP ${name}"
    return 0
  fi
  log "START ${name}"
  "$@"
  mark_done "${name}"
  log "DONE ${name}"
}

download_file() {
  local url="$1"
  local out="$2"
  local done="${out}.done"
  if [[ -s "${out}" && -s "${done}" ]]; then
    log "SKIP download $(basename "${out}")"
    return 0
  fi
  mkdir -p "$(dirname "${out}")"
  wget -c -O "${out}" "${url}"
  touch "${done}"
}

check_inputs() {
  require_file "${GENOME}"
  require_file "${GENOME}.fai"
  require_file "${JASPAR2026_MOTIFS}"
  require_exec "${SAMTOOLS}"
  require_exec "${BEDTOOLS}"
  require_exec "${ATAC_CORRECT}"
  require_exec "${CALL_FOOTPRINTS}"
  require_exec "${DIFF_FOOTPRINTS}"
  require_exec "${NORMALIZE_BIGWIG}"
}

prepare_reference() {
  awk 'BEGIN{OFS="\t"}{print $1,$2}' "${GENOME}.fai" > "${CHROM_SIZES}"
  if [[ ! -s "${BLACKLIST}" ]]; then
    log "Downloading hg38 blacklist"
    wget -O "${BLACKLIST}.gz" "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
    gzip -dc "${BLACKLIST}.gz" > "${BLACKLIST}"
  fi
}

download_encode_inputs() {
  local sample condition experiment accession url
  for entry in "${SAMPLES[@]}"; do
    read -r sample condition experiment accession url <<< "${entry}"
    download_file "${url}" "${RAW_DIR}/bam/${sample}.${accession}.bam"
  done
  for entry in "${PEAK_FILES[@]}"; do
    read -r condition experiment accession url <<< "${entry}"
    download_file "${url}" "${RAW_DIR}/peaks/${condition}.${experiment}.${accession}.bed.gz"
  done
}

index_bams() {
  local sample condition experiment accession url bam
  for entry in "${SAMPLES[@]}"; do
    read -r sample condition experiment accession url <<< "${entry}"
    bam="${RAW_DIR}/bam/${sample}.${accession}.bam"
    run_step "${sample}.bam_index" "${SAMTOOLS}" index -@ "${THREADS}" "${bam}"
  done
}

merge_peaks() {
  gzip -dc "${RAW_DIR}/peaks/"*.bed.gz |
    awk 'BEGIN{OFS="\t"} $1 !~ /_/ && $1 != "chrM" && $1 != "MT" {print $1,$2,$3}' |
    "${BEDTOOLS}" sort -i - |
    "${BEDTOOLS}" merge -i - > "${MERGED_PEAKS}"
  "${BEDTOOLS}" makewindows -b "${MERGED_PEAKS}" -w 50 > "${PEAK_BINS}"
  awk 'BEGIN{OFS="\t"}{print $0,"peak_"NR}' "${MERGED_PEAKS}" > "${PEAK_DIR}/merged_peaks_named.bed"
}

write_metadata() {
  {
    echo -e "sample\tcondition\texperiment_accession\tfile_accession\tlocal_bam"
    local sample condition experiment accession url
    for entry in "${SAMPLES[@]}"; do
      read -r sample condition experiment accession url <<< "${entry}"
      echo -e "${sample}\t${condition}\t${experiment}\t${accession}\t${RAW_DIR}/bam/${sample}.${accession}.bam"
    done
  } > "${RUN_MANIFEST}"
  {
    echo -e "tool\tpath"
    echo -e "samtools\t${SAMTOOLS}"
    echo -e "bedtools\t${BEDTOOLS}"
    echo -e "atac-correct\t${ATAC_CORRECT}"
    echo -e "call-footprints\t${CALL_FOOTPRINTS}"
    echo -e "normalize-bigwig\t${NORMALIZE_BIGWIG}"
    echo -e "diff-footprints\t${DIFF_FOOTPRINTS}"
    echo -e "JASPAR2026 vertebrates\t${JASPAR2026_MOTIFS}"
  } > "${OUT_DIR}/software_paths.tsv"
  {
    echo -e "analysis\tparameters"
    echo -e "encode_inputs\tK562 ENCSR868FGK vs HepG2 ENCSR291GJU; three released GRCh38 alignment BAM biological replicates per condition"
    echo -e "merged_peaks\tbedtools sort and merge across released ENCODE IDR thresholded peak BEDs; exclude '_' contigs and chrM/MT"
    echo -e "atac-correct\t--bam ENCODE alignment BAM --genome hg38.fa --peaks merged_peaks.bed --blacklist hg38-blacklist.v2.bed --outdir sample"
    echo -e "call-footprints\t--signal sample_corrected.bw --regions merged_peaks.bed --score footprint"
    echo -e "normalize-bigwig\t--background merged_peaks.50bp_bins.bed --method background-scale --stat q95 --target median"
    echo -e "diff-footprints_sample_quantile\tJASPAR2026 CORE vertebrates non-redundant; --normalization sample-quantile --cond-names K562 K562 K562 HepG2 HepG2 HepG2 --aggregate-site-set bound --plot-aggregate sig --aggregate-flank 100"
    echo -e "diff-footprints_none\tSame inputs with --normalization none; aggregate tracks use unnormalized corrected bigWigs"
  } > "${OUT_DIR}/analysis_parameters.tsv"
}

corrected_bw_for_sample() {
  local sample="$1"
  find "${FP_DIR}/atac_correct/${sample}" -maxdepth 1 -name '*_corrected.bw' | sort | head -1
}

run_fp_tools_for_sample() {
  local sample="$1"
  local accession="$2"
  local bam="${RAW_DIR}/bam/${sample}.${accession}.bam"
  local atac_dir="${FP_DIR}/atac_correct/${sample}"
  local footprint_bw="${FP_DIR}/footprints/${sample}.footprints.bw"
  mkdir -p "${atac_dir}" "${FP_DIR}/footprints"

  run_step "${sample}.atac_correct" "${ATAC_CORRECT}" \
    --bam "${bam}" \
    --genome "${GENOME}" \
    --peaks "${MERGED_PEAKS}" \
    --blacklist "${BLACKLIST}" \
    --outdir "${atac_dir}"

  local corrected_bw
  corrected_bw="$(corrected_bw_for_sample "${sample}")"
  require_file "${corrected_bw}"

  run_step "${sample}.call_footprints" "${CALL_FOOTPRINTS}" \
    --signal "${corrected_bw}" \
    --regions "${MERGED_PEAKS}" \
    --output "${footprint_bw}" \
    --score footprint
}

normalize_corrected_bigwigs() {
  local norm_dir="${FP_DIR}/normalized_corrected_bigwigs/peak_q95"
  mkdir -p "${norm_dir}"
  run_step "normalize_corrected_bigwigs.peak_q95" "${NORMALIZE_BIGWIG}" \
    --bigwigs \
      "$(corrected_bw_for_sample K562_rep1)" \
      "$(corrected_bw_for_sample K562_rep2)" \
      "$(corrected_bw_for_sample K562_rep3)" \
      "$(corrected_bw_for_sample HepG2_rep1)" \
      "$(corrected_bw_for_sample HepG2_rep2)" \
      "$(corrected_bw_for_sample HepG2_rep3)" \
    --background "${PEAK_BINS}" \
    --method background-scale \
    --stat q95 \
    --target median \
    --chrom-sizes "${CHROM_SIZES}" \
    --outdir "${norm_dir}"
}

normalized_bw_for_sample() {
  local sample="$1"
  find "${FP_DIR}/normalized_corrected_bigwigs/peak_q95" -maxdepth 1 -name "*${sample}*.background_scale_q95.bw" | sort | head -1
}

run_diff_footprints() {
  local normalization="$1"
  local outdir="${FP_DIR}/diff_footprints_jaspar2026_vertebrates_norm_${normalization//-/_}"
  mkdir -p "${outdir}"
  local aggregate_args=()
  if [[ "${normalization}" == "sample-quantile" ]]; then
    aggregate_args=(
      "$(normalized_bw_for_sample K562_rep1)"
      "$(normalized_bw_for_sample K562_rep2)"
      "$(normalized_bw_for_sample K562_rep3)"
      "$(normalized_bw_for_sample HepG2_rep1)"
      "$(normalized_bw_for_sample HepG2_rep2)"
      "$(normalized_bw_for_sample HepG2_rep3)"
    )
  else
    aggregate_args=(
      "$(corrected_bw_for_sample K562_rep1)"
      "$(corrected_bw_for_sample K562_rep2)"
      "$(corrected_bw_for_sample K562_rep3)"
      "$(corrected_bw_for_sample HepG2_rep1)"
      "$(corrected_bw_for_sample HepG2_rep2)"
      "$(corrected_bw_for_sample HepG2_rep3)"
    )
  fi
  run_step "diff_footprints.${normalization}" "${DIFF_FOOTPRINTS}" \
    --motifs "${JASPAR2026_MOTIFS}" \
    --signals \
      "${FP_DIR}/footprints/K562_rep1.footprints.bw" \
      "${FP_DIR}/footprints/K562_rep2.footprints.bw" \
      "${FP_DIR}/footprints/K562_rep3.footprints.bw" \
      "${FP_DIR}/footprints/HepG2_rep1.footprints.bw" \
      "${FP_DIR}/footprints/HepG2_rep2.footprints.bw" \
      "${FP_DIR}/footprints/HepG2_rep3.footprints.bw" \
    --genome "${GENOME}" \
    --peaks "${MERGED_PEAKS}" \
    --outdir "${outdir}" \
    --prefix diff_footprints \
    --cond-names K562 K562 K562 HepG2 HepG2 HepG2 \
    --normalization "${normalization}" \
    --replicate-report auto \
    --aggregate-signals "${aggregate_args[@]}" \
    --aggregate-normalization none \
    --aggregate-site-set bound \
    --plot-aggregate sig \
    --aggregate-flank 100 \
    --skip-excel \
    --cores "${THREADS}"
}

main() {
  log "ENCODE K562 vs HepG2 workflow started"
  check_inputs
  write_metadata
  run_step "reference.prepare" prepare_reference
  run_step "encode.download" download_encode_inputs
  index_bams
  run_step "peaks.merge" merge_peaks
  local sample condition experiment accession url
  for entry in "${SAMPLES[@]}"; do
    read -r sample condition experiment accession url <<< "${entry}"
    run_fp_tools_for_sample "${sample}" "${accession}"
  done
  normalize_corrected_bigwigs
  run_diff_footprints sample-quantile
  run_diff_footprints none
  log "ENCODE K562 vs HepG2 workflow finished"
  log "Primary report: ${FP_DIR}/diff_footprints_jaspar2026_vertebrates_norm_sample_quantile/diff_footprints_K562_HepG2.html"
  log "Baseline report: ${FP_DIR}/diff_footprints_jaspar2026_vertebrates_norm_none/diff_footprints_K562_HepG2.html"
}

main "$@"
