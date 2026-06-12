#!/usr/bin/env bash
set -Eeuo pipefail
on_error() {
  local rc=$?
  printf '[%s] ERROR line %s rc=%s: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${BASH_LINENO[0]:-${LINENO}}" "$rc" "$BASH_COMMAND" >&2
}
trap on_error ERR
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="${ROOT_DIR}/data/public/raw/buenrostro_atac_replicates"
GENOME="${ROOT_DIR}/data/public/raw/genome/hg38.fa"
OUT_DIR="${ROOT_DIR}/data/public/processed/buenrostro_atac_replicates"
ENV_DIR="${FP_TOOLS_ATAC_ENV:-/home/exouser/miniforge3/envs/fp-tools-atac}"
FP_TOOLS_ENV="${FP_TOOLS_ENV:-${ROOT_DIR}/.venv}"
THREADS="${THREADS:-24}"
BOWTIE2_THREADS="${BOWTIE2_THREADS:-${THREADS}}"
MACS2_QVALUE="${MACS2_QVALUE:-0.01}"

BIN_DIR="${ENV_DIR}/bin"
FASTP="${BIN_DIR}/fastp"
BOWTIE2="${BIN_DIR}/bowtie2"
BOWTIE2_BUILD="${BIN_DIR}/bowtie2-build"
SAMTOOLS="${BIN_DIR}/samtools"
BEDTOOLS="${BIN_DIR}/bedtools"
if [[ -x "${BIN_DIR}/macs3" ]]; then
  MACS_CALLPEAK="${BIN_DIR}/macs3"
else
  MACS_CALLPEAK="${BIN_DIR}/macs2"
fi
ATAC_CORRECT="${FP_TOOLS_ENV}/bin/atac-correct"
CALL_FOOTPRINTS="${FP_TOOLS_ENV}/bin/call-footprints"
DIFF_FOOTPRINTS="${FP_TOOLS_ENV}/bin/diff-footprints"
JASPAR2026_MOTIFS="${ROOT_DIR}/data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

TRIM_DIR="${OUT_DIR}/trimmed_fastq"
BAM_DIR="${OUT_DIR}/bam"
PEAK_DIR="${OUT_DIR}/peaks"
FP_DIR="${OUT_DIR}/fp_tools"
REF_DIR="${OUT_DIR}/reference"
LOG_DIR="${OUT_DIR}/logs"
STATUS_DIR="${OUT_DIR}/status"
INDEX_PREFIX="${REF_DIR}/bowtie2/hg38"
BLACKLIST="${REF_DIR}/hg38-blacklist.v2.bed"
CHROM_SIZES="${REF_DIR}/hg38.chrom.sizes"
MERGED_PEAKS="${PEAK_DIR}/merged_peaks.bed"
RUN_MANIFEST="${OUT_DIR}/run_manifest.tsv"

mkdir -p "${TRIM_DIR}" "${BAM_DIR}" "${PEAK_DIR}" "${FP_DIR}" "${REF_DIR}/bowtie2" "${LOG_DIR}" "${STATUS_DIR}"

SAMPLES=(
  "Bcell_rep1 Bcell SRR891268"
  "Bcell_rep2 Bcell SRR891269"
  "Tcell_rep1 Tcell SRR891275"
  "Tcell_rep2 Tcell SRR891276"
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

bowtie2_index_complete() {
  local suffix
  if [[ -s "${INDEX_PREFIX}.1.bt2l" ]]; then
    for suffix in .1.bt2l .2.bt2l .3.bt2l .4.bt2l .rev.1.bt2l .rev.2.bt2l; do
      [[ -s "${INDEX_PREFIX}${suffix}" ]] || return 1
    done
    return 0
  fi
  for suffix in .1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2; do
    [[ -s "${INDEX_PREFIX}${suffix}" ]] || return 1
  done
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

check_inputs() {
  require_file "${GENOME}"
  require_file "${GENOME}.fai"
  require_file "${JASPAR2026_MOTIFS}"
  require_exec "${FASTP}"
  require_exec "${BOWTIE2}"
  require_exec "${BOWTIE2_BUILD}"
  require_exec "${SAMTOOLS}"
  require_exec "${BEDTOOLS}"
  require_exec "${MACS_CALLPEAK}"
  require_exec "${ATAC_CORRECT}"
  require_exec "${CALL_FOOTPRINTS}"
  require_exec "${DIFF_FOOTPRINTS}"
  for entry in "${SAMPLES[@]}"; do
    read -r _sample _condition run <<< "${entry}"
    require_file "${RAW_DIR}/${run}_1.fastq.gz"
    require_file "${RAW_DIR}/${run}_2.fastq.gz"
  done
}

write_versions() {
  {
    echo -e "tool	path"
    echo -e "fastp	${FASTP}"
    echo -e "bowtie2	${BOWTIE2}"
    echo -e "samtools	${SAMTOOLS}"
    echo -e "bedtools	${BEDTOOLS}"
    echo -e "macs-callpeak	${MACS_CALLPEAK}"
    echo -e "atac-correct	${ATAC_CORRECT}"
    echo -e "call-footprints	${CALL_FOOTPRINTS}"
    echo -e "diff-footprints	${DIFF_FOOTPRINTS}"
    echo -e "JASPAR2026 vertebrates	${JASPAR2026_MOTIFS}"
  } > "${OUT_DIR}/software_paths.tsv"

  {
    echo -e "tool	version"
    echo -e "fastp	$(${FASTP} --version 2>&1 | awk '{print $2}')"
    echo -e "bowtie2	$(${BOWTIE2} --version 2>&1 | awk '/version/{print $NF; exit}')"
    echo -e "samtools	$(${SAMTOOLS} --version 2>&1 | awk 'NR==1{print $2}')"
    echo -e "htslib	$(${SAMTOOLS} --version 2>&1 | awk '/Using htslib/{print $3; exit}')"
    echo -e "bedtools	$(${BEDTOOLS} --version 2>&1 | awk '{print $2}')"
    echo -e "macs-callpeak	$(${MACS_CALLPEAK} --version 2>&1 | awk '{print $2}')"
    "${FP_TOOLS_ENV}/bin/python" - <<'PYVERS'
import sys
import fp_tools
mods = ['numpy', 'pandas', 'scipy', 'matplotlib', 'pyBigWig', 'pysam', 'Bio']
print(f"python	{sys.version.split()[0]}")
print(f"fp-tools	{fp_tools.__version__}")
for name in mods:
    try:
        mod = __import__(name)
        print(f"{name}	{getattr(mod, '__version__', 'unknown')}")
    except Exception:
        print(f"{name}	unavailable")
PYVERS
    echo -e "JASPAR	2026 CORE vertebrates non-redundant"
    echo -e "Cell Ranger ARC	2.0.0 (10x PBMC source processing)"
  } > "${OUT_DIR}/software_versions.tsv"
}

write_parameters() {
  {
    echo -e "analysis	parameters"
    echo -e "fastp	--detect_adapter_for_pe --thread ${THREADS}"
    echo -e "bowtie2	--very-sensitive -X 2000 -p ${BOWTIE2_THREADS}"
    echo -e "samtools_filter	view -b -f 2 -F 2828 -q 30; fixmate -m; markdup -r"
    echo -e "blacklist_filter	exclude chrM/MT and hg38-blacklist.v2 regions"
    echo -e "macs_callpeak	callpeak -f BAMPE -g hs --keep-dup all -q ${MACS2_QVALUE}"
    echo -e "merged_peaks	bedtools sort and merge across four replicate narrowPeak files; exclude '_' contigs and chrM/MT"
    echo -e "atac-correct	--peaks merged_peaks.bed --blacklist hg38-blacklist.v2.bed --cores ${THREADS}; defaults: extend=100, k_flank=12, read_shift=4,-5, bg_shift=100, window=100, score_mat=DWM"
    echo -e "call-footprints	--score footprint --regions merged_peaks.bed --cores ${THREADS}; defaults: fp_min=20, fp_max=50, flank_min=10, flank_max=30, smooth=1"
    echo -e "diff-footprints_none	JASPAR2026 CORE vertebrates non-redundant; --cond-names Bcell Bcell Tcell Tcell --normalization none --replicate-report auto --aggregate-signals corrected cut-site bigWigs --plot-aggregate top --plot-aggregate-top-n 5 --skip-excel --cores ${THREADS}"
    echo -e "diff-footprints_sample_quantile	JASPAR2026 CORE vertebrates non-redundant; --cond-names Bcell Bcell Tcell Tcell --normalization sample-quantile --replicate-report auto --aggregate-signals corrected cut-site bigWigs --plot-aggregate top --plot-aggregate-top-n 5 --skip-excel --cores ${THREADS}"
    echo -e "pseudobulk_pbmc	10x PBMC Multiome fragments grouped by broad immune labels; min_cells=300, min_fragments=50000, CPM-normalized cut-site bigWigs; chr1-chr22 and chrX for figures"
  } > "${OUT_DIR}/analysis_parameters.tsv"
}


prepare_reference() {
  awk 'BEGIN{OFS="\t"}{print $1,$2}' "${GENOME}.fai" > "${CHROM_SIZES}"
  if [[ ! -s "${BLACKLIST}" ]]; then
    log "Downloading hg38 blacklist"
    wget -O "${BLACKLIST}.gz" "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
    gzip -dc "${BLACKLIST}.gz" > "${BLACKLIST}"
  fi
  if ! bowtie2_index_complete; then
    rm -f "${INDEX_PREFIX}".*.bt2 "${INDEX_PREFIX}".*.bt2l "${INDEX_PREFIX}".*.tmp
    "${BOWTIE2_BUILD}" --threads "${THREADS}" "${GENOME}" "${INDEX_PREFIX}"
  fi
}

process_sample() {
  local sample="$1"
  local condition="$2"
  local run="$3"
  local r1="${RAW_DIR}/${run}_1.fastq.gz"
  local r2="${RAW_DIR}/${run}_2.fastq.gz"
  local trim1="${TRIM_DIR}/${sample}_R1.trim.fastq.gz"
  local trim2="${TRIM_DIR}/${sample}_R2.trim.fastq.gz"
  local name_bam="${BAM_DIR}/${sample}.name_sorted.bam"
  local fixmate_bam="${BAM_DIR}/${sample}.fixmate.bam"
  local coord_bam="${BAM_DIR}/${sample}.coord_sorted.bam"
  local dedup_bam="${BAM_DIR}/${sample}.dedup.bam"
  local nomito_bam="${BAM_DIR}/${sample}.dedup.no_mito.bam"
  local final_bam="${BAM_DIR}/${sample}.filtered.bam"
  local keep_chroms="${BAM_DIR}/${sample}.keep_chroms.txt"
  local macs_name="${sample}"

  run_step "${sample}.fastp" "${FASTP}" \
    --in1 "${r1}" --in2 "${r2}" \
    --out1 "${trim1}" --out2 "${trim2}" \
    --detect_adapter_for_pe \
    --thread "${THREADS}" \
    --html "${LOG_DIR}/${sample}.fastp.html" \
    --json "${LOG_DIR}/${sample}.fastp.json"

  run_step "${sample}.align_fixmate_markdup" bash -lc "
    set -euo pipefail
    '${BOWTIE2}' --very-sensitive -X 2000 -p '${BOWTIE2_THREADS}' -x '${INDEX_PREFIX}' -1 '${trim1}' -2 '${trim2}' 2> '${LOG_DIR}/${sample}.bowtie2.log' |
      '${SAMTOOLS}' view -@ '${THREADS}' -b -f 2 -F 2828 -q 30 - |
      '${SAMTOOLS}' sort -@ '${THREADS}' -n -o '${name_bam}' -
    '${SAMTOOLS}' fixmate -@ '${THREADS}' -m '${name_bam}' '${fixmate_bam}'
    '${SAMTOOLS}' sort -@ '${THREADS}' -o '${coord_bam}' '${fixmate_bam}'
    '${SAMTOOLS}' markdup -@ '${THREADS}' -r '${coord_bam}' '${dedup_bam}'
    '${SAMTOOLS}' index -@ '${THREADS}' '${dedup_bam}'
    '${SAMTOOLS}' flagstat -@ '${THREADS}' '${dedup_bam}' > '${LOG_DIR}/${sample}.dedup.flagstat.txt'
    rm -f '${name_bam}' '${fixmate_bam}' '${coord_bam}'
  "

  run_step "${sample}.filter_blacklist" bash -lc "
    set -euo pipefail
    '${SAMTOOLS}' idxstats '${dedup_bam}' | awk '\$1 != \"*\" && \$1 != \"chrM\" && \$1 != \"MT\" {print \$1}' > '${keep_chroms}'
    '${SAMTOOLS}' view -@ '${THREADS}' -b '${dedup_bam}' \$(cat '${keep_chroms}') > '${nomito_bam}'
    '${BEDTOOLS}' intersect -v -abam '${nomito_bam}' -b '${BLACKLIST}' |
      '${SAMTOOLS}' sort -@ '${THREADS}' -o '${final_bam}' -
    '${SAMTOOLS}' index -@ '${THREADS}' '${final_bam}'
    '${SAMTOOLS}' flagstat -@ '${THREADS}' '${final_bam}' > '${LOG_DIR}/${sample}.filtered.flagstat.txt'
    '${SAMTOOLS}' idxstats '${final_bam}' > '${LOG_DIR}/${sample}.filtered.idxstats.txt'
  "

  run_step "${sample}.macs2" "${MACS_CALLPEAK}" callpeak \
    -t "${final_bam}" \
    -f BAMPE \
    -g hs \
    -n "${macs_name}" \
    --outdir "${PEAK_DIR}" \
    --keep-dup all \
    -q "${MACS2_QVALUE}"
}

merge_peaks() {
  local all_peaks="${PEAK_DIR}/all_replicate_peaks.bed"
  cat "${PEAK_DIR}"/*_peaks.narrowPeak |
    awk 'BEGIN{OFS="\t"} $1 !~ /_/ && $1 != "chrM" && $1 != "MT" {print $1,$2,$3}' |
    "${BEDTOOLS}" sort -i - |
    "${BEDTOOLS}" merge -i - > "${MERGED_PEAKS}"
  cp "${MERGED_PEAKS}" "${all_peaks}"
  awk 'BEGIN{OFS="\t"}{print $0,"peak_"NR}' "${MERGED_PEAKS}" > "${PEAK_DIR}/merged_peaks_named.bed"
}

run_fp_tools_for_sample() {
  local sample="$1"
  local bam="${BAM_DIR}/${sample}.filtered.bam"
  local atac_dir="${FP_DIR}/atac_correct/${sample}"
  local footprint_bw="${FP_DIR}/footprints/${sample}.footprints.bw"
  mkdir -p "${atac_dir}" "${FP_DIR}/footprints"

  run_step "${sample}.atac_correct" "${ATAC_CORRECT}" \
    --bam "${bam}" \
    --genome "${GENOME}" \
    --peaks "${MERGED_PEAKS}" \
    --blacklist "${BLACKLIST}" \
    --outdir "${atac_dir}" \
    --cores "${THREADS}"

  local corrected_bw
  corrected_bw="$(find "${atac_dir}" -maxdepth 1 -name '*_corrected.bw' | head -1)"
  require_file "${corrected_bw}"

  run_step "${sample}.call_footprints" "${CALL_FOOTPRINTS}" \
    --signal "${corrected_bw}" \
    --regions "${MERGED_PEAKS}" \
    --output "${footprint_bw}" \
    --score footprint \
    --cores "${THREADS}"
}

run_diff_footprints() {
  local normalization="$1"
  local outdir="${FP_DIR}/diff_footprints_jaspar2026_vertebrates_norm_${normalization//-/_}"
  mkdir -p "${outdir}"
  run_step "diff_footprints.${normalization}" "${DIFF_FOOTPRINTS}" \
    --motifs "${JASPAR2026_MOTIFS}" \
    --signals \
      "${FP_DIR}/footprints/Bcell_rep1.footprints.bw" \
      "${FP_DIR}/footprints/Bcell_rep2.footprints.bw" \
      "${FP_DIR}/footprints/Tcell_rep1.footprints.bw" \
      "${FP_DIR}/footprints/Tcell_rep2.footprints.bw" \
    --genome "${GENOME}" \
    --peaks "${MERGED_PEAKS}" \
    --outdir "${outdir}" \
    --cond-names Bcell Bcell Tcell Tcell \
    --normalization "${normalization}" \
    --replicate-report auto \
    --aggregate-signals \
      "${FP_DIR}/atac_correct/Bcell_rep1/Bcell_rep1.filtered_corrected.bw" \
      "${FP_DIR}/atac_correct/Bcell_rep2/Bcell_rep2.filtered_corrected.bw" \
      "${FP_DIR}/atac_correct/Tcell_rep1/Tcell_rep1.filtered_corrected.bw" \
      "${FP_DIR}/atac_correct/Tcell_rep2/Tcell_rep2.filtered_corrected.bw" \
    --plot-aggregate top \
    --plot-aggregate-top-n 5 \
    --aggregate-flank 100 \
    --skip-excel \
    --cores "${THREADS}"
}

write_manifest() {
  {
    echo -e "sample\tcondition\trun\tfastq1\tfastq2\tfiltered_bam\tpeaks\tfootprint_bigwig"
    for entry in "${SAMPLES[@]}"; do
      read -r sample condition run <<< "${entry}"
      echo -e "${sample}\t${condition}\t${run}\t${RAW_DIR}/${run}_1.fastq.gz\t${RAW_DIR}/${run}_2.fastq.gz\t${BAM_DIR}/${sample}.filtered.bam\t${PEAK_DIR}/${sample}_peaks.narrowPeak\t${FP_DIR}/footprints/${sample}.footprints.bw"
    done
  } > "${RUN_MANIFEST}"
}

main() {
  log "Workflow started"
  check_inputs
  write_versions
  write_parameters
  run_step "reference.prepare" prepare_reference
  for entry in "${SAMPLES[@]}"; do
    read -r sample condition run <<< "${entry}"
    process_sample "${sample}" "${condition}" "${run}"
  done
  run_step "peaks.merge" merge_peaks
  for entry in "${SAMPLES[@]}"; do
    read -r sample _condition _run <<< "${entry}"
    run_fp_tools_for_sample "${sample}"
  done
  run_diff_footprints none
  run_diff_footprints sample-quantile
  write_manifest
  log "Workflow finished"
}

main "$@"
