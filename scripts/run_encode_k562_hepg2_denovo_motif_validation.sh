#!/usr/bin/env bash
set -Eeuo pipefail

on_error() {
  local rc=$?
  printf '[%s] ERROR line %s rc=%s: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${BASH_LINENO[0]:-${LINENO}}" "$rc" "$BASH_COMMAND" >&2
}
trap on_error ERR

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOME="${ROOT_DIR}/data/public/raw/genome/hg38.fa"
JASPAR2026="${ROOT_DIR}/data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
BASE_DIR="${ROOT_DIR}/data/public/processed/encode_k562_hepg2_atac_replicates"
PEAKS="${BASE_DIR}/peaks/merged_peaks.bed"
FP_DIR="${BASE_DIR}/fp_tools"
OUT_DIR="${OUT_DIR:-${FP_DIR}/denovo_motif_validation}"

FP_TOOLS_ENV="${FP_TOOLS_ENV:-${ROOT_DIR}/.venv}"
CALL_FOOTPRINTS="${FP_TOOLS_ENV}/bin/call-footprints"
MOTIF_DISCOVERY="${FP_TOOLS_ENV}/bin/motif-discovery"
DIFF_FOOTPRINTS="${FP_TOOLS_ENV}/bin/diff-footprints"
PYTHON="${FP_TOOLS_ENV}/bin/python"
STREME="${STREME:-/home/exouser/miniforge3/envs/fp-tools-atac/bin/streme}"
TOMTOM="${TOMTOM:-/home/exouser/miniforge3/envs/fp-tools-atac/bin/tomtom}"

TOP_N_PER_REPLICATE="${TOP_N_PER_REPLICATE:-5000}"
TOP_N_PER_CONDITION="${TOP_N_PER_CONDITION:-8000}"
CALL_WIDTH="${CALL_WIDTH:-50}"
MIN_DISTANCE="${MIN_DISTANCE:-40}"
FLANK="${FLANK:-75}"
STREME_NMOTIFS="${STREME_NMOTIFS-8}"
STREME_THRESH="${STREME_THRESH:-0.05}"
STREME_PATIENCE="${STREME_PATIENCE:-10}"
THREADS="${THREADS:-$(nproc)}"
DIFF_PLOT_AGGREGATE="${DIFF_PLOT_AGGREGATE:-sig}"
DIFF_REUSE_EXISTING="${DIFF_REUSE_EXISTING:-0}"

CONDITIONS=(K562 HepG2)
SAMPLES=(K562_rep1 K562_rep2 K562_rep3 HepG2_rep1 HepG2_rep2 HepG2_rep3)

mkdir -p "${OUT_DIR}/status" "${OUT_DIR}/candidate_calls" "${OUT_DIR}/candidate_fastas" "${OUT_DIR}/motifs" "${OUT_DIR}/diff_footprints"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
require_file() { [[ -s "$1" ]] || { log "ERROR missing file: $1"; exit 1; }; }
require_exec() { [[ -x "$1" ]] || { log "ERROR missing executable: $1"; exit 1; }; }
done_marker() { [[ -s "${OUT_DIR}/status/$1.done" ]]; }
mark_done() { date '+%Y-%m-%d %H:%M:%S' > "${OUT_DIR}/status/$1.done"; }

run_step() {
  local name="$1"
  shift
  if done_marker "$name"; then
    log "SKIP $name"
    return 0
  fi
  log "START $name"
  "$@"
  mark_done "$name"
  log "DONE $name"
}

corrected_bw_for_sample() {
  local sample="$1"
  local path
  path="$(ls "${FP_DIR}/atac_correct/${sample}/"*"_corrected.bw" 2>/dev/null | sort | head -1 || true)"
  [[ -n "${path}" ]] || { log "ERROR no corrected bigWig found for ${sample}"; exit 1; }
  printf '%s\n' "${path}"
}

check_inputs() {
  require_file "$GENOME"
  require_file "$PEAKS"
  require_file "$JASPAR2026"
  for sample in "${SAMPLES[@]}"; do
    require_file "${FP_DIR}/footprints/${sample}.footprints.bw"
    require_file "$(corrected_bw_for_sample "$sample")"
  done
  require_exec "$CALL_FOOTPRINTS"
  require_exec "$MOTIF_DISCOVERY"
  require_exec "$DIFF_FOOTPRINTS"
  require_exec "$PYTHON"
  require_exec "$STREME"
  require_exec "$TOMTOM"
}

write_versions() {
  {
    echo -e "tool\tversion"
    "$STREME" --version 2>&1 | awk 'NR==1{print "STREME\t"$0}'
    "$TOMTOM" --version 2>&1 | awk 'NR==1{print "Tomtom\t"$0}'
    "$PYTHON" -c 'import sys, fp_tools; print(f"python\t{sys.version.split()[0]}"); print(f"fp-tools\t{fp_tools.__version__}")'
    echo -e "JASPAR\t2026 CORE vertebrates non-redundant"
  } > "${OUT_DIR}/software_versions.tsv"
  {
    echo -e "analysis_step\tparameters"
    echo -e "candidate_generation\tcall-footprints --score footprint --output-bed --top-n ${TOP_N_PER_REPLICATE} --call-width ${CALL_WIDTH} --min-distance ${MIN_DISTANCE} on each corrected ENCODE replicate"
    echo -e "condition_candidates\tmerge replicate candidates by condition; retain top ${TOP_N_PER_CONDITION} non-overlapping centers per condition"
    echo -e "sequence_export\tmotif-discovery --method streme candidate-centered FASTA with flank +/-${FLANK} bp"
    if [[ -n "${STREME_NMOTIFS}" && "${STREME_NMOTIFS}" != "0" ]]; then
      echo -e "streme\tK562 primary vs HepG2 control and HepG2 primary vs K562 control; --dna --nmotifs ${STREME_NMOTIFS}"
    else
      echo -e "streme\tK562 primary vs HepG2 control and HepG2 primary vs K562 control; --dna --thresh ${STREME_THRESH} --patience ${STREME_PATIENCE}"
    fi
    echo -e "tomtom\tdiscovered motifs compared to JASPAR2026 CORE vertebrates non-redundant motifs"
    echo -e "diff_footprints\tsample-quantile normalization, replicate report auto, corrected aggregate signals, plot-aggregate ${DIFF_PLOT_AGGREGATE}, reuse-existing-results ${DIFF_REUSE_EXISTING}"
    echo -e "restricted_database\tJASPAR2026 with common erythroid/K562 and liver/HepG2 motif families removed for sensitivity rescue demo"
  } > "${OUT_DIR}/analysis_parameters.tsv"
}

convert_jaspar_for_tomtom() {
  "$PYTHON" - "$JASPAR2026" "${OUT_DIR}/motifs/jaspar2026_vertebrates_for_tomtom.meme" <<'PYJASPAR'
from pathlib import Path
import sys
from fp_tools.utils.motifs import MotifList

MotifList().from_file(sys.argv[1]).to_file(sys.argv[2], fmt="meme")
print(f"Wrote Tomtom-compatible MEME motif database to {Path(sys.argv[2])}")
PYJASPAR
}

call_candidates_for_sample() {
  local sample="$1"
  local corrected score_bw bed
  corrected="$(corrected_bw_for_sample "$sample")"
  score_bw="${OUT_DIR}/candidate_calls/${sample}.candidate_scores.bw"
  bed="${OUT_DIR}/candidate_calls/${sample}.candidate_footprints.bed"
  "$CALL_FOOTPRINTS" \
    --signal "$corrected" \
    --regions "$PEAKS" \
    --output "$score_bw" \
    --score footprint \
    --output-bed "$bed" \
    --top-n "$TOP_N_PER_REPLICATE" \
    --call-width "$CALL_WIDTH" \
    --min-distance "$MIN_DISTANCE"
}

merge_condition_candidates() {
  "$PYTHON" - "$OUT_DIR" "$TOP_N_PER_CONDITION" "$MIN_DISTANCE" <<'PYMERGE'
from pathlib import Path
import sys

out = Path(sys.argv[1])
top_n = int(sys.argv[2])
min_distance = int(sys.argv[3])
condition_map = {
    "K562": ["K562_rep1", "K562_rep2", "K562_rep3"],
    "HepG2": ["HepG2_rep1", "HepG2_rep2", "HepG2_rep3"],
}
for condition, samples in condition_map.items():
    rows = []
    for sample in samples:
        path = out / "candidate_calls" / f"{sample}.candidate_footprints.bed"
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                score = float(fields[4]) if len(fields) > 4 and fields[4] not in {"", "."} else 0.0
                center = int(fields[7]) if len(fields) > 7 else (int(fields[1]) + int(fields[2])) // 2
                rows.append((fields[0], int(fields[1]), int(fields[2]), score, center, sample))
    rows.sort(key=lambda row: (-row[3], row[0], row[4]))
    kept = []
    centers_by_chrom = {}
    for chrom, start, end, score, center, sample in rows:
        used = centers_by_chrom.setdefault(chrom, [])
        if any(abs(center - other) < min_distance for other in used):
            continue
        used.append(center)
        kept.append((chrom, start, end, f"{condition}_denovo_candidate_{len(kept) + 1}", score, ".", sample, center))
        if len(kept) >= top_n:
            break
    out_path = out / "candidate_calls" / f"{condition}.merged_candidates.bed"
    with out_path.open("w") as handle:
        handle.write("#chrom\tstart\tend\tname\tscore\tstrand\tsource_sample\tcenter\n")
        for chrom, start, end, name, score, strand, sample, center in kept:
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{score:.6g}\t{strand}\t{sample}\t{center}\n")
    print(f"Wrote {len(kept)} candidates to {out_path}")
PYMERGE
}

export_condition_fastas() {
  "$PYTHON" - "$OUT_DIR" "$GENOME" "$FLANK" <<'PYFASTA'
from pathlib import Path
import sys
from fp_tools.tools.motif_discovery import export_candidate_fasta

out = Path(sys.argv[1])
genome = Path(sys.argv[2])
flank = int(sys.argv[3])
for condition in ["K562", "HepG2"]:
    candidates = out / "candidate_calls" / f"{condition}.merged_candidates.bed"
    fasta = out / "candidate_fastas" / f"{condition}.candidate_sequences.fa"
    written = export_candidate_fasta(candidates, genome, fasta, flank=flank)
    print(f"Wrote {written} sequences to {fasta}")
PYFASTA
}

prepare_streme_plan() {
  local condition="$1"
  local control="$2"
  local outdir="${OUT_DIR}/motifs/${condition}_vs_${control}_streme"
  local extra_args=(--dna)
  if [[ -n "${STREME_NMOTIFS}" && "${STREME_NMOTIFS}" != "0" ]]; then
    extra_args+=(--nmotifs "$STREME_NMOTIFS")
  else
    extra_args+=(--thresh "$STREME_THRESH" --patience "$STREME_PATIENCE")
  fi
  extra_args+=(--n "${OUT_DIR}/candidate_fastas/${control}.candidate_sequences.fa")
  "$MOTIF_DISCOVERY" \
    --fasta "${OUT_DIR}/candidate_fastas/${condition}.candidate_sequences.fa" \
    --outdir "$outdir" \
    --method streme \
    --known-motifs "${OUT_DIR}/motifs/jaspar2026_vertebrates_for_tomtom.meme" \
    --extra-args "${extra_args[@]}"
}

run_streme_script() {
  local condition="$1"
  local control="$2"
  local outdir="${OUT_DIR}/motifs/${condition}_vs_${control}_streme"
  PATH="$(dirname "$STREME"):$(dirname "$TOMTOM"):$(dirname "$MOTIF_DISCOVERY"):$PATH" bash "${outdir}/run_motif_discovery.sh"
}

build_motif_sets() {
  "$PYTHON" - "$OUT_DIR" "$JASPAR2026" <<'PYSETS'
from pathlib import Path
import re
import sys
from fp_tools.utils.motifs import MotifList

out = Path(sys.argv[1])
jaspar = Path(sys.argv[2])
directions = [
    ("K562_denovo", out / "motifs" / "K562_vs_HepG2_streme" / "streme" / "streme.txt"),
    ("HepG2_denovo", out / "motifs" / "HepG2_vs_K562_streme" / "streme" / "streme.txt"),
]
merged = MotifList()
for label, path in directions:
    motifs = MotifList().from_file(str(path))
    for idx, motif in enumerate(motifs, start=1):
        motif.id = f"{label}_{idx}_{motif.id}"
        motif.name = f"{label}_{idx}"
        merged.append(motif)
denovo_path = out / "motifs" / "encode_k562_hepg2_denovo_streme.meme"
merged.to_file(str(denovo_path), fmt="meme")

exclude_re = re.compile(
    r"(GATA|TAL|KLF|NFE2|BACH|HNF|FOXA|CEBP|ONECUT|NR[0-9A-Z]*|PPAR|RXR|"
    r"JUN|FOS|ATF|MAF|MEF2|TEAD|AP-1|EBOX)",
    re.I,
)
restricted_lines = []
keep = True
with jaspar.open() as handle:
    for line in handle:
        if line.startswith(">"):
            keep = exclude_re.search(line) is None
        if keep:
            restricted_lines.append(line)
restricted_path = out / "motifs" / "jaspar2026_vertebrates_restricted_k562_hepg2_sensitivity.jaspar"
restricted_path.write_text("".join(restricted_lines))

full = MotifList().from_file(str(jaspar))
for motif in merged:
    full.append(motif)
full_plus_path = out / "motifs" / "jaspar2026_plus_denovo_streme.meme"
full.to_file(str(full_plus_path), fmt="meme")

restricted = MotifList().from_file(str(restricted_path))
for motif in merged:
    restricted.append(motif)
restricted_plus_path = out / "motifs" / "restricted_jaspar_plus_denovo_streme.meme"
restricted.to_file(str(restricted_plus_path), fmt="meme")

with (out / "motifs" / "motif_set_summary.tsv").open("w") as handle:
    handle.write("motif_set\tn_motifs\tpath\n")
    handle.write(f"de_novo_only\t{len(merged)}\t{denovo_path}\n")
    handle.write(f"jaspar2026_full\t{len(full) - len(merged)}\t{jaspar}\n")
    handle.write(f"jaspar2026_plus_denovo\t{len(full)}\t{full_plus_path}\n")
    handle.write(f"jaspar2026_restricted\t{len(restricted) - len(merged)}\t{restricted_path}\n")
    handle.write(f"jaspar2026_restricted_plus_denovo\t{len(restricted)}\t{restricted_plus_path}\n")
print("Wrote motif sets under", out / "motifs")
PYSETS
}

run_diff_for_set() {
  local label="$1"
  local motif_file="$2"
  local outdir="${OUT_DIR}/diff_footprints/${label}"
  local extra_args=()
  mkdir -p "$outdir"
  if [[ "${DIFF_REUSE_EXISTING}" == "1" && -s "${outdir}/diff_footprints_results.txt" ]]; then
    extra_args+=(--reuse-existing-results)
  fi
  "$DIFF_FOOTPRINTS" \
    --motifs "$motif_file" \
    --signals \
      "${FP_DIR}/footprints/K562_rep1.footprints.bw" \
      "${FP_DIR}/footprints/K562_rep2.footprints.bw" \
      "${FP_DIR}/footprints/K562_rep3.footprints.bw" \
      "${FP_DIR}/footprints/HepG2_rep1.footprints.bw" \
      "${FP_DIR}/footprints/HepG2_rep2.footprints.bw" \
      "${FP_DIR}/footprints/HepG2_rep3.footprints.bw" \
    --genome "$GENOME" \
    --peaks "$PEAKS" \
    --outdir "$outdir" \
    --prefix diff_footprints \
    --cond-names K562 K562 K562 HepG2 HepG2 HepG2 \
    --normalization sample-quantile \
    --replicate-report auto \
    --aggregate-signals \
      "$(corrected_bw_for_sample K562_rep1)" \
      "$(corrected_bw_for_sample K562_rep2)" \
      "$(corrected_bw_for_sample K562_rep3)" \
      "$(corrected_bw_for_sample HepG2_rep1)" \
      "$(corrected_bw_for_sample HepG2_rep2)" \
      "$(corrected_bw_for_sample HepG2_rep3)" \
    --aggregate-normalization none \
    --aggregate-site-set bound \
    --plot-aggregate "$DIFF_PLOT_AGGREGATE" \
    --aggregate-flank 100 \
    --skip-excel \
    --cores "$THREADS" \
    "${extra_args[@]}"
}

main() {
  log "ENCODE K562 vs HepG2 de novo motif validation started"
  check_inputs
  write_versions
  for sample in "${SAMPLES[@]}"; do
    run_step "candidates.${sample}" call_candidates_for_sample "$sample"
  done
  run_step "candidates.merge_conditions" merge_condition_candidates
  run_step "candidate_fastas.export" export_condition_fastas
  run_step "motifs.convert_jaspar_for_tomtom" convert_jaspar_for_tomtom
  run_step "motif_plan.K562_vs_HepG2" prepare_streme_plan K562 HepG2
  run_step "motif_plan.HepG2_vs_K562" prepare_streme_plan HepG2 K562
  run_step "streme.K562_vs_HepG2" run_streme_script K562 HepG2
  run_step "streme.HepG2_vs_K562" run_streme_script HepG2 K562
  run_step "motif_sets.build" build_motif_sets
  run_step "diff.denovo_only" run_diff_for_set denovo_only "${OUT_DIR}/motifs/encode_k562_hepg2_denovo_streme.meme"
  run_step "diff.jaspar_plus_denovo" run_diff_for_set jaspar2026_plus_denovo "${OUT_DIR}/motifs/jaspar2026_plus_denovo_streme.meme"
  run_step "diff.restricted_jaspar" run_diff_for_set restricted_jaspar "${OUT_DIR}/motifs/jaspar2026_vertebrates_restricted_k562_hepg2_sensitivity.jaspar"
  run_step "diff.restricted_plus_denovo" run_diff_for_set restricted_jaspar_plus_denovo "${OUT_DIR}/motifs/restricted_jaspar_plus_denovo_streme.meme"
  log "ENCODE K562 vs HepG2 de novo motif validation finished: ${OUT_DIR}"
}

main "$@"
