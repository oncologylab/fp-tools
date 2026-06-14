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
BASE_DIR="${ROOT_DIR}/data/public/processed/buenrostro_atac_replicates"
PEAKS="${BASE_DIR}/peaks/merged_peaks.bed"
FP_DIR="${BASE_DIR}/fp_tools"
OUT_DIR="${FP_DIR}/denovo_motif_validation"
FP_TOOLS_ENV="${FP_TOOLS_ENV:-${ROOT_DIR}/.venv}"
CALL_FOOTPRINTS="${FP_TOOLS_ENV}/bin/call-footprints"
MOTIF_DISCOVERY="${FP_TOOLS_ENV}/bin/motif-discovery"
DIFF_FOOTPRINTS="${FP_TOOLS_ENV}/bin/diff-footprints"
PYTHON="${FP_TOOLS_ENV}/bin/python"
STREME="${STREME:-streme}"
TOMTOM="${TOMTOM:-tomtom}"
TOP_N_PER_REPLICATE="${TOP_N_PER_REPLICATE:-5000}"
TOP_N_PER_CONDITION="${TOP_N_PER_CONDITION:-8000}"
CALL_WIDTH="${CALL_WIDTH:-50}"
MIN_DISTANCE="${MIN_DISTANCE:-40}"
FLANK="${FLANK:-75}"
STREME_NMOTIFS="${STREME_NMOTIFS:-8}"

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
require_file() { [[ -s "$1" ]] || { log "ERROR missing file: $1"; exit 1; }; }
require_exec() { command -v "$1" >/dev/null 2>&1 || { log "ERROR missing executable on PATH: $1"; exit 1; }; }
done_marker() { [[ -s "${OUT_DIR}/status/$1.done" ]]; }
mark_done() { date '+%Y-%m-%d %H:%M:%S' > "${OUT_DIR}/status/$1.done"; }
run_step() {
  local name="$1"; shift
  if done_marker "$name"; then
    log "SKIP $name"
    return 0
  fi
  log "START $name"
  "$@"
  mark_done "$name"
  log "DONE $name"
}

SAMPLES=(Bcell_rep1 Bcell_rep2 Tcell_rep1 Tcell_rep2)
mkdir -p "${OUT_DIR}/status" "${OUT_DIR}/candidate_calls" "${OUT_DIR}/candidate_fastas" "${OUT_DIR}/motifs" "${OUT_DIR}/diff_footprints"

check_inputs() {
  require_file "$GENOME"
  require_file "$PEAKS"
  require_file "$JASPAR2026"
  for sample in "${SAMPLES[@]}"; do
    require_file "${FP_DIR}/footprints/${sample}.footprints.bw"
    require_file "${FP_DIR}/atac_correct/${sample}/${sample}.filtered_corrected.bw"
  done
  [[ -x "$CALL_FOOTPRINTS" ]] || { log "ERROR missing $CALL_FOOTPRINTS"; exit 1; }
  [[ -x "$MOTIF_DISCOVERY" ]] || { log "ERROR missing $MOTIF_DISCOVERY"; exit 1; }
  [[ -x "$DIFF_FOOTPRINTS" ]] || { log "ERROR missing $DIFF_FOOTPRINTS"; exit 1; }
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
    echo -e "candidate_generation\tcall-footprints --score footprint --output-bed --top-n ${TOP_N_PER_REPLICATE} --call-width ${CALL_WIDTH} --min-distance ${MIN_DISTANCE} on each corrected cut-site replicate"
    echo -e "condition_candidates\tmerge replicate candidates by condition; retain top ${TOP_N_PER_CONDITION} non-overlapping centers per condition"
    echo -e "sequence_export\tmotif-discovery --method streme candidate-centered FASTA with flank +/-${FLANK} bp"
    echo -e "streme\tBcell primary vs Tcell control and Tcell primary vs Bcell control; --dna --nmotifs ${STREME_NMOTIFS}"
    echo -e "tomtom\tdiscovered motifs compared to JASPAR2026 CORE vertebrates non-redundant motifs"
    echo -e "diff_footprints\tsample-quantile normalization, replicate report auto, aggregate signals, plot-aggregate sig"
    echo -e "restricted_database\tJASPAR2026 with common immune/AP-1/IRF/ETS/T-cell/B-cell marker motif families removed for sensitivity rescue demo"
  } > "${OUT_DIR}/analysis_parameters.tsv"
}

convert_jaspar_for_tomtom() {
  "$PYTHON" - "$JASPAR2026" "${OUT_DIR}/motifs/jaspar2026_vertebrates_for_tomtom.meme" <<'PYJASPAR'
from pathlib import Path
import sys
from fp_tools.utils.motifs import MotifList
source = Path(sys.argv[1])
output = Path(sys.argv[2])
MotifList().from_file(str(source)).to_file(str(output), fmt='meme')
print(f'Wrote Tomtom-compatible MEME motif database to {output}')
PYJASPAR
}

call_candidates_for_sample() {
  local sample="$1"
  local corrected="${FP_DIR}/atac_correct/${sample}/${sample}.filtered_corrected.bw"
  local score_bw="${OUT_DIR}/candidate_calls/${sample}.candidate_scores.bw"
  local bed="${OUT_DIR}/candidate_calls/${sample}.candidate_footprints.bed"
  "$CALL_FOOTPRINTS" --signal "$corrected" --regions "$PEAKS" --output "$score_bw" --score footprint --output-bed "$bed" --top-n "$TOP_N_PER_REPLICATE" --call-width "$CALL_WIDTH" --min-distance "$MIN_DISTANCE"
}

merge_condition_candidates() {
  "$PYTHON" - "$OUT_DIR" "$TOP_N_PER_CONDITION" "$MIN_DISTANCE" <<'PYMERGE'
from pathlib import Path
import sys
out = Path(sys.argv[1]); top_n = int(sys.argv[2]); min_distance = int(sys.argv[3])
condition_map = {'Bcell': ['Bcell_rep1', 'Bcell_rep2'], 'Tcell': ['Tcell_rep1', 'Tcell_rep2']}
for condition, samples in condition_map.items():
    rows = []
    for sample in samples:
        path = out / 'candidate_calls' / f'{sample}.candidate_footprints.bed'
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                score = float(fields[4]) if len(fields) > 4 and fields[4] not in {'', '.'} else 0.0
                center = int(fields[7]) if len(fields) > 7 else (int(fields[1]) + int(fields[2])) // 2
                rows.append((fields[0], int(fields[1]), int(fields[2]), score, center, sample))
    rows.sort(key=lambda row: (-row[3], row[0], row[4]))
    kept = []; centers_by_chrom = {}
    for chrom, start, end, score, center, sample in rows:
        used = centers_by_chrom.setdefault(chrom, [])
        if any(abs(center - other) < min_distance for other in used):
            continue
        used.append(center)
        kept.append((chrom, start, end, f'{condition}_denovo_candidate_{len(kept)+1}', score, '.', sample, center))
        if len(kept) >= top_n:
            break
    out_path = out / 'candidate_calls' / f'{condition}.merged_candidates.bed'
    with out_path.open('w') as handle:
        handle.write('#chrom\tstart\tend\tname\tscore\tstrand\tsource_sample\tcenter\n')
        for row in kept:
            chrom, start, end, name, score, strand, sample, center = row
            handle.write(f'{chrom}\t{start}\t{end}\t{name}\t{score:.6g}\t{strand}\t{sample}\t{center}\n')
    print(f'Wrote {len(kept)} candidates to {out_path}')
PYMERGE
}

export_condition_fastas() {
  "$PYTHON" - "$OUT_DIR" "$GENOME" "$FLANK" <<'PYFASTA'
from pathlib import Path
import sys
from fp_tools.tools.motif_discovery import export_candidate_fasta
out = Path(sys.argv[1]); genome = Path(sys.argv[2]); flank = int(sys.argv[3])
for condition in ['Bcell', 'Tcell']:
    candidates = out / 'candidate_calls' / f'{condition}.merged_candidates.bed'
    fasta = out / 'candidate_fastas' / f'{condition}.candidate_sequences.fa'
    written = export_candidate_fasta(candidates, genome, fasta, flank=flank)
    print(f'Wrote {written} sequences to {fasta}')
PYFASTA
}

prepare_streme_plan() {
  local condition="$1"
  local control="$2"
  local outdir="${OUT_DIR}/motifs/${condition}_vs_${control}_streme"
  "$MOTIF_DISCOVERY" --fasta "${OUT_DIR}/candidate_fastas/${condition}.candidate_sequences.fa" --outdir "$outdir" --method streme --known-motifs "${OUT_DIR}/motifs/jaspar2026_vertebrates_for_tomtom.meme" --extra-args --dna --nmotifs "$STREME_NMOTIFS" --n "${OUT_DIR}/candidate_fastas/${control}.candidate_sequences.fa"
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
import re, sys
from fp_tools.utils.motifs import MotifList
out = Path(sys.argv[1]); jaspar = Path(sys.argv[2])
streme_paths = [out/'motifs'/'Bcell_vs_Tcell_streme'/'streme'/'streme.txt', out/'motifs'/'Tcell_vs_Bcell_streme'/'streme'/'streme.txt']
merged = MotifList()
for label, path in zip(['Bcell_denovo', 'Tcell_denovo'], streme_paths):
    motifs = MotifList().from_file(str(path))
    for idx, motif in enumerate(motifs, start=1):
        motif.id = f'{label}_{idx}_{motif.id}'
        motif.name = f'{label}_{idx}'
        merged.append(motif)
merged.to_file(str(out/'motifs'/'buenrostro_denovo_streme.meme'), fmt='meme')
exclude_re = re.compile(r'(BACH|BATF|FOS|JUN|IRF|ATF|TCF|LEF|GATA|RUNX|ETS|ELF|ELK|ETV|BCL11|EBF|PAX|POU2F2|SPIB|SPI1|ROR|TBX21|EOMES|STAT|TOX|ZBTB7B)', re.I)
restricted_lines = []; keep = True
with jaspar.open() as handle:
    for line in handle:
        if line.startswith('>'):
            keep = exclude_re.search(line) is None
        if keep:
            restricted_lines.append(line)
restricted_path = out/'motifs'/'jaspar2026_vertebrates_restricted_sensitivity.jaspar'
restricted_path.write_text(''.join(restricted_lines))
full = MotifList().from_file(str(jaspar))
for motif in merged: full.append(motif)
full.to_file(str(out/'motifs'/'jaspar2026_plus_denovo_streme.meme'), fmt='meme')
restricted = MotifList().from_file(str(restricted_path))
for motif in merged: restricted.append(motif)
restricted.to_file(str(out/'motifs'/'restricted_jaspar_plus_denovo_streme.meme'), fmt='meme')
with (out/'motifs'/'motif_set_summary.tsv').open('w') as handle:
    handle.write('motif_set\tn_motifs\tpath\n')
    handle.write(f'de_novo_only\t{len(merged)}\t{out/"motifs"/"buenrostro_denovo_streme.meme"}\n')
    handle.write(f'jaspar2026_full\t{len(full)-len(merged)}\t{jaspar}\n')
    handle.write(f'jaspar2026_plus_denovo\t{len(full)}\t{out/"motifs"/"jaspar2026_plus_denovo_streme.meme"}\n')
    handle.write(f'jaspar2026_restricted\t{len(restricted)-len(merged)}\t{restricted_path}\n')
    handle.write(f'jaspar2026_restricted_plus_denovo\t{len(restricted)}\t{out/"motifs"/"restricted_jaspar_plus_denovo_streme.meme"}\n')
print('Wrote motif sets under', out/'motifs')
PYSETS
}

run_diff_for_set() {
  local label="$1"
  local motif_file="$2"
  local outdir="${OUT_DIR}/diff_footprints/${label}"
  mkdir -p "$outdir"
  "$DIFF_FOOTPRINTS" --motifs "$motif_file" --signals "${FP_DIR}/footprints/Bcell_rep1.footprints.bw" "${FP_DIR}/footprints/Bcell_rep2.footprints.bw" "${FP_DIR}/footprints/Tcell_rep1.footprints.bw" "${FP_DIR}/footprints/Tcell_rep2.footprints.bw" --genome "$GENOME" --peaks "$PEAKS" --outdir "$outdir" --cond-names Bcell Bcell Tcell Tcell --normalization sample-quantile --replicate-report auto --aggregate-signals "${FP_DIR}/atac_correct/Bcell_rep1/Bcell_rep1.filtered_corrected.bw" "${FP_DIR}/atac_correct/Bcell_rep2/Bcell_rep2.filtered_corrected.bw" "${FP_DIR}/atac_correct/Tcell_rep1/Tcell_rep1.filtered_corrected.bw" "${FP_DIR}/atac_correct/Tcell_rep2/Tcell_rep2.filtered_corrected.bw" --plot-aggregate sig --aggregate-flank 100 --skip-excel
}

main() {
  check_inputs
  write_versions
  for sample in "${SAMPLES[@]}"; do run_step "candidates.${sample}" call_candidates_for_sample "$sample"; done
  run_step "candidates.merge_conditions" merge_condition_candidates
  run_step "candidate_fastas.export" export_condition_fastas
  run_step "motifs.convert_jaspar_for_tomtom" convert_jaspar_for_tomtom
  run_step "motif_plan.Bcell_vs_Tcell" prepare_streme_plan Bcell Tcell
  run_step "motif_plan.Tcell_vs_Bcell" prepare_streme_plan Tcell Bcell
  run_step "streme.Bcell_vs_Tcell" run_streme_script Bcell Tcell
  run_step "streme.Tcell_vs_Bcell" run_streme_script Tcell Bcell
  run_step "motif_sets.build" build_motif_sets
  run_step "diff.denovo_only" run_diff_for_set denovo_only "${OUT_DIR}/motifs/buenrostro_denovo_streme.meme"
  run_step "diff.jaspar_plus_denovo" run_diff_for_set jaspar2026_plus_denovo "${OUT_DIR}/motifs/jaspar2026_plus_denovo_streme.meme"
  run_step "diff.restricted_jaspar" run_diff_for_set restricted_jaspar "${OUT_DIR}/motifs/jaspar2026_vertebrates_restricted_sensitivity.jaspar"
  run_step "diff.restricted_plus_denovo" run_diff_for_set restricted_jaspar_plus_denovo "${OUT_DIR}/motifs/restricted_jaspar_plus_denovo_streme.meme"
  log "De novo motif validation finished: ${OUT_DIR}"
}

main "$@"
