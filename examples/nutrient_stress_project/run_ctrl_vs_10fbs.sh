#!/usr/bin/env bash
set -euo pipefail

# Portable nutrient-stress fp-tools workflow.
#
# Before running on a new server:
# 1. Install fp-tools:
#      python -m pip install fp-tools-bio==0.1.18
#
# 2. Prepare raw data under RAW:
#      RAW/ATAC_Nutrients_hg38_*.txt
#      RAW/<sample_id>/<sample_id>.hg38.filtered.bam
#      RAW/<sample_id>/<sample_id>.hg38.rp10m.narrowpeaks.bed
#
# 3. Prepare references under REF_ROOT:
#      GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
#      GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.fai
#      ENCFF419RSJ_hg38.bed
#
# 4. Prepare metadata files:
#      PROJECT/metadata/samples.tsv
#        columns: sample, condition, bam, peaks
#        use one row per sample; remove TGFB rows; rename 10_FBS to 10_FBS_Ctrl.
#
#      PROJECT/metadata/comparisons.tsv
#        columns: comparison, cond1, cond2
#        compare each nutrient condition against 10_FBS_Ctrl.
#
# 5. Sanity-check inputs without running analysis:
#      CHECK_ONLY=1 bash run_ctrl_vs_10fbs.sh
#
# 6. Run the full workflow:
#      CORES=16 bash run_ctrl_vs_10fbs.sh

ROOT="${ROOT:-/path/to/fp-tools-or-workspace}"
RAW="${RAW:-$ROOT/data/public/raw/nutrient_project}"
PROJECT="${PROJECT:-$ROOT/data/public/processed/nutrient_project_ctrl_vs_10fbs}"
REF_ROOT="${REF_ROOT:-$RAW/references}"
CORES="${CORES:-$(nproc)}"
VERSION="${FP_TOOLS_VERSION:-0.1.18}"
MOTIF_DB="${MOTIF_DB:-jaspar2026_vertebrates}"

if [[ -d "$ROOT/.venv/bin" ]]; then
  export PATH="$ROOT/.venv/bin:$PATH"
fi

GENOME="${GENOME:-$REF_ROOT/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta}"
BLACKLIST="${BLACKLIST:-$REF_ROOT/ENCFF419RSJ_hg38.bed}"
SAMPLES="$PROJECT/metadata/samples.tsv"
COMPARISONS="$PROJECT/metadata/comparisons.tsv"
PEAKS="$PROJECT/peaks/merged_peaks_filtered.bed"

python - <<PY
import fp_tools
expected = "$VERSION"
found = fp_tools.__version__
if found != expected:
    raise SystemExit(f"Expected fp-tools-bio {expected}, but found {found}. Install with: python -m pip install fp-tools-bio=={expected}")
print(f"Using fp-tools-bio {found}")
PY

for f in "$GENOME" "$GENOME.fai" "$BLACKLIST" "$SAMPLES" "$COMPARISONS"; do
  if [[ ! -s "$f" ]]; then
    echo "Missing required file: $f" >&2
    exit 1
  fi
done

python - <<PY
from pathlib import Path
import pandas as pd

samples = pd.read_csv("$SAMPLES", sep="\t")
comparisons = pd.read_csv("$COMPARISONS", sep="\t")
missing = {"sample", "condition", "bam", "peaks"} - set(samples.columns)
if missing:
    raise SystemExit(f"samples.tsv missing columns: {sorted(missing)}")
missing = {"comparison", "cond1", "cond2"} - set(comparisons.columns)
if missing:
    raise SystemExit(f"comparisons.tsv missing columns: {sorted(missing)}")
for col in ("bam", "peaks"):
    bad = [p for p in samples[col].astype(str) if not Path(p).exists()]
    if bad:
        raise SystemExit(f"Missing {col} files:\\n" + "\\n".join(bad[:20]))
conds = set(samples["condition"].astype(str))
bad = sorted((set(comparisons["cond1"].astype(str)) | set(comparisons["cond2"].astype(str))) - conds)
if bad:
    raise SystemExit(f"Comparison conditions not present in samples.tsv: {bad}")
print(f"Input check passed: {len(samples)} samples, {len(comparisons)} comparisons")
PY

if [[ "${CHECK_ONLY:-0}" == "1" ]]; then
  echo "CHECK_ONLY=1: inputs passed; workflow commands were not run."
  exit 0
fi

COMMON=(--sample-table "$SAMPLES" --outdir "$PROJECT")

# 1. Bias-correct all BAM files and merge all sample peak BEDs.
atac-correct \
  "${COMMON[@]}" \
  --genome "$GENOME" \
  --blacklist "$BLACKLIST" \
  --cores "$CORES"

# 2. q95-scale corrected cut-site tracks over the shared peak universe.
normalize-bigwig \
  "${COMMON[@]}" \
  --background "$PEAKS" \
  --method background-scale \
  --stat q95 \
  --target median \
  --workers "$CORES"

# 3. Call footprint score tracks.
call-footprints \
  "${COMMON[@]}" \
  --regions "$PEAKS" \
  --cores "$CORES"

# 4. Match JASPAR 2026 motifs for each sample.
match-motifs \
  "${COMMON[@]}" \
  --genome "$GENOME" \
  --peaks "$PEAKS" \
  --motif-db "$MOTIF_DB" \
  --cores "$CORES"

# 5. Compare each nutrient-stress condition against 10_FBS_Ctrl.
diff-footprints \
  "${COMMON[@]}" \
  --comparison-table "$COMPARISONS" \
  --genome "$GENOME" \
  --peaks "$PEAKS" \
  --motif-db "$MOTIF_DB" \
  --cores "$CORES"

# 6. Build one review page for all differential-footprint reports.
review-multi-comparisons \
  --outdir "$PROJECT" \
  --title "Nutrient stress vs 10_FBS_Ctrl"

echo "Done. Project output: $PROJECT"
