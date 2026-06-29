#!/usr/bin/env bash
set -euo pipefail

# Portable nutrient-stress fp-tools workflow.
#
# Data preparation on a new server:
# 1. Install fp-tools:
#      python -m pip install "fp-tools-bio==${FP_TOOLS_VERSION:-0.1.12}"
# 2. Put the raw project folder at RAW, or edit RAW below. The raw folder should contain:
#      ATAC_Nutrients_hg38_*.txt
#      <sample_id>/<sample_id>.hg38.filtered.bam
#      <sample_id>/<sample_id>.hg38.rp10m.narrowpeaks.bed
# 3. Put references at REF_ROOT, or edit REF_ROOT below. Required files are:
#      GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
#      GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.fai
#      ENCFF419RSJ_hg38.bed
# 4. The script prepares clean metadata automatically:
#      metadata/samples.tsv with columns: sample, condition, bam, peaks
#      metadata/comparisons.tsv with columns: comparison, cond1, cond2
#    It keeps ATAC high-confidence Ctrl rows, removes TGFB rows, renames 10_FBS
#    to 10_FBS_Ctrl, and compares every non-control condition against 10_FBS_Ctrl.
# 5. To use manually curated metadata instead, place those two TSVs under
#    metadata/ and run with PREPARE_METADATA=0.
#
# Run:
#      CHECK_ONLY=1 bash run_ctrl_vs_10fbs.sh
#      CORES=16 bash run_ctrl_vs_10fbs.sh

ROOT="${ROOT:-/path/to/fp-tools-or-workspace}"
RAW="${RAW:-$ROOT/data/public/raw/nutrient_project}"
PROJECT="${PROJECT:-$ROOT/data/public/processed/nutrient_project_ctrl_vs_10fbs}"
REF_ROOT="${REF_ROOT:-$RAW/references}"
CORES="${CORES:-$(nproc)}"
VERSION="${FP_TOOLS_VERSION:-0.1.12}"
MOTIF_DB="${MOTIF_DB:-jaspar2026_vertebrates}"
PREPARE_METADATA="${PREPARE_METADATA:-1}"

if [[ -d "$ROOT/.venv/bin" ]]; then
  export PATH="$ROOT/.venv/bin:$PATH"
fi

GENOME="${GENOME:-$REF_ROOT/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta}"
BLACKLIST="${BLACKLIST:-$REF_ROOT/ENCFF419RSJ_hg38.bed}"
SAMPLES="$PROJECT/metadata/samples.tsv"
COMPARISONS="$PROJECT/metadata/comparisons.tsv"
PEAKS="$PROJECT/peaks/merged_peaks_filtered.bed"

require_file() {
  if [[ ! -s "$1" ]]; then
    echo "Missing required file: $1" >&2
    exit 1
  fi
}

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing command: $1" >&2
    echo "Install fp-tools with: python -m pip install fp-tools-bio==$VERSION" >&2
    exit 1
  fi
}

require_command python
for cmd in atac-correct normalize-bigwig call-footprints match-motifs diff-footprints review-multi-comparisons; do
  require_command "$cmd"
done
require_file "$GENOME"
require_file "$GENOME.fai"
require_file "$BLACKLIST"

python - <<PY_VERSION
import fp_tools
expected = "$VERSION"
found = fp_tools.__version__
if found != expected:
    raise SystemExit(f"Expected fp-tools-bio {expected}, but found {found}. Install with: python -m pip install fp-tools-bio=={expected}")
print(f"Using fp-tools-bio {found}")
PY_VERSION

mkdir -p "$PROJECT/metadata"

if [[ "$PREPARE_METADATA" == "1" ]]; then
  python - <<PY_METADATA
from pathlib import Path
import pandas as pd
import re

raw = Path("$RAW")
project = Path("$PROJECT")
meta_files = sorted(raw.glob("ATAC_Nutrients_hg38_*.txt"))
if not meta_files:
    raise SystemExit(f"No ATAC_Nutrients_hg38_*.txt file found in {raw}")
df = pd.read_csv(meta_files[0], sep="\t")
required = {"ID", "Antibody", "Sample", "Condition", "Confidence"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"Metadata file is missing columns: {sorted(missing)}")

def clean_condition(value):
    value = str(value).strip()
    value = re.sub(r"^DMEM_", "", value)
    value = value.replace("uM", "")
    value = re.sub(r"_+", "_", value).strip("_")
    if value == "10_FBS":
        return "10_FBS_Ctrl"
    if value.startswith("10_FBS_"):
        return value[len("10_FBS_"):]
    return value

keep = (
    df["Antibody"].astype(str).eq("ATAC")
    & df["Condition"].astype(str).str.upper().ne("TGFB")
    & df["Confidence"].astype(str).str.lower().eq("high")
)
df = df.loc[keep].copy()
df["condition"] = df["Sample"].map(clean_condition)
rows, errors = [], []
for _, row in df.iterrows():
    sample = str(row["ID"]).strip()
    sample_dir = raw / sample
    bam = sample_dir / f"{sample}.hg38.filtered.bam"
    peaks = sample_dir / f"{sample}.hg38.rp10m.narrowpeaks.bed"
    if not bam.exists():
        errors.append(str(bam))
    if not peaks.exists():
        errors.append(str(peaks))
    rows.append({"sample": sample, "condition": row["condition"], "bam": str(bam), "peaks": str(peaks)})
if errors:
    raise SystemExit("Missing input files:\n" + "\n".join(errors[:20]))

samples = pd.DataFrame(rows).sort_values(["condition", "sample"])
if "10_FBS_Ctrl" not in set(samples["condition"]):
    raise SystemExit("No 10_FBS_Ctrl control condition after cleaning metadata")
comparisons = pd.DataFrame([
    {"comparison": f"{cond}_vs_10_FBS_Ctrl", "cond1": cond, "cond2": "10_FBS_Ctrl"}
    for cond in sorted(c for c in samples["condition"].unique() if c != "10_FBS_Ctrl")
])
metadata = project / "metadata"
metadata.mkdir(parents=True, exist_ok=True)
samples.to_csv(metadata / "samples.tsv", sep="\t", index=False)
comparisons.to_csv(metadata / "comparisons.tsv", sep="\t", index=False)
print(f"Wrote {metadata / 'samples.tsv'} with {len(samples)} samples")
print(f"Wrote {metadata / 'comparisons.tsv'} with {len(comparisons)} comparisons")
PY_METADATA
fi

require_file "$SAMPLES"
require_file "$COMPARISONS"

python - <<PY_CHECK
from pathlib import Path
import pandas as pd
samples = pd.read_csv("$SAMPLES", sep="\t")
comparisons = pd.read_csv("$COMPARISONS", sep="\t")
for name, cols, df in [
    ("samples.tsv", {"sample", "condition", "bam", "peaks"}, samples),
    ("comparisons.tsv", {"comparison", "cond1", "cond2"}, comparisons),
]:
    missing = cols - set(df.columns)
    if missing:
        raise SystemExit(f"{name} is missing columns: {sorted(missing)}")
for col in ["bam", "peaks"]:
    missing = [p for p in samples[col].astype(str) if not Path(p).exists()]
    if missing:
        raise SystemExit(f"Missing {col} files:\n" + "\n".join(missing[:20]))
conds = set(samples["condition"].astype(str))
bad = sorted((set(comparisons["cond1"].astype(str)) | set(comparisons["cond2"].astype(str))) - conds)
if bad:
    raise SystemExit(f"Comparison conditions not present in samples.tsv: {bad}")
print(f"Input check passed: {len(samples)} samples, {len(comparisons)} comparisons, {len(conds)} conditions")
PY_CHECK

if [[ "${CHECK_ONLY:-0}" == "1" ]]; then
  echo "CHECK_ONLY=1: metadata and input checks passed; workflow commands were not run."
  exit 0
fi

COMMON=(--sample-table "$SAMPLES" --outdir "$PROJECT")

atac-correct "${COMMON[@]}" --genome "$GENOME" --blacklist "$BLACKLIST" --cores "$CORES"
normalize-bigwig "${COMMON[@]}" --background "$PEAKS" --method background-scale --stat q95 --target median --workers "$CORES"
call-footprints "${COMMON[@]}" --regions "$PEAKS" --cores "$CORES"
match-motifs "${COMMON[@]}" --genome "$GENOME" --peaks "$PEAKS" --motif-db "$MOTIF_DB" --cores "$CORES"
diff-footprints "${COMMON[@]}" --comparison-table "$COMPARISONS" --genome "$GENOME" --peaks "$PEAKS" --motif-db "$MOTIF_DB" --cores "$CORES"
review-multi-comparisons --outdir "$PROJECT" --title "Nutrient stress vs 10_FBS_Ctrl"

echo "Done. Project output: $PROJECT"
