# Nutrient Stress Project Template

This folder contains a portable shell template for running a nutrient-stress
ATAC footprinting project against `10_FBS_Ctrl`.

On a new server, copy `run_ctrl_vs_10fbs.sh` into the project output folder,
edit `RAW`, `PROJECT`, and `REF_ROOT`, then run:

```bash
CHECK_ONLY=1 bash run_ctrl_vs_10fbs.sh
CORES=16 bash run_ctrl_vs_10fbs.sh
```

The script prepares clean `metadata/samples.tsv` and
`metadata/comparisons.tsv` from an `ATAC_Nutrients_hg38_*.txt` table by
removing TGFB rows, keeping high-confidence ATAC Ctrl rows, normalizing condition
names, and comparing every non-control condition against `10_FBS_Ctrl`.
