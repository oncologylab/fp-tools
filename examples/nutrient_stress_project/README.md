# Nutrient Stress Project Template

Use this shell template to compare nutrient-stress ATAC-seq samples with
`10_FBS_Ctrl`. It starts from aligned BAM files and peak BEDs.

1. Install `fp-tools-bio==0.2.3` in the Python environment used by the script.
2. Copy `run_ctrl_vs_10fbs.sh` to your project folder and set `RAW`, `PROJECT`, and `REF_ROOT` to your input, output, and reference directories.
3. Provide the coordinate-sorted BAM/BAI and peak BED files, plus the hg38 FASTA, FASTA index, and blacklist named in the script's setup comments.
4. Create `PROJECT/metadata/samples.tsv` with `sample`, `condition`, `bam`, and `peaks` columns. Use one row per biological sample and file paths that are valid from the launch directory.
5. Create `PROJECT/metadata/comparisons.tsv` with `comparison`, `cond1`, and `cond2` columns. Put each nutrient condition in `cond1` and `10_FBS_Ctrl` in `cond2`.

For the original `ATAC_Nutrients_hg38_*.txt` metadata, the template's design
excludes TGFB samples and names the 10% FBS control `10_FBS_Ctrl`. Prepare these
tables before running the script; it does not create them automatically.

Check the required columns, condition labels, and input file paths, then run
the analysis:

```bash
CHECK_ONLY=1 bash run_ctrl_vs_10fbs.sh
CORES=16 bash run_ctrl_vs_10fbs.sh
```

Results are written under `PROJECT/samples/`, `PROJECT/comparisons/`, and
`PROJECT/reports/`. Open the combined comparison review under `PROJECT/reports/`
to explore the results.
