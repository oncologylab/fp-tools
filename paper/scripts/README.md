# Paper Scripts

This directory will hold reproducible figure and table builders for the planned BioMedInformatics manuscript.

Outputs should be generated from benchmark tables, not edited by hand:

- `paper/figures/*.pdf` and `*.svg` for manuscript submission.
- `paper/figures/*.png` for GitHub previews.
- `paper/tables/*.tsv` and `*.csv` for source data.

## Available Builders

- `plot_benchmark_panels.py`: render AUROC/AUPRC/recall/calibration benchmark summaries from metrics TSV files.
- `plot_multiscale_npz.py`: render PDF/SVG/PNG multi-panel summaries from `score-footprints --output-multiscale-npz` sidecars.
- `prepare_biomedinformatics_template.py`: create a manuscript workspace from the BioMedInformatics/MDPI templates.
