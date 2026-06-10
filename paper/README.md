# Paper Workspace

This directory holds the BioMedInformatics manuscript assets and reproducible paper-output scaffolding from `DEV_PLAN.md`.

## Local templates

- `MDPI_template_ACS.zip`: MDPI LaTeX template archive. This is the preferred source for the BioMedInformatics manuscript scaffold.
- `biomedinformatics-template.dot`: Microsoft Word template fallback.
- `BioMedInformatics`: local copy of journal instructions for authors.

## Manuscript and generated outputs

- `BioMedInformatics_manuscript/`: tracked BioMedInformatics LaTeX manuscript workspace, including `main.tex`, selected figures, availability notes, and the preview PDF.
- `figures/`: generated PDF/SVG/PNG manuscript figures outside the journal workspace. Large generated figures are ignored by default.
- `tables/`: generated TSV/CSV source tables. Large generated tables are ignored by default.
- `scripts/`: helper scripts for figures and manuscript preparation.

## Prepare or refresh the LaTeX manuscript scaffold

```bash
python paper/scripts/prepare_biomedinformatics_template.py   --template-zip paper/MDPI_template_ACS.zip   --outdir paper/BioMedInformatics_manuscript
```

The scaffold extraction does not run benchmarks or generate figures; it only prepares or refreshes the template workspace. Do not overwrite manuscript edits unless the refresh is intentional.
