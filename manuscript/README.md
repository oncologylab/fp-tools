# Manuscript Workspace

This directory is the canonical flat BioMedInformatics manuscript workspace for `fp-tools`.

## Layout

- `main.tex`, `references.bib`, and `main.pdf`: manuscript source, bibliography, and compiled preview.
- `Definitions/`: MDPI/BioMedInformatics LaTeX class, bibliography style, and logo assets.
- `figures/`: curated manuscript figures and source TSVs that are intentionally versioned.
- `tables/`: curated manuscript source tables, including software versions and analysis parameters.
- `scripts/`: reproducible figure and table builders.
- `templates/`: original journal templates and author instructions.
- `code_availability.md`, `data_availability.md`, and `cover_letter.md`: submission support text.

Generated LaTeX auxiliary files and large exploratory figure/table outputs are ignored by default. Keep only manuscript-ready source assets in this directory.

## Compile

```bash
cd manuscript && latexmk -pdf -shell-escape -interaction=nonstopmode main.tex
```

## Refresh the LaTeX scaffold

```bash
python manuscript/scripts/prepare_biomedinformatics_template.py   --template-zip manuscript/templates/MDPI_template_ACS.zip   --outdir manuscript
```

Do not refresh with `--force` unless intentionally replacing local manuscript edits.
