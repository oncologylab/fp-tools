#!/usr/bin/env python3
"""Prepare a BioMedInformatics/MDPI LaTeX manuscript workspace."""

from __future__ import annotations

import argparse
import shutil
import zipfile
from pathlib import Path


REQUIRED_TEMPLATE_FILES = [
    "template.tex",
    "template.pdf",
    "Definitions/mdpi.cls",
    "Definitions/mdpi.bst",
]


def validate_template_zip(path: str | Path) -> list[str]:
    archive = Path(path)
    with zipfile.ZipFile(archive) as zf:
        names = set(zf.namelist())
    missing = [name for name in REQUIRED_TEMPLATE_FILES if name not in names]
    if missing:
        raise ValueError(f"Template archive is missing required files: {', '.join(missing)}")
    return sorted(names)


def prepare_template(template_zip: str | Path, outdir: str | Path, *, force: bool = False) -> Path:
    template_zip = Path(template_zip)
    outdir = Path(outdir)
    validate_template_zip(template_zip)
    if outdir.exists():
        if not force:
            raise FileExistsError(f"{outdir} already exists; use --force to replace it")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(template_zip) as zf:
        zf.extractall(outdir)
    main_tex = outdir / "main.tex"
    if not main_tex.exists() and (outdir / "template.tex").exists():
        shutil.copy2(outdir / "template.tex", main_tex)
    for subdir in ("figures", "tables", "supplement"):
        (outdir / subdir).mkdir(exist_ok=True)
    for filename, heading in [
        ("cover_letter.md", "# Cover Letter"),
        ("data_availability.md", "# Data Availability Statement"),
        ("code_availability.md", "# Code Availability Statement"),
    ]:
        path = outdir / filename
        if not path.exists():
            path.write_text(heading + "\n\nTODO\n", encoding="utf-8")
    return outdir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-zip", default="paper/MDPI_template_ACS.zip")
    parser.add_argument("--outdir", default="paper/manuscript")
    parser.add_argument("--force", action="store_true", help="Replace an existing output directory.")
    args = parser.parse_args()
    outdir = prepare_template(args.template_zip, args.outdir, force=args.force)
    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
