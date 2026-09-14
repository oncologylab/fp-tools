"""Shared input schema for built-in single-cell signature analysis."""

import gzip
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {"barcode", "cell_type", "snap_cell_type", "umap_1", "umap_2"}


def read_signature_annotations(path, *, header_only=False):
    """Read TSV/CSV annotations, or validate just their header for preflight."""
    path = Path(path).expanduser()
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8-sig") as stream:
            first = stream.readline()
        separator = "," if path.suffix.lower() == ".csv" or first.count(",") > first.count("\t") else "\t"
        annotations = pd.read_csv(path, sep=separator, nrows=0 if header_only else None)
    except (OSError, UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError(f"Could not read annotation table: {exc}") from exc
    missing = REQUIRED_COLUMNS.difference(annotations.columns)
    if missing:
        raise ValueError(f"Annotation table is missing required columns: {', '.join(sorted(missing))}")
    return annotations
