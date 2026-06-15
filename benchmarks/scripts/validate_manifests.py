#!/usr/bin/env python
"""Validate benchmark manifest TSV schemas."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


FULL_MANIFEST_COLUMNS = [
    "source",
    "benchmark_tier",
    "cell_type",
    "donor",
    "tf",
    "assay",
    "experiment_accession",
    "file_accession",
    "assembly",
    "output_type",
    "file_format",
    "url",
    "checksum",
    "status",
    "local_path",
    "split",
    "notes",
]

COMPACT_MANIFEST_SCHEMAS = {
    "a549_tasks.tsv": ["cell", "tf", "assay", "experiment", "file_accession", "url"],
    "10x_pbmc_pseudobulk.tsv": ["dataset", "asset", "url"],
}


def read_header(path: Path) -> list[str]:
    return list(pd.read_csv(path, sep="\t", nrows=0).columns)


def validate_manifest(path: Path, compact_root: Path | None = None) -> list[str]:
    errors: list[str] = []
    header = read_header(path)
    if compact_root is not None and compact_root in path.parents:
        expected = COMPACT_MANIFEST_SCHEMAS.get(path.name)
        if expected is None:
            errors.append(f"{path}: no compact schema is registered")
        elif header != expected:
            errors.append(f"{path}: expected compact columns {expected}, found {header}")
        return errors
    missing = [column for column in FULL_MANIFEST_COLUMNS if column not in header]
    extra = [column for column in header if column not in FULL_MANIFEST_COLUMNS]
    if missing:
        errors.append(f"{path}: missing full-manifest columns {missing}")
    if extra:
        errors.append(f"{path}: unexpected full-manifest columns {extra}")
    return errors


def validate_manifests(manifest_dir: Path) -> list[str]:
    compact_root = manifest_dir / "compact"
    errors: list[str] = []
    for path in sorted(manifest_dir.rglob("*.tsv")):
        errors.extend(validate_manifest(path, compact_root=compact_root))
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-dir", type=Path, default=Path("benchmarks/manifests"))
    args = parser.parse_args(argv)
    errors = validate_manifests(args.manifest_dir)
    if errors:
        for error in errors:
            print(error)
        return 1
    print(f"Validated benchmark manifests under {args.manifest_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
