#!/usr/bin/env python3
"""Export a leakage-safe site table from a functional-profile artifact.

The output preserves site order while removing every column whose name can
carry ChIP or label information.  Filters are allowed only on non-label
columns, which makes the helper suitable for constructing chromosome-heldout
naked-DNA control panels from evaluation artifacts.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

import pandas as pd


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_label_column(column: str) -> bool:
    lowered = str(column).lower()
    return "label" in lowered or "chip" in lowered


def parse_filter(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("filters must use COLUMN=VALUE")
    column, expected = value.split("=", 1)
    if not column or expected == "":
        raise argparse.ArgumentTypeError("filters must use COLUMN=VALUE")
    if is_label_column(column):
        raise argparse.ArgumentTypeError(f"label-derived filters are forbidden: {column}")
    return column, expected


def export_sites(
    artifact: Path,
    output: Path,
    filters: Sequence[tuple[str, str]],
) -> tuple[Path, Path]:
    document = json.loads(artifact.read_text(encoding="utf-8"))
    sites_path = Path(document["sites"])
    if file_sha256(sites_path) != document["sites_sha256"]:
        raise ValueError(f"site checksum mismatch: {sites_path}")
    sites = pd.read_csv(sites_path, sep="\t")
    selected = pd.Series(True, index=sites.index)
    for column, expected in filters:
        if is_label_column(column):
            raise ValueError(f"label-derived filters are forbidden: {column}")
        if column not in sites:
            raise ValueError(f"artifact sites are missing filter column: {column}")
        selected &= sites[column].astype(str).eq(expected)
    sites = sites.loc[selected].reset_index(drop=True)
    if sites.empty:
        rendered = ", ".join(f"{column}={value}" for column, value in filters)
        raise ValueError(f"filters selected no sites: {rendered}")
    dropped = [column for column in sites.columns if is_label_column(column)]
    sites = sites.drop(columns=dropped)
    if any(is_label_column(column) for column in sites.columns):
        raise AssertionError("label-bearing columns remain after export")

    output.parent.mkdir(parents=True, exist_ok=True)
    sites.to_csv(output, sep="\t", index=False)
    manifest_path = Path(str(output) + ".manifest.json")
    manifest = {
        "schema": "fp-tools-label-free-profile-sites-v1",
        "labels_used": False,
        "source_artifact": str(artifact),
        "source_artifact_sha256": file_sha256(artifact),
        "source_sites": str(sites_path),
        "source_sites_sha256": str(document["sites_sha256"]),
        "filters": [
            {"column": column, "equals": expected} for column, expected in filters
        ],
        "dropped_columns": dropped,
        "rows": int(len(sites)),
        "columns": list(sites.columns),
        "output": str(output),
        "output_sha256": file_sha256(output),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return output, manifest_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument(
        "--where",
        action="append",
        type=parse_filter,
        default=[],
        metavar="COLUMN=VALUE",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    paths = export_sites(args.artifact, args.output, args.where)
    print("\n".join(str(path) for path in paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
