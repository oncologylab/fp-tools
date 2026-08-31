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
    *,
    bed_output: Path | None = None,
    flank: int = 100,
) -> tuple[Path, Path]:
    if flank < 0:
        raise ValueError("flank must be non-negative")
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
    bed_record = None
    if bed_output is not None:
        required = {"TFBS_chr", "TFBS_start", "TFBS_end"}
        missing = required.difference(sites.columns)
        if missing:
            raise ValueError(
                "cannot write BED; sites are missing columns: "
                + ", ".join(sorted(missing))
            )
        center = (
            sites["TFBS_start"].to_numpy(dtype=int)
            + sites["TFBS_end"].to_numpy(dtype=int)
        ) // 2
        bed = pd.DataFrame(
            {
                "chromosome": sites["TFBS_chr"].astype(str),
                "start": (center - flank).clip(min=0),
                "end": center + flank + 1,
            }
        )
        bed_output.parent.mkdir(parents=True, exist_ok=True)
        bed.to_csv(bed_output, sep="\t", index=False, header=False)
        bed_record = {
            "path": str(bed_output),
            "sha256": file_sha256(bed_output),
            "flank": int(flank),
            "rows": int(len(bed)),
        }
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
        "bed_output": bed_record,
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
    parser.add_argument("--bed-output", type=Path)
    parser.add_argument("--flank", type=int, default=100)
    args = parser.parse_args(argv)
    paths = export_sites(
        args.artifact,
        args.output,
        args.where,
        bed_output=args.bed_output,
        flank=args.flank,
    )
    print("\n".join(str(path) for path in paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
