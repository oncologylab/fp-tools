#!/usr/bin/env python3
"""Create a row-filtered copy of a label-free functional-profile artifact.

This helper is intended for leakage-safe benchmark preparation, for example
extracting development chromosomes before fitting a detector that will be
applied to a negative-control dataset.  It preserves every profile channel,
recomputes site hashes, and records the exact parent artifact and filters.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import site_hashes  # noqa: E402
from pool_functional_profile_artifacts import load_artifact  # noqa: E402


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_filter(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("filters must use COLUMN=VALUE")
    column, expected = value.split("=", 1)
    if not column or expected == "":
        raise argparse.ArgumentTypeError("filters must use COLUMN=VALUE")
    return column, expected


def subset_artifact(
    artifact: Path,
    output_prefix: Path,
    filters: Sequence[tuple[str, str]],
) -> tuple[Path, Path, Path]:
    if not filters:
        raise ValueError("at least one row filter is required")
    document, sites, arrays = load_artifact(artifact)
    selected = np.ones(len(sites), dtype=bool)
    for column, expected in filters:
        if "chip" in column.lower() or "label" in column.lower():
            raise ValueError(f"label-derived filters are forbidden: {column}")
        if column not in sites:
            raise ValueError(f"artifact sites are missing filter column: {column}")
        selected &= sites[column].astype(str).eq(expected).to_numpy()
    if not selected.any():
        rendered = ", ".join(f"{column}={value}" for column, value in filters)
        raise ValueError(f"filters selected no sites: {rendered}")

    subset_sites = sites.loc[selected].reset_index(drop=True)
    subset_arrays: dict[str, np.ndarray] = {}
    for key, values in arrays.items():
        values = np.asarray(values)
        if values.ndim and values.shape[0] == len(sites):
            subset_arrays[key] = values[selected]
        else:
            subset_arrays[key] = values
    subset_arrays["site_hash"] = site_hashes(subset_sites)

    npz_path = Path(str(output_prefix) + ".npz")
    json_path = Path(str(output_prefix) + ".json")
    sites_path = Path(str(output_prefix) + ".sites.tsv.gz")
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz_path, **subset_arrays)
    subset_sites.to_csv(sites_path, sep="\t", index=False)

    metadata = {
        "labels_used": False,
        "parent_artifact": str(artifact),
        "parent_artifact_sha256": file_sha256(artifact),
        "parent_profiles_sha256": str(document["profiles_sha256"]),
        "parent_sites_sha256": str(document["sites_sha256"]),
        "row_filters": [
            {"column": column, "equals": expected} for column, expected in filters
        ],
    }
    output = {
        "schema": str(document["schema"]),
        "profiles_npz": str(npz_path),
        "profiles_sha256": file_sha256(npz_path),
        "sites": str(sites_path),
        "sites_sha256": file_sha256(sites_path),
        "sites_total": int(len(subset_sites)),
        "sites_valid": int(np.asarray(subset_arrays["valid"], dtype=bool).sum()),
        "metadata": metadata,
    }
    json_path.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return npz_path, json_path, sites_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument(
        "--where",
        action="append",
        type=parse_filter,
        required=True,
        metavar="COLUMN=VALUE",
    )
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args(argv)
    paths = subset_artifact(args.artifact, args.out_prefix, args.where)
    print("\n".join(str(path) for path in paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
