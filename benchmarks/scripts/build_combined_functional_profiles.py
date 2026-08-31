#!/usr/bin/env python3
"""Build orientation-aligned observed/expected profiles from DWM bigWigs.

The artifact is the combined-signal counterpart of
``build_strand_functional_profiles.py``.  It is intentionally label-free and
is used to apply frozen DWM reference detectors to naked-DNA and holdout data.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import site_hashes  # noqa: E402
from search_tf_footprint_models import extract_profiles  # noqa: E402
from fp_tools.tools.functional_footprints import deviance_profiles, orient_profiles  # noqa: E402


SCHEMA = "fp-tools-combined-functional-profiles-v1"
REQUIRED_COLUMNS = {
    "TFBS_chr",
    "TFBS_start",
    "TFBS_end",
    "TFBS_strand",
}


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_unlabeled_sites(sites: pd.DataFrame, path: str | Path) -> None:
    missing = REQUIRED_COLUMNS.difference(sites.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    forbidden = [
        column
        for column in sites.columns
        if "chip" in column.lower() or "label" in column.lower()
    ]
    if forbidden:
        raise ValueError(
            f"label-free profile construction refuses columns: {', '.join(forbidden)}"
        )


def extract_combined_profiles(
    sites: pd.DataFrame,
    observed_track: str | Path,
    expected_track: str | Path,
    flank: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    observed, observed_valid = extract_profiles(sites, observed_track, flank)
    expected, expected_valid = extract_profiles(sites, expected_track, flank)
    strands = sites["TFBS_strand"].astype(str).to_numpy()
    observed = orient_profiles(np.nan_to_num(observed), strands)
    expected = orient_profiles(np.nan_to_num(expected), strands)
    return observed, expected, np.asarray(observed_valid & expected_valid, dtype=bool)


def write_artifact(
    prefix: str | Path,
    sites: pd.DataFrame,
    observed: np.ndarray,
    expected: np.ndarray,
    valid: np.ndarray,
    *,
    dispersion: float,
    metadata: dict,
) -> tuple[Path, Path, Path]:
    prefix = Path(prefix)
    npz_path = Path(str(prefix) + ".npz")
    json_path = Path(str(prefix) + ".json")
    sites_path = Path(str(prefix) + ".sites.tsv.gz")
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        npz_path,
        observed=np.asarray(observed, dtype=np.float32),
        expected=np.asarray(expected, dtype=np.float32),
        combined_residual=np.asarray(
            deviance_profiles(observed, expected, dispersion), dtype=np.float32
        ),
        valid=np.asarray(valid, dtype=bool),
        site_hash=site_hashes(sites),
    )
    sites.to_csv(sites_path, sep="\t", index=False)
    document = {
        "schema": SCHEMA,
        "profiles_npz": str(npz_path),
        "profiles_sha256": file_sha256(npz_path),
        "sites": str(sites_path),
        "sites_sha256": file_sha256(sites_path),
        "sites_total": int(len(sites)),
        "sites_valid": int(np.sum(valid)),
        "metadata": metadata,
    }
    json_path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return npz_path, json_path, sites_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--cell")
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--dispersion", type=float, default=0.0)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args(argv)
    sites = pd.read_csv(args.sites, sep="\t").reset_index(drop=True)
    validate_unlabeled_sites(sites, args.sites)
    if args.cell is not None:
        if "cell" not in sites:
            raise SystemExit("--cell requires a cell column in --sites")
        sites = sites[sites["cell"].astype(str).eq(args.cell)].reset_index(drop=True)
        if sites.empty:
            raise SystemExit(f"--sites contains no rows for cell {args.cell}")
    observed, expected, valid = extract_combined_profiles(
        sites, args.observed, args.expected, args.flank
    )
    write_artifact(
        args.out_prefix,
        sites,
        observed,
        expected,
        valid,
        dispersion=args.dispersion,
        metadata={
            "observed": str(args.observed),
            "observed_sha256": file_sha256(args.observed),
            "expected": str(args.expected),
            "expected_sha256": file_sha256(args.expected),
            "flank": int(args.flank),
            "dispersion": float(args.dispersion),
            "cell": args.cell,
            "labels_used": False,
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
