#!/usr/bin/env python3
"""Pool label-free functional profile artifacts across biological replicates.

Strand-aware raw-count artifacts can be depth weighted (``sum``) or
library-equalized before averaging.  Combined DWM artifacts normally use an
equal replicate mean because the input bigWigs are already depth normalized.
All residual channels are recomputed after pooling rather than averaged.
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

from build_combined_functional_profiles import write_artifact as write_combined  # noqa: E402
from build_strand_functional_profiles import write_profiles as write_strand  # noqa: E402
from fp_tools.tools.functional_footprints import construct_strand_functional_profiles  # noqa: E402


SUPPORTED_SCHEMAS = {
    "fp-tools-combined-functional-profiles-v1",
    "fp-tools-strand-functional-profiles-v1",
}


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve(document_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_file():
        return path
    candidate = document_path.parent / path
    if candidate.is_file():
        return candidate
    raise FileNotFoundError(path)


def load_artifact(path: Path) -> tuple[dict, pd.DataFrame, dict[str, np.ndarray]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") not in SUPPORTED_SCHEMAS:
        raise ValueError(f"unsupported artifact schema in {path}: {document.get('schema')}")
    if bool(document.get("metadata", {}).get("labels_used", False)):
        raise ValueError(f"{path} reports that labels were used")
    sites_path = _resolve(path, str(document["sites"]))
    profiles_path = _resolve(path, str(document["profiles_npz"]))
    if file_sha256(sites_path) != str(document["sites_sha256"]):
        raise ValueError(f"site checksum mismatch for {path}")
    if file_sha256(profiles_path) != str(document["profiles_sha256"]):
        raise ValueError(f"profile checksum mismatch for {path}")
    sites = pd.read_csv(sites_path, sep="\t").reset_index(drop=True)
    forbidden = [
        column for column in sites if "chip" in column.lower() or "label" in column.lower()
    ]
    if forbidden:
        raise ValueError(f"{path} contains forbidden label columns: {', '.join(forbidden)}")
    with np.load(profiles_path, allow_pickle=False) as payload:
        arrays = {key: payload[key] for key in payload.files}
    return document, sites, arrays


def observed_total(schema: str, arrays: dict[str, np.ndarray]) -> float:
    valid = np.asarray(arrays["valid"], dtype=bool)
    if schema == "fp-tools-combined-functional-profiles-v1":
        return float(np.asarray(arrays["observed"])[valid].sum())
    return float(
        np.asarray(arrays["plus_observed"])[valid].sum()
        + np.asarray(arrays["minus_observed"])[valid].sum()
    )


def pooling_weights(
    schema: str,
    arrays: list[dict[str, np.ndarray]],
    mode: str,
) -> np.ndarray:
    if mode == "sum":
        return np.ones(len(arrays), dtype=float)
    if mode == "mean":
        return np.full(len(arrays), 1.0 / len(arrays), dtype=float)
    totals = np.asarray([observed_total(schema, item) for item in arrays], dtype=float)
    if np.any(totals <= 0):
        raise ValueError("library-equalized pooling requires positive observed totals")
    target = float(np.median(totals))
    # Each replicate is first scaled to the median library and then averaged.
    return (target / totals) / len(arrays)


def validate_alignment(
    documents: list[dict],
    sites: list[pd.DataFrame],
    arrays: list[dict[str, np.ndarray]],
) -> str:
    schemas = {str(document["schema"]) for document in documents}
    if len(schemas) != 1:
        raise ValueError("all input artifacts must use the same schema")
    schema = schemas.pop()
    reference = arrays[0]["site_hash"]
    for index in range(1, len(arrays)):
        if not np.array_equal(reference, arrays[index]["site_hash"]):
            raise ValueError("replicate artifact site hashes or order differ")
        if list(sites[0].columns) != list(sites[index].columns):
            raise ValueError("replicate artifact site columns differ")
    return schema


def weighted_sum(
    arrays: list[dict[str, np.ndarray]], key: str, weights: np.ndarray
) -> np.ndarray:
    return np.sum(
        [float(weight) * np.asarray(item[key], dtype=np.float64) for weight, item in zip(weights, arrays)],
        axis=0,
    )


def pool_artifacts(
    inputs: Sequence[Path],
    output_prefix: Path,
    *,
    mode: str,
    dispersion: float,
) -> tuple[Path, Path, Path]:
    if len(inputs) < 2:
        raise ValueError("pooling requires at least two replicate artifacts")
    loaded = [load_artifact(path) for path in inputs]
    documents = [item[0] for item in loaded]
    sites = [item[1] for item in loaded]
    arrays = [item[2] for item in loaded]
    schema = validate_alignment(documents, sites, arrays)
    weights = pooling_weights(schema, arrays, mode)
    valid = np.logical_and.reduce([np.asarray(item["valid"], dtype=bool) for item in arrays])
    metadata = {
        "labels_used": False,
        "pooling_mode": mode,
        "replicate_weights": weights.tolist(),
        "replicates": [
            {
                "artifact": str(path),
                "artifact_sha256": file_sha256(path),
                "profiles_sha256": str(document["profiles_sha256"]),
                "sites_sha256": str(document["sites_sha256"]),
                "observed_total": observed_total(schema, item),
            }
            for path, document, item in zip(inputs, documents, arrays)
        ],
        "dispersion": float(dispersion),
    }
    if schema == "fp-tools-combined-functional-profiles-v1":
        observed = weighted_sum(arrays, "observed", weights)
        expected = weighted_sum(arrays, "expected", weights)
        return write_combined(
            output_prefix,
            sites[0],
            observed,
            expected,
            valid,
            dispersion=dispersion,
            metadata=metadata,
        )
    plus_observed = weighted_sum(arrays, "plus_observed", weights)
    minus_observed = weighted_sum(arrays, "minus_observed", weights)
    plus_expected = weighted_sum(arrays, "plus_expected", weights)
    minus_expected = weighted_sum(arrays, "minus_expected", weights)
    profiles = construct_strand_functional_profiles(
        plus_observed,
        minus_observed,
        plus_expected,
        minus_expected,
        sites[0]["TFBS_strand"].astype(str),
        dispersion=dispersion,
    )
    return write_strand(output_prefix, sites[0], profiles, valid, metadata)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", action="append", type=Path, required=True)
    parser.add_argument("--mode", choices=("sum", "mean", "library-equalized-mean"), required=True)
    parser.add_argument("--dispersion", type=float, default=0.0)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args(argv)
    if len(args.artifact) < 2:
        raise SystemExit("provide at least two --artifact inputs")
    if args.dispersion < 0:
        raise SystemExit("--dispersion must be non-negative")
    paths = pool_artifacts(
        args.artifact,
        args.out_prefix,
        mode=args.mode,
        dispersion=args.dispersion,
    )
    print("\n".join(str(path) for path in paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
