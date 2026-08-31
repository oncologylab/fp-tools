#!/usr/bin/env python3
"""Build a checksummed geometric coefficient ensemble from compatible fits."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from fp_tools.tools.parametric_bias import (  # noqa: E402
    ConditionalSequenceBiasModel,
    ensemble_sequence_bias_models,
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def shared_provenance(paths: list[Path]) -> dict:
    documents = [json.loads(path.with_suffix(".json").read_text(encoding="utf-8")) for path in paths]
    keys = (
        "training_source",
        "configuration",
        "read_shift",
        "training_depth_cuts_per_sample",
        "l2",
    )
    shared = {}
    for key in keys:
        values = [document.get("metadata", {}).get(key) for document in documents]
        encoded = {json.dumps(value, sort_keys=True) for value in values}
        if len(encoded) != 1:
            raise ValueError(f"ensemble members disagree on {key}: {values}")
        shared[key] = values[0]
    return shared


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, action="append", required=True)
    parser.add_argument("--weight", type=float, action="append")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--label", required=True)
    args = parser.parse_args(argv)
    paths = [path if path.suffix == ".npz" else Path(str(path) + ".npz") for path in args.model]
    if len(paths) < 2:
        raise SystemExit("at least two --model values are required")
    if args.weight is not None and len(args.weight) != len(paths):
        raise SystemExit("--weight must be supplied once per --model")
    provenance = shared_provenance(paths)
    members = [ConditionalSequenceBiasModel.load(path) for path in paths]
    ensemble = ensemble_sequence_bias_models(members, args.weight)
    npz_path, json_path = ensemble.save(
        args.out,
        metadata={
            **provenance,
            "ensemble_label": args.label,
            "member_models": [
                {
                    "path": str(path),
                    "npz_sha256": file_sha256(path),
                    "json_sha256": file_sha256(path.with_suffix(".json")),
                }
                for path in paths
            ],
            "ensemble_size": len(paths),
            "ensemble_weights": (
                args.weight if args.weight is not None else [1.0 / len(paths)] * len(paths)
            ),
        },
    )
    print(npz_path)
    print(json_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
