#!/usr/bin/env python3
"""Create label-free train/test site plans for locked holdout profiling.

All test-chromosome sites are retained.  Training sites are deterministically
sampled to the study's preregistered per-TF ceiling.  Validation chromosomes
are intentionally excluded because no holdout model choice is permitted.
"""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


SCHEMA = "fp-tools-locked-holdout-profile-sites-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*values: object, seed: int) -> int:
    digest = blake2b(digest_size=8)
    digest.update(str(seed).encode())
    for value in values:
        digest.update(b"\0")
        digest.update(str(value).encode())
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def validate_label_free(frame: pd.DataFrame, path: Path) -> None:
    forbidden = [
        column for column in frame if "chip" in column.lower() or "label" in column.lower()
    ]
    if forbidden:
        raise ValueError(f"{path} contains forbidden label columns: {', '.join(forbidden)}")
    required = {
        "cell",
        "tf",
        "motif_id",
        "motif_family",
        "TFBS_chr",
        "TFBS_start",
        "TFBS_end",
        "TFBS_strand",
        "motif_score",
        "peak_start",
        "peak_end",
        "chromosome_split",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} lacks columns: {', '.join(sorted(missing))}")
    observed = set(frame["chromosome_split"].astype(str))
    if not observed.issubset({"train", "validation", "test"}):
        raise ValueError(f"{path} contains unexpected chromosome splits: {sorted(observed)}")


def select_profile_sites(
    frame: pd.DataFrame,
    *,
    maximum_train_per_tf: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if maximum_train_per_tf < 1:
        raise ValueError("maximum_train_per_tf must be positive")
    rows = []
    counts = []
    for (cell, tf), group in frame.groupby(["cell", "tf"], sort=True):
        train = group[group["chromosome_split"].astype(str).eq("train")]
        test = group[group["chromosome_split"].astype(str).eq("test")]
        if len(train) > maximum_train_per_tf:
            rng = np.random.default_rng(stable_seed(cell, tf, "train", seed=seed))
            indexes = np.sort(
                rng.choice(train.index.to_numpy(dtype=int), size=maximum_train_per_tf, replace=False)
            )
            train = frame.loc[indexes]
        selected = pd.concat([train, test], ignore_index=False)
        rows.append(selected)
        counts.append(
            {
                "cell": str(cell),
                "tf": str(tf),
                "source_train": int(
                    group["chromosome_split"].astype(str).eq("train").sum()
                ),
                "selected_train": int(len(train)),
                "source_validation": int(
                    group["chromosome_split"].astype(str).eq("validation").sum()
                ),
                "selected_validation": 0,
                "source_test": int(
                    group["chromosome_split"].astype(str).eq("test").sum()
                ),
                "selected_test": int(len(test)),
            }
        )
    output = pd.concat(rows, ignore_index=False).sort_values(
        ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"],
        kind="mergesort",
    ).reset_index(drop=True)
    return output, pd.DataFrame(counts)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--counts-out", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    parser.add_argument("--maximum-train-per-tf", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument(
        "--routes",
        type=Path,
        help="Optional frozen route table used to retain only one bias configuration's tasks.",
    )
    parser.add_argument(
        "--bias-configuration",
        help="Bias configuration selected from --routes (for example MT_SELMA10_4m4).",
    )
    args = parser.parse_args(argv)
    source = pd.read_csv(args.sites, sep="\t")
    validate_label_free(source, args.sites)
    full_source_sites = len(source)
    route_metadata = None
    if (args.routes is None) != (args.bias_configuration is None):
        raise SystemExit("--routes and --bias-configuration must be provided together")
    if args.routes is not None:
        routes = pd.read_csv(args.routes, sep="\t")
        required = {"cell", "tf", "bias_configuration"}
        missing = required.difference(routes.columns)
        if missing:
            raise SystemExit("route table lacks columns: " + ", ".join(sorted(missing)))
        selected_routes = routes[
            routes["bias_configuration"].astype(str).eq(args.bias_configuration)
        ]
        keys = set(zip(selected_routes["cell"].astype(str), selected_routes["tf"].astype(str)))
        source = source[
            [
                (str(cell), str(tf)) in keys
                for cell, tf in zip(source["cell"], source["tf"])
            ]
        ].copy()
        if source.empty:
            raise SystemExit(
                f"no site tasks use bias configuration {args.bias_configuration}"
            )
        route_metadata = {
            "path": str(args.routes),
            "sha256": file_sha256(args.routes),
            "bias_configuration": str(args.bias_configuration),
            "tasks": [f"{cell}/{tf}" for cell, tf in sorted(keys)],
        }
    selected, counts = select_profile_sites(
        source,
        maximum_train_per_tf=args.maximum_train_per_tf,
        seed=args.seed,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.counts_out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(
        args.out,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    counts.to_csv(args.counts_out, sep="\t", index=False)
    document = {
        "schema": SCHEMA,
        "locked_holdout_labels_read": False,
        "selection_uses_labels": False,
        "source": str(args.sites),
        "source_sha256": file_sha256(args.sites),
        "maximum_train_per_tf": int(args.maximum_train_per_tf),
        "seed": int(args.seed),
        "validation_sites_retained": 0,
        "all_test_sites_retained": bool(
            counts["source_test"].eq(counts["selected_test"]).all()
        ),
        "source_sites": int(len(source)),
        "unfiltered_source_sites": int(full_source_sites),
        "selected_sites": int(len(selected)),
        "route_filter": route_metadata,
        "output": str(args.out),
        "output_sha256": file_sha256(args.out),
        "counts": str(args.counts_out),
        "counts_sha256": file_sha256(args.counts_out),
    }
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(counts.to_string(index=False))
    print(f"selected {len(selected):,} of {len(source):,} label-free sites")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
