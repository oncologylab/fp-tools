#!/usr/bin/env python3
"""Create deterministic, nested fragment-level BAM depth subsets.

Both mates receive the same hash decision because selection is based on the
query name.  Reusing a seed at increasing target depths creates nested subsets,
which is essential for paired correction/scoring ablations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def fragment_uniform(query_name: str, seed: int) -> float:
    payload = f"{seed}\0{query_name}".encode("utf-8")
    value = int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big")
    return value / float(2**64)


def usable_alignment(read) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
    )


def fragment_representative(read) -> bool:
    return usable_alignment(read) and (not read.is_paired or read.is_read1)


def count_fragments(path: Path) -> int:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - environment-dependent message
        raise RuntimeError("pysam is required to downsample BAM files") from exc
    count = 0
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            count += int(fragment_representative(read))
    return count


def downsample_bam(
    source: Path,
    output: Path,
    target_fragments: int,
    seed: int,
    create_index: bool = True,
    available_fragments: int | None = None,
) -> dict[str, object]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - environment-dependent message
        raise RuntimeError("pysam is required to downsample BAM files") from exc
    if target_fragments <= 0:
        raise ValueError("target_fragments must be positive")
    if available_fragments is not None and available_fragments <= 0:
        raise ValueError("available_fragments must be positive when supplied")
    available = (
        int(available_fragments)
        if available_fragments is not None
        else count_fragments(source)
    )
    if available == 0:
        raise ValueError(f"{source} contains no usable fragments")
    fraction = min(1.0, target_fragments / available)
    output.parent.mkdir(parents=True, exist_ok=True)
    selected = 0
    with pysam.AlignmentFile(str(source), "rb") as source_bam:
        with pysam.AlignmentFile(str(output), "wb", template=source_bam) as output_bam:
            for read in source_bam.fetch(until_eof=True):
                if not usable_alignment(read):
                    continue
                if fragment_uniform(read.query_name, seed) < fraction:
                    output_bam.write(read)
                    selected += int(fragment_representative(read))
    if create_index:
        pysam.index(str(output))
    return {
        "source_bam": str(source),
        "output_bam": str(output),
        "seed": int(seed),
        "available_fragments": int(available),
        "available_fragments_source": "provided" if available_fragments is not None else "counted",
        "target_fragments": int(target_fragments),
        "sampling_fraction": float(fraction),
        "selected_fragments": int(selected),
        "selection_rule": "blake2b(seed,query_name)<sampling_fraction",
        "nested_with_same_seed": True,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--target-fragments", type=int, required=True)
    parser.add_argument(
        "--available-fragments",
        type=int,
        help="Validated usable-fragment count; avoids a separate full BAM counting pass.",
    )
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--metadata-out", type=Path)
    parser.add_argument("--no-index", action="store_true")
    args = parser.parse_args(argv)

    metadata = downsample_bam(
        args.bam,
        args.out,
        args.target_fragments,
        args.seed,
        create_index=not args.no_index,
        available_fragments=args.available_fragments,
    )
    metadata_out = args.metadata_out or args.out.with_suffix(".downsampling.json")
    metadata_out.parent.mkdir(parents=True, exist_ok=True)
    metadata_out.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
