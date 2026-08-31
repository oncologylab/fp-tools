#!/usr/bin/env python3
"""Build strand-aware motif-centered profiles for functional research models.

The builder retains forward and reverse observed/expected cuts, swaps strands
correctly when orienting reverse motifs, and emits shared and antisymmetric
signed-deviance channels.  It is benchmark infrastructure only; production
commands and current DWM outputs are unchanged.
"""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from fp_tools.tools.functional_footprints import (  # noqa: E402
    construct_strand_functional_profiles,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    ConditionalSequenceBiasModel,
    cut_position_from_alignment,
    strand_log_bias,
)
from fp_tools.utils.fasta import open_fasta  # noqa: E402


SCHEMA = "fp-tools-strand-functional-profiles-v1"
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


def parse_shift(value: str) -> tuple[int, int]:
    fields = value.replace("/", ",").split(",")
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("read shift must use FORWARD,REVERSE")
    try:
        return int(fields[0]), int(fields[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("read shift must contain two integers") from exc


def site_hashes(sites: pd.DataFrame) -> np.ndarray:
    output = []
    for row in sites.itertuples(index=False):
        digest = blake2b(digest_size=8)
        for value in (
            row.TFBS_chr,
            int(row.TFBS_start),
            int(row.TFBS_end),
            row.TFBS_strand,
            getattr(row, "tf", ""),
        ):
            digest.update(str(value).encode())
            digest.update(b"\0")
        output.append(int.from_bytes(digest.digest(), "little"))
    return np.asarray(output, dtype=np.uint64)


def _usable_read(read, minimum_mapq: int, keep_duplicates: bool) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or (read.is_duplicate and not keep_duplicates)
        or int(read.mapping_quality) < minimum_mapq
    )


def _merged_ranges(starts: np.ndarray, ends: np.ndarray) -> list[tuple[int, int]]:
    ranges: list[list[int]] = []
    for start, end in zip(starts, ends):
        if ranges and int(start) <= ranges[-1][1]:
            ranges[-1][1] = max(ranges[-1][1], int(end))
        else:
            ranges.append([int(start), int(end)])
    return [(start, end) for start, end in ranges]


def extract_strand_cut_profiles(
    sites: pd.DataFrame,
    bam_path: str | Path,
    *,
    flank: int,
    read_shift: tuple[int, int],
    minimum_mapq: int = 30,
    keep_duplicates: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract genomic plus/minus cuts with overlap-aware BAM traversal."""

    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - research extraction is Linux-only
        raise RuntimeError("pysam is required to build strand profiles") from exc
    missing = REQUIRED_COLUMNS.difference(sites.columns)
    if missing:
        raise ValueError("sites are missing columns: " + ", ".join(sorted(missing)))
    width = flank * 2 + 1
    plus = np.zeros((len(sites), width), dtype=np.uint32)
    minus = np.zeros_like(plus)
    valid = np.zeros(len(sites), dtype=bool)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        references = dict(zip(bam.references, bam.lengths))
        for chromosome, group in sites.groupby("TFBS_chr", sort=True):
            chromosome = str(chromosome)
            if chromosome not in references:
                continue
            rows = group.index.to_numpy(dtype=int)
            centers = (
                group["TFBS_start"].to_numpy(dtype=int)
                + group["TFBS_end"].to_numpy(dtype=int)
            ) // 2
            window_starts = centers - flank
            window_ends = centers + flank + 1
            in_bounds = (window_starts >= 0) & (window_ends <= int(references[chromosome]))
            valid[rows[in_bounds]] = True
            rows = rows[in_bounds]
            window_starts = window_starts[in_bounds]
            window_ends = window_ends[in_bounds]
            if len(rows) == 0:
                continue
            order = np.argsort(window_starts, kind="mergesort")
            rows = rows[order]
            window_starts = window_starts[order]
            window_ends = window_ends[order]
            for region_start, region_end in _merged_ranges(window_starts, window_ends):
                for read in bam.fetch(chromosome, region_start, region_end):
                    if not _usable_read(read, minimum_mapq, keep_duplicates):
                        continue
                    cut = cut_position_from_alignment(read, read_shift)
                    first = int(np.searchsorted(window_starts, cut - width + 1, side="left"))
                    last = int(np.searchsorted(window_starts, cut, side="right"))
                    for position in range(first, last):
                        if cut >= window_ends[position]:
                            continue
                        offset = cut - int(window_starts[position])
                        if 0 <= offset < width:
                            target = minus if read.is_reverse else plus
                            target[int(rows[position]), offset] += 1
    return plus, minus, valid


def _expected_from_scores(observed: np.ndarray, scores: np.ndarray, valid: np.ndarray) -> np.ndarray:
    output = np.zeros_like(observed, dtype=np.float64)
    total = float(np.sum(observed))
    if total <= 0 or not np.any(valid):
        return output
    centered = np.where(valid, scores - np.max(scores[valid]), -np.inf)
    propensity = np.exp(centered)
    propensity /= max(float(propensity.sum()), np.finfo(float).tiny)
    return total * propensity


def predict_strand_expected_profiles(
    sites: pd.DataFrame,
    plus_observed: np.ndarray,
    minus_observed: np.ndarray,
    model: ConditionalSequenceBiasModel,
    genome: str | Path,
    *,
    flank: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Predict forward/reverse expectations while preserving strand totals."""

    width = flank * 2 + 1
    if plus_observed.shape != (len(sites), width) or minus_observed.shape != plus_observed.shape:
        raise ValueError("observed strand profiles do not match sites and flank")
    plus_expected = np.zeros_like(plus_observed, dtype=np.float64)
    minus_expected = np.zeros_like(minus_observed, dtype=np.float64)
    valid_rows = np.zeros(len(sites), dtype=bool)
    margin = max(41, model.feature_spec.context_length // 2 + 1)
    positions = margin + np.arange(width)
    with open_fasta(genome) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for index, row in enumerate(sites.itertuples(index=False)):
            chromosome = str(row.TFBS_chr)
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            start = center - flank - margin
            end = center + flank + margin + 1
            if chromosome not in lengths or start < 0 or end > lengths[chromosome]:
                continue
            sequence = fasta.fetch(chromosome, start, end).upper()
            plus_scores, minus_scores, valid = strand_log_bias(model, sequence, positions)
            if not valid.any():
                continue
            plus_expected[index] = _expected_from_scores(plus_observed[index], plus_scores, valid)
            minus_expected[index] = _expected_from_scores(minus_observed[index], minus_scores, valid)
            valid_rows[index] = True
    return plus_expected, minus_expected, valid_rows


def write_profiles(
    prefix: str | Path,
    sites: pd.DataFrame,
    profiles,
    valid: np.ndarray,
    metadata: dict,
) -> tuple[Path, Path, Path]:
    prefix = Path(prefix)
    npz_path = prefix.with_suffix(".npz")
    json_path = prefix.with_suffix(".json")
    sites_path = prefix.with_suffix(".sites.tsv.gz")
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        "plus_observed": profiles.plus_observed,
        "minus_observed": profiles.minus_observed,
        "plus_expected": profiles.plus_expected,
        "minus_expected": profiles.minus_expected,
        "combined_residual": profiles.combined_residual,
        "shared_strand_residual": profiles.shared_strand_residual,
        "antisymmetric_strand_residual": profiles.antisymmetric_strand_residual,
        "valid": np.asarray(valid, dtype=bool),
        "site_hash": site_hashes(sites),
    }
    np.savez_compressed(npz_path, **arrays)
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
    json_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return npz_path, json_path, sites_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--bias-model", type=Path, required=True)
    parser.add_argument("--read-shift", type=parse_shift, default=(4, -5))
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--minimum-mapq", type=int, default=30)
    parser.add_argument("--keep-duplicates", action="store_true")
    parser.add_argument("--dispersion", type=float, default=0.0)
    parser.add_argument("--out-prefix", type=Path, required=True)
    args = parser.parse_args(argv)
    sites = pd.read_csv(args.sites, sep="\t").reset_index(drop=True)
    plus, minus, bam_valid = extract_strand_cut_profiles(
        sites,
        args.bam,
        flank=args.flank,
        read_shift=args.read_shift,
        minimum_mapq=args.minimum_mapq,
        keep_duplicates=args.keep_duplicates,
    )
    model = ConditionalSequenceBiasModel.load(args.bias_model)
    plus_expected, minus_expected, sequence_valid = predict_strand_expected_profiles(
        sites,
        plus,
        minus,
        model,
        args.genome,
        flank=args.flank,
    )
    profiles = construct_strand_functional_profiles(
        plus,
        minus,
        plus_expected,
        minus_expected,
        sites["TFBS_strand"].astype(str),
        dispersion=args.dispersion,
    )
    write_profiles(
        args.out_prefix,
        sites,
        profiles,
        bam_valid & sequence_valid,
        {
            "bam": str(args.bam),
            "bam_sha256": file_sha256(args.bam),
            "genome": str(args.genome),
            "genome_sha256": file_sha256(args.genome),
            "bias_model": str(args.bias_model),
            "bias_model_sha256": file_sha256(args.bias_model),
            "read_shift": list(args.read_shift),
            "flank": int(args.flank),
            "minimum_mapq": int(args.minimum_mapq),
            "keep_duplicates": bool(args.keep_duplicates),
            "dispersion": float(args.dispersion),
            "labels_used": False,
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
