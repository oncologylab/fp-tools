#!/usr/bin/env python
"""Build motif-centric TFBS feature tables from fp-tools candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

from fp_tools.tools.motif_discovery import CandidateSite, load_genome_fasta, read_candidate_sites


def _read_overlap_intervals(path: str | Path | None) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    if path is None:
        return intervals
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            intervals.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return intervals


def _overlaps(site: CandidateSite, intervals: dict[str, list[tuple[int, int]]], min_overlap_bp: int = 1) -> bool:
    for start, end in intervals.get(site.chrom, []):
        if max(0, min(site.end, end) - max(site.start, start)) >= min_overlap_bp:
            return True
    return False


def _gc_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    valid = [base for base in sequence if base in {"A", "C", "G", "T"}]
    if not valid:
        return float("nan")
    return float(sum(base in {"G", "C"} for base in valid) / len(valid))


def _signal_features_from_bw(site: CandidateSite, bw: pyBigWig.pyBigWig, label: str) -> dict[str, float]:
    chroms = bw.chroms()
    if site.chrom not in chroms:
        return {f"{label}_mean": float("nan"), f"{label}_max": float("nan")}
    start = max(0, site.start)
    end = min(int(chroms[site.chrom]), site.end)
    values = np.asarray(bw.values(site.chrom, start, end), dtype=float) if end > start else np.array([], dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {f"{label}_mean": float("nan"), f"{label}_max": float("nan")}
    return {f"{label}_mean": float(np.mean(finite)), f"{label}_max": float(np.max(finite))}


def build_feature_table(
    candidates: str | Path,
    output: str | Path,
    genome: str | Path | None = None,
    signals: list[str | Path] | None = None,
    signal_labels: list[str] | None = None,
    labels_bed: str | Path | None = None,
    min_label_overlap_bp: int = 1,
) -> pd.DataFrame:
    """Build a tabular feature matrix for supervised TFBS training or prediction."""

    sites = read_candidate_sites(candidates)
    genome_records = load_genome_fasta(genome) if genome is not None else {}
    label_intervals = _read_overlap_intervals(labels_bed)
    signal_paths = list(signals or [])
    if signal_labels is None:
        signal_labels = [Path(path).stem.replace(".", "_") for path in signal_paths]
    if len(signal_labels) != len(signal_paths):
        raise ValueError("--signal-labels must match --signals length")

    rows = []
    signal_handles = [(pyBigWig.open(str(path)), label) for path, label in zip(signal_paths, signal_labels)]
    try:
        for site in sites:
            sequence = genome_records.get(site.chrom, "")[site.start:site.end] if genome_records else ""
            try:
                candidate_score = float(site.score)
            except ValueError:
                candidate_score = float("nan")
            row = {
                "site_id": site.name,
                "chrom": site.chrom,
                "start": site.start,
                "end": site.end,
                "candidate_score": candidate_score,
                "length": int(site.end - site.start),
                "gc": _gc_fraction(sequence) if genome_records else float("nan"),
            }
            if labels_bed is not None:
                row["label"] = int(_overlaps(site, label_intervals, min_overlap_bp=min_label_overlap_bp))
            for bw, label in signal_handles:
                row.update(_signal_features_from_bw(site, bw, label))
            rows.append(row)
    finally:
        for bw, _label in signal_handles:
            bw.close()

    frame = pd.DataFrame(rows)
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, sep="\t", index=False)
    return frame


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build supervised TFBS feature TSVs from candidate BED rows.")
    parser.add_argument("--candidates", required=True, help="BED or fp-tools candidate table.")
    parser.add_argument("--out", required=True, help="Output feature TSV.")
    parser.add_argument("--genome", default=None, help="Optional genome FASTA/FASTA.GZ for GC content.")
    parser.add_argument("--signals", nargs="*", default=[], help="Optional bigWig signal tracks to summarize over candidates.")
    parser.add_argument("--signal-labels", nargs="*", default=None, help="Labels for --signals; defaults to file stems.")
    parser.add_argument("--labels-bed", default=None, help="Optional BED intervals treated as positive labels by overlap.")
    parser.add_argument("--min-label-overlap-bp", type=int, default=1, help="Minimum overlap bp for assigning label=1.")
    args = parser.parse_args(argv)

    frame = build_feature_table(
        args.candidates,
        args.out,
        genome=args.genome,
        signals=args.signals,
        signal_labels=args.signal_labels,
        labels_bed=args.labels_bed,
        min_label_overlap_bp=args.min_label_overlap_bp,
    )
    print(f"Wrote {len(frame)} feature rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
