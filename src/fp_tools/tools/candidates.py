#!/usr/bin/env python
"""Motif-free footprint candidate generation."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyBigWig


@dataclass(frozen=True)
class Candidate:
    chrom: str
    start: int
    end: int
    name: str
    score: float
    strand: str = "."
    generator: str = "signal-local-max"
    rank: int = 0
    source_region: str = ""
    candidate_score: float | None = None
    motif_id: str = "."
    motif_score: float | None = None
    motif_source: str = "."

    def to_bed_row(self) -> list[str]:
        candidate_score = self.score if self.candidate_score is None else self.candidate_score
        motif_score = "" if self.motif_score is None else f"{self.motif_score:.6f}"
        return [
            self.chrom,
            str(self.start),
            str(self.end),
            self.name,
            f"{self.score:.6f}",
            self.strand,
            self.generator,
            str(self.rank),
            self.source_region,
            f"{candidate_score:.6f}",
            self.motif_id,
            motif_score,
            self.motif_source,
        ]


def call_candidates_from_array(
    values: np.ndarray,
    chrom: str,
    region_start: int,
    region_end: int,
    candidate_width: int = 20,
    window: int = 10,
    min_score: float | None = None,
    top_n: int | None = 5,
    generator: str = "signal-local-max",
) -> list[Candidate]:
    """Call local score maxima from one genomic interval."""

    arr = np.asarray(values, dtype=float)
    arr = np.nan_to_num(arr, nan=-np.inf)
    if arr.size == 0:
        return []

    candidates: list[Candidate] = []
    half_width = max(1, int(candidate_width) // 2)
    local_window = max(1, int(window))
    source = f"{chrom}:{region_start}-{region_end}"

    for idx, score in enumerate(arr):
        if not np.isfinite(score):
            continue
        if min_score is not None and score < min_score:
            continue
        left = max(0, idx - local_window)
        right = min(arr.size, idx + local_window + 1)
        if score < np.max(arr[left:right]):
            continue
        # Keep the first base of a plateau to avoid duplicate calls.
        if idx > 0 and arr[idx - 1] == score:
            continue

        center = region_start + idx
        start = max(region_start, center - half_width)
        end = min(region_end, center + half_width)
        candidates.append(
            Candidate(
                chrom=chrom,
                start=int(start),
                end=int(max(start + 1, end)),
                name=f"candidate_{len(candidates) + 1}",
                score=float(score),
                generator=generator,
                source_region=source,
            )
        )

    candidates.sort(key=lambda candidate: (-candidate.score, candidate.start))
    if top_n is not None:
        candidates = candidates[: int(top_n)]
    return [
        Candidate(
            chrom=c.chrom,
            start=c.start,
            end=c.end,
            name=f"{generator}_{rank}",
            score=c.score,
            strand=c.strand,
            generator=c.generator,
            rank=rank,
            source_region=c.source_region,
            candidate_score=c.score,
        )
        for rank, c in enumerate(candidates, start=1)
    ]


def _read_bed_regions(path: str | Path) -> list[tuple[str, int, int]]:
    regions = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            regions.append((fields[0], int(fields[1]), int(fields[2])))
    return regions


def _overlaps_any_peak(chrom: str, start: int, end: int, peaks: list[tuple[str, int, int]]) -> bool:
    return any(chrom == peak_chrom and min(end, peak_end) > max(start, peak_start) for peak_chrom, peak_start, peak_end in peaks)


def _summarize_bigwig_interval(bw: pyBigWig.pyBigWig, chrom: str, start: int, end: int, signal_math: str) -> float:
    values = np.asarray(bw.values(chrom, start, end), dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0
    if signal_math == "mean":
        return float(np.mean(finite))
    return float(np.max(finite))


def read_motif_site_candidates(
    motif_sites: list[str | Path],
    peaks: list[tuple[str, int, int]],
    bw: pyBigWig.pyBigWig,
    chrom_sizes: dict[str, int],
    signal_math: str = "max",
    generator: str = "motif-relaxed",
) -> list[Candidate]:
    """Read lower-threshold motif BEDs and score their intervals from the signal bigWig."""

    candidates: list[Candidate] = []
    for motif_path in motif_sites:
        source = str(motif_path)
        with Path(motif_path).open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                if chrom not in chrom_sizes:
                    continue
                start = max(0, int(fields[1]))
                end = min(int(chrom_sizes[chrom]), int(fields[2]))
                if end <= start or not _overlaps_any_peak(chrom, start, end, peaks):
                    continue
                motif_id = fields[3] if len(fields) > 3 and fields[3] else f"motif_site_{len(candidates) + 1}"
                try:
                    motif_score = float(fields[4]) if len(fields) > 4 and fields[4] != "." else None
                except ValueError:
                    motif_score = None
                strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-", "."} else "."
                signal_score = _summarize_bigwig_interval(bw, chrom, start, end, signal_math)
                candidates.append(
                    Candidate(
                        chrom=chrom,
                        start=start,
                        end=end,
                        name=f"{generator}_{len(candidates) + 1}",
                        score=signal_score,
                        strand=strand,
                        generator=generator,
                        rank=len(candidates) + 1,
                        source_region=f"{chrom}:{start}-{end}",
                        candidate_score=signal_score,
                        motif_id=motif_id,
                        motif_score=motif_score,
                        motif_source=source,
                    )
                )
    return candidates


def generate_candidates(
    signal: str | Path,
    peaks: str | Path,
    output: str | Path,
    candidate_width: int = 20,
    window: int = 10,
    min_score: float | None = None,
    top_n_per_region: int | None = 5,
    generator: str = "motif-free",
    motif_sites: list[str | Path] | None = None,
    motif_signal_math: str = "max",
    motif_generator: str = "motif-relaxed",
) -> list[Candidate]:
    """Generate motif-free signal candidates and optional motif-relaxed candidates."""

    regions = _read_bed_regions(peaks)
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_candidates: list[Candidate] = []

    with pyBigWig.open(str(signal)) as bw:
        chrom_sizes = bw.chroms()
        for chrom, start, end in regions:
            if chrom not in chrom_sizes:
                continue
            start = max(0, start)
            end = min(int(chrom_sizes[chrom]), end)
            if end <= start:
                continue
            values = np.asarray(bw.values(chrom, start, end), dtype=float)
            candidates = call_candidates_from_array(
                values,
                chrom,
                start,
                end,
                candidate_width=candidate_width,
                window=window,
                min_score=min_score,
                top_n=top_n_per_region,
                generator=generator,
            )
            all_candidates.extend(candidates)
        if motif_sites:
            all_candidates.extend(
                read_motif_site_candidates(
                    motif_sites,
                    regions,
                    bw,
                    chrom_sizes,
                    signal_math=motif_signal_math,
                    generator=motif_generator,
                )
            )

    all_candidates.sort(key=lambda c: (c.chrom, c.start, -c.score))
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "#chrom\tstart\tend\tname\tscore\tstrand\tgenerator\trank\tsource_region\t"
            "candidate_score\tmotif_id\tmotif_score\tmotif_source\n"
        )
        for candidate in all_candidates:
            handle.write("\t".join(candidate.to_bed_row()) + "\n")
    return all_candidates


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate motif-free footprint candidates from score bigWig signal.")
    parser.add_argument("--signal", required=True, help="Footprint or multiscale summary bigWig.")
    parser.add_argument("--peaks", required=True, help="Accessible peak/region BED file.")
    parser.add_argument("--out", required=True, help="Output BED-like candidate table.")
    parser.add_argument("--candidate-width", type=int, default=20, help="Candidate interval width in bp (default: 20).")
    parser.add_argument("--window", type=int, default=10, help="Local maximum window in bp (default: 10).")
    parser.add_argument("--min-score", type=float, default=None, help="Minimum score required for a candidate.")
    parser.add_argument("--top-n-per-region", type=int, default=5, help="Maximum candidates per input region (default: 5).")
    parser.add_argument("--generator", default="motif-free", help="Generator label written to output rows.")
    parser.add_argument("--motif-sites", nargs="*", default=None, help="Optional lower-threshold motif BED(s) to merge as motif-relaxed candidates.")
    parser.add_argument("--motif-signal-math", choices=["max", "mean"], default="max", help="Signal summary used to score motif-relaxed intervals (default: max).")
    parser.add_argument("--motif-generator", default="motif-relaxed", help="Generator label for rows from --motif-sites.")
    args = parser.parse_args(argv)

    candidates = generate_candidates(
        args.signal,
        args.peaks,
        args.out,
        candidate_width=args.candidate_width,
        window=args.window,
        min_score=args.min_score,
        top_n_per_region=args.top_n_per_region,
        generator=args.generator,
        motif_sites=args.motif_sites,
        motif_signal_math=args.motif_signal_math,
        motif_generator=args.motif_generator,
    )
    print(f"Wrote {len(candidates)} candidates to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
