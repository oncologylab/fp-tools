#!/usr/bin/env python3
"""Footprint-occupancy scoring from Tn5 cut sites in an ATAC BAM.

Builds a per-base Tn5 insertion (cut-site) count track for one chromosome from a
coordinate-sorted ATAC BAM, then for each peak finds the best PWM match and scores
the footprint as the cut-site depletion at the motif site relative to its flanks.
Unlike coverage tracks, cut-site density shows the canonical footprint dip at bound
sites, so this is the substrate footprinting actually needs. Writes a table with
motif, footprint, and combined (standardized) scores per peak.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent.parent / "src"))

from fp_tools.tools.variants import read_pwm_motifs, reverse_complement, score_pwm_window  # noqa: E402


def build_cutsites(bam_path: str, chrom: str, chrom_len: int, shift_fwd: int = 4, shift_rev: int = -5) -> np.ndarray:
    """Return a per-base Tn5 insertion-count array for one chromosome."""

    counts = np.zeros(chrom_len, dtype=np.int32)
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        if read.is_reverse:
            pos = read.reference_end - 1 + shift_rev
        else:
            pos = read.reference_start + shift_fwd
        if 0 <= pos < chrom_len:
            counts[pos] += 1
    bam.close()
    return counts


def best_match(sequence: str, motif) -> tuple[float, int]:
    width = len(motif.probabilities)
    seq = sequence.upper()
    best, best_off = float("-inf"), 0
    for idx in range(0, max(0, len(seq) - width + 1)):
        window = seq[idx:idx + width]
        s = max(score_pwm_window(window, motif), score_pwm_window(reverse_complement(window), motif))
        if s > best:
            best, best_off = s, idx
    return best, best_off


def footprint_score(counts: np.ndarray, lo: int, hi: int, flank: int) -> float:
    center = counts[lo:hi]
    left = counts[max(0, lo - flank):lo]
    right = counts[hi:hi + flank]
    flanks = np.concatenate([left, right])
    if not len(center) or not len(flanks):
        return float("nan")
    # Positive when the motif site is depleted relative to its flanks (a footprint).
    return float(flanks.mean() - center.mean())


def score_peaks(peaks, genome, motif_file, bam, output, chrom="chr4", motif_index=0, flank=50) -> int:
    motif = read_pwm_motifs(motif_file)[motif_index]
    width = len(motif.probabilities)
    fasta = pysam.FastaFile(str(genome))
    chrom_len = fasta.get_reference_length(chrom)
    counts = build_cutsites(str(bam), chrom, chrom_len)

    records = []
    with Path(peaks).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not f[1].isdigit() or f[0] != chrom:
                continue
            start, end = int(f[1]), int(f[2])
            name = f[3] if len(f) > 3 and f[3] not in (".", "") else f"peak_{len(records)+1}"
            seq = fasta.fetch(chrom, start, end)
            mscore, off = best_match(seq, motif)
            if mscore == float("-inf"):
                continue
            lo = start + off
            fos = footprint_score(counts, lo, lo + width, flank)
            records.append((chrom, start, end, name, mscore, fos))
    fasta.close()

    motif_scores = np.array([r[4] for r in records], float)
    fos_scores = np.array([r[5] for r in records], float)

    def z(a):
        s = np.nanstd(a)
        return (a - np.nanmean(a)) / s if s > 0 else np.zeros_like(a)

    combined = z(motif_scores) + z(np.nan_to_num(fos_scores, nan=np.nanmin(fos_scores)))

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("#chrom\tstart\tend\tname\tmotif_score\tfootprint_score\tcombined_score\n")
        for (c, s, e, name, ms, fos), comb in zip(records, combined):
            handle.write(f"{c}\t{s}\t{e}\t{name}\t{ms:.5f}\t{0.0 if math.isnan(fos) else fos:.5f}\t{comb:.5f}\n")
    return len(records)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--motif", required=True)
    parser.add_argument("--bam", required=True, help="Coordinate-sorted, indexed ATAC BAM.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--chrom", default="chr4")
    parser.add_argument("--motif-index", type=int, default=0)
    parser.add_argument("--flank", type=int, default=50)
    args = parser.parse_args(argv)
    n = score_peaks(args.peaks, args.genome, args.motif, args.bam, args.out, args.chrom, args.motif_index, args.flank)
    print(f"scored {n} peaks to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
