#!/usr/bin/env python3
"""Score peaks by motif match and by ATAC footprint depletion at the motif site.

For each peak this finds the best PWM match position (best two-strand log-odds),
then measures the ATAC signal depletion at that site relative to its flanks (a
footprint-occupancy score, FOS). It writes a table with the motif score, the
footprint score, and a combined standardized score, so each can be benchmarked
against TF-binding labels. The combined score tests whether footprint evidence
adds discrimination beyond motif presence alone.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np
import pyBigWig
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent.parent / "src"))

from fp_tools.tools.variants import read_pwm_motifs, reverse_complement, score_pwm_window  # noqa: E402


def best_match(sequence: str, motif) -> tuple[float, int]:
    """Return (best log-odds score, offset of motif start within sequence)."""

    width = len(motif.probabilities)
    seq = sequence.upper()
    best, best_off = float("-inf"), 0
    for idx in range(0, max(0, len(seq) - width + 1)):
        window = seq[idx:idx + width]
        fwd = score_pwm_window(window, motif)
        rev = score_pwm_window(reverse_complement(window), motif)
        s = max(fwd, rev)
        if s > best:
            best, best_off = s, idx
    return best, best_off


def footprint_score(signal: np.ndarray, center_lo: int, center_hi: int, flank: int) -> float:
    """Footprint-occupancy score: mean flank signal minus mean center signal."""

    n = len(signal)
    center = signal[max(0, center_lo):min(n, center_hi)]
    left = signal[max(0, center_lo - flank):max(0, center_lo)]
    right = signal[min(n, center_hi):min(n, center_hi + flank)]
    flanks = np.concatenate([left, right]) if (len(left) + len(right)) else np.array([])
    if not len(center) or not len(flanks):
        return float("nan")
    return float(np.nanmean(flanks) - np.nanmean(center))


def score_peaks(
    peaks: str | Path,
    genome: str | Path,
    motif_file: str | Path,
    signal_bw: str | Path,
    output: str | Path,
    motif_index: int = 0,
    flank: int = 50,
) -> int:
    """Write per-peak motif, footprint, and combined scores."""

    motif = read_pwm_motifs(motif_file)[motif_index]
    width = len(motif.probabilities)
    fasta = pysam.FastaFile(str(genome))
    bw = pyBigWig.open(str(signal_bw))
    chroms = set(fasta.references) & set(bw.chroms().keys())

    records = []
    with Path(peaks).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not f[1].isdigit():
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if chrom not in chroms:
                continue
            name = f[3] if len(f) > 3 and f[3] not in (".", "") else f"peak_{len(records)+1}"
            seq = fasta.fetch(chrom, start, end)
            mscore, off = best_match(seq, motif)
            if mscore == float("-inf"):
                continue
            # Genomic window around the motif site for the footprint measurement.
            site_lo = start + off
            site_hi = site_lo + width
            win_lo = max(0, site_lo - flank)
            win_hi = site_hi + flank
            values = np.nan_to_num(np.array(bw.values(chrom, win_lo, win_hi), dtype=float), nan=0.0)
            fos = footprint_score(values, site_lo - win_lo, site_hi - win_lo, flank)
            records.append((chrom, start, end, name, mscore, fos))
    fasta.close(); bw.close()

    motif_scores = np.array([r[4] for r in records], dtype=float)
    fos_scores = np.array([r[5] for r in records], dtype=float)

    def z(a):
        m, s = np.nanmean(a), np.nanstd(a)
        return (a - m) / s if s > 0 else np.zeros_like(a)

    combined = z(motif_scores) + z(np.nan_to_num(fos_scores, nan=np.nanmin(fos_scores)))

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("#chrom\tstart\tend\tname\tmotif_score\tfootprint_score\tcombined_score\n")
        for (chrom, start, end, name, ms, fos), comb in zip(records, combined):
            fos_out = fos if not math.isnan(fos) else 0.0
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{ms:.5f}\t{fos_out:.5f}\t{comb:.5f}\n")
    return len(records)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--motif", required=True)
    parser.add_argument("--signal", required=True, help="ATAC signal bigWig.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--motif-index", type=int, default=0)
    parser.add_argument("--flank", type=int, default=50)
    args = parser.parse_args(argv)
    n = score_peaks(args.peaks, args.genome, args.motif, args.signal, args.out, args.motif_index, args.flank)
    print(f"scored {n} peaks to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
