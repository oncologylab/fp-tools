#!/usr/bin/env python3
"""Score BED-like peaks by their best PWM match over a genome FASTA.

Writes a scored BED (``#chrom start end name score``) where ``score`` is the best
log2-odds PWM match (both strands) anywhere inside each peak. This provides a
motif-aware predictor that can be benchmarked against an accessibility-only
baseline with the same label-overlap and metric scripts.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
SRC = SCRIPT_DIR.parent.parent / "src"
sys.path.insert(0, str(SRC))

from fp_tools.tools.variants import best_pwm_score, read_pwm_motifs  # noqa: E402


def read_peaks(path: str | Path, chroms: set[str] | None = None):
    """Yield (chrom, start, end, name) from a BED, optionally filtered by chrom."""

    with Path(path).open(encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[1].isdigit():
                continue
            chrom = fields[0]
            if chroms is not None and chrom not in chroms:
                continue
            name = fields[3] if len(fields) > 3 and fields[3] not in (".", "") else f"peak_{idx}"
            yield chrom, int(fields[1]), int(fields[2]), name


def score_peaks_with_pwm(
    peaks: str | Path,
    genome: str | Path,
    motif_file: str | Path,
    output: str | Path,
    motif_index: int = 0,
    max_width: int = 2000,
) -> int:
    """Score each peak by its best PWM match and write a scored BED."""

    motifs = read_pwm_motifs(motif_file)
    if not motifs:
        raise ValueError(f"No motifs parsed from {motif_file}")
    motif = motifs[motif_index]
    fasta = pysam.FastaFile(str(genome))
    available = set(fasta.references)

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with out.open("w", encoding="utf-8") as handle:
        handle.write("#chrom\tstart\tend\tname\tscore\n")
        for chrom, start, end, name in read_peaks(peaks, chroms=available):
            # Cap very wide peaks around their center to bound runtime.
            if end - start > max_width:
                center = (start + end) // 2
                fetch_start = max(0, center - max_width // 2)
                fetch_end = fetch_start + max_width
            else:
                fetch_start, fetch_end = start, end
            sequence = fasta.fetch(chrom, fetch_start, fetch_end)
            score = best_pwm_score(sequence, motif)
            if score != score:  # NaN (sequence shorter than motif)
                continue
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{score:.5f}\n")
            written += 1
    fasta.close()
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", required=True, help="BED-like peak intervals.")
    parser.add_argument("--genome", required=True, help="Indexed genome FASTA (.fa/.fa.gz with .fai).")
    parser.add_argument("--motif", required=True, help="JASPAR or MEME motif file.")
    parser.add_argument("--out", required=True, help="Output scored BED.")
    parser.add_argument("--motif-index", type=int, default=0, help="Which motif in the file to use.")
    parser.add_argument("--max-width", type=int, default=2000, help="Cap scanned width per peak (bp).")
    args = parser.parse_args(argv)

    n = score_peaks_with_pwm(
        args.peaks, args.genome, args.motif, args.out,
        motif_index=args.motif_index, max_width=args.max_width,
    )
    print(f"scored {n} peaks to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
