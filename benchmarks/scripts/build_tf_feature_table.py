#!/usr/bin/env python3
"""Build a per-peak TF-binding feature table for method comparison.

For each accessible peak this assembles the features used by competing scoring
strategies: raw accessibility magnitude, best PWM motif log-odds, GC content,
and (optionally) a Tn5 cut-site footprint-occupancy score at the motif site,
plus the binding label from matched ChIP overlap. The resulting table lets a
single accessibility/motif/footprint baseline and an integrated supervised model
be benchmarked on identical inputs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR.parent.parent / "src"))

from footprint_from_bam import best_match, build_cutsites, footprint_score  # noqa: E402
from build_label_overlap_benchmark import read_label_intervals, overlap_bp  # noqa: E402
from fp_tools.tools.variants import read_pwm_motifs  # noqa: E402


def read_accessibility(path: str | Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[1].isdigit():
                try:
                    scores[f[3]] = float(f[4])
                except ValueError:
                    continue
    return scores


def gc_content(seq: str) -> float:
    seq = seq.upper()
    n = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / n if n else 0.0


def build_table(
    peaks: str | Path,
    genome: str | Path,
    motif_file: str | Path,
    labels_bed: str | Path,
    output: str | Path,
    bam: str | Path | None = None,
    chrom: str = "chr4",
    motif_index: int = 0,
    flank: int = 50,
) -> int:
    motif = read_pwm_motifs(motif_file)[motif_index]
    width = len(motif.probabilities)
    fasta = pysam.FastaFile(str(genome))
    access = read_accessibility(peaks)
    labels = read_label_intervals(labels_bed)
    counts = build_cutsites(str(bam), chrom, fasta.get_reference_length(chrom)) if bam else None

    rows = []
    with Path(peaks).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or not f[1].isdigit() or f[0] != chrom:
                continue
            name = f[3]
            start, end = int(f[1]), int(f[2])
            seq = fasta.fetch(chrom, start, end)
            mscore, off = best_match(seq, motif)
            if mscore == float("-inf"):
                continue
            fos = 0.0
            if counts is not None:
                lo = start + off
                fos_val = footprint_score(counts, lo, lo + width, flank)
                fos = 0.0 if fos_val != fos_val else fos_val
            label = int(overlap_bp(f[0], start, end, labels) >= 1)
            rows.append((f[0], start, end, name, access.get(name, 0.0), mscore, gc_content(seq), fos, label))
    fasta.close()

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        handle.write("chrom\tstart\tend\tname\taccessibility\tmotif\tgc\tfootprint\tlabel\n")
        for r in rows:
            handle.write("\t".join(str(x) for x in r) + "\n")
    return len(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--motif", required=True)
    parser.add_argument("--labels-bed", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--bam", help="Optional ATAC BAM for cut-site footprint feature.")
    parser.add_argument("--chrom", default="chr4")
    parser.add_argument("--motif-index", type=int, default=0)
    parser.add_argument("--flank", type=int, default=50)
    args = parser.parse_args(argv)
    n = build_table(args.peaks, args.genome, args.motif, args.labels_bed, args.out,
                    bam=args.bam, chrom=args.chrom, motif_index=args.motif_index, flank=args.flank)
    print(f"wrote {n} feature rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
