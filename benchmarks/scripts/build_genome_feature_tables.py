#!/usr/bin/env python3
"""Build genome-wide per-peak TF-binding feature tables for many TFs of one cell.

This is the genome-scale generalization of ``build_tf_feature_table.py``. It walks
the accessible-peak universe of a single cell line across a chromosome whitelist
and, for each requested transcription factor, emits a feature table with columns
``chrom start end name accessibility motif gc footprint label`` -- the identical
schema the evaluator consumes, so chromosome-held-out evaluation is possible.

Two things make the genome-wide pass tractable:

* The PWM scan is vectorized with NumPy (a sliding-window log-odds scan over both
  strands) instead of the per-offset Python loop in ``footprint_from_bam.best_match``.
  It reproduces ``fp_tools.tools.variants.score_pwm_window`` exactly: each position
  contributes ``log2(p_base / 0.25)`` and any non-ACGT base makes the window -inf.
* Peaks are processed one chromosome at a time, so the Tn5 cut-site array for a
  chromosome is built once from the BAM and reused across every TF for that cell,
  then freed before the next chromosome.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import sys

import numpy as np
import pysam

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR.parent.parent / "src"))

from footprint_from_bam import build_cutsites, footprint_score  # noqa: E402
from build_label_overlap_benchmark import overlap_bp  # noqa: E402
from fp_tools.tools.variants import read_pwm_motifs  # noqa: E402

BASE_TO_CODE = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    BASE_TO_CODE[ord(_b)] = _i
    BASE_TO_CODE[ord(_b.lower())] = _i
# Complement code mapping A<->T (0<->3), C<->G (1<->2).
COMPLEMENT_CODE = np.array([3, 2, 1, 0], dtype=np.int8)

DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def log_odds_matrices(motif) -> tuple[np.ndarray, np.ndarray]:
    """Return (forward, reverse-complement) (width, 4) log2(p/0.25) matrices."""
    width = len(motif.probabilities)
    fwd = np.full((width, 4), -np.inf, dtype=np.float64)
    with np.errstate(divide="ignore"):
        for j, probs in enumerate(motif.probabilities):
            for b, base in enumerate("ACGT"):
                p = probs.get(base, 0.0)
                fwd[j, b] = np.log2(p / 0.25) if p > 0 else -np.inf
    # Reverse complement: position width-1-j, complemented base.
    rc = fwd[::-1][:, COMPLEMENT_CODE]
    return fwd, rc


def scan_best(codes: np.ndarray, fwd: np.ndarray, rc: np.ndarray) -> tuple[float, int]:
    """Best two-strand PWM score and forward offset over a coded sequence."""
    width = fwd.shape[0]
    n = codes.shape[0] - width + 1
    if n <= 0:
        return float("-inf"), 0
    invalid = codes < 0
    safe = np.where(invalid, 0, codes)
    fwd_s = np.zeros(n, dtype=np.float64)
    rc_s = np.zeros(n, dtype=np.float64)
    bad = np.zeros(n, dtype=bool)
    for j in range(width):
        cj = safe[j:j + n]
        fwd_s += fwd[j][cj]
        rc_s += rc[j][cj]
        bad |= invalid[j:j + n]
    fwd_s[bad] = -np.inf
    rc_s[bad] = -np.inf
    combined = np.maximum(fwd_s, rc_s)
    off = int(np.argmax(combined))
    return float(combined[off]), off


def gc_content(seq: str) -> float:
    seq = seq.upper()
    n = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / n if n else 0.0


def read_peaks_by_chrom(path: str | Path, chroms: set[str]) -> dict[str, list[tuple[int, int, float]]]:
    """Parse an ENCODE narrowPeak (gz) into {chrom: [(start, end, signalValue)]}.

    Accessibility is the narrowPeak signalValue (column 7); peaks off the whitelist
    are dropped.
    """
    out: dict[str, list[tuple[int, int, float]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] not in chroms or not f[1].isdigit():
                continue
            signal = float(f[6]) if len(f) >= 7 else float(f[4]) if len(f) >= 5 else 0.0
            out.setdefault(f[0], []).append((int(f[1]), int(f[2]), signal))
    for chrom in out:
        out[chrom].sort()
    return out


def read_labels_by_chrom(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    """Parse an ENCODE narrowPeak (gz or plain) into {chrom: [(start, end)]} labels."""
    out: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not f[1].isdigit():
                continue
            out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return out


def build_for_cell(
    atac_narrowpeak: str | Path,
    genome: str | Path,
    tasks: list[tuple[str, str, str]],
    outdir: str | Path,
    bam: str | Path | None = None,
    chroms: list[str] | None = None,
    flank: int = 50,
    motif_index: int = 0,
) -> dict[str, int]:
    chroms = chroms or DEFAULT_CHROMS
    chrom_set = set(chroms)
    fasta = pysam.FastaFile(str(genome))
    available = set(fasta.references)
    chroms = [c for c in chroms if c in available]

    peaks_by_chrom = read_peaks_by_chrom(atac_narrowpeak, chrom_set)

    # Per-TF motif matrices and label intervals.
    motifs = {}
    labels = {}
    handles = {}
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for tf, motif_file, chip in tasks:
        motif = read_pwm_motifs(motif_file)[motif_index]
        motifs[tf] = (log_odds_matrices(motif), len(motif.probabilities))
        labels[tf] = read_labels_by_chrom(chip)
        h = (outdir / f"{tf}.tsv").open("w", encoding="utf-8")
        h.write("chrom\tstart\tend\tname\taccessibility\tmotif\tgc\tfootprint\tlabel\n")
        handles[tf] = h

    written = {tf: 0 for tf, _, _ in tasks}
    for chrom in chroms:
        chrom_peaks = peaks_by_chrom.get(chrom, [])
        if not chrom_peaks:
            continue
        counts = None
        if bam is not None:
            counts = build_cutsites(str(bam), chrom, fasta.get_reference_length(chrom))
        for start, end, signal in chrom_peaks:
            seq = fasta.fetch(chrom, start, end)
            codes = BASE_TO_CODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
            gc = gc_content(seq)
            name = f"{chrom}:{start}-{end}"
            for tf, _, _ in tasks:
                (fwd, rc), width = motifs[tf]
                mscore, off = scan_best(codes, fwd, rc)
                if mscore == float("-inf"):
                    continue
                fos = 0.0
                if counts is not None:
                    lo = start + off
                    val = footprint_score(counts, lo, lo + width, flank)
                    fos = 0.0 if val != val else val
                label = int(overlap_bp(chrom, start, end, labels[tf]) >= 1)
                handles[tf].write(
                    f"{chrom}\t{start}\t{end}\t{name}\t{signal}\t{mscore}\t{gc}\t{fos}\t{label}\n"
                )
                written[tf] += 1
        del counts
    fasta.close()
    for h in handles.values():
        h.close()
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--atac-narrowpeak", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--bam", help="Optional ATAC BAM for cut-site footprint feature.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--chroms", nargs="*", default=DEFAULT_CHROMS)
    parser.add_argument("--flank", type=int, default=50)
    parser.add_argument("--motif-index", type=int, default=0)
    parser.add_argument(
        "--task", action="append", required=True, dest="tasks",
        help="TF:motif_file:chip_narrowpeak triple; repeatable.",
    )
    args = parser.parse_args(argv)

    tasks = []
    for spec in args.tasks:
        tf, motif_file, chip = spec.split(":", 2)
        tasks.append((tf, motif_file, chip))

    written = build_for_cell(
        args.atac_narrowpeak, args.genome, tasks, args.outdir,
        bam=args.bam, chroms=args.chroms, flank=args.flank, motif_index=args.motif_index,
    )
    for tf, n in written.items():
        print(f"{args.cell} {tf}: wrote {n} rows -> {Path(args.outdir) / (tf + '.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
