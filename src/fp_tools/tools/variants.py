#!/usr/bin/env python
"""Footprint-aware variant scoring scaffold."""

from __future__ import annotations

import argparse
import math
import pickle
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from fp_tools.tools.motif_discovery import load_genome_fasta
from fp_tools.tools.tfbs_model import MODEL_VERSION


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    start: int
    end: int
    name: str
    ref: str
    alt: str


@dataclass(frozen=True)
class PwmMotif:
    motif_id: str
    name: str
    probabilities: list[dict[str, float]]


def read_variants(path: str | Path) -> list[VariantRecord]:
    """Read BED-like variants with columns chrom, start, end, name, ref, alt."""

    variants: list[VariantRecord] = []
    with Path(path).open(encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError("Variants must have at least 6 BED-like columns: chrom start end name ref alt")
            variants.append(VariantRecord(fields[0], int(fields[1]), int(fields[2]), fields[3] or f"variant_{idx}", fields[4].upper(), fields[5].upper()))
    return variants


def read_score_intervals(path: str | Path | None) -> dict[str, list[tuple[int, int, str, float]]]:
    """Read BED-like scored intervals keyed by chromosome."""

    intervals: dict[str, list[tuple[int, int, str, float]]] = {}
    if path is None:
        return intervals
    with Path(path).open(encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            name = fields[3] if len(fields) > 3 and fields[3] else f"interval_{idx}"
            try:
                score = float(fields[4]) if len(fields) > 4 else 0.0
            except ValueError:
                score = 0.0
            intervals.setdefault(fields[0], []).append((int(fields[1]), int(fields[2]), name, score))
    return intervals


def best_overlap(variant: VariantRecord, intervals: dict[str, list[tuple[int, int, str, float]]]) -> tuple[str, float, int]:
    """Return best overlapping interval name, score, and overlap bp."""

    best_name = "."
    best_score = 0.0
    best_bp = 0
    for start, end, name, score in intervals.get(variant.chrom, []):
        overlap = max(0, min(variant.end, end) - max(variant.start, start))
        if overlap > 0 and (score > best_score or best_name == "."):
            best_name = name
            best_score = float(score)
            best_bp = int(overlap)
    return best_name, best_score, best_bp


def _normalize_columns(rows: dict[str, list[float]], pseudocount: float = 0.1) -> list[dict[str, float]]:
    width = min(len(values) for values in rows.values())
    probabilities = []
    for idx in range(width):
        values = {base: float(rows[base][idx]) + pseudocount for base in "ACGT"}
        total = sum(values.values())
        probabilities.append({base: values[base] / total for base in "ACGT"})
    return probabilities


def read_pwm_motifs(path: str | Path) -> list[PwmMotif]:
    """Read simple JASPAR or MEME/PFM motifs for ref/alt PWM delta scoring."""

    text = Path(path).read_text(encoding="utf-8")
    motifs: list[PwmMotif] = []
    if "letter-probability matrix" in text:
        current_id = ""
        current_name = ""
        rows: list[list[float]] = []
        remaining = 0
        for raw_line in text.splitlines() + ["MOTIF __flush__"]:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("MOTIF"):
                if current_id and rows:
                    probabilities = [{base: max(float(row[idx]), 1e-9) for idx, base in enumerate("ACGT")} for row in rows if len(row) >= 4]
                    motifs.append(PwmMotif(current_id, current_name, probabilities))
                parts = line.split()
                current_id = parts[1] if len(parts) > 1 else "motif"
                current_name = parts[2] if len(parts) > 2 else current_id
                rows = []
                remaining = 0
                continue
            if line.startswith("letter-probability matrix"):
                match = re.search(r"\bw\s*=\s*(\d+)", line)
                remaining = int(match.group(1)) if match else 0
                continue
            if current_id and (remaining > 0 or line[0].isdigit()):
                values = [float(value) for value in line.split()[:4]]
                if len(values) == 4:
                    rows.append(values)
                    remaining = max(0, remaining - 1)
        return motifs

    motif_id = ""
    motif_name = ""
    rows: dict[str, list[float]] = {}
    for raw_line in text.splitlines() + [">__flush__"]:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if motif_id and set(rows) >= set("ACGT"):
                motifs.append(PwmMotif(motif_id, motif_name, _normalize_columns(rows)))
            parts = line[1:].split()
            motif_id = parts[0] if parts else "motif"
            motif_name = parts[1] if len(parts) > 1 else motif_id
            rows = {}
            continue
        base = line[0].upper()
        if base in "ACGT":
            cleaned = line.replace("[", " ").replace("]", " ")
            rows[base] = [float(value) for value in cleaned.split()[1:]]
    return motifs


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def score_pwm_window(sequence: str, motif: PwmMotif) -> float:
    score = 0.0
    for base, probabilities in zip(sequence.upper(), motif.probabilities):
        if base not in probabilities:
            return float("-inf")
        score += math.log2(probabilities[base] / 0.25)
    return float(score)


def best_pwm_score(sequence: str, motif: PwmMotif) -> float:
    width = len(motif.probabilities)
    sequence = sequence.upper()
    if width == 0 or len(sequence) < width:
        return float("nan")
    best = float("-inf")
    for idx in range(0, len(sequence) - width + 1):
        window = sequence[idx:idx + width]
        best = max(best, score_pwm_window(window, motif), score_pwm_window(reverse_complement(window), motif))
    return float(best)


def best_motif_delta(
    variant: VariantRecord,
    chrom_seq: str,
    motifs: list[PwmMotif],
    flank: int,
) -> dict[str, float | str]:
    """Return the motif with the largest absolute ref/alt PWM score delta."""

    if not motifs or not chrom_seq:
        return {}
    ref_context, alt_context, _, _ = variant_sequence_context(variant, chrom_seq, flank=flank)
    best: dict[str, float | str] | None = None
    for motif in motifs:
        ref_score = best_pwm_score(ref_context, motif)
        alt_score = best_pwm_score(alt_context, motif)
        if math.isnan(ref_score) or math.isnan(alt_score):
            continue
        delta = alt_score - ref_score
        if best is None or abs(delta) > abs(float(best["motif_delta_score"])):
            best = {
                "best_motif_id": motif.motif_id,
                "best_motif_name": motif.name,
                "ref_motif_score": ref_score,
                "alt_motif_score": alt_score,
                "motif_delta_score": delta,
                "motif_delta_direction": "gain" if delta > 0 else "loss" if delta < 0 else "unchanged",
            }
    return best or {}


def _model_feature_value(row: dict[str, object], feature: str, allele: str) -> object:
    if feature in {"motif_score", "pwm_score"}:
        return row.get(f"{allele}_motif_score", float("nan"))
    if feature in {"gc", "sequence_gc"}:
        return row.get(f"{allele}_gc", float("nan"))
    if feature in {"length", "site_length"}:
        context = row.get(f"{allele}_context", "")
        return len(str(context)) if context else row.get("end", 0) - row.get("start", 0)
    return row.get(feature, float("nan"))


def model_delta_features(row: dict[str, object], model_bundle: dict[str, object] | None) -> dict[str, float]:
    """Predict ref/alt TFBS probabilities from a saved fp-tools tabular model."""

    if model_bundle is None:
        return {}
    if model_bundle.get("version") != MODEL_VERSION:
        raise ValueError(f"Unsupported model version: {model_bundle.get('version')}")
    feature_columns = list(model_bundle["feature_columns"])  # type: ignore[index]
    ref_values = {feature: _model_feature_value(row, feature, "ref") for feature in feature_columns}
    alt_values = {feature: _model_feature_value(row, feature, "alt") for feature in feature_columns}
    model = model_bundle["model"]
    probabilities = model.predict_proba(pd.DataFrame([ref_values, alt_values], columns=feature_columns))[:, 1]  # type: ignore[attr-defined]
    ref_probability = float(probabilities[0])
    alt_probability = float(probabilities[1])
    return {
        "ref_model_probability": ref_probability,
        "alt_model_probability": alt_probability,
        "model_delta_probability": alt_probability - ref_probability,
    }



def gc_fraction(sequence: str) -> float:
    """Return GC fraction over non-N bases."""

    sequence = sequence.upper()
    informative = [base for base in sequence if base in {"A", "C", "G", "T"}]
    if not informative:
        return 0.0
    return float(sum(base in {"G", "C"} for base in informative) / len(informative))


def kmer_set(sequence: str, k: int) -> set[str]:
    """Return canonical exact k-mers without ambiguous bases."""

    k = int(k)
    if k <= 0 or len(sequence) < k:
        return set()
    sequence = sequence.upper()
    return {sequence[idx:idx + k] for idx in range(0, len(sequence) - k + 1) if set(sequence[idx:idx + k]) <= {"A", "C", "G", "T"}}


def variant_sequence_context(variant: VariantRecord, chrom_seq: str, flank: int = 0) -> tuple[str, str, int, int]:
    """Return reference and alternate sequence contexts around a variant."""

    flank = max(0, int(flank))
    context_start = max(0, variant.start - flank)
    context_end = min(len(chrom_seq), variant.end + flank)
    ref_context = chrom_seq[context_start:context_end].upper()
    rel_start = max(0, variant.start - context_start)
    rel_end = max(rel_start, variant.end - context_start)
    alt_context = (ref_context[:rel_start] + variant.alt + ref_context[rel_end:]).upper()
    return ref_context, alt_context, context_start, context_end


def sequence_delta_features(variant: VariantRecord, chrom_seq: str, flank: int = 0, kmer_size: int = 3) -> dict[str, float | int | str]:
    """Compute deterministic ref/alt sequence delta features for a variant."""

    ref_context, alt_context, context_start, context_end = variant_sequence_context(variant, chrom_seq, flank=flank)
    ref_gc = gc_fraction(ref_context)
    alt_gc = gc_fraction(alt_context)
    ref_kmers = kmer_set(ref_context, kmer_size)
    alt_kmers = kmer_set(alt_context, kmer_size)
    union = ref_kmers | alt_kmers
    shared = ref_kmers & alt_kmers
    return {
        "context_start": context_start,
        "context_end": context_end,
        "ref_context": ref_context,
        "alt_context": alt_context,
        "ref_gc": ref_gc,
        "alt_gc": alt_gc,
        "delta_gc": alt_gc - ref_gc,
        "allele_length_delta": len(variant.alt) - len(variant.ref),
        "kmer_size": int(kmer_size),
        "ref_kmers": len(ref_kmers),
        "alt_kmers": len(alt_kmers),
        "shared_kmers": len(shared),
        "lost_kmers": len(ref_kmers - alt_kmers),
        "gained_kmers": len(alt_kmers - ref_kmers),
        "kmer_jaccard": float(len(shared) / len(union)) if union else 1.0,
    }


def score_variants(
    variants_path: str | Path,
    genome_fasta: str | Path,
    output: str | Path,
    candidate_scores: str | Path | None = None,
    sequence_flank: int = 0,
    kmer_size: int = 3,
    motifs: list[str | Path] | None = None,
    motif_flank: int = 30,
    tfbs_model: str | Path | None = None,
) -> pd.DataFrame:
    """Score variants with allele checks, overlaps, sequence deltas, and optional PWM deltas."""

    variants = read_variants(variants_path)
    genome = load_genome_fasta(genome_fasta)
    intervals = read_score_intervals(candidate_scores)
    pwm_motifs: list[PwmMotif] = []
    for motif_path in motifs or []:
        pwm_motifs.extend(read_pwm_motifs(motif_path))
    model_bundle = None
    if tfbs_model is not None:
        with Path(tfbs_model).open("rb") as handle:
            model_bundle = pickle.load(handle)
    rows = []
    for variant in variants:
        chrom_seq = genome.get(variant.chrom, "")
        genome_ref = chrom_seq[variant.start:variant.end].upper() if chrom_seq else ""
        ref_matches = genome_ref == variant.ref
        interval_name, interval_score, overlap_bp = best_overlap(variant, intervals)
        sequence_features = sequence_delta_features(variant, chrom_seq, flank=sequence_flank, kmer_size=kmer_size) if chrom_seq else {}
        row = {
            "chrom": variant.chrom,
            "start": variant.start,
            "end": variant.end,
            "name": variant.name,
            "ref": variant.ref,
            "alt": variant.alt,
            "genome_ref": genome_ref,
            "ref_matches_genome": ref_matches,
            "candidate_name": interval_name,
            "candidate_score": interval_score,
            "candidate_overlap_bp": overlap_bp,
            "overlaps_candidate": overlap_bp > 0,
            "effect_hint": "candidate_overlap" if overlap_bp > 0 else "no_candidate_overlap",
        }
        row.update(sequence_features)
        row.update(best_motif_delta(variant, chrom_seq, pwm_motifs, flank=motif_flank))
        row.update(model_delta_features(row, model_bundle))
        rows.append(row)
    frame = pd.DataFrame(rows)
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, sep="\t", index=False)
    return frame


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Annotate variants with genome allele checks and footprint/candidate overlaps.")
    parser.add_argument("--variants", required=True, help="BED-like variants: chrom start end name ref alt.")
    parser.add_argument("--genome", required=True, help="Genome FASTA, optionally gzipped.")
    parser.add_argument("--out", required=True, help="Output TSV.")
    parser.add_argument("--candidate-scores", default=None, help="Optional BED-like scored candidates or footprint intervals.")
    parser.add_argument("--sequence-flank", type=int, default=0, help="Flanking bases on each side for ref/alt sequence-context delta features.")
    parser.add_argument("--kmer-size", type=int, default=3, help="K-mer size for exact ref/alt disruption features.")
    parser.add_argument("--motifs", nargs="*", default=[], help="Optional JASPAR/MEME motif files for best ref/alt PWM delta scoring.")
    parser.add_argument("--motif-flank", type=int, default=30, help="Flanking bases on each side for motif ref/alt delta scoring.")
    parser.add_argument("--tfbs-model", help="Optional fp-tools tabular TFBS model pickle for ref/alt probability deltas.")
    args = parser.parse_args(argv)

    frame = score_variants(
        args.variants,
        args.genome,
        args.out,
        candidate_scores=args.candidate_scores,
        sequence_flank=args.sequence_flank,
        kmer_size=args.kmer_size,
        motifs=args.motifs,
        motif_flank=args.motif_flank,
        tfbs_model=args.tfbs_model,
    )
    print(f"Wrote {len(frame)} scored variants to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
