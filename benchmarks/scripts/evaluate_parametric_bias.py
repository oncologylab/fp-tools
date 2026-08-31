#!/usr/bin/env python3
"""Fit and evaluate fast conditional Tn5 sequence-bias models.

This is a research benchmark, not a production command.  It extracts
label-free enzyme-control windows, compares the legacy +4/-5 and aligned
+4/-4 cut conventions, and fits sample, pooled, and pooled-plus-adapted
conditional log-linear models.  ChIP labels and motif identities are never
used for fitting.

The output models use the package's checksummed NPZ + JSON format.  Test
chromosomes are deliberately unavailable here: model selection is based only
on the development and validation chromosome partitions in the locked study
specification.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import blake2b, sha256
import heapq
import json
from pathlib import Path
import sys
from time import perf_counter
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.special import rel_entr

try:
    import resource
except ImportError:  # pragma: no cover - Windows imports benchmark helpers only
    resource = None

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from fp_tools.tools.parametric_bias import (  # noqa: E402
    BiasFeatureSpec,
    ConditionalSequenceBiasModel,
    cut_position_from_alignment,
    encode_sequence,
    reverse_complement_contexts,
)
from fp_tools.utils.fasta import open_fasta  # noqa: E402
from fp_tools.utils.intervals import IntervalIndex  # noqa: E402


DATASET_SCHEMA = "fp-tools-control-windows-v1"
MODEL_NAMES = ("selma10", "loglinear81")
CONFIGURATIONS = ("sample_specific", "cross_cell_pooled", "pooled_plus_sample_adaptation")


def stable_u64(*values: object, seed: int = 2026) -> int:
    """Return a deterministic unsigned hash suitable for ranked sampling."""

    digest = blake2b(digest_size=8)
    digest.update(str(seed).encode())
    for value in values:
        digest.update(b"\0")
        digest.update(str(value).encode())
    return int.from_bytes(digest.digest(), "little")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cache_identity(paths: Iterable[str | Path | None], *parameters: object) -> str:
    """Key caches by input identity and every extraction-affecting parameter."""

    values: list[object] = list(parameters)
    for raw_path in paths:
        if raw_path is None:
            values.extend((None, None, None))
            continue
        path = Path(raw_path).expanduser().resolve()
        stat = path.stat()
        values.extend((str(path), int(stat.st_size), int(stat.st_mtime_ns)))
    return f"{stable_u64(*values):016x}"


def parse_name_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("samples must use NAME=PATH")
    name, raw_path = value.split("=", 1)
    if not name or not raw_path:
        raise argparse.ArgumentTypeError("samples must use NAME=PATH")
    return name, Path(raw_path)


def parse_shift(value: str) -> tuple[int, int]:
    fields = value.replace("/", ",").split(",")
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("read shifts must use FORWARD,REVERSE")
    try:
        return int(fields[0]), int(fields[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("read shifts must contain two integers") from exc


def parse_depth(value: str) -> int | None:
    if value.lower() == "full":
        return None
    try:
        depth = int(value.replace(",", ""))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("training depths must be integers or 'full'") from exc
    if depth < 1:
        raise argparse.ArgumentTypeError("training depths must be positive")
    return depth


def cut_position(read, shift: tuple[int, int]) -> int:
    """Return the zero-based cut coordinate used by ``OneRead.get_cutsite``.

    Soft-clipped bases are included to preserve the current fp-tools cut
    convention exactly.  Subtracting one from OneRead's one-based coordinate
    gives these expressions.
    """

    return cut_position_from_alignment(read, shift)


def usable_read(read, *, keep_duplicates: bool = False, minimum_mapq: int = 1) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or (read.is_duplicate and not keep_duplicates)
        or int(read.mapping_quality) < minimum_mapq
    )


@dataclass
class ControlWindowDataset:
    """Sequence and strand-specific cut counts for fixed genomic windows."""

    sample: str
    split: str
    source: str
    shift: tuple[int, int]
    window_size: int
    margin: int
    chromosomes: np.ndarray
    starts: np.ndarray
    sequences: np.ndarray
    forward_counts: np.ndarray
    reverse_counts: np.ndarray
    gc_fraction: np.ndarray

    def __post_init__(self) -> None:
        n_windows = len(self.starts)
        width = int(self.window_size)
        expected_sequence_length = width + 2 * int(self.margin)
        if self.forward_counts.shape != (n_windows, width):
            raise ValueError("forward count matrix does not match dataset dimensions")
        if self.reverse_counts.shape != (n_windows, width):
            raise ValueError("reverse count matrix does not match dataset dimensions")
        if len(self.chromosomes) != n_windows or len(self.sequences) != n_windows:
            raise ValueError("window metadata does not match count matrices")
        if any(len(str(value)) != expected_sequence_length for value in self.sequences):
            raise ValueError("control-window sequences have an unexpected length")

    @property
    def cuts(self) -> int:
        return int(self.forward_counts.sum() + self.reverse_counts.sum())

    def subset(self, indices: Sequence[int] | np.ndarray, split: str | None = None) -> "ControlWindowDataset":
        selected = np.asarray(indices, dtype=int)
        return ControlWindowDataset(
            sample=self.sample,
            split=self.split if split is None else split,
            source=self.source,
            shift=self.shift,
            window_size=self.window_size,
            margin=self.margin,
            chromosomes=self.chromosomes[selected],
            starts=self.starts[selected],
            sequences=self.sequences[selected],
            forward_counts=self.forward_counts[selected],
            reverse_counts=self.reverse_counts[selected],
            gc_fraction=self.gc_fraction[selected],
        )

    def model_arrays(self, spec: BiasFeatureSpec) -> tuple[np.ndarray, np.ndarray]:
        """Build forward and reverse-complemented candidate-cut contexts."""

        context_length = int(spec.context_length)
        left = context_length // 2
        first = int(self.margin) - left
        if first < 0:
            raise ValueError("dataset margin is too small for the requested context")
        encoded = np.stack([encode_sequence(str(sequence)) for sequence in self.sequences])
        windows = np.lib.stride_tricks.sliding_window_view(encoded, context_length, axis=1)
        forward = np.asarray(windows[:, first:first + self.window_size], dtype=np.uint8)
        reverse = reverse_complement_contexts(forward)
        contexts = np.concatenate([forward, reverse], axis=0)
        counts = np.concatenate([self.forward_counts, self.reverse_counts], axis=0).astype(float)
        valid = np.all(contexts < 4, axis=2)
        counts[~valid] = 0.0
        keep = (counts.sum(axis=1) > 0) & valid.any(axis=1)
        return contexts[keep], counts[keep]

    def save(self, path: str | Path, metadata: dict | None = None) -> tuple[Path, Path]:
        npz_path = Path(path)
        if npz_path.suffix != ".npz":
            npz_path = npz_path.with_suffix(".npz")
        json_path = npz_path.with_suffix(".json")
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz_path,
            chromosomes=np.asarray(self.chromosomes, dtype=str),
            starts=np.asarray(self.starts, dtype=np.int64),
            sequences=np.asarray(self.sequences, dtype=str),
            forward_counts=np.asarray(self.forward_counts),
            reverse_counts=np.asarray(self.reverse_counts),
            gc_fraction=np.asarray(self.gc_fraction, dtype=np.float32),
        )
        document = {
            "schema": DATASET_SCHEMA,
            "sample": self.sample,
            "split": self.split,
            "source": self.source,
            "shift": list(self.shift),
            "window_size": int(self.window_size),
            "margin": int(self.margin),
            "windows": int(len(self.starts)),
            "cuts": self.cuts,
            "npz_sha256": file_sha256(npz_path),
            "metadata": dict(metadata or {}),
        }
        json_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return npz_path, json_path

    @classmethod
    def load(cls, path: str | Path) -> "ControlWindowDataset":
        npz_path = Path(path)
        if npz_path.suffix != ".npz":
            npz_path = npz_path.with_suffix(".npz")
        document = json.loads(npz_path.with_suffix(".json").read_text(encoding="utf-8"))
        if document.get("schema") != DATASET_SCHEMA:
            raise ValueError("unsupported control-window dataset schema")
        if document.get("npz_sha256") != file_sha256(npz_path):
            raise ValueError("control-window checksum does not match its metadata")
        with np.load(npz_path, allow_pickle=False) as arrays:
            return cls(
                sample=str(document["sample"]),
                split=str(document["split"]),
                source=str(document["source"]),
                shift=tuple(int(value) for value in document["shift"]),
                window_size=int(document["window_size"]),
                margin=int(document["margin"]),
                chromosomes=np.asarray(arrays["chromosomes"], dtype=str),
                starts=np.asarray(arrays["starts"], dtype=np.int64),
                sequences=np.asarray(arrays["sequences"], dtype=str),
                forward_counts=np.asarray(arrays["forward_counts"]),
                reverse_counts=np.asarray(arrays["reverse_counts"]),
                gc_fraction=np.asarray(arrays["gc_fraction"], dtype=np.float32),
            )


def gc_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    valid = sum(sequence.count(base) for base in "ACGT")
    return float((sequence.count("G") + sequence.count("C")) / valid) if valid else np.nan


def gc_matched_indices(
    candidate_gc: np.ndarray,
    target_gc: np.ndarray,
    maximum: int,
    *,
    bins: int = 20,
    seed: int = 2026,
) -> np.ndarray:
    """Select candidates to approximate a target GC distribution."""

    candidate = np.asarray(candidate_gc, dtype=float)
    target = np.asarray(target_gc, dtype=float)
    valid_candidate = np.flatnonzero(np.isfinite(candidate))
    target = target[np.isfinite(target)]
    maximum = min(int(maximum), len(valid_candidate))
    if maximum <= 0:
        return np.asarray([], dtype=int)
    hashes = np.asarray([stable_u64("gc", int(index), seed=seed) for index in valid_candidate], dtype=np.uint64)
    if len(target) == 0:
        return np.sort(valid_candidate[np.argsort(hashes)[:maximum]])

    edges = np.linspace(0.0, 1.0, int(bins) + 1)
    target_bin = np.clip(np.digitize(target, edges[1:-1]), 0, bins - 1)
    candidate_bin = np.clip(np.digitize(candidate[valid_candidate], edges[1:-1]), 0, bins - 1)
    target_counts = np.bincount(target_bin, minlength=bins).astype(float)
    desired = target_counts / max(target_counts.sum(), 1.0) * maximum
    allocation = np.floor(desired).astype(int)
    remainder_order = np.argsort(-(desired - allocation))
    for bin_index in remainder_order[: maximum - int(allocation.sum())]:
        allocation[bin_index] += 1

    selected: list[int] = []
    selected_set: set[int] = set()
    for bin_index in range(bins):
        rows = np.flatnonzero(candidate_bin == bin_index)
        ranked = rows[np.argsort(hashes[rows])]
        for local_index in ranked[: allocation[bin_index]]:
            value = int(valid_candidate[local_index])
            selected.append(value)
            selected_set.add(value)
    if len(selected) < maximum:
        ranked_all = valid_candidate[np.argsort(hashes)]
        for value in ranked_all:
            value = int(value)
            if value not in selected_set:
                selected.append(value)
                selected_set.add(value)
                if len(selected) == maximum:
                    break
    return np.sort(np.asarray(selected, dtype=int))


def _ranked_candidate_windows(
    bam_path: Path,
    chromosomes: Sequence[str],
    shift: tuple[int, int],
    window_size: int,
    margin: int,
    maximum: int,
    exclusions: Sequence[IntervalIndex],
    chromosome_lengths: dict[str, int],
    *,
    seed: int,
    keep_duplicates: bool,
    minimum_mapq: int,
) -> list[tuple[str, int]]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - benchmark is Linux-only
        raise RuntimeError("pysam is required to extract bias-control windows") from exc

    references: set[str]
    selected: list[tuple[int, str, int]] = []
    per_chromosome = max(1, int(np.ceil(maximum / max(len(chromosomes), 1))))
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        references = set(bam.references)
        for chromosome in chromosomes:
            if chromosome not in references or chromosome not in chromosome_lengths:
                continue
            seen: set[int] = set()
            heap: list[tuple[int, int]] = []
            length = chromosome_lengths[chromosome]
            for read in bam.fetch(chromosome):
                if not usable_read(read, keep_duplicates=keep_duplicates, minimum_mapq=minimum_mapq):
                    continue
                cut = cut_position(read, shift)
                start = (cut // window_size) * window_size
                if start in seen:
                    continue
                seen.add(start)
                if start < margin or start + window_size + margin > length:
                    continue
                if any(index.overlaps(chromosome, start - margin, start + window_size + margin) for index in exclusions):
                    continue
                rank = stable_u64(chromosome, start, shift, seed=seed)
                item = (-rank, start)
                if len(heap) < per_chromosome:
                    heapq.heappush(heap, item)
                elif item > heap[0]:
                    heapq.heapreplace(heap, item)
            selected.extend((-negative_rank, chromosome, start) for negative_rank, start in heap)
    selected.sort()
    return [(chromosome, start) for _rank, chromosome, start in selected[:maximum]]


def _count_selected_windows(
    bam_path: Path,
    windows: Sequence[tuple[str, int]],
    shift: tuple[int, int],
    window_size: int,
    *,
    keep_duplicates: bool,
    minimum_mapq: int,
) -> tuple[np.ndarray, np.ndarray]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pysam is required to count bias-control cuts") from exc

    forward = np.zeros((len(windows), window_size), dtype=np.uint32)
    reverse = np.zeros_like(forward)
    by_chromosome: dict[str, dict[int, int]] = {}
    for index, (chromosome, start) in enumerate(windows):
        by_chromosome.setdefault(chromosome, {})[int(start)] = index
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for chromosome, start_to_row in by_chromosome.items():
            for read in bam.fetch(chromosome):
                if not usable_read(read, keep_duplicates=keep_duplicates, minimum_mapq=minimum_mapq):
                    continue
                cut = cut_position(read, shift)
                start = (cut // window_size) * window_size
                row = start_to_row.get(start)
                if row is None:
                    continue
                offset = cut - start
                if 0 <= offset < window_size:
                    (reverse if read.is_reverse else forward)[row, offset] += 1
    return forward, reverse


def extract_control_windows(
    sample: str,
    split: str,
    source: str,
    bam_path: str | Path,
    genome: str | Path,
    chromosomes: Sequence[str],
    shift: tuple[int, int],
    *,
    window_size: int,
    margin: int,
    maximum_windows: int,
    candidate_factor: int = 4,
    exclusions: Sequence[IntervalIndex] = (),
    target_gc: np.ndarray | None = None,
    low_signal_quantile: float | None = 0.75,
    seed: int = 2026,
    keep_duplicates: bool = False,
    minimum_mapq: int = 1,
) -> ControlWindowDataset:
    """Extract deterministic cut-containing control windows from one BAM."""

    bam_path = Path(bam_path)
    with open_fasta(genome) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        candidates = _ranked_candidate_windows(
            bam_path,
            chromosomes,
            shift,
            window_size,
            margin,
            max(maximum_windows, maximum_windows * int(candidate_factor)),
            exclusions,
            lengths,
            seed=seed,
            keep_duplicates=keep_duplicates,
            minimum_mapq=minimum_mapq,
        )
        if not candidates:
            raise ValueError(f"no eligible {source} windows were found for {sample} {split}")
        forward, reverse = _count_selected_windows(
            bam_path,
            candidates,
            shift,
            window_size,
            keep_duplicates=keep_duplicates,
            minimum_mapq=minimum_mapq,
        )
        sequences = np.asarray(
            [fasta.fetch(chromosome, start - margin, start + window_size + margin).upper() for chromosome, start in candidates],
            dtype=str,
        )
    core_sequences = [sequence[margin:margin + window_size] for sequence in sequences]
    gc_values = np.asarray([gc_fraction(sequence) for sequence in core_sequences], dtype=float)
    total = forward.sum(axis=1) + reverse.sum(axis=1)
    keep = (total > 0) & np.isfinite(gc_values)
    if low_signal_quantile is not None and np.any(keep):
        threshold = float(np.quantile(total[keep], low_signal_quantile))
        keep &= total <= max(threshold, 1.0)
    eligible = np.flatnonzero(keep)
    if len(eligible) == 0:
        raise ValueError(f"no nonzero low-signal windows remained for {sample} {split}")
    matched_local = gc_matched_indices(
        gc_values[eligible],
        np.asarray([] if target_gc is None else target_gc),
        maximum_windows,
        seed=stable_u64(sample, split, shift, seed=seed) % (2**32 - 1),
    )
    chosen = eligible[matched_local]
    return ControlWindowDataset(
        sample=sample,
        split=split,
        source=source,
        shift=shift,
        window_size=window_size,
        margin=margin,
        chromosomes=np.asarray([candidates[index][0] for index in chosen], dtype=str),
        starts=np.asarray([candidates[index][1] for index in chosen], dtype=np.int64),
        sequences=sequences[chosen],
        forward_counts=forward[chosen],
        reverse_counts=reverse[chosen],
        gc_fraction=gc_values[chosen],
    )


def sample_peak_gc(
    genome: str | Path,
    peaks: str | Path,
    chromosomes: Sequence[str],
    *,
    window_size: int,
    maximum: int = 5000,
    seed: int = 2026,
) -> np.ndarray:
    """Sample the GC distribution of development peaks without occupancy labels."""

    chromosome_set = set(chromosomes)
    heap: list[tuple[int, str, int, int]] = []
    with Path(peaks).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3 or fields[0] not in chromosome_set:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            rank = stable_u64(fields[0], start, end, "peak-gc", seed=seed)
            item = (-rank, fields[0], start, end)
            if len(heap) < maximum:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
    values: list[float] = []
    with open_fasta(genome) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for _rank, chromosome, start, end in heap:
            center = (start + end) // 2
            left = center - window_size // 2
            right = left + window_size
            if left < 0 or right > lengths.get(chromosome, 0):
                continue
            value = gc_fraction(fasta.fetch(chromosome, left, right))
            if np.isfinite(value):
                values.append(value)
    return np.asarray(values, dtype=float)


def split_mitochondrial_dataset(
    dataset: ControlWindowDataset,
    training_windows: int,
    validation_windows: int,
    *,
    seed: int,
) -> tuple[ControlWindowDataset, ControlWindowDataset]:
    """Hash-split chrM windows when chromosome splitting is impossible."""

    order = np.argsort(
        np.asarray(
            [stable_u64(chromosome, int(start), "mt-split", seed=seed) for chromosome, start in zip(dataset.chromosomes, dataset.starts)],
            dtype=np.uint64,
        )
    )
    requested = min(len(order), int(training_windows) + int(validation_windows))
    order = order[:requested]
    validation_n = min(int(validation_windows), max(1, int(round(requested * 0.2))))
    training_n = min(int(training_windows), requested - validation_n)
    if training_n < 1 or validation_n < 1:
        raise ValueError("mitochondrial control requires at least two eligible windows")
    return (
        dataset.subset(order[:training_n], split="train"),
        dataset.subset(order[training_n:training_n + validation_n], split="validation"),
    )


def thin_counts(counts: np.ndarray, target_cuts: int | None, *, seed: int) -> np.ndarray:
    observed = np.asarray(counts, dtype=float)
    total = float(observed.sum())
    if target_cuts is None or target_cuts >= total:
        return observed.copy()
    probability = float(target_cuts / total)
    rng = np.random.default_rng(seed)
    return rng.binomial(np.rint(observed).astype(np.int64), probability).astype(float)


def conditional_metrics(
    model: ConditionalSequenceBiasModel,
    contexts: np.ndarray,
    counts: np.ndarray,
) -> dict[str, float | int]:
    """Calculate held-out likelihood, deviance, calibration, and JSD."""

    observed = np.asarray(counts, dtype=float)
    valid = np.all(contexts < 4, axis=2)
    observed = np.where(valid, observed, 0.0)
    totals = observed.sum(axis=1)
    keep = (totals > 0) & valid.any(axis=1)
    contexts = contexts[keep]
    observed = observed[keep]
    totals = totals[keep]
    probabilities = model.probabilities(contexts)
    expected = probabilities * totals[:, None]
    nll = model.conditional_nll(contexts, observed)
    valid_counts = valid[keep].sum(axis=1)
    null_nll = float(np.sum(totals * np.log(valid_counts)) / np.sum(totals))
    with np.errstate(divide="ignore", invalid="ignore"):
        saturated = np.where(observed > 0, observed * np.log(observed / totals[:, None]), 0.0)
        fitted = np.where(observed > 0, observed * np.log(np.maximum(probabilities, 1e-300)), 0.0)
    deviance = float(2.0 * np.sum(saturated - fitted) / np.sum(totals))

    aggregate_observed = observed.sum(axis=0)
    aggregate_expected = expected.sum(axis=0)
    aggregate_observed /= max(float(aggregate_observed.sum()), 1.0)
    aggregate_expected /= max(float(aggregate_expected.sum()), 1.0)
    midpoint = 0.5 * (aggregate_observed + aggregate_expected)
    jsd = 0.5 * (
        np.sum(rel_entr(aggregate_observed, midpoint))
        + np.sum(rel_entr(aggregate_expected, midpoint))
    )

    positive = expected > 0
    ratio = np.divide(observed, expected, out=np.zeros_like(observed), where=positive)
    probability_flat = probabilities[positive]
    ratio_flat = ratio[positive]
    if len(probability_flat) >= 10:
        deciles = np.quantile(probability_flat, np.linspace(0.0, 1.0, 11))
        decile = np.clip(np.digitize(probability_flat, np.unique(deciles)[1:-1]), 0, 9)
        calibration_error = float(
            np.mean(
                [abs(np.mean(ratio_flat[decile == group]) - 1.0) for group in np.unique(decile)]
            )
        )
    else:
        calibration_error = np.nan
    return {
        "windows": int(len(contexts)),
        "cuts": int(np.sum(observed)),
        "conditional_nll": nll,
        "null_nll": null_nll,
        "nll_gain": null_nll - nll,
        "multinomial_deviance_per_cut": deviance,
        "aggregate_jsd": float(jsd),
        "calibration_error": calibration_error,
    }


def aggregate_cut_motif(
    contexts: np.ndarray,
    counts: np.ndarray,
    *,
    sample: str,
    split: str,
    shift: tuple[int, int],
    model: str,
) -> pd.DataFrame:
    """Return the strand-aligned sequence motif around observed cuts."""

    flat_contexts = contexts.reshape(-1, contexts.shape[-1])
    flat_counts = counts.reshape(-1)
    base_counts = np.zeros((contexts.shape[-1], 4), dtype=float)
    for position in range(contexts.shape[-1]):
        codes = flat_contexts[:, position]
        valid = codes < 4
        base_counts[position] = np.bincount(
            codes[valid], weights=flat_counts[valid], minlength=4
        )[:4]
    probabilities = base_counts / np.maximum(base_counts.sum(axis=1, keepdims=True), 1.0)
    center = contexts.shape[-1] // 2
    rows = []
    for position in range(contexts.shape[-1]):
        for base, probability in zip("ACGT", probabilities[position]):
            rows.append(
                {
                    "sample": sample,
                    "split": split,
                    "shift_forward": int(shift[0]),
                    "shift_reverse": int(shift[1]),
                    "model": model,
                    "relative_position": int(position - center),
                    "base": base,
                    "probability": float(probability),
                    "weighted_cuts": float(base_counts[position].sum()),
                }
            )
    return pd.DataFrame(rows)


def feature_spec(name: str) -> BiasFeatureSpec:
    if name == "selma10":
        return BiasFeatureSpec.selma10()
    if name == "loglinear81":
        return BiasFeatureSpec.loglinear81()
    raise ValueError(f"unknown parametric bias model: {name}")


def _fit_one(
    spec: BiasFeatureSpec,
    contexts: np.ndarray,
    counts: np.ndarray,
    *,
    l2: float,
    epochs: int,
    batch_windows: int,
    seed: int,
    prior: ConditionalSequenceBiasModel | None = None,
    prior_strength: float = 0.0,
) -> tuple[ConditionalSequenceBiasModel, float, float]:
    before_memory = (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource is not None else 0.0
    )
    started = perf_counter()
    model = ConditionalSequenceBiasModel(spec).fit(
        contexts,
        counts,
        epochs=epochs,
        batch_windows=batch_windows,
        learning_rate=0.03,
        l2=l2,
        prior=prior,
        prior_strength=prior_strength,
        seed=seed,
        patience=15,
    )
    runtime = perf_counter() - started
    after_memory = (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource is not None else before_memory
    )
    memory_mb = max(0.0, float(after_memory - before_memory) / 1024.0)
    return model, runtime, memory_mb


def _safe_token(value: object) -> str:
    return str(value).replace("-", "m").replace(".", "p").replace("/", "_")


def evaluate_models(
    datasets: dict[tuple[tuple[int, int], str, str], ControlWindowDataset],
    outdir: Path,
    *,
    models: Sequence[str],
    l2_values: Sequence[float],
    depths: Sequence[int | None],
    seeds: Sequence[int],
    epochs: int,
    batch_windows: int,
    adaptation_strength: float,
    source: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Fit the full sample/pooled/adapted development factorial."""

    metric_rows: list[dict] = []
    artifact_rows: list[dict] = []
    motif_frames: list[pd.DataFrame] = []
    shifts = sorted({key[0] for key in datasets})
    samples = sorted({key[1] for key in datasets})
    model_dir = outdir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)

    for shift in shifts:
        for model_name in models:
            spec = feature_spec(model_name)
            train_arrays = {
                sample: datasets[(shift, sample, "train")].model_arrays(spec)
                for sample in samples
            }
            validation_arrays = {
                sample: datasets[(shift, sample, "validation")].model_arrays(spec)
                for sample in samples
            }
            for sample in samples:
                contexts, counts = train_arrays[sample]
                motif_frames.append(
                    aggregate_cut_motif(
                        contexts,
                        counts,
                        sample=sample,
                        split="train",
                        shift=shift,
                        model=model_name,
                    )
                )
            for l2 in l2_values:
                for depth in depths:
                    depth_name = "full" if depth is None else str(depth)
                    for seed in seeds:
                        thinned = {
                            sample: thin_counts(
                                train_arrays[sample][1],
                                depth,
                                seed=stable_u64(sample, shift, model_name, l2, depth_name, seed=seed) % (2**32 - 1),
                            )
                            for sample in samples
                        }
                        pooled_contexts = np.concatenate([train_arrays[sample][0] for sample in samples])
                        pooled_counts = np.concatenate([thinned[sample] for sample in samples])
                        pooled, pooled_runtime, pooled_memory = _fit_one(
                            spec,
                            pooled_contexts,
                            pooled_counts,
                            l2=l2,
                            epochs=epochs,
                            batch_windows=batch_windows,
                            seed=stable_u64("pooled", shift, model_name, l2, depth_name, seed=seed) % (2**32 - 1),
                        )
                        pooled_stem = model_dir / (
                            f"pooled.{source}.{model_name}.shift_{shift[0]}_{shift[1]}."
                            f"l2_{_safe_token(l2)}.depth_{depth_name}.seed_{seed}"
                        )
                        pooled_npz, _ = pooled.save(
                            pooled_stem,
                            metadata={
                                "training_source": source,
                                "configuration": "cross_cell_pooled",
                                "read_shift": list(shift),
                                "training_depth_cuts_per_sample": depth_name,
                                "samples": samples,
                            },
                        )
                        pooled_size = pooled_npz.stat().st_size / (1024 * 1024)
                        artifact_rows.append(
                            {
                                "configuration": "cross_cell_pooled",
                                "sample": "pooled",
                                "shift_forward": shift[0],
                                "shift_reverse": shift[1],
                                "model": model_name,
                                "l2": l2,
                                "training_depth": depth_name,
                                "seed": seed,
                                "model_npz": str(pooled_npz),
                                "model_json": str(pooled_npz.with_suffix(".json")),
                                "runtime_seconds": pooled_runtime,
                                "peak_memory_increment_mb": pooled_memory,
                                "model_size_mb": pooled_size,
                            }
                        )

                        for sample in samples:
                            sample_contexts, _sample_counts = train_arrays[sample]
                            sample_counts = thinned[sample]
                            sample_seed = stable_u64(sample, shift, model_name, l2, depth_name, seed=seed) % (2**32 - 1)
                            independent, independent_runtime, independent_memory = _fit_one(
                                spec,
                                sample_contexts,
                                sample_counts,
                                l2=l2,
                                epochs=epochs,
                                batch_windows=batch_windows,
                                seed=sample_seed,
                            )
                            adapted, adapted_runtime, adapted_memory = _fit_one(
                                spec,
                                sample_contexts,
                                sample_counts,
                                l2=l2,
                                epochs=epochs,
                                batch_windows=batch_windows,
                                seed=sample_seed,
                                prior=pooled,
                                prior_strength=adaptation_strength,
                            )
                            configurations = (
                                ("sample_specific", independent, independent_runtime, independent_memory),
                                ("cross_cell_pooled", pooled, pooled_runtime, pooled_memory),
                                ("pooled_plus_sample_adaptation", adapted, adapted_runtime, adapted_memory),
                            )
                            for configuration, fitted, runtime, memory_mb in configurations:
                                if configuration == "cross_cell_pooled":
                                    model_npz = pooled_npz
                                    model_size = pooled_size
                                else:
                                    stem = model_dir / (
                                        f"{sample}.{source}.{configuration}.{model_name}."
                                        f"shift_{shift[0]}_{shift[1]}.l2_{_safe_token(l2)}."
                                        f"depth_{depth_name}.seed_{seed}"
                                    )
                                    model_npz, _ = fitted.save(
                                        stem,
                                        metadata={
                                            "training_source": source,
                                            "configuration": configuration,
                                            "read_shift": list(shift),
                                            "training_depth_cuts": depth_name,
                                            "sample": sample,
                                            "pooled_prior": str(pooled_npz) if configuration.endswith("adaptation") else None,
                                        },
                                    )
                                    model_size = model_npz.stat().st_size / (1024 * 1024)
                                    artifact_rows.append(
                                        {
                                            "configuration": configuration,
                                            "sample": sample,
                                            "shift_forward": shift[0],
                                            "shift_reverse": shift[1],
                                            "model": model_name,
                                            "l2": l2,
                                            "training_depth": depth_name,
                                            "seed": seed,
                                            "model_npz": str(model_npz),
                                            "model_json": str(model_npz.with_suffix(".json")),
                                            "runtime_seconds": runtime,
                                            "peak_memory_increment_mb": memory_mb,
                                            "model_size_mb": model_size,
                                        }
                                    )
                                for split in ("train", "validation"):
                                    contexts, counts = (
                                        (sample_contexts, sample_counts)
                                        if split == "train"
                                        else validation_arrays[sample]
                                    )
                                    metrics = conditional_metrics(fitted, contexts, counts)
                                    metric_rows.append(
                                        {
                                            "source": source,
                                            "sample": sample,
                                            "split": split,
                                            "shift_forward": shift[0],
                                            "shift_reverse": shift[1],
                                            "model": model_name,
                                            "configuration": configuration,
                                            "l2": l2,
                                            "training_depth": depth_name,
                                            "seed": seed,
                                            "runtime_seconds": runtime,
                                            "peak_memory_increment_mb": memory_mb,
                                            "model_size_mb": model_size,
                                            "model_npz": str(model_npz),
                                            **metrics,
                                        }
                                    )
    metrics = pd.DataFrame(metric_rows)
    artifacts = pd.DataFrame(artifact_rows).drop_duplicates(subset=["model_npz"])
    motifs = pd.concat(motif_frames, ignore_index=True).drop_duplicates(
        subset=["sample", "split", "shift_forward", "shift_reverse", "model", "relative_position", "base"]
    )
    return metrics, artifacts, motifs


def select_bias_configurations(metrics: pd.DataFrame, maximum_model_size_mb: float = 25.0) -> pd.DataFrame:
    """Rank development configurations using validation likelihood only."""

    validation = metrics[(metrics["split"] == "validation") & (metrics["model_size_mb"] <= maximum_model_size_mb)].copy()
    group_columns = [
        "source",
        "shift_forward",
        "shift_reverse",
        "model",
        "configuration",
        "l2",
        "training_depth",
        "seed",
    ]
    ranked = (
        validation.groupby(group_columns, as_index=False)
        .agg(
            mean_conditional_nll=("conditional_nll", "mean"),
            mean_nll_gain=("nll_gain", "mean"),
            mean_deviance=("multinomial_deviance_per_cut", "mean"),
            mean_calibration_error=("calibration_error", "mean"),
            mean_runtime_seconds=("runtime_seconds", "mean"),
            maximum_model_size_mb=("model_size_mb", "max"),
            evaluated_samples=("sample", "nunique"),
        )
        .sort_values(
            ["mean_conditional_nll", "mean_calibration_error", "mean_runtime_seconds"],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )
    ranked.insert(0, "rank", np.arange(1, len(ranked) + 1))
    ranked["passed_control_likelihood"] = ranked["mean_nll_gain"] > 0
    passing_order = ranked["passed_control_likelihood"].cumsum()
    ranked["retained_for_functional_screen"] = (
        ranked["passed_control_likelihood"] & (passing_order <= 2)
    )
    return ranked


def _cache_path(cache_dir: Path, sample: str, source: str, split: str, shift: tuple[int, int]) -> Path:
    return cache_dir / f"{sample}.{source}.{split}.shift_{shift[0]}_{shift[1]}.npz"


def build_or_load_datasets(args: argparse.Namespace, study: dict) -> tuple[dict, pd.DataFrame]:
    datasets: dict[tuple[tuple[int, int], str, str], ControlWindowDataset] = {}
    manifest_rows: list[dict] = []
    cache_dir = Path(args.cache_dir or (args.outdir / "control_window_cache"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    exclusions: list[IntervalIndex] = []
    if args.source == "gc_matched_low_signal_nonpeak":
        if args.peaks is None:
            raise ValueError("--peaks is required for GC-matched nonpeak controls")
        exclusions.append(IntervalIndex.from_bed(args.peaks))
        if args.blacklist is not None:
            exclusions.append(IntervalIndex.from_bed(args.blacklist))

    train_chromosomes = list(study["chromosome_split"]["train"])
    validation_chromosomes = list(study["chromosome_split"]["validation"])
    if args.source == "mitochondrial":
        train_gc = validation_gc = np.asarray([])
    elif args.peaks is not None:
        train_gc = sample_peak_gc(
            args.genome,
            args.peaks,
            train_chromosomes,
            window_size=args.window_size,
            seed=args.seed,
        )
        validation_gc = sample_peak_gc(
            args.genome,
            args.peaks,
            validation_chromosomes,
            window_size=args.window_size,
            seed=args.seed,
        )
    else:
        train_gc = validation_gc = np.asarray([])

    for shift in args.read_shifts:
        for sample, bam_path in args.sample:
            identity = cache_identity(
                (bam_path, args.genome, args.peaks, args.blacklist),
                args.source,
                shift,
                args.window_size,
                args.margin,
                args.train_windows,
                args.validation_windows,
                args.candidate_factor,
                args.low_signal_quantile,
                args.minimum_mapq,
                args.keep_duplicates,
                args.seed,
            )
            sample_cache_dir = cache_dir / f"{sample}.{identity}"
            if args.source == "mitochondrial":
                combined_cache = _cache_path(sample_cache_dir, sample, args.source, "combined", shift)
                if combined_cache.is_file() and combined_cache.with_suffix(".json").is_file():
                    combined = ControlWindowDataset.load(combined_cache)
                else:
                    combined = extract_control_windows(
                        sample,
                        "combined",
                        args.source,
                        bam_path,
                        args.genome,
                        args.mitochondrial_chromosomes,
                        shift,
                        window_size=args.window_size,
                        margin=args.margin,
                        maximum_windows=args.train_windows + args.validation_windows,
                        candidate_factor=1,
                        exclusions=(),
                        target_gc=None,
                        low_signal_quantile=None,
                        seed=args.seed,
                        keep_duplicates=args.keep_duplicates,
                        minimum_mapq=args.minimum_mapq,
                    )
                    combined.save(combined_cache, {"bam": str(bam_path), "genome": str(args.genome)})
                train, validation = split_mitochondrial_dataset(
                    combined,
                    args.train_windows,
                    args.validation_windows,
                    seed=args.seed,
                )
                for dataset in (train, validation):
                    cache_path = _cache_path(sample_cache_dir, sample, args.source, dataset.split, shift)
                    dataset.save(cache_path, {"derived_from": str(combined_cache)})
                    datasets[(shift, sample, dataset.split)] = dataset
            else:
                for split, chromosomes, maximum, target_gc in (
                    ("train", train_chromosomes, args.train_windows, train_gc),
                    ("validation", validation_chromosomes, args.validation_windows, validation_gc),
                ):
                    cache_path = _cache_path(sample_cache_dir, sample, args.source, split, shift)
                    if cache_path.is_file() and cache_path.with_suffix(".json").is_file():
                        dataset = ControlWindowDataset.load(cache_path)
                    else:
                        dataset = extract_control_windows(
                            sample,
                            split,
                            args.source,
                            bam_path,
                            args.genome,
                            chromosomes,
                            shift,
                            window_size=args.window_size,
                            margin=args.margin,
                            maximum_windows=maximum,
                            candidate_factor=args.candidate_factor,
                            exclusions=exclusions,
                            target_gc=target_gc,
                            low_signal_quantile=(
                                args.low_signal_quantile
                                if args.source == "gc_matched_low_signal_nonpeak"
                                else None
                            ),
                            seed=args.seed,
                            keep_duplicates=args.keep_duplicates,
                            minimum_mapq=args.minimum_mapq,
                        )
                        dataset.save(
                            cache_path,
                            {
                                "bam": str(bam_path),
                                "genome": str(args.genome),
                                "peaks": str(args.peaks) if args.peaks else None,
                                "blacklist": str(args.blacklist) if args.blacklist else None,
                            },
                        )
                    datasets[(shift, sample, split)] = dataset

            for split in ("train", "validation"):
                dataset = datasets[(shift, sample, split)]
                manifest_rows.append(
                    {
                        "sample": sample,
                        "source": args.source,
                        "split": split,
                        "shift_forward": shift[0],
                        "shift_reverse": shift[1],
                        "windows": len(dataset.starts),
                        "cuts": dataset.cuts,
                        "mean_gc": float(np.mean(dataset.gc_fraction)),
                        "cache_npz": str(_cache_path(sample_cache_dir, sample, args.source, split, shift)),
                        "cache_identity": identity,
                        "bam": str(bam_path),
                    }
                )
    return datasets, pd.DataFrame(manifest_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--sample", type=parse_name_path, action="append", required=True, help="NAME=coordinate-sorted-BAM")
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--peaks", type=Path)
    parser.add_argument("--blacklist", type=Path)
    parser.add_argument(
        "--source",
        choices=("gc_matched_low_signal_nonpeak", "mitochondrial", "naked_dna"),
        default="gc_matched_low_signal_nonpeak",
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--read-shift", dest="read_shifts", type=parse_shift, action="append")
    parser.add_argument("--model", dest="models", choices=MODEL_NAMES, action="append")
    parser.add_argument("--l2", dest="l2_values", type=float, action="append")
    parser.add_argument("--training-depth", dest="depths", type=parse_depth, action="append")
    parser.add_argument("--fit-seed", dest="seeds", type=int, action="append")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--window-size", type=int, default=200)
    parser.add_argument("--margin", type=int, default=41)
    parser.add_argument("--train-windows", type=int, default=1200)
    parser.add_argument("--validation-windows", type=int, default=400)
    parser.add_argument("--candidate-factor", type=int, default=4)
    parser.add_argument("--low-signal-quantile", type=float, default=0.75)
    parser.add_argument("--minimum-mapq", type=int, default=30)
    parser.add_argument("--keep-duplicates", action="store_true")
    parser.add_argument("--mitochondrial-chromosome", dest="mitochondrial_chromosomes", action="append")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-windows", type=int, default=64)
    parser.add_argument("--adaptation-strength", type=float, default=0.1)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.read_shifts = args.read_shifts or [(4, -5), (4, -4)]
    args.models = args.models or list(MODEL_NAMES)
    args.l2_values = args.l2_values or [0.001]
    args.depths = args.depths or [None]
    args.seeds = args.seeds or [args.seed]
    args.mitochondrial_chromosomes = args.mitochondrial_chromosomes or ["chrM", "MT"]
    if args.window_size < 20 or args.margin < 41:
        raise ValueError("window size must be >=20 and margin must be >=41")
    if not 0 < args.low_signal_quantile <= 1:
        raise ValueError("low-signal quantile must be in (0, 1]")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    if study.get("status") != "development_locked_holdout_unscored":
        raise ValueError("parametric model selection requires the locked, unscored study")
    if any(chromosome in study["chromosome_split"]["test"] for chromosome in study["chromosome_split"]["train"] + study["chromosome_split"]["validation"]):
        raise ValueError("study chromosome partitions overlap the locked test set")

    args.outdir.mkdir(parents=True, exist_ok=True)
    datasets, window_manifest = build_or_load_datasets(args, study)
    metrics, artifacts, motifs = evaluate_models(
        datasets,
        args.outdir,
        models=args.models,
        l2_values=args.l2_values,
        depths=args.depths,
        seeds=args.seeds,
        epochs=args.epochs,
        batch_windows=args.batch_windows,
        adaptation_strength=args.adaptation_strength,
        source=args.source,
    )
    selection = select_bias_configurations(
        metrics,
        maximum_model_size_mb=float(study["promotion_gates"]["maximum_model_size_mb"]),
    )
    window_manifest.to_csv(args.outdir / "control_windows.tsv", sep="\t", index=False)
    metrics.to_csv(args.outdir / "bias_model_metrics.tsv", sep="\t", index=False)
    artifacts.to_csv(args.outdir / "bias_model_artifacts.tsv", sep="\t", index=False)
    motifs.to_csv(args.outdir / "strand_aligned_cut_motifs.tsv.gz", sep="\t", index=False)
    selection.to_csv(args.outdir / "bias_model_selection.tsv", sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-parametric-bias-benchmark-v1",
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "source": args.source,
        "genome": str(args.genome),
        "peaks": str(args.peaks) if args.peaks else None,
        "blacklist": str(args.blacklist) if args.blacklist else None,
        "samples": [{"name": name, "bam": str(path)} for name, path in args.sample],
        "read_shifts": [list(value) for value in args.read_shifts],
        "models": args.models,
        "l2_values": args.l2_values,
        "training_depths": ["full" if value is None else value for value in args.depths],
        "seeds": args.seeds,
        "test_chromosomes_scored": False,
        "retained_configurations": selection[selection["retained_for_functional_screen"]].to_dict("records"),
        "outputs": {
            name: {"path": str(args.outdir / name), "sha256": file_sha256(args.outdir / name)}
            for name in (
                "control_windows.tsv",
                "bias_model_metrics.tsv",
                "bias_model_artifacts.tsv",
                "strand_aligned_cut_motifs.tsv.gz",
                "bias_model_selection.tsv",
            )
        },
    }
    (args.outdir / "parametric_bias_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(selection.head(10).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
