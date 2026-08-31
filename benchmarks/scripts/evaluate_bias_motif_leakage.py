#!/usr/bin/env python3
"""Measure whether a control-trained bias model produces TF-motif responses.

The benchmark never reads ChIP labels.  It compares the model's log-propensity
profile at JASPAR motif sites with a nearby sequence control, reports broad
center/shoulder response and uncertainty, and measures cross-model/cross-cell
concordance.  A detected response is a review flag rather than automatic proof
of leakage because genuine Tn5 sequence preference can overlap a TF motif.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
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

from fp_tools.tools.functional_footprints import profile_descriptors  # noqa: E402
from fp_tools.tools.parametric_bias import (  # noqa: E402
    ConditionalSequenceBiasModel,
    encode_sequence,
    reverse_complement_contexts,
)
from fp_tools.utils.fasta import open_fasta  # noqa: E402


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*values: object, seed: int = 2026) -> int:
    digest = blake2b(digest_size=8)
    digest.update(str(seed).encode())
    for value in values:
        digest.update(b"\0")
        digest.update(str(value).encode())
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def parse_cell_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("motif sites must use CELL=PATH")
    cell, path = value.split("=", 1)
    return cell, Path(path)


def parse_model(value: str) -> tuple[tuple[str, str], Path]:
    if "=" not in value or ":" not in value.split("=", 1)[0]:
        raise argparse.ArgumentTypeError("models must use CELL:LABEL=MODEL.npz")
    identity, raw_path = value.split("=", 1)
    cell, label = identity.split(":", 1)
    if not cell or not label or not raw_path:
        raise argparse.ArgumentTypeError("models must use CELL:LABEL=MODEL.npz")
    return (cell, label), Path(raw_path)


def select_unlabeled_motif_sites(
    source: pd.DataFrame,
    motif_id: str,
    train_chromosomes: set[str],
    maximum: int,
    *,
    seed: int,
) -> pd.DataFrame:
    if any("chip" in column.lower() or "label" in column.lower() for column in source.columns):
        raise ValueError("motif-leakage inputs must not contain ChIP or label columns")
    required = {"motif", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError("motif sites are missing columns: " + ", ".join(sorted(missing)))
    selected = source[
        source["motif"].astype(str).str.contains(str(motif_id), regex=False)
        & source["TFBS_chr"].astype(str).isin(train_chromosomes)
    ].copy()
    if len(selected) > maximum:
        hashes = pd.util.hash_pandas_object(
            selected[["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"]],
            index=False,
        ).to_numpy(dtype=np.uint64)
        order = np.argsort(hashes ^ np.uint64(seed), kind="mergesort")[:maximum]
        selected = selected.iloc[order]
    return selected.reset_index(drop=True)


def fetch_motif_and_local_control_sequences(
    sites: pd.DataFrame,
    genome: str | Path,
    *,
    flank: int,
    margin: int,
    control_offset: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    width = flank * 2 + 1
    sequence_length = width + 2 * margin
    motif_sequences = np.full(len(sites), "N" * sequence_length, dtype=f"U{sequence_length}")
    control_sequences = motif_sequences.copy()
    valid = np.zeros(len(sites), dtype=bool)
    with open_fasta(genome) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for index, row in enumerate(sites.itertuples(index=False)):
            chromosome = str(row.TFBS_chr)
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            direction = 1 if stable_seed(chromosome, center, seed=seed) % 2 else -1
            control_center = center + direction * control_offset
            motif_start = center - flank - margin
            control_start = control_center - flank - margin
            if (
                chromosome not in lengths
                or motif_start < 0
                or control_start < 0
                or motif_start + sequence_length > lengths[chromosome]
                or control_start + sequence_length > lengths[chromosome]
            ):
                direction *= -1
                control_center = center + direction * control_offset
                control_start = control_center - flank - margin
            if (
                chromosome not in lengths
                or motif_start < 0
                or control_start < 0
                or motif_start + sequence_length > lengths[chromosome]
                or control_start + sequence_length > lengths[chromosome]
            ):
                continue
            motif_sequence = fasta.fetch(chromosome, motif_start, motif_start + sequence_length).upper()
            control_sequence = fasta.fetch(chromosome, control_start, control_start + sequence_length).upper()
            if len(motif_sequence) != sequence_length or len(control_sequence) != sequence_length:
                continue
            motif_sequences[index] = motif_sequence
            control_sequences[index] = control_sequence
            valid[index] = all(base in "ACGT" for base in motif_sequence + control_sequence)
    return motif_sequences, control_sequences, valid


def score_sequence_profiles(
    model: ConditionalSequenceBiasModel,
    sequences: np.ndarray,
    *,
    width: int,
    margin: int,
    batch_size: int = 128,
) -> np.ndarray:
    """Vectorize combined-strand log-bias scoring across many sequences."""

    context_length = model.feature_spec.context_length
    left = context_length // 2
    first = margin - left
    if first < 0:
        raise ValueError("margin is too small for the model context")
    output = np.full((len(sequences), width), np.nan, dtype=np.float32)
    for start in range(0, len(sequences), batch_size):
        batch = sequences[start:start + batch_size]
        encoded = np.stack([encode_sequence(str(sequence)) for sequence in batch])
        windows = np.lib.stride_tricks.sliding_window_view(encoded, context_length, axis=1)
        forward = np.asarray(windows[:, first:first + width], dtype=np.uint8)
        reverse = reverse_complement_contexts(forward)
        forward_score = model.log_scores(forward)
        reverse_score = model.log_scores(reverse)
        maximum = np.maximum(forward_score, reverse_score)
        combined = maximum + np.log(
            0.5 * np.exp(forward_score - maximum) + 0.5 * np.exp(reverse_score - maximum)
        )
        combined[~np.isfinite(combined)] = np.nan
        output[start:start + len(batch)] = combined.astype(np.float32)
    return output


def orient_score_profiles(profiles: np.ndarray, strands: Sequence[str]) -> np.ndarray:
    output = np.asarray(profiles, dtype=float).copy()
    reverse = np.asarray([str(value) == "-" for value in strands], dtype=bool)
    output[reverse] = output[reverse, ::-1]
    return output


def center_flank_effect(profiles: np.ndarray, positions: np.ndarray) -> np.ndarray:
    center = np.abs(positions) <= 15
    flanks = (np.abs(positions) >= 40) & (np.abs(positions) <= 80)
    return np.mean(profiles[:, flanks], axis=1) - np.mean(profiles[:, center], axis=1)


def summarize_response(
    motif_profiles: np.ndarray,
    control_profiles: np.ndarray,
    positions: np.ndarray,
    *,
    bootstraps: int,
    seed: int,
    review_threshold: float,
) -> tuple[dict[str, float | bool], pd.DataFrame]:
    """Summarize motif-minus-local-control model response and uncertainty."""

    motif = np.asarray(motif_profiles, dtype=float)
    control = np.asarray(control_profiles, dtype=float)
    if motif.shape != control.shape or motif.ndim != 2:
        raise ValueError("motif and control profiles must be equal two-dimensional arrays")
    outer = np.abs(positions) >= 70
    motif -= np.mean(motif[:, outer], axis=1, keepdims=True)
    control -= np.mean(control[:, outer], axis=1, keepdims=True)
    response = motif - control
    effect = center_flank_effect(response, positions)
    rng = np.random.default_rng(seed)
    sampled_effect = np.empty(bootstraps, dtype=float)
    sampled_profiles = np.empty((bootstraps, response.shape[1]), dtype=np.float32)
    for index in range(bootstraps):
        selected = rng.integers(0, len(response), size=len(response))
        sampled_effect[index] = np.mean(effect[selected])
        sampled_profiles[index] = np.mean(response[selected], axis=0)
    effect_lower, effect_upper = np.quantile(sampled_effect, [0.025, 0.975])
    mean_response = np.mean(response, axis=0)
    lower, upper = np.quantile(sampled_profiles, [0.025, 0.975], axis=0)
    descriptors = profile_descriptors(mean_response, positions)
    effect_mean = float(np.mean(effect))
    detected = bool(
        abs(effect_mean) >= review_threshold
        and (effect_lower > 0 or effect_upper < 0)
    )
    summary = {
        "sites": int(len(response)),
        "center_flank_log_bias_effect": effect_mean,
        "effect_lower_95": float(effect_lower),
        "effect_upper_95": float(effect_upper),
        "maximum_absolute_response": float(np.max(np.abs(mean_response))),
        "potential_motif_response_requires_review": detected,
        **asdict(descriptors),
    }
    curves = pd.DataFrame(
        {
            "position": positions.astype(int),
            "motif_mean_log_bias": np.mean(motif, axis=0),
            "local_control_mean_log_bias": np.mean(control, axis=0),
            "response": mean_response,
            "response_lower_95": lower,
            "response_upper_95": upper,
        }
    )
    return summary, curves


def cross_model_concordance(curves: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tf, group in curves.groupby("tf", sort=True):
        profiles = {
            (str(cell), str(label)): values.sort_values("position")["response"].to_numpy(dtype=float)
            for (cell, label), values in group.groupby(["cell", "model_label"], sort=True)
        }
        identities = sorted(profiles)
        for first_index, first in enumerate(identities):
            for second in identities[first_index + 1:]:
                rows.append(
                    {
                        "tf": tf,
                        "cell_a": first[0],
                        "model_a": first[1],
                        "cell_b": second[0],
                        "model_b": second[1],
                        "response_correlation": float(np.corrcoef(profiles[first], profiles[second])[0, 1]),
                    }
                )
    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--motif-sites", type=parse_cell_path, action="append", required=True)
    parser.add_argument("--model", type=parse_model, action="append", required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--control-offset", type=int, default=250)
    parser.add_argument("--maximum-sites-per-tf", type=int, default=1000)
    parser.add_argument("--bootstraps", type=int, default=500)
    parser.add_argument("--review-threshold", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    train_chromosomes = set(study["chromosome_split"]["train"])
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == "development"]
    sources = {cell: pd.read_csv(path, sep="\t") for cell, path in args.motif_sites}
    models = dict(args.model)
    summary_rows = []
    curve_frames = []
    positions = np.arange(-args.flank, args.flank + 1, dtype=float)
    width = len(positions)
    for (cell, label), model_path in sorted(models.items()):
        if cell not in sources:
            raise ValueError(f"no motif-site source was supplied for {cell}")
        model = ConditionalSequenceBiasModel.load(model_path)
        margin = max(41, model.feature_spec.context_length // 2 + 1)
        for task in tasks[tasks["cell"] == cell].itertuples(index=False):
            sites = select_unlabeled_motif_sites(
                sources[cell],
                str(task.motif_id),
                train_chromosomes,
                args.maximum_sites_per_tf,
                seed=stable_seed(cell, task.tf, label, seed=args.seed),
            )
            if len(sites) < 50:
                continue
            motif_sequences, control_sequences, valid = fetch_motif_and_local_control_sequences(
                sites,
                args.genome,
                flank=args.flank,
                margin=margin,
                control_offset=args.control_offset,
                seed=args.seed,
            )
            sites = sites.loc[valid].reset_index(drop=True)
            motif_sequences = motif_sequences[valid]
            control_sequences = control_sequences[valid]
            if len(sites) < 50:
                continue
            motif_profiles = score_sequence_profiles(
                model, motif_sequences, width=width, margin=margin
            )
            control_profiles = score_sequence_profiles(
                model, control_sequences, width=width, margin=margin
            )
            strands = sites["TFBS_strand"].astype(str).tolist()
            motif_profiles = orient_score_profiles(motif_profiles, strands)
            control_profiles = orient_score_profiles(control_profiles, strands)
            finite = np.isfinite(motif_profiles).all(axis=1) & np.isfinite(control_profiles).all(axis=1)
            if np.sum(finite) < 50:
                continue
            summary, curves = summarize_response(
                motif_profiles[finite],
                control_profiles[finite],
                positions,
                bootstraps=args.bootstraps,
                seed=stable_seed(cell, task.tf, label, "bootstrap", seed=args.seed),
                review_threshold=args.review_threshold,
            )
            summary_rows.append(
                {
                    "cell": cell,
                    "tf": str(task.tf),
                    "motif_family": str(task.motif_family),
                    "motif_id": str(task.motif_id),
                    "model_label": label,
                    "model_path": str(model_path),
                    **summary,
                }
            )
            curves.insert(0, "model_label", label)
            curves.insert(0, "motif_family", str(task.motif_family))
            curves.insert(0, "tf", str(task.tf))
            curves.insert(0, "cell", cell)
            curve_frames.append(curves)
    summary = pd.DataFrame(summary_rows)
    curves = pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame()
    concordance = cross_model_concordance(curves) if len(curves) else pd.DataFrame()
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "bias_motif_response_summary.tsv"
    curves_path = args.outdir / "bias_motif_response_curves.tsv.gz"
    concordance_path = args.outdir / "bias_motif_response_concordance.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    concordance.to_csv(concordance_path, sep="\t", index=False)
    document = {
        "schema": "fp-tools-bias-motif-leakage-v1",
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "genome": str(args.genome),
        "models": [
            {"cell": cell, "label": label, "path": str(path), "sha256": file_sha256(path)}
            for (cell, label), path in args.model
        ],
        "chIP_labels_read": False,
        "review_threshold_log_propensity": args.review_threshold,
        "important_interpretation": (
            "A flagged motif response requires review; sequence-only Tn5 preference is not by itself occupancy leakage."
        ),
        "outputs": {
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "curves": {"path": str(curves_path), "sha256": file_sha256(curves_path)},
            "concordance": {"path": str(concordance_path), "sha256": file_sha256(concordance_path)},
        },
    }
    (args.outdir / "bias_motif_response_manifest.json").write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
