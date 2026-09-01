#!/usr/bin/env python3
"""Evaluate a frozen TF-specific geometry on naked-DNA motif sites."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from apply_frozen_dual_null_thresholds import wilson_interval  # noqa: E402
from compare_frozen_tf_candidates import score_centers  # noqa: E402
from search_tf_footprint_models import (  # noqa: E402
    candidate_from_row,
    extract_profiles,
    score_candidate,
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_cell_path(value: str) -> tuple[str, Path]:
    fields = value.split(",", 1)
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError("value must use CELL,PATH")
    return fields[0], Path(fields[1])


def parse_method_path(value: str) -> tuple[str, Path]:
    fields = value.split(",", 1)
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError("value must use METHOD,PATH")
    return fields[0], Path(fields[1])


def frozen_upper_tail_threshold(scores: np.ndarray, alpha: float) -> float:
    """Choose a deterministic threshold with at most alpha training calls."""

    values = np.asarray(scores, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        raise ValueError("cannot freeze a threshold from zero finite scores")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    maximum_calls = int(np.floor(alpha * len(values)))
    if maximum_calls == 0:
        return float(np.nextafter(np.max(values), np.inf))
    descending = np.sort(values)[::-1]
    threshold = float(descending[maximum_calls - 1])
    if int(np.sum(values >= threshold)) > maximum_calls:
        threshold = float(np.nextafter(threshold, np.inf))
    return threshold


def _rate_row(
    *,
    scores: np.ndarray,
    valid: np.ndarray,
    threshold: float,
    prefix: str,
) -> tuple[dict[str, float | int], np.ndarray]:
    valid = np.asarray(valid, dtype=bool) & np.isfinite(scores)
    calls = valid & (scores >= threshold)
    n_valid = int(valid.sum())
    n_calls = int(calls.sum())
    lower, upper = wilson_interval(n_calls, n_valid)
    return (
        {
            f"{prefix}_threshold": float(threshold),
            f"{prefix}_valid": n_valid,
            f"{prefix}_calls": n_calls,
            f"{prefix}_false_positive_rate": (
                n_calls / n_valid if n_valid else np.nan
            ),
            f"{prefix}_false_positive_rate_lower_95": lower,
            f"{prefix}_false_positive_rate_upper_95": upper,
        },
        calls,
    )


def evaluate(
    *,
    development_sites: pd.DataFrame,
    winners: pd.DataFrame,
    development_baselines: pd.DataFrame,
    profile_cache: Path,
    naked_sites: dict[str, pd.DataFrame],
    naked_corrected: dict[str, Path],
    naked_legacy: Path,
    tf: str,
    split: str,
    alpha: float,
    flank: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    site_frames: list[pd.DataFrame] = []
    cells = sorted(set(winners.loc[winners["tf"].astype(str).eq(tf), "cell"]))
    for cell in cells:
        cell_sites = development_sites[
            development_sites["cell"].astype(str).eq(cell)
        ].reset_index(drop=True)
        winner_rows = winners[
            winners["cell"].astype(str).eq(cell)
            & winners["tf"].astype(str).eq(tf)
        ]
        baseline_rows = development_baselines[
            development_baselines["cell"].astype(str).eq(cell)
        ]
        if len(winner_rows) != 1 or len(baseline_rows) != 1:
            raise ValueError(f"expected one winner and baseline for {cell}/{tf}")
        winner = winner_rows.iloc[0]
        candidate = candidate_from_row(winner)
        development_indexes = np.flatnonzero(
            cell_sites["tf"].astype(str).to_numpy() == tf
        )
        development_indexes = development_indexes[
            cell_sites.iloc[development_indexes]["chromosome_split"]
            .astype(str)
            .eq(split)
            .to_numpy()
        ]
        candidate_profiles = np.load(
            profile_cache / f"{cell}.{candidate.correction}.flank{flank}.npz"
        )["profiles"][development_indexes]
        candidate_development = score_candidate(candidate_profiles, candidate)
        baseline_development = score_centers(
            cell_sites, Path(baseline_rows.iloc[0].signal)
        )[development_indexes]
        labels = cell_sites.iloc[development_indexes]["chip_label"].to_numpy(
            dtype=int
        )
        negative = labels == 0
        candidate_threshold = frozen_upper_tail_threshold(
            candidate_development[negative], alpha
        )
        baseline_threshold = frozen_upper_tail_threshold(
            baseline_development[negative], alpha
        )

        if cell not in naked_sites:
            raise ValueError(f"no naked-DNA motif sites supplied for {cell}")
        if candidate.correction not in naked_corrected:
            raise ValueError(
                f"no naked-DNA {candidate.correction} corrected signal supplied"
            )
        sites = naked_sites[cell]
        sites = sites[sites["tf"].astype(str).eq(tf)].reset_index(drop=True)
        profiles, profile_valid = extract_profiles(
            sites, naked_corrected[candidate.correction], flank
        )
        candidate_scores = score_candidate(profiles, candidate)
        legacy_scores = score_centers(sites, naked_legacy)
        candidate_valid = profile_valid & np.isfinite(candidate_scores)
        legacy_valid = np.isfinite(legacy_scores)
        common_valid = (
            candidate_valid
            & legacy_valid
        )
        candidate_rate, candidate_calls = _rate_row(
            scores=candidate_scores,
            valid=candidate_valid,
            threshold=candidate_threshold,
            prefix="candidate",
        )
        baseline_rate, baseline_calls = _rate_row(
            scores=legacy_scores,
            valid=legacy_valid,
            threshold=baseline_threshold,
            prefix="legacy",
        )
        paired_candidate_rate, _ = _rate_row(
            scores=candidate_scores,
            valid=common_valid,
            threshold=candidate_threshold,
            prefix="paired_candidate",
        )
        paired_baseline_rate, _ = _rate_row(
            scores=legacy_scores,
            valid=common_valid,
            threshold=baseline_threshold,
            prefix="paired_legacy",
        )
        summary_rows.append(
            {
                "cell": cell,
                "tf": tf,
                "candidate_id": candidate.identifier,
                "candidate_correction": candidate.correction,
                "development_split": split,
                "development_negative_sites": int(negative.sum()),
                "target_false_positive_rate": float(alpha),
                **candidate_rate,
                **baseline_rate,
                **paired_candidate_rate,
                **paired_baseline_rate,
                "false_positive_rate_increase": (
                    paired_candidate_rate[
                        "paired_candidate_false_positive_rate"
                    ]
                    - paired_baseline_rate[
                        "paired_legacy_false_positive_rate"
                    ]
                ),
            }
        )
        site_frames.append(
            pd.DataFrame(
                {
                    "cell": cell,
                    "tf": tf,
                    "TFBS_chr": sites["TFBS_chr"],
                    "TFBS_start": sites["TFBS_start"],
                    "TFBS_end": sites["TFBS_end"],
                    "candidate_score": candidate_scores,
                    "legacy_score": legacy_scores,
                    "candidate_valid": candidate_valid,
                    "legacy_valid": legacy_valid,
                    "paired_valid": common_valid,
                    "candidate_call": candidate_calls,
                    "legacy_call": baseline_calls,
                }
            )
        )
    return pd.DataFrame(summary_rows), pd.concat(site_frames, ignore_index=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--development-baselines", type=Path, required=True)
    parser.add_argument("--profile-cache", type=Path, required=True)
    parser.add_argument(
        "--naked-sites",
        action="append",
        type=parse_cell_path,
        required=True,
        metavar="CELL,TSV",
    )
    parser.add_argument(
        "--naked-corrected",
        action="append",
        type=parse_method_path,
        required=True,
        metavar="METHOD,BIGWIG",
    )
    parser.add_argument("--naked-legacy", type=Path, required=True)
    parser.add_argument("--tf", default="CTCF")
    parser.add_argument(
        "--candidate-id",
        help="Optional exact candidate identifier to select from --winners.",
    )
    parser.add_argument("--winner-stage", help="Optional winner-table stage filter.")
    parser.add_argument(
        "--cell",
        action="append",
        help="Optional target cell to evaluate; repeat for multiple cells.",
    )
    parser.add_argument("--split", default="validation")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--maximum-rate-increase", type=float, default=0.01)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    naked_site_paths = dict(args.naked_sites)
    naked_corrected_paths = dict(args.naked_corrected)

    winners = pd.read_csv(args.winners, sep="\t")
    if args.cell:
        winners = winners[
            winners["cell"].astype(str).isin([str(cell) for cell in args.cell])
        ].reset_index(drop=True)
    if args.candidate_id:
        winners = winners[
            winners["candidate"].astype(str).eq(args.candidate_id)
        ].reset_index(drop=True)
    if args.winner_stage:
        winners = winners[
            winners["stage"].astype(str).eq(args.winner_stage)
        ].reset_index(drop=True)
    summary, scores = evaluate(
        development_sites=pd.read_csv(args.development_sites, sep="\t"),
        winners=winners,
        development_baselines=pd.read_csv(
            args.development_baselines, sep="\t"
        ),
        profile_cache=args.profile_cache,
        naked_sites={
            cell: pd.read_csv(path, sep="\t")
            for cell, path in naked_site_paths.items()
        },
        naked_corrected=naked_corrected_paths,
        naked_legacy=args.naked_legacy,
        tf=args.tf,
        split=args.split,
        alpha=args.alpha,
        flank=args.flank,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / f"{args.tf}_naked_dna_geometry_summary.tsv"
    scores_path = args.outdir / f"{args.tf}_naked_dna_geometry_scores.tsv.gz"
    summary.to_csv(summary_path, sep="\t", index=False)
    scores.to_csv(scores_path, sep="\t", index=False)
    passes = bool(
        (summary["candidate_false_positive_rate"] <= args.alpha).all()
        and (
            summary["false_positive_rate_increase"]
            <= args.maximum_rate_increase
        ).all()
    )
    manifest = {
        "schema": "fp-tools-tf-geometry-naked-dna-v1",
        "tf": args.tf,
        "cells": sorted(summary["cell"].astype(str).tolist()),
        "threshold_source": f"ChIP-negative {args.split} sites",
        "target_false_positive_rate": float(args.alpha),
        "maximum_rate_increase": float(args.maximum_rate_increase),
        "passes_point_estimate_gate": passes,
        "inputs": {
            "development_sites": {
                "path": str(args.development_sites),
                "sha256": file_sha256(args.development_sites),
            },
            "winners": {
                "path": str(args.winners),
                "sha256": file_sha256(args.winners),
            },
            "development_baselines": {
                "path": str(args.development_baselines),
                "sha256": file_sha256(args.development_baselines),
            },
            "naked_legacy": {
                "path": str(args.naked_legacy),
                "sha256": file_sha256(args.naked_legacy),
            },
            "naked_sites": {
                cell: {"path": str(path), "sha256": file_sha256(path)}
                for cell, path in sorted(naked_site_paths.items())
            },
            "naked_corrected": {
                method: {"path": str(path), "sha256": file_sha256(path)}
                for method, path in sorted(naked_corrected_paths.items())
            },
        },
        "outputs": {
            summary_path.name: file_sha256(summary_path),
            scores_path.name: file_sha256(scores_path),
        },
    }
    manifest_path = args.outdir / f"{args.tf}_naked_dna_geometry_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0 if passes else 2


if __name__ == "__main__":
    raise SystemExit(main())
