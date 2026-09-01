#!/usr/bin/env python3
"""Test a frozen TF geometry in external cells without retuning it."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from fp_tools.utils.signals import footprint_score_array_fast  # noqa: E402
from render_tf_before_after_report import crossfit_covariate_residuals  # noqa: E402
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


def parse_signal(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("signal must use CELL,REPLICATE,PATH")
    return fields[0], fields[1], Path(fields[2])


def conventional_profile_scores(profiles: np.ndarray, center: int) -> np.ndarray:
    """Return the conventional footprint score at each profile center."""

    values = np.asarray(profiles, dtype=float)
    if values.ndim != 2 or not 0 <= center < values.shape[1]:
        raise ValueError("profiles must be a matrix with a valid center")
    if not np.isfinite(values).all():
        raise ValueError("conventional scoring requires finite profiles")
    return np.asarray(
        [
            footprint_score_array_fast(row, 10, 30, 20, 50)[center]
            for row in values
        ],
        dtype=float,
    )


def paired_metrics(
    labels: np.ndarray, conventional: np.ndarray, candidate: np.ndarray
) -> dict[str, float]:
    labels = np.asarray(labels, dtype=int)
    conventional = np.asarray(conventional, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    conventional_auroc = float(roc_auc_score(labels, conventional))
    candidate_auroc = float(roc_auc_score(labels, candidate))
    conventional_auprc = float(average_precision_score(labels, conventional))
    candidate_auprc = float(average_precision_score(labels, candidate))
    return {
        "conventional_auroc": conventional_auroc,
        "candidate_auroc": candidate_auroc,
        "delta_auroc": candidate_auroc - conventional_auroc,
        "conventional_auprc": conventional_auprc,
        "candidate_auprc": candidate_auprc,
        "delta_auprc": candidate_auprc - conventional_auprc,
    }


def chromosome_block_bootstrap(
    labels: np.ndarray,
    conventional: np.ndarray,
    candidate: np.ndarray,
    chromosomes: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> dict[str, float | int]:
    """Bootstrap paired metric gains by resampling chromosome blocks."""

    labels = np.asarray(labels, dtype=int)
    conventional = np.asarray(conventional, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    chromosomes = np.asarray(chromosomes)
    blocks = np.unique(chromosomes)
    if len(blocks) < 2:
        raise ValueError("chromosome bootstrap requires at least two blocks")
    indexes = {block: np.flatnonzero(chromosomes == block) for block in blocks}
    rng = np.random.default_rng(seed)
    auroc: list[float] = []
    auprc: list[float] = []
    for _ in range(iterations):
        sampled = rng.choice(blocks, size=len(blocks), replace=True)
        selected = np.concatenate([indexes[block] for block in sampled])
        if len(np.unique(labels[selected])) != 2:
            continue
        metrics = paired_metrics(
            labels[selected], conventional[selected], candidate[selected]
        )
        auroc.append(metrics["delta_auroc"])
        auprc.append(metrics["delta_auprc"])
    if not auroc:
        raise ValueError("no valid chromosome bootstrap samples")
    return {
        "bootstrap_successful": len(auroc),
        "delta_auroc_ci_low": float(np.quantile(auroc, 0.025)),
        "delta_auroc_ci_high": float(np.quantile(auroc, 0.975)),
        "delta_auprc_ci_low": float(np.quantile(auprc, 0.025)),
        "delta_auprc_ci_high": float(np.quantile(auprc, 0.975)),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--source-cell", required=True)
    parser.add_argument(
        "--signal",
        action="append",
        type=parse_signal,
        required=True,
        metavar="CELL,REPLICATE,BIGWIG",
    )
    parser.add_argument("--tf", default="CTCF")
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    sites_all = pd.read_csv(args.evaluation_sites, sep="\t")
    winners = pd.read_csv(args.winners, sep="\t")
    winner_rows = winners[
        winners["cell"].astype(str).eq(args.source_cell)
        & winners["tf"].astype(str).eq(args.tf)
    ]
    if len(winner_rows) != 1:
        parser.error(
            f"expected one winner for {args.source_cell}/{args.tf}, "
            f"found {len(winner_rows)}"
        )
    candidate = candidate_from_row(winner_rows.iloc[0])
    if candidate.correction != "DWM":
        parser.error("external transfer currently requires a frozen DWM geometry")

    signals: dict[str, list[tuple[str, Path]]] = {}
    for cell, replicate, path in args.signal:
        signals.setdefault(cell, []).append((replicate, path))

    summary_rows: list[dict[str, object]] = []
    replicate_rows: list[dict[str, object]] = []
    site_frames: list[pd.DataFrame] = []
    profile_outputs: list[Path] = []
    args.outdir.mkdir(parents=True, exist_ok=True)
    covariate_columns = [
        "motif_score",
        "log_accessibility",
        "gc_fraction",
        "peak_position_abs",
    ]
    for cell in sorted(signals):
        sites = sites_all[
            sites_all["cell"].astype(str).eq(cell)
            & sites_all["tf"].astype(str).eq(args.tf)
        ].reset_index(drop=True)
        if sites.empty:
            raise ValueError(f"no evaluation sites for {cell}/{args.tf}")
        profiles: list[np.ndarray] = []
        validity: list[np.ndarray] = []
        replicate_names: list[str] = []
        for replicate, signal in sorted(signals[cell]):
            profile, valid = extract_profiles(sites, signal, args.flank)
            profiles.append(profile)
            validity.append(valid)
            replicate_names.append(replicate)
        stacked = np.stack(profiles)
        valid = np.logical_and.reduce(validity)
        valid &= np.isfinite(stacked).all(axis=(0, 2))
        covariates = sites[covariate_columns].to_numpy(dtype=float)
        valid &= np.isfinite(covariates).all(axis=1)
        pooled = np.mean(stacked, axis=0)
        labels = sites["label"].to_numpy(dtype=int)
        candidate_scores = score_candidate(pooled, candidate)
        conventional_scores = conventional_profile_scores(pooled, args.flank)
        valid &= np.isfinite(candidate_scores) & np.isfinite(conventional_scores)
        if len(np.unique(labels[valid])) != 2:
            raise ValueError(f"{cell}/{args.tf} does not contain two valid classes")

        report_sites = sites.loc[valid].reset_index(drop=True)
        report_labels = labels[valid]
        report_candidate = candidate_scores[valid]
        report_conventional = conventional_scores[valid]
        report_covariates = covariates[valid]
        chromosomes = report_sites["TFBS_chr"].astype(str).to_numpy()
        raw_metrics = paired_metrics(
            report_labels, report_conventional, report_candidate
        )
        residual_conventional = crossfit_covariate_residuals(
            report_conventional, report_covariates, chromosomes
        )
        residual_candidate = crossfit_covariate_residuals(
            report_candidate, report_covariates, chromosomes
        )
        residual_metrics = {
            f"residual_{key}": value
            for key, value in paired_metrics(
                report_labels, residual_conventional, residual_candidate
            ).items()
        }
        bootstrap = chromosome_block_bootstrap(
            report_labels,
            report_conventional,
            report_candidate,
            chromosomes,
            iterations=args.bootstrap,
            seed=args.seed,
        )

        for replicate, profile in zip(replicate_names, profiles):
            replicate_candidate = score_candidate(profile, candidate)[valid]
            replicate_conventional = conventional_profile_scores(
                profile[valid], args.flank
            )
            replicate_rows.append(
                {
                    "cell": cell,
                    "tf": args.tf,
                    "replicate": replicate,
                    **paired_metrics(
                        report_labels,
                        replicate_conventional,
                        replicate_candidate,
                    ),
                }
            )
        cell_replicates = [
            row for row in replicate_rows if row["cell"] == cell
        ]
        summary_rows.append(
            {
                "cell": cell,
                "tf": args.tf,
                "source_cell": args.source_cell,
                "candidate_id": candidate.identifier,
                "replicates": len(profiles),
                "sites": int(valid.sum()),
                "positive_sites": int(report_labels.sum()),
                **raw_metrics,
                **residual_metrics,
                **bootstrap,
                "replicate_min_delta_auroc": min(
                    row["delta_auroc"] for row in cell_replicates
                ),
                "replicate_min_delta_auprc": min(
                    row["delta_auprc"] for row in cell_replicates
                ),
            }
        )
        site_frame = report_sites.copy()
        site_frame["conventional_score"] = report_conventional
        site_frame["candidate_score"] = report_candidate
        site_frame["residual_conventional_score"] = residual_conventional
        site_frame["residual_candidate_score"] = residual_candidate
        site_frames.append(site_frame)
        profile_path = args.outdir / f"{cell}.{args.tf}.pooled_profiles.npz"
        np.savez_compressed(
            profile_path,
            profiles=pooled[valid].astype(np.float32),
            labels=report_labels.astype(np.int8),
        )
        profile_outputs.append(profile_path)

    summary = pd.DataFrame(summary_rows)
    replicate_metrics = pd.DataFrame(replicate_rows)
    site_scores = pd.concat(site_frames, ignore_index=True)
    summary_path = args.outdir / f"{args.tf}_external_transfer_summary.tsv"
    replicate_path = args.outdir / f"{args.tf}_external_transfer_replicates.tsv"
    scores_path = args.outdir / f"{args.tf}_external_transfer_scores.tsv.gz"
    summary.to_csv(summary_path, sep="\t", index=False)
    replicate_metrics.to_csv(replicate_path, sep="\t", index=False)
    site_scores.to_csv(scores_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-tf-geometry-external-transfer-v1",
        "tf": args.tf,
        "source_cell": args.source_cell,
        "candidate_id": candidate.identifier,
        "labels_used_for_selection": False,
        "analysis_status": "posthoc_no_retuning",
        "inputs": {
            "evaluation_sites": {
                "path": str(args.evaluation_sites),
                "sha256": file_sha256(args.evaluation_sites),
            },
            "winners": {
                "path": str(args.winners),
                "sha256": file_sha256(args.winners),
            },
            "signals": [
                {
                    "cell": cell,
                    "replicate": replicate,
                    "path": str(path),
                    "sha256": file_sha256(path),
                }
                for cell, replicate, path in args.signal
            ],
        },
        "outputs": {
            path.name: file_sha256(path)
            for path in [
                summary_path,
                replicate_path,
                scores_path,
                *profile_outputs,
            ]
        },
    }
    manifest_path = args.outdir / f"{args.tf}_external_transfer_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(replicate_metrics.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
