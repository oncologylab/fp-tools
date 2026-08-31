#!/usr/bin/env python3
"""Score frozen detectors on motif-misaligned, label-free ATAC profiles.

Naked DNA controls enzyme and sequence artifacts, but it does not reproduce
the broad accessibility and nucleosomal structure present in cellular ATAC.
This helper creates a second null by cyclically shifting each motif-oriented
profile away from its true motif center.  The shift preserves each site's
coverage and cut distribution while destroying motif-centered alignment.

Only label-free development artifacts are accepted.  The output can be
combined with naked-DNA scores to calibrate a detector against both nulls.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Sequence

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import stable_seed  # noqa: E402
from evaluate_naked_dna_functional_policy import (  # noqa: E402
    _candidate_lookup,
    file_sha256,
    fit_dwm_detector,
    fit_strand_detector,
    load_dwm_training_source,
    load_strand_artifact,
    parse_dwm_training_artifact,
    parse_training_artifact,
    predict_detector,
)


DEFAULT_SHIFTS = (-45, -35, 35, 45)


def cyclic_shift_profiles(
    profiles: dict[str, np.ndarray], shift: int
) -> dict[str, np.ndarray]:
    """Shift every position-indexed profile without changing its row totals."""

    if shift == 0:
        raise ValueError("zero is not a valid motif-misalignment shift")
    shifted: dict[str, np.ndarray] = {}
    for name, values in profiles.items():
        array = np.asarray(values)
        if array.ndim != 2:
            raise ValueError(f"profile array {name} must be two-dimensional")
        shifted[name] = np.roll(array, int(shift), axis=1)
    return shifted


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def evaluate(
    *,
    study: dict[str, Any],
    policy: pd.DataFrame,
    strand_training_paths: dict[tuple[str, str], Path],
    dwm_training_paths: dict[str, Path],
    dwm_base_run: Path | None,
    shifts: Sequence[int],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, str]]]:
    shifts = tuple(int(value) for value in shifts)
    if not shifts or 0 in shifts or len(set(shifts)) != len(shifts):
        raise ValueError("shifts must be unique, nonzero integers")
    promoted = policy[policy["passes_development_gates"].astype(bool)].copy()
    routes = promoted.set_index("motif_family")
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[
        tasks["split"].eq("development")
        & tasks["motif_family"].astype(str).isin(routes.index.astype(str))
    ].copy()
    flank = int(study["profile_flank_bp"])
    positions = np.arange(-flank, flank + 1, dtype=float)
    strand_candidates, dwm_candidates = _candidate_lookup()
    training_cache: dict[tuple[str, str], tuple] = {}
    dwm_cache: dict[str, tuple] = {}
    score_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    inputs: list[dict[str, str]] = []

    for task in tasks.sort_values(["cell", "tf"]).itertuples(index=False):
        cell = str(task.cell)
        tf = str(task.tf)
        motif_family = str(task.motif_family)
        route = routes.loc[motif_family]
        bias = str(route["candidate_bias_configuration"])
        candidate_id = str(route["candidate_id"])
        reference_id = str(route["reference_candidate_id"])
        if candidate_id not in strand_candidates:
            raise ValueError(f"unknown frozen strand candidate: {candidate_id}")
        if reference_id not in dwm_candidates:
            raise ValueError(f"unknown frozen DWM candidate: {reference_id}")

        key = (bias, cell)
        if key not in training_cache:
            path = strand_training_paths[key]
            training_cache[key] = load_strand_artifact(path, cell)
            inputs.append(
                {"purpose": "strand_training", "path": str(path), "sha256": file_sha256(path)}
            )
        sites, strand_profiles, strand_valid, strand_hashes, _document = training_cache[key]
        if not np.all(strand_valid):
            sites = sites.loc[strand_valid].reset_index(drop=True)
            strand_profiles = {
                name: values[strand_valid] for name, values in strand_profiles.items()
            }
            strand_hashes = strand_hashes[strand_valid]
        candidate = fit_strand_detector(
            strand_candidates[candidate_id],
            sites,
            strand_profiles,
            tf=tf,
            motif_family=motif_family,
            positions=positions,
            seed=stable_seed(cell, tf, candidate_id, "shift-null", seed=seed),
        )

        if cell not in dwm_cache:
            dwm_cache[cell] = load_dwm_training_source(
                cell,
                artifact_paths=dwm_training_paths,
                base_run=dwm_base_run,
                flank=flank,
            )
            for record in dwm_cache[cell][4]:
                inputs.append({"purpose": "dwm_training", **record})
        dwm_sites, dwm_profiles, dwm_valid, dwm_hashes, _records = dwm_cache[cell]
        dwm_sites = dwm_sites.loc[dwm_valid].reset_index(drop=True)
        dwm_profiles = {name: values[dwm_valid] for name, values in dwm_profiles.items()}
        dwm_hashes = dwm_hashes[dwm_valid]
        if not np.array_equal(strand_hashes, dwm_hashes):
            raise ValueError(f"strand and DWM training sites differ for {cell}")
        reference = fit_dwm_detector(
            dwm_candidates[reference_id],
            dwm_sites,
            dwm_profiles,
            cell=cell,
            tf=tf,
            motif_family=motif_family,
            positions=positions,
            seed=stable_seed(cell, tf, reference_id, "shift-null", seed=seed),
        )

        tf_mask = sites["tf"].astype(str).eq(tf).to_numpy()
        for shift in shifts:
            for method, fitted, source_profiles, method_id, method_bias in (
                (
                    "frozen_policy_candidate",
                    candidate,
                    strand_profiles,
                    candidate_id,
                    bias,
                ),
                (
                    "frozen_dwm_reference",
                    reference,
                    dwm_profiles,
                    reference_id,
                    "DWM",
                ),
            ):
                shifted = cyclic_shift_profiles(source_profiles, shift)
                probabilities, total_signal, _residual = predict_detector(
                    fitted, shifted, sites
                )
                valid = tf_mask & np.isfinite(probabilities)
                informative = valid & (total_signal > 0)
                selected = np.flatnonzero(tf_mask)
                metadata = {
                    "cell": cell,
                    "tf": tf,
                    "motif_family": motif_family,
                    "method": method,
                    "candidate_id": method_id,
                    "bias_configuration": method_bias,
                    "null_source": "motif_misaligned_atac",
                    "shift_bp": int(shift),
                }
                score_parts.append(
                    pd.DataFrame(
                        {
                            **metadata,
                            "site_hash": strand_hashes[selected],
                            "TFBS_chr": sites.loc[tf_mask, "TFBS_chr"].to_numpy(),
                            "TFBS_start": sites.loc[tf_mask, "TFBS_start"].to_numpy(),
                            "TFBS_end": sites.loc[tf_mask, "TFBS_end"].to_numpy(),
                            "total_signal": total_signal[selected],
                            "binding_probability": probabilities[selected],
                            "valid": valid[selected],
                            "informative": informative[selected],
                        }
                    )
                )
                summary_rows.append(
                    {
                        **metadata,
                        "sites": int(tf_mask.sum()),
                        "valid_sites": int(valid.sum()),
                        "informative_sites": int(informative.sum()),
                        "mean_probability": (
                            float(np.mean(probabilities[informative]))
                            if informative.any()
                            else np.nan
                        ),
                        "q95_probability": (
                            float(np.quantile(probabilities[informative], 0.95))
                            if informative.any()
                            else np.nan
                        ),
                        "q975_probability": (
                            float(np.quantile(probabilities[informative], 0.975))
                            if informative.any()
                            else np.nan
                        ),
                    }
                )
    scores = pd.concat(score_parts, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    unique_inputs = list({record["path"]: record for record in inputs}.values())
    return scores, summary, unique_inputs


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"),
    )
    parser.add_argument(
        "--policy",
        type=Path,
        default=Path("benchmarks/manifests/compact/functional_detector_policy_v1.tsv"),
    )
    parser.add_argument("--dwm-base-run", type=Path)
    parser.add_argument(
        "--dwm-training-artifact",
        action="append",
        type=parse_dwm_training_artifact,
        default=[],
        metavar="CELL,JSON",
    )
    parser.add_argument(
        "--training-artifact",
        action="append",
        type=parse_training_artifact,
        default=[],
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument("--shifts", nargs="+", type=int, default=list(DEFAULT_SHIFTS))
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    policy = pd.read_csv(args.policy, sep="\t")
    strand_paths = {(model, cell): path for model, cell, path in args.training_artifact}
    dwm_paths = {cell: path for cell, path in args.dwm_training_artifact}
    scores, summary, inputs = evaluate(
        study=study,
        policy=policy,
        strand_training_paths=strand_paths,
        dwm_training_paths=dwm_paths,
        dwm_base_run=args.dwm_base_run,
        shifts=args.shifts,
        seed=args.seed,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    scores_path = args.outdir / "shifted_atac_null_scores.tsv.gz"
    summary_path = args.outdir / "shifted_atac_null_summary.tsv"
    scores.to_csv(scores_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-shifted-atac-null-policy-v1",
        "labels_used": False,
        "null": "cyclic motif-profile misalignment",
        "shifts_bp": [int(value) for value in args.shifts],
        "seed": int(args.seed),
        "study": {"path": str(args.study), "sha256": _sha256(args.study)},
        "policy": {"path": str(args.policy), "sha256": _sha256(args.policy)},
        "inputs": inputs,
        "outputs": {
            "scores": {"path": str(scores_path), "sha256": _sha256(scores_path)},
            "summary": {"path": str(summary_path), "sha256": _sha256(summary_path)},
        },
    }
    (args.outdir / "shifted_atac_null_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
