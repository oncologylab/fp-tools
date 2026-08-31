#!/usr/bin/env python3
"""Evaluate shape-only TF/family functional templates on development chromosomes.

This supervised diagnostic never reads locked test or holdout labels. It asks
whether a footprint shape learned on train chromosomes, another cell, another
TF in the same motif family, or other TFs can discriminate ChIP occupancy on
validation chromosomes without motif-score or accessibility covariates.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import (  # noqa: E402
    RESIDUAL_MODES,
    _evaluation_profiles,
    binary_metrics,
    chromosome_split,
    residual_score,
    stable_seed,
    validate_sites,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    FunctionalTemplateDetector,
    normalize_functional_profiles,
    profile_descriptors,
    site_accessibility_background,
    standardized_functional_separation,
)
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


TRAINING_SCOPES = (
    "same_cell_ceiling",
    "cross_cell_tf",
    "family_leave_tf_out",
    "global_leave_tf_out",
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selection_score(metrics: dict[str, float | int]) -> float:
    prevalence = float(metrics["prevalence"])
    adjusted_ap = (float(metrics["auprc"]) - prevalence) / max(1.0 - prevalence, 1e-8)
    return float(metrics["auroc"]) + adjusted_ap


def balanced_training_indexes(
    sites: pd.DataFrame,
    indexes: np.ndarray,
    *,
    maximum_per_tf_class: int,
    seed: int,
) -> np.ndarray:
    """Deterministically cap each TF/class so pooled templates are not dominated."""

    selected = sites.iloc[np.asarray(indexes, dtype=int)].copy()
    selected["_original_index"] = np.asarray(indexes, dtype=int)
    parts = []
    for (tf, label), group in selected.groupby(["tf", "chip_label"], sort=True):
        if len(group) > maximum_per_tf_class:
            group = group.sample(
                maximum_per_tf_class,
                random_state=stable_seed(tf, label, seed=seed),
                replace=False,
            )
        parts.append(group["_original_index"].to_numpy(dtype=int))
    return np.sort(np.concatenate(parts)) if parts else np.array([], dtype=int)


def training_indexes(
    sites: pd.DataFrame,
    *,
    target_cell: str,
    target_tf: str,
    target_family: str,
    scope: str,
) -> np.ndarray:
    train = sites["chromosome_split"].astype(str).to_numpy() == "train"
    cells = sites["cell"].astype(str).to_numpy()
    tfs = sites["tf"].astype(str).to_numpy()
    families = sites["motif_family"].astype(str).to_numpy()
    if scope == "same_cell_ceiling":
        mask = train & (cells == target_cell) & (tfs == target_tf)
    elif scope == "cross_cell_tf":
        mask = train & (cells != target_cell) & (tfs == target_tf)
    elif scope == "family_leave_tf_out":
        mask = train & (families == target_family) & (tfs != target_tf)
    elif scope == "global_leave_tf_out":
        mask = train & (tfs != target_tf)
    else:
        raise ValueError(f"unknown training scope: {scope}")
    return np.flatnonzero(mask)


def residual_profiles(
    observed: np.ndarray,
    expected: np.ndarray,
    positions: np.ndarray,
    *,
    residual_mode: str,
    background: str,
    dispersion: float,
) -> np.ndarray:
    adjusted = site_accessibility_background(
        observed,
        expected,
        positions,
        method=background,
    )
    residual, _score = residual_score(observed, adjusted, residual_mode, dispersion)
    return residual


def _fit_and_score(
    *,
    sites: pd.DataFrame,
    profiles_by_cell: dict[str, np.ndarray],
    target_cell: str,
    target_tf: str,
    target_family: str,
    correction: str,
    residual_mode: str,
    background: str,
    smoother: str,
    window_limit: float,
    scope: str,
    positions: np.ndarray,
    minimum_sites_per_class: int,
    maximum_train_per_tf_class: int,
    seed: int,
) -> tuple[dict[str, Any], pd.DataFrame | None]:
    started = perf_counter()
    validation = np.flatnonzero(
        (sites["cell"].astype(str).to_numpy() == target_cell)
        & (sites["tf"].astype(str).to_numpy() == target_tf)
        & (sites["chromosome_split"].astype(str).to_numpy() == "validation")
    )
    train = training_indexes(
        sites,
        target_cell=target_cell,
        target_tf=target_tf,
        target_family=target_family,
        scope=scope,
    )
    train = balanced_training_indexes(
        sites,
        train,
        maximum_per_tf_class=maximum_train_per_tf_class,
        seed=stable_seed(
            target_cell,
            target_tf,
            correction,
            residual_mode,
            background,
            scope,
            seed=seed,
        ),
    )
    train_labels = sites.iloc[train]["chip_label"].to_numpy(dtype=int)
    validation_labels = sites.iloc[validation]["chip_label"].to_numpy(dtype=int)
    base = {
        "cell": target_cell,
        "tf": target_tf,
        "motif_family": target_family,
        "correction": correction,
        "residual": residual_mode,
        "background": background,
        "smoother": smoother,
        "window_limit": float(window_limit),
        "training_scope": scope,
        "training_labels_used": True,
        "motif_or_accessibility_features_used": False,
        "train_sites": int(len(train)),
        "validation_sites": int(len(validation)),
        "train_positive_sites": int(np.sum(train_labels == 1)),
        "train_negative_sites": int(np.sum(train_labels == 0)),
        "validation_positive_sites": int(np.sum(validation_labels == 1)),
        "validation_negative_sites": int(np.sum(validation_labels == 0)),
    }
    if (
        min(np.sum(train_labels == 0), np.sum(train_labels == 1)) < minimum_sites_per_class
        or min(np.sum(validation_labels == 0), np.sum(validation_labels == 1))
        < minimum_sites_per_class
    ):
        return {**base, "status": "insufficient_sites"}, None

    # Rows from different cells live in separate profile arrays. Construct the
    # pooled train matrix in global-table order without copying unrelated rows.
    train_parts = []
    train_label_parts = []
    train_site_parts = []
    for cell in sorted(sites.iloc[train]["cell"].astype(str).unique()):
        local_global = train[sites.iloc[train]["cell"].astype(str).to_numpy() == cell]
        cell_rows = np.flatnonzero(sites["cell"].astype(str).to_numpy() == cell)
        local_lookup = pd.Series(np.arange(len(cell_rows)), index=cell_rows)
        local = local_lookup.loc[local_global].to_numpy(dtype=int)
        train_parts.append(profiles_by_cell[cell][local])
        train_label_parts.append(sites.iloc[local_global]["chip_label"].to_numpy(dtype=int))
        train_site_parts.append(sites.iloc[local_global])
    train_profiles = np.vstack(train_parts)
    train_labels = np.concatenate(train_label_parts)
    train_sites_frame = pd.concat(train_site_parts, ignore_index=True)
    target_rows = np.flatnonzero(sites["cell"].astype(str).to_numpy() == target_cell)
    target_lookup = pd.Series(np.arange(len(target_rows)), index=target_rows)
    local_validation = target_lookup.loc[validation].to_numpy(dtype=int)
    validation_profiles = profiles_by_cell[target_cell][local_validation]
    model = FunctionalTemplateDetector(
        positions,
        smoother=smoother,
        window_limit=window_limit,
    ).fit(train_profiles, train_labels)
    scores = model.predict_proba(validation_profiles)
    metrics = binary_metrics(validation_labels, scores)
    normalized_validation = normalize_functional_profiles(validation_profiles, positions)
    positive_mean = np.mean(normalized_validation[validation_labels == 1], axis=0)
    negative_mean = np.mean(normalized_validation[validation_labels == 0], axis=0)
    descriptors = asdict(profile_descriptors(model.footprint_template_, positions))
    row = {
        **base,
        "status": "ok",
        "training_cells": ",".join(sorted(train_sites_frame["cell"].astype(str).unique())),
        "training_tfs": ",".join(sorted(train_sites_frame["tf"].astype(str).unique())),
        "effective_train_sites": model.effective_sites_,
        "fit_seconds": perf_counter() - started,
        "selection_score": selection_score(metrics),
        "functional_separation": standardized_functional_separation(
            validation_profiles,
            validation_labels,
            positions,
        ),
        **descriptors,
        **metrics,
    }
    curve = pd.DataFrame(
        {
            "cell": target_cell,
            "tf": target_tf,
            "motif_family": target_family,
            "correction": correction,
            "residual": residual_mode,
            "background": background,
            "smoother": smoother,
            "window_limit": float(window_limit),
            "training_scope": scope,
            "position": positions.astype(int),
            "template": model.footprint_template_,
            "template_standard_error": model.template_standard_error_,
            "positive_mean": positive_mean,
            "negative_mean": negative_mean,
            "positive_minus_negative": positive_mean - negative_mean,
        }
    )
    return row, curve


def evaluate(
    study: dict,
    sites: pd.DataFrame,
    tracks: pd.DataFrame,
    *,
    corrections: tuple[str, ...],
    residual_modes: tuple[str, ...],
    backgrounds: tuple[str, ...],
    smoothers: tuple[str, ...],
    windows: tuple[float, ...],
    cache_dir: Path,
    genome: Path | None,
    flank: int,
    minimum_sites_per_class: int,
    maximum_train_per_tf_class: int,
    seed: int,
    workers: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    positions = np.arange(-flank, flank + 1, dtype=float)
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == "development"]
    cell_sites = {
        cell: sites[sites["cell"].astype(str) == cell].reset_index(drop=True)
        for cell in sorted(tasks["cell"].unique())
    }
    profile_sets: dict[tuple[str, str, str, str], np.ndarray] = {}
    dispersions: dict[tuple[str, str], float] = {}
    for correction in corrections:
        for cell, selected_sites in cell_sites.items():
            observed, expected = _evaluation_profiles(
                selected_sites,
                tracks,
                cell,
                correction,
                cache_dir,
                "development",
                flank,
                genome,
            )
            dispersion = estimate_nb_dispersion(observed, expected)
            dispersions[(cell, correction)] = dispersion
            for background in backgrounds:
                for residual_mode in residual_modes:
                    profile_sets[(cell, correction, background, residual_mode)] = residual_profiles(
                        observed,
                        expected,
                        positions,
                        residual_mode=residual_mode,
                        background=background,
                        dispersion=dispersion,
                    )

    rows = []
    curves = []
    futures = {}
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=workers) as executor:
        for correction in corrections:
            for background in backgrounds:
                for residual_mode in residual_modes:
                    profiles_by_cell = {
                        cell: profile_sets[(cell, correction, background, residual_mode)]
                        for cell in cell_sites
                    }
                    for task in tasks.itertuples(index=False):
                        for scope in TRAINING_SCOPES:
                            for smoother in smoothers:
                                for window in windows:
                                    future = executor.submit(
                                        _fit_and_score,
                                        sites=sites,
                                        profiles_by_cell=profiles_by_cell,
                                        target_cell=str(task.cell),
                                        target_tf=str(task.tf),
                                        target_family=str(task.motif_family),
                                        correction=correction,
                                        residual_mode=residual_mode,
                                        background=background,
                                        smoother=smoother,
                                        window_limit=window,
                                        scope=scope,
                                        positions=positions,
                                        minimum_sites_per_class=minimum_sites_per_class,
                                        maximum_train_per_tf_class=maximum_train_per_tf_class,
                                        seed=seed,
                                    )
                                    futures[future] = (
                                        str(task.cell),
                                        str(task.tf),
                                        correction,
                                        residual_mode,
                                        background,
                                        smoother,
                                        window,
                                        scope,
                                    )
        for future in as_completed(futures):
            identity = futures[future]
            try:
                row, curve = future.result()
            except Exception as error:
                cell, tf, correction, residual_mode, background, smoother, window, scope = identity
                row = {
                    "cell": cell,
                    "tf": tf,
                    "correction": correction,
                    "residual": residual_mode,
                    "background": background,
                    "smoother": smoother,
                    "window_limit": window,
                    "training_scope": scope,
                    "status": "error",
                    "error": f"{type(error).__name__}: {error}",
                }
                curve = None
            rows.append(row)
            if curve is not None:
                curves.append(curve)
    metrics = pd.DataFrame(rows)
    profile_curves = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    return metrics, profile_curves


def select_winners(metrics: pd.DataFrame) -> pd.DataFrame:
    passing = metrics[metrics["status"] == "ok"].copy()
    if passing.empty:
        return passing
    keys = ["cell", "tf", "correction", "training_scope"]
    return (
        passing.sort_values(
            keys + ["selection_score", "auprc", "auroc"],
            ascending=[True] * len(keys) + [False, False, False],
            kind="mergesort",
        )
        .groupby(keys, sort=True, as_index=False)
        .head(1)
        .reset_index(drop=True)
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-run", type=Path, required=True)
    parser.add_argument("--tracks", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--corrections", nargs="+")
    parser.add_argument("--residuals", nargs="+", choices=RESIDUAL_MODES, default=list(RESIDUAL_MODES))
    parser.add_argument(
        "--backgrounds",
        nargs="+",
        choices=("none", "linear", "quadratic", "gp-long"),
        default=["none", "linear", "gp-long"],
    )
    parser.add_argument("--smoothers", nargs="+", choices=("spline", "gp"), default=["spline", "gp"])
    parser.add_argument("--windows", type=float, nargs="+", default=[30.0, 50.0, 80.0])
    parser.add_argument("--minimum-sites-per-class", type=int, default=100)
    parser.add_argument("--maximum-train-per-tf-class", type=int, default=2000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args(argv)
    if args.workers < 1 or args.minimum_sites_per_class < 2:
        raise SystemExit("workers and minimum sites must be positive")

    base_manifest_path = args.base_run / "functional_benchmark_manifest.json"
    base_manifest = json.loads(base_manifest_path.read_text(encoding="utf-8"))
    if base_manifest.get("test_unlocked"):
        raise SystemExit("template transfer refuses a base run that opened test labels")
    study_path = Path(base_manifest["study"])
    sites_path = Path(base_manifest["development_sites"])
    tracks_path = args.tracks or Path(base_manifest["tracks"])
    genome = Path(base_manifest["genome"]) if base_manifest.get("genome") else None
    study = json.loads(study_path.read_text(encoding="utf-8"))
    sites = validate_sites(pd.read_csv(sites_path, sep="\t"), sites_path)
    sites["chromosome_split"] = sites["TFBS_chr"].map(
        lambda chromosome: chromosome_split(str(chromosome), study)
    )
    sites = sites[sites["chromosome_split"].isin(["train", "validation"])].reset_index(drop=True)
    tracks = pd.read_csv(tracks_path, sep="\t")
    corrections = tuple(args.corrections or base_manifest["corrections"])
    unknown = set(corrections).difference(tracks["model"].astype(str))
    if unknown:
        raise SystemExit("unknown corrections: " + ", ".join(sorted(unknown)))
    seed = int(args.seed if args.seed is not None else base_manifest["seed"])
    flank = int(base_manifest["flank"])
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics, curves = evaluate(
        study,
        sites,
        tracks,
        corrections=corrections,
        residual_modes=tuple(args.residuals),
        backgrounds=tuple(args.backgrounds),
        smoothers=tuple(args.smoothers),
        windows=tuple(args.windows),
        cache_dir=args.base_run / "profile_cache",
        genome=genome,
        flank=flank,
        minimum_sites_per_class=args.minimum_sites_per_class,
        maximum_train_per_tf_class=args.maximum_train_per_tf_class,
        seed=seed,
        workers=args.workers,
    )
    winners = select_winners(metrics)
    metrics_path = args.outdir / "functional_template_transfer_metrics.tsv.gz"
    curves_path = args.outdir / "functional_template_transfer_profiles.tsv.gz"
    winners_path = args.outdir / "functional_template_transfer_winners.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    winners.to_csv(winners_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-functional-template-transfer-v1",
        "locked_test_labels_read": False,
        "selection_split": "validation",
        "training_labels_used": True,
        "motif_or_accessibility_features_used": False,
        "base_manifest": str(base_manifest_path),
        "base_manifest_sha256": file_sha256(base_manifest_path),
        "study": str(study_path),
        "study_sha256": file_sha256(study_path),
        "development_sites": str(sites_path),
        "development_sites_sha256": file_sha256(sites_path),
        "tracks": str(tracks_path),
        "tracks_sha256": file_sha256(tracks_path),
        "corrections": list(corrections),
        "residuals": args.residuals,
        "backgrounds": args.backgrounds,
        "smoothers": args.smoothers,
        "windows": args.windows,
        "training_scopes": list(TRAINING_SCOPES),
        "minimum_sites_per_class": args.minimum_sites_per_class,
        "maximum_train_per_tf_class": args.maximum_train_per_tf_class,
        "seed": seed,
        "workers": args.workers,
        "metrics_rows": int(len(metrics)),
        "winner_rows": int(len(winners)),
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {"path": str(curves_path), "sha256": file_sha256(curves_path)},
            "winners": {"path": str(winners_path), "sha256": file_sha256(winners_path)},
        },
    }
    (args.outdir / "functional_template_transfer_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    columns = [
        "cell",
        "tf",
        "training_scope",
        "correction",
        "residual",
        "background",
        "smoother",
        "window_limit",
        "auroc",
        "auprc",
    ]
    print(winners[columns].to_string(index=False) if len(winners) else "no eligible winners")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
