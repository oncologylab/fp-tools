#!/usr/bin/env python3
"""Evaluate pinned ChromBPNet predictions on frozen fp-tools motif sites.

The benchmark places every signal on identical finite motif sites and uses the
same conventional center-versus-shoulder geometry before comparing more
specialized detectors.  The frozen parametric and ChromBPNet bias predictions
are treated as expected cut profiles and locally scaled to the observed total.
ChromBPNet's regulatory and bias-free predictions are total-matched before
geometry scoring so count-head calibration cannot create an artificial
advantage.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_frozen_functional_policy import (  # noqa: E402
    aggregate_curve,
    load_test_sites,
    preflight_test_artifact,
)
from evaluate_functional_footprints import binary_metrics  # noqa: E402
from evaluate_parametric_factorization import (  # noqa: E402
    block_bootstrap_delta,
    geometry_score,
    residual_score,
)
from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from fp_tools.tools.functional_footprints import (  # noqa: E402
    normalize_functional_profiles,
    orient_profiles,
    profile_descriptors,
    standardized_functional_separation,
)
from fp_tools.utils import bigwig as pyBigWig  # noqa: E402


SCHEMA = "fp-tools-chrombpnet-frozen-comparison-v1"
INPUT_SCHEMA = "fp-tools-chrombpnet-comparison-inputs-v1"
REFERENCE_SCHEMA = "fp-tools-chrombpnet-reference-run-v1"
PINNED_COMMIT = "09938fdb4397ec0006510e5251e48920a505d4de"


def parse_cell_path(value: str) -> tuple[str, Path]:
    fields = value.split("=", 1)
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError("value must use CELL=PATH")
    return fields[0], Path(fields[1])


def unique_cell_paths(values: Iterable[tuple[str, Path]], label: str) -> dict[str, Path]:
    output: dict[str, Path] = {}
    for cell, path in values:
        if cell in output:
            raise ValueError(f"duplicate {label} for {cell}")
        output[cell] = path
    return output


def canonical_id(document: dict) -> str:
    encoded = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return sha256(encoded.encode()).hexdigest()


def immutable_json(path: Path, document: dict) -> None:
    rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if path.exists() and path.read_text(encoding="utf-8") != rendered:
        raise ValueError(f"immutable artifact differs: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding="utf-8")


def validate_bigwig(path: Path) -> dict[str, int]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(path)
    handle = pyBigWig.open(str(path))
    try:
        chromosomes = handle.chroms()
        header = handle.header()
    finally:
        handle.close()
    covered = int(header.get("nBasesCovered", 0))
    if not chromosomes or covered <= 0:
        raise ValueError(f"bigWig has no nonempty coverage: {path}")
    return {"chromosomes": len(chromosomes), "covered_bases": covered}


def validate_reference_manifest(path: Path) -> dict:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != REFERENCE_SCHEMA:
        raise ValueError(f"unsupported ChromBPNet reference manifest: {path}")
    if document.get("source_commit") != PINNED_COMMIT:
        raise ValueError(f"ChromBPNet source is not pinned in {path}")
    if document.get("completed") is not True:
        raise ValueError(f"ChromBPNet stage is incomplete: {path}")
    for record in document.get("outputs", []):
        output = Path(record["path"])
        if not output.is_file() or file_sha256(output) != record["sha256"]:
            raise ValueError(f"ChromBPNet output checksum mismatch: {output}")
    return document


def extract_prediction_profiles(
    sites: pd.DataFrame,
    signal: Path,
    flank: int,
    *,
    require_dense: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract and motif-orient a prediction track.

    Dense mode distinguishes unpredicted ChromBPNet regions from genuine zero
    predictions.  Sparse mode is used for the conventional DWM track, where an
    absent interval denotes zero expected cuts.
    """

    width = 2 * flank + 1
    profiles = np.full((len(sites), width), np.nan, dtype=np.float64)
    valid = np.zeros(len(sites), dtype=bool)
    handle = pyBigWig.open(str(signal))
    try:
        chromosomes = handle.chroms()
        for index, row in enumerate(sites.itertuples(index=False)):
            chromosome = str(row.TFBS_chr)
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            start = center - flank
            end = center + flank + 1
            if (
                chromosome not in chromosomes
                or start < 0
                or end > int(chromosomes[chromosome])
            ):
                continue
            values = np.asarray(
                handle.values(chromosome, start, end, numpy=True), dtype=float
            )
            if len(values) != width:
                continue
            finite = np.isfinite(values)
            if require_dense and not np.all(finite):
                continue
            profiles[index] = np.nan_to_num(values, nan=0.0)
            valid[index] = True
    finally:
        handle.close()
    profiles = orient_profiles(profiles, sites["TFBS_strand"].astype(str))
    return profiles, valid


def scale_to_observed_total(observed: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    """Match prediction shape to each observed profile total."""

    observed = np.asarray(observed, dtype=float)
    prediction = np.maximum(np.asarray(prediction, dtype=float), 0.0)
    if observed.ndim != 2 or prediction.shape != observed.shape:
        raise ValueError("observed and prediction profiles must have equal 2D shape")
    observed_total = observed.sum(axis=1)
    predicted_total = prediction.sum(axis=1)
    return np.divide(
        prediction * observed_total[:, None],
        predicted_total[:, None],
        out=np.zeros_like(prediction),
        where=predicted_total[:, None] > 0,
    )


def score_methods(
    observed: np.ndarray,
    dwm_expected: np.ndarray,
    parametric_expected: np.ndarray,
    bias_prediction: np.ndarray,
    regulatory_prediction: np.ndarray,
    nobias_prediction: np.ndarray,
    positions: np.ndarray,
    *,
    dispersion: float,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return conventional-geometry scores and profiles for every arm."""

    dwm_expected = scale_to_observed_total(observed, dwm_expected)
    parametric_expected = scale_to_observed_total(observed, parametric_expected)
    bias_prediction = scale_to_observed_total(observed, bias_prediction)
    regulatory_prediction = scale_to_observed_total(observed, regulatory_prediction)
    nobias_prediction = scale_to_observed_total(observed, nobias_prediction)
    dwm_score, dwm_residual = residual_score(
        observed, dwm_expected, positions, "deviance", dispersion
    )
    parametric_score, parametric_residual = residual_score(
        observed, parametric_expected, positions, "deviance", dispersion
    )
    bias_score, bias_residual = residual_score(
        observed, bias_prediction, positions, "deviance", dispersion
    )
    return {
        "raw_geometry": (geometry_score(observed, positions), observed),
        "DWM_conventional_geometry": (
            dwm_score,
            normalize_functional_profiles(dwm_residual, positions),
        ),
        "frozen_parametric_bias_conventional_geometry": (
            parametric_score,
            normalize_functional_profiles(parametric_residual, positions),
        ),
        "ChromBPNet_bias_conventional_geometry": (
            bias_score,
            normalize_functional_profiles(bias_residual, positions),
        ),
        "ChromBPNet_regulatory_geometry": (
            geometry_score(regulatory_prediction, positions),
            regulatory_prediction,
        ),
        "ChromBPNet_nobias_geometry": (
            geometry_score(nobias_prediction, positions),
            nobias_prediction,
        ),
    }


def metric_row(
    sites: pd.DataFrame,
    score: np.ndarray,
    profiles: np.ndarray,
    positions: np.ndarray,
    *,
    method: str,
    minimum_sites_per_class: int,
) -> dict[str, object]:
    labels = sites["chip_label"].to_numpy(dtype=int)
    metrics = binary_metrics(labels, score)
    positive = int(np.sum(labels == 1))
    negative = int(np.sum(labels == 0))
    first = sites.iloc[0]
    return {
        "cell": str(first["cell"]),
        "tf": str(first["tf"]),
        "motif_family": str(first["motif_family"]),
        "method": method,
        "split": "test",
        "status": (
            "eligible"
            if min(positive, negative) >= minimum_sites_per_class
            else "underpowered"
        ),
        "n_sites": int(len(sites)),
        "n_positive": positive,
        "n_negative": negative,
        "functional_separation": standardized_functional_separation(
            profiles, labels, positions
        ),
        **{
            key: metrics[key]
            for key in ("auroc", "auprc", "brier", "prevalence")
        },
    }


def add_baseline_deltas(metrics: pd.DataFrame) -> pd.DataFrame:
    keys = ["cell", "tf"]
    dwm = metrics[metrics["method"].eq("DWM_conventional_geometry")][
        keys + ["auroc", "auprc", "functional_separation"]
    ].rename(
        columns={
            "auroc": "dwm_auroc",
            "auprc": "dwm_auprc",
            "functional_separation": "dwm_functional_separation",
        }
    )
    raw = metrics[metrics["method"].eq("raw_geometry")][
        keys + ["auroc", "auprc", "functional_separation"]
    ].rename(
        columns={
            "auroc": "raw_auroc",
            "auprc": "raw_auprc",
            "functional_separation": "raw_functional_separation",
        }
    )
    output = metrics.merge(dwm, on=keys, validate="many_to_one")
    output = output.merge(raw, on=keys, validate="many_to_one")
    output["auroc_gain_over_dwm"] = output["auroc"] - output["dwm_auroc"]
    output["relative_auprc_gain_over_dwm"] = (
        output["auprc"] - output["dwm_auprc"]
    ) / output["dwm_auprc"].clip(lower=1e-8)
    output["functional_separation_relative_change_over_dwm"] = (
        output["functional_separation"]
        / output["dwm_functional_separation"].clip(lower=1e-8)
        - 1.0
    )
    output["auroc_gain_over_raw"] = output["auroc"] - output["raw_auroc"]
    output["relative_auprc_gain_over_raw"] = (
        output["auprc"] - output["raw_auprc"]
    ) / output["raw_auprc"].clip(lower=1e-8)
    return output


def input_record(path: Path, purpose: str) -> dict[str, object]:
    return {
        "purpose": purpose,
        "path": str(path),
        "bytes": int(path.stat().st_size),
        "sha256": file_sha256(path),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument(
        "--test-artifact", type=parse_cell_path, action="append", required=True
    )
    parser.add_argument(
        "--dwm-expected", type=parse_cell_path, action="append", required=True
    )
    parser.add_argument(
        "--bias-prediction", type=parse_cell_path, action="append", required=True
    )
    parser.add_argument(
        "--regulatory-prediction", type=parse_cell_path, action="append", required=True
    )
    parser.add_argument(
        "--nobias-prediction", type=parse_cell_path, action="append", required=True
    )
    parser.add_argument(
        "--reference-manifest", type=Path, action="append", default=[]
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimum-sites-per-class", type=int, default=200)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--aggregate-bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--dispersion", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=2026)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.minimum_sites_per_class < 1:
        raise ValueError("minimum sites per class must be positive")
    if args.bootstrap_iterations < 1 or args.aggregate_bootstrap_iterations < 0:
        raise ValueError("bootstrap iteration counts are invalid")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    artifacts = unique_cell_paths(args.test_artifact, "test artifact")
    tracks = {
        "dwm": unique_cell_paths(args.dwm_expected, "DWM track"),
        "bias": unique_cell_paths(args.bias_prediction, "bias prediction"),
        "regulatory": unique_cell_paths(
            args.regulatory_prediction, "regulatory prediction"
        ),
        "nobias": unique_cell_paths(args.nobias_prediction, "no-bias prediction"),
    }
    cells = set(artifacts)
    if not cells or any(set(values) != cells for values in tracks.values()):
        raise ValueError("all test artifacts and prediction arms must have identical cells")
    reference_documents = [
        validate_reference_manifest(path) for path in args.reference_manifest
    ]
    args.outdir.mkdir(parents=True, exist_ok=True)

    loaded: dict[str, tuple[pd.DataFrame, dict[str, np.ndarray], dict]] = {}
    inputs = [input_record(args.study, "study")]
    for cell in sorted(cells):
        document, arrays = preflight_test_artifact(
            artifacts[cell], expected_cell=cell
        )
        sites = load_test_sites(document, arrays, study)
        loaded[cell] = (sites, arrays, document)
        inputs.extend(
            [
                input_record(artifacts[cell], f"{cell}-test-manifest"),
                input_record(Path(document["profiles_npz"]), f"{cell}-test-profiles"),
                input_record(Path(document["sites"]), f"{cell}-test-sites"),
            ]
        )
        for label, paths in tracks.items():
            validate_bigwig(paths[cell])
            inputs.append(input_record(paths[cell], f"{cell}-{label}"))
    for path in args.reference_manifest:
        inputs.append(input_record(path, "ChromBPNet-stage-manifest"))
    freeze = {
        "schema": INPUT_SCHEMA,
        "cells": sorted(cells),
        "pinned_chrombpnet_commit": PINNED_COMMIT,
        "reference_stages": [
            {"stage": document["stage"], "source_commit": document["source_commit"]}
            for document in reference_documents
        ],
        "common_finite_support": True,
        "prediction_totals_matched_to_observed": True,
        "frozen_parametric_expected_from_test_artifact": True,
        "head_to_head_parametric_vs_chrombpnet_bias": True,
        "test_labels_previously_opened_for_other_models": True,
        "external_reference_used_for_policy_selection": False,
        "minimum_sites_per_class": args.minimum_sites_per_class,
        "bootstrap_iterations": args.bootstrap_iterations,
        "aggregate_bootstrap_iterations": args.aggregate_bootstrap_iterations,
        "dispersion": args.dispersion,
        "seed": args.seed,
        "inputs": inputs,
    }
    freeze["input_id"] = canonical_id(freeze)
    freeze_path = args.outdir / "chrombpnet_comparison_inputs.freeze.json"
    immutable_json(freeze_path, freeze)

    metric_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []
    score_frames: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    for cell in sorted(cells):
        sites, arrays, document = loaded[cell]
        observed = np.asarray(arrays["plus_observed"] + arrays["minus_observed"])
        parametric_expected = np.asarray(
            arrays["plus_expected"] + arrays["minus_expected"]
        )
        flank = int(document["metadata"]["flank"])
        positions = np.arange(-flank, flank + 1, dtype=float)
        extracted: dict[str, np.ndarray] = {}
        valid = np.asarray(arrays["valid"], dtype=bool).copy()
        for label, paths in tracks.items():
            values, method_valid = extract_prediction_profiles(
                sites,
                paths[cell],
                flank,
                require_dense=label != "dwm",
            )
            extracted[label] = values
            valid &= method_valid
        for tf, task in sites.groupby("tf", sort=True):
            indexes = task.index.to_numpy(dtype=int)
            common = indexes[valid[indexes]]
            coverage_rows.append(
                {
                    "cell": cell,
                    "tf": str(tf),
                    "motif_family": str(task.iloc[0]["motif_family"]),
                    "total_sites": int(len(indexes)),
                    "common_finite_sites": int(len(common)),
                    "common_finite_fraction": float(len(common) / max(len(indexes), 1)),
                    "common_positive_sites": int(
                        np.sum(sites.iloc[common]["chip_label"].to_numpy(dtype=int) == 1)
                    ),
                    "common_negative_sites": int(
                        np.sum(sites.iloc[common]["chip_label"].to_numpy(dtype=int) == 0)
                    ),
                }
            )
            if not len(common):
                continue
            task_sites = sites.iloc[common].reset_index(drop=True)
            methods = score_methods(
                observed[common],
                extracted["dwm"][common],
                parametric_expected[common],
                extracted["bias"][common],
                extracted["regulatory"][common],
                extracted["nobias"][common],
                positions,
                dispersion=args.dispersion,
            )
            labels = task_sites["chip_label"].to_numpy(dtype=int)
            if np.unique(labels).size != 2:
                continue
            score_frame = task_sites[
                [
                    "cell",
                    "tf",
                    "motif_family",
                    "TFBS_chr",
                    "TFBS_start",
                    "TFBS_end",
                    "TFBS_strand",
                    "site_id",
                    "chip_label",
                ]
            ].copy()
            for method, (score, profiles) in methods.items():
                score_frame[method] = score
                metric_rows.append(
                    metric_row(
                        task_sites,
                        score,
                        profiles,
                        positions,
                        method=method,
                        minimum_sites_per_class=args.minimum_sites_per_class,
                    )
                )
                curve = aggregate_curve(
                    profiles,
                    labels,
                    task_sites["TFBS_chr"].astype(str).to_numpy(),
                    iterations=args.aggregate_bootstrap_iterations,
                    seed=args.seed,
                )
                descriptors = asdict(
                    profile_descriptors(curve["difference"], positions)
                )
                for offset, position in enumerate(positions.astype(int)):
                    profile_rows.append(
                        {
                            "cell": cell,
                            "tf": str(tf),
                            "motif_family": str(task.iloc[0]["motif_family"]),
                            "method": method,
                            "position": int(position),
                            "positive_mean": curve["positive_mean"][offset],
                            "negative_mean": curve["negative_mean"][offset],
                            "positive_minus_negative": curve["difference"][offset],
                            "lower_95": curve["lower_95"][offset],
                            "upper_95": curve["upper_95"][offset],
                            **descriptors,
                        }
                    )
            score_frames.append(score_frame)
            if min(np.sum(labels == 1), np.sum(labels == 0)) < args.minimum_sites_per_class:
                continue
            comparison_targets = {
                "DWM_conventional_geometry": [
                    method for method in methods if method != "DWM_conventional_geometry"
                ],
                "raw_geometry": [
                    method for method in methods if method != "raw_geometry"
                ],
                "ChromBPNet_bias_conventional_geometry": [
                    "frozen_parametric_bias_conventional_geometry"
                ],
            }
            for baseline_name, target_methods in comparison_targets.items():
                baseline_score = methods[baseline_name][0]
                for method in target_methods:
                    score = methods[method][0]
                    bootstrap_rows.append(
                        {
                            "cell": cell,
                            "tf": str(tf),
                            "motif_family": str(task.iloc[0]["motif_family"]),
                            "method": method,
                            "baseline": baseline_name,
                            **block_bootstrap_delta(
                                task_sites,
                                score,
                                baseline_score,
                                iterations=args.bootstrap_iterations,
                                seed=args.seed,
                            ),
                        }
                    )

    metrics = add_baseline_deltas(pd.DataFrame(metric_rows))
    bootstrap = pd.DataFrame(bootstrap_rows)
    profiles = pd.DataFrame(profile_rows)
    site_scores = pd.concat(score_frames, ignore_index=True)
    coverage = pd.DataFrame(coverage_rows)
    paths = {
        "metrics": args.outdir / "chrombpnet_comparison_metrics.tsv",
        "bootstrap": args.outdir / "chrombpnet_comparison_bootstrap.tsv",
        "profiles": args.outdir / "chrombpnet_comparison_profiles.tsv.gz",
        "site_scores": args.outdir / "chrombpnet_comparison_site_scores.tsv.gz",
        "coverage": args.outdir / "chrombpnet_comparison_coverage.tsv",
    }
    metrics.to_csv(paths["metrics"], sep="\t", index=False)
    bootstrap.to_csv(paths["bootstrap"], sep="\t", index=False)
    profiles.to_csv(paths["profiles"], sep="\t", index=False)
    site_scores.to_csv(paths["site_scores"], sep="\t", index=False)
    coverage.to_csv(paths["coverage"], sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "input_id": freeze["input_id"],
        "pinned_chrombpnet_commit": PINNED_COMMIT,
        "external_reference_used_for_policy_selection": False,
        "outputs": {
            name: {"path": str(path), "sha256": file_sha256(path)}
            for name, path in paths.items()
        },
    }
    manifest["result_id"] = canonical_id(manifest)
    immutable_json(args.outdir / "chrombpnet_comparison_manifest.json", manifest)
    print(metrics.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
