#!/usr/bin/env python3
"""Fit and evaluate the frozen parametric bias factorization by chromosome.

Tune mode fits without ChIP labels on training chromosomes and scores only
chr16--18.  It then freezes the winning residual.  Test mode requires that
freeze and opens only chr19--22/X without refitting.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_parametric_bias import ControlWindowDataset  # noqa: E402
from fp_tools.tools.functional_footprints import (  # noqa: E402
    orient_profiles,
    standardized_functional_separation,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    ConditionalSequenceBiasModel,
    calibrated_residuals,
)
from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenBiasStrengthCalibrator,
    FrozenParametricFactorization,
    expected_profile_counts,
)


RESIDUALS = ("difference", "pearson", "deviance", "log-ratio", "nb-center-flank")
RESIDUAL_TIE_ORDER = (
    "deviance",
    "pearson",
    "difference",
    "log-ratio",
    "nb-center-flank",
)


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_name_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("values must use CELL=PATH")
    name, raw_path = value.split("=", 1)
    if not name or not raw_path:
        raise argparse.ArgumentTypeError("values must use CELL=PATH")
    return name, Path(raw_path)


def artifact_paths(prefix: Path) -> tuple[Path, Path, Path]:
    if prefix.suffix in {".npz", ".json"}:
        prefix = prefix.with_suffix("")
    return (
        Path(str(prefix) + ".npz"),
        Path(str(prefix) + ".json"),
        Path(str(prefix) + ".sites.tsv.gz"),
    )


def load_profiles(
    prefix: Path, *, require_log_bias: bool
) -> tuple[dict[str, np.ndarray], pd.DataFrame, dict]:
    npz_path, json_path, sites_path = artifact_paths(prefix)
    document = json.loads(json_path.read_text(encoding="utf-8"))
    if document.get("profiles_sha256") != sha256_file(npz_path):
        raise ValueError(f"profile checksum mismatch: {npz_path}")
    if document.get("sites_sha256") != sha256_file(sites_path):
        raise ValueError(f"site-table checksum mismatch: {sites_path}")
    with np.load(npz_path, allow_pickle=False) as arrays:
        values = {key: np.asarray(arrays[key]) for key in arrays.files}
    required = {
        "plus_observed",
        "minus_observed",
        "plus_expected",
        "minus_expected",
        "valid",
        "site_hash",
    }
    if require_log_bias:
        required.add("combined_log_bias")
    missing = required.difference(values)
    if missing:
        raise ValueError(
            f"profile artifact is missing arrays: {', '.join(sorted(missing))}"
        )
    sites = pd.read_csv(sites_path, sep="\t")
    if len(sites) != len(values["site_hash"]):
        raise ValueError("profile arrays and site table have different lengths")
    return values, sites, document


def load_dwm_baseline(prefix: Path) -> tuple[dict[str, np.ndarray], list[Path]]:
    """Load a verified conventional-DWM expected profile artifact.

    The depth-matrix cache is a single NPZ with an embedded source identity.
    Full strand artifacts are accepted only when their metadata explicitly
    identifies DWM.  This prevents older parametric profiles from being
    silently relabeled as the conventional baseline.
    """

    if prefix.suffix == ".npz" and prefix.is_file():
        with np.load(prefix, allow_pickle=False) as arrays:
            required = {"profiles", "valid", "site_hash", "signal_identity"}
            missing = required.difference(arrays.files)
            if missing:
                raise ValueError(
                    "DWM baseline cache is missing arrays: "
                    + ", ".join(sorted(missing))
                )
            identity = str(np.asarray(arrays["signal_identity"]).item())
            if "fp_tools_dwm" not in identity.lower():
                raise ValueError(
                    "baseline cache does not identify a conventional DWM track"
                )
            output = {
                "expected": np.asarray(arrays["profiles"], dtype=float),
                "valid": np.asarray(arrays["valid"], dtype=bool),
                "site_hash": np.asarray(arrays["site_hash"], dtype=np.uint64),
                "orientation_aligned": np.asarray(False),
            }
        return output, [prefix]

    arrays, _sites, document = load_profiles(prefix, require_log_bias=False)
    metadata = document.get("metadata", {})
    correction = str(
        metadata.get("correction", metadata.get("bias_method", ""))
    ).lower()
    bias_model = str(metadata.get("bias_model", "")).lower()
    if (
        correction != "dwm"
        and "atacbias" not in bias_model
        and "fp_tools_dwm" not in bias_model
    ):
        raise ValueError(
            "baseline strand artifact is not explicitly identified as conventional DWM"
        )
    return (
        {
            "expected": arrays["plus_expected"] + arrays["minus_expected"],
            "valid": arrays["valid"].astype(bool),
            "site_hash": arrays["site_hash"],
            "orientation_aligned": np.asarray(True),
        },
        list(artifact_paths(prefix)),
    )


def load_pwm_baseline(prefix: Path) -> tuple[dict[str, np.ndarray], list[Path]]:
    """Load a cached expected profile explicitly sourced from the PWM arm."""

    if prefix.suffix != ".npz" or not prefix.is_file():
        raise ValueError("PWM baseline must be a direct expected-profile NPZ cache")
    with np.load(prefix, allow_pickle=False) as arrays:
        required = {"profiles", "valid", "site_hash", "signal_identity"}
        missing = required.difference(arrays.files)
        if missing:
            raise ValueError(
                "PWM baseline cache is missing arrays: " + ", ".join(sorted(missing))
            )
        identity = str(np.asarray(arrays["signal_identity"]).item())
        if "fp_tools_pwm" not in identity.lower():
            raise ValueError(
                "baseline cache does not identify a conventional PWM track"
            )
        output = {
            "expected": np.asarray(arrays["profiles"], dtype=float),
            "valid": np.asarray(arrays["valid"], dtype=bool),
            "site_hash": np.asarray(arrays["site_hash"], dtype=np.uint64),
            "orientation_aligned": np.asarray(False),
        }
    return output, [prefix]


def align_baseline(
    candidate: dict[str, np.ndarray],
    baseline: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    lookup = {int(value): index for index, value in enumerate(baseline["site_hash"])}
    try:
        order = np.asarray(
            [lookup[int(value)] for value in candidate["site_hash"]], dtype=int
        )
    except KeyError as exc:
        raise ValueError(
            f"baseline is missing candidate site hash {exc.args[0]}"
        ) from exc
    return baseline["expected"][order], baseline["valid"][order]


def orient_aligned_baseline(
    expected: np.ndarray,
    baseline: dict[str, np.ndarray],
    sites: pd.DataFrame,
) -> np.ndarray:
    """Orient direct genomic caches after aligning them to candidate rows."""

    if bool(np.asarray(baseline["orientation_aligned"]).item()):
        return np.asarray(expected, dtype=float)
    return orient_profiles(expected, sites["TFBS_strand"].astype(str))


def geometry_score(profiles: np.ndarray, positions: np.ndarray) -> np.ndarray:
    values = np.asarray(profiles, dtype=np.float64)
    center = np.abs(positions) <= 7
    shoulders = (np.abs(positions) >= 12) & (np.abs(positions) <= 40)
    if values.ndim != 2 or values.shape[1] != len(positions):
        raise ValueError("profiles must match geometry positions")
    return np.mean(values[:, shoulders], axis=1) - np.mean(values[:, center], axis=1)


def expected_calibration_error(
    labels: np.ndarray, probabilities: np.ndarray, bins: int = 10
) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    value = 0.0
    for index in range(bins):
        member = (probabilities >= edges[index]) & (
            probabilities <= edges[index + 1]
            if index + 1 == bins
            else probabilities < edges[index + 1]
        )
        if np.any(member):
            value += float(np.mean(member)) * abs(
                float(np.mean(labels[member])) - float(np.mean(probabilities[member]))
            )
    return float(value)


def metric_row(
    sites: pd.DataFrame,
    score: np.ndarray,
    profiles: np.ndarray,
    *,
    cell: str,
    tf: str,
    method: str,
    split: str,
    positions: np.ndarray,
    probability: bool = False,
) -> dict[str, object] | None:
    selected = (sites["cell"].astype(str) == cell) & (sites["tf"].astype(str) == tf)
    labels = sites.loc[selected, "chip_label"].to_numpy(dtype=int)
    if len(labels) < 2 or np.unique(labels).size != 2:
        return None
    values = np.asarray(score[selected], dtype=np.float64)
    tf_profiles = np.asarray(profiles[selected], dtype=np.float64)
    row = sites.loc[selected].iloc[0]
    return {
        "cell": cell,
        "tf": tf,
        "motif_id": str(row["motif"]),
        "motif_family": str(row["motif_family"]),
        "role": str(row["role"]),
        "split": split,
        "method": method,
        "n_sites": int(len(labels)),
        "n_positive": int(np.sum(labels == 1)),
        "n_negative": int(np.sum(labels == 0)),
        "auroc": float(roc_auc_score(labels, values)),
        "auprc": float(average_precision_score(labels, values)),
        "brier": float(brier_score_loss(labels, np.clip(values, 0, 1)))
        if probability
        else np.nan,
        "calibration_error": expected_calibration_error(labels, np.clip(values, 0, 1))
        if probability
        else np.nan,
        "functional_separation": standardized_functional_separation(
            tf_profiles, labels, positions
        ),
    }


def block_bootstrap_delta(
    sites: pd.DataFrame,
    candidate: np.ndarray,
    baseline: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> dict[str, float | int]:
    chromosomes = np.asarray(sorted(sites["TFBS_chr"].astype(str).unique()))
    chrom = sites["TFBS_chr"].astype(str).to_numpy()
    labels = sites["chip_label"].to_numpy(dtype=int)
    rng = np.random.default_rng(seed)
    auroc, relative_ap = [], []
    for _ in range(iterations):
        sampled = rng.choice(chromosomes, size=len(chromosomes), replace=True)
        indexes = np.concatenate([np.flatnonzero(chrom == value) for value in sampled])
        if np.unique(labels[indexes]).size != 2:
            continue
        base_auc = roc_auc_score(labels[indexes], baseline[indexes])
        candidate_auc = roc_auc_score(labels[indexes], candidate[indexes])
        base_ap = average_precision_score(labels[indexes], baseline[indexes])
        candidate_ap = average_precision_score(labels[indexes], candidate[indexes])
        auroc.append(candidate_auc - base_auc)
        relative_ap.append((candidate_ap - base_ap) / max(base_ap, 1e-8))
    output: dict[str, float | int] = {"bootstrap_successful": len(auroc)}
    for name, values in (("auroc_gain", auroc), ("relative_auprc_gain", relative_ap)):
        array = np.asarray(values, dtype=float)
        output[f"{name}_lower_95"] = (
            float(np.quantile(array, 0.025)) if len(array) else np.nan
        )
        output[f"{name}_upper_95"] = (
            float(np.quantile(array, 0.975)) if len(array) else np.nan
        )
        output[f"{name}_probability_positive"] = (
            float(np.mean(array > 0)) if len(array) else np.nan
        )
    return output


def select_label_free_sites(sites: pd.DataFrame, maximum_per_tf: int) -> np.ndarray:
    selected: list[int] = []
    for (_cell, _tf), group in sites.groupby(["cell", "tf"], sort=True):
        indexes = group.index.to_numpy(dtype=int)
        if len(indexes) > maximum_per_tf:
            hashes = group["site_hash"].to_numpy(dtype=np.uint64)
            indexes = indexes[np.argsort(hashes, kind="mergesort")[:maximum_per_tf]]
        selected.extend(indexes.tolist())
    return np.asarray(sorted(selected), dtype=int)


def calibration_from_controls(
    cells: list[str],
    controls: dict[str, Path],
    bias_models: dict[str, Path],
) -> FrozenBiasStrengthCalibrator:
    all_counts, all_scores, all_samples = [], [], []
    for cell in cells:
        dataset = ControlWindowDataset.load(controls[cell])
        if dataset.split != "train":
            raise ValueError(
                f"calibration control for {cell} is not a training-split artifact"
            )
        model = ConditionalSequenceBiasModel.load(bias_models[cell])
        contexts, counts = dataset.model_arrays(model.feature_spec)
        all_counts.append(counts)
        all_scores.append(model.log_scores(contexts))
        all_samples.extend([cell] * len(counts))
    return FrozenBiasStrengthCalibrator().fit(
        np.concatenate(all_counts),
        np.concatenate(all_scores),
        all_samples,
    )


def collect_datasets(
    candidates: dict[str, Path],
    baselines: dict[str, Path],
    pwm_baselines: dict[str, Path] | None = None,
) -> tuple[dict[str, np.ndarray], pd.DataFrame, list[Path]]:
    arrays: dict[str, list[np.ndarray]] = {}
    site_frames = []
    baseline_inputs: list[Path] = []
    for cell in sorted(candidates):
        candidate, sites, _document = load_profiles(
            candidates[cell], require_log_bias=True
        )
        baseline, baseline_paths = load_dwm_baseline(baselines[cell])
        baseline_inputs.extend(baseline_paths)
        baseline_expected, baseline_valid = align_baseline(candidate, baseline)
        baseline_expected = orient_aligned_baseline(
            baseline_expected, baseline, sites
        )
        pwm_expected = None
        pwm_valid = np.ones(len(sites), dtype=bool)
        if pwm_baselines is not None:
            pwm, pwm_paths = load_pwm_baseline(pwm_baselines[cell])
            baseline_inputs.extend(pwm_paths)
            pwm_expected, pwm_valid = align_baseline(candidate, pwm)
            pwm_expected = orient_aligned_baseline(pwm_expected, pwm, sites)
        sites = sites.copy()
        sites["cell"] = cell
        sites["site_hash"] = candidate["site_hash"]
        arrays.setdefault("counts", []).append(
            candidate["plus_observed"] + candidate["minus_observed"]
        )
        arrays.setdefault("parametric_expected", []).append(
            candidate["plus_expected"] + candidate["minus_expected"]
        )
        arrays.setdefault("baseline_expected", []).append(baseline_expected)
        if pwm_expected is not None:
            arrays.setdefault("pwm_expected", []).append(pwm_expected)
        arrays.setdefault("log_bias", []).append(candidate["combined_log_bias"])
        arrays.setdefault("valid", []).append(
            candidate["valid"].astype(bool)
            & baseline_valid
            & pwm_valid
            & np.isfinite(candidate["combined_log_bias"]).all(axis=1)
        )
        site_frames.append(sites)
    return (
        {key: np.concatenate(values) for key, values in arrays.items()},
        pd.concat(site_frames, ignore_index=True),
        baseline_inputs,
    )


def residual_score(
    counts: np.ndarray,
    expected: np.ndarray,
    positions: np.ndarray,
    residual: str,
    dispersion: float,
) -> tuple[np.ndarray, np.ndarray]:
    if residual == "nb-center-flank":
        from fp_tools.tools.parametric_bias import center_flank_likelihood_score

        score = center_flank_likelihood_score(counts, expected, dispersion=dispersion)
        profiles = calibrated_residuals(
            counts, expected, "deviance", dispersion=dispersion
        )
        return score, profiles
    profiles = calibrated_residuals(counts, expected, residual, dispersion=dispersion)
    return geometry_score(profiles, positions), profiles


def score_methods(
    counts: np.ndarray,
    parametric_expected: np.ndarray,
    baseline_expected: np.ndarray,
    pwm_expected: np.ndarray | None,
    log_bias: np.ndarray,
    strengths: np.ndarray,
    factor_result,
    positions: np.ndarray,
    dispersion: float,
    residuals: Iterable[str],
) -> dict[str, tuple[np.ndarray, np.ndarray, bool]]:
    methods: dict[str, tuple[np.ndarray, np.ndarray, bool]] = {}
    raw_score = geometry_score(counts, positions)
    methods["raw"] = (raw_score, counts, False)
    baseline_score, baseline_profiles = residual_score(
        counts, baseline_expected, positions, "deviance", dispersion
    )
    methods["DWM"] = (baseline_score, baseline_profiles, False)
    if pwm_expected is not None:
        pwm_score, pwm_profiles = residual_score(
            counts, pwm_expected, positions, "deviance", dispersion
        )
        methods["PWM"] = (pwm_score, pwm_profiles, False)
    lambda_expected = expected_profile_counts(counts, strengths[:, None] * log_bias)
    for residual in residuals:
        direct_score, direct_profiles = residual_score(
            counts, parametric_expected, positions, residual, dispersion
        )
        methods[f"parametric_direct_{residual}"] = (
            direct_score,
            direct_profiles,
            False,
        )
        lambda_score, lambda_profiles = residual_score(
            counts, lambda_expected, positions, residual, dispersion
        )
        methods[f"parametric_lambda_{residual}"] = (
            lambda_score,
            lambda_profiles,
            False,
        )
        factor_score, factor_profiles = residual_score(
            counts,
            factor_result.expected_unbound,
            positions,
            residual,
            dispersion,
        )
        methods[f"factorized_residual_{residual}"] = (
            factor_score,
            factor_profiles,
            False,
        )
    methods["factorized_posterior"] = (
        factor_result.posterior_bound,
        calibrated_residuals(
            counts, factor_result.expected_unbound, "deviance", dispersion=dispersion
        ),
        True,
    )
    return methods


def select_residual(metrics: pd.DataFrame) -> tuple[str, pd.DataFrame]:
    baseline = metrics[metrics["method"] == "DWM"][
        ["cell", "tf", "auroc", "auprc"]
    ].rename(columns={"auroc": "baseline_auroc", "auprc": "baseline_auprc"})
    candidates = metrics[
        metrics["method"].str.startswith("factorized_residual_")
    ].merge(
        baseline,
        on=["cell", "tf"],
        validate="many_to_one",
    )
    candidates["auroc_gain"] = candidates["auroc"] - candidates["baseline_auroc"]
    candidates["relative_auprc_gain"] = (
        candidates["auprc"] - candidates["baseline_auprc"]
    ) / candidates["baseline_auprc"].clip(lower=1e-8)
    candidates["residual"] = candidates["method"].str.removeprefix(
        "factorized_residual_"
    )
    ranking = candidates[candidates["role"].astype(str) == "difficult"].copy()
    if ranking.empty:
        raise RuntimeError(
            "residual selection requires at least one difficult development task"
        )
    summaries = ranking.groupby("residual", as_index=False).agg(
        mean_relative_auprc_gain=("relative_auprc_gain", "mean"),
        mean_auroc_gain=("auroc_gain", "mean"),
        task_sd_relative_auprc_gain=("relative_auprc_gain", "std"),
        task_count=("tf", "size"),
    )
    ctcf = (
        candidates[candidates["tf"] == "CTCF"].groupby("residual")["auroc_gain"].min()
    )
    summaries["minimum_ctcf_auroc_gain"] = summaries["residual"].map(ctcf)
    summaries["passes_ctcf_gate"] = summaries["minimum_ctcf_auroc_gain"] >= -0.02
    summaries["standard_error"] = summaries["task_sd_relative_auprc_gain"].fillna(
        0.0
    ) / np.sqrt(summaries["task_count"].clip(lower=1))
    eligible = summaries[summaries["passes_ctcf_gate"]].copy()
    if eligible.empty:
        raise RuntimeError("every residual failed the frozen CTCF loss gate")
    best = eligible.sort_values(
        "mean_relative_auprc_gain", ascending=False, kind="mergesort"
    ).iloc[0]
    threshold = float(best["mean_relative_auprc_gain"] - best["standard_error"])
    within = eligible[eligible["mean_relative_auprc_gain"] >= threshold]
    tie_rank = {name: index for index, name in enumerate(RESIDUAL_TIE_ORDER)}
    selected = min(within["residual"].astype(str), key=lambda value: tie_rank[value])
    summaries["within_one_se"] = summaries["mean_relative_auprc_gain"] >= threshold
    summaries["selected"] = summaries["residual"] == selected
    return selected, summaries.sort_values(
        ["selected", "mean_relative_auprc_gain"],
        ascending=[False, False],
        kind="mergesort",
    )


def evaluate_split(
    arrays: dict[str, np.ndarray],
    sites: pd.DataFrame,
    indexes: np.ndarray,
    model: FrozenParametricFactorization,
    *,
    split_name: str,
    residuals: Iterable[str],
    bootstrap_iterations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, tuple[np.ndarray, np.ndarray, bool]]]:
    selected_sites = sites.iloc[indexes].reset_index(drop=True)
    counts = arrays["counts"][indexes]
    log_bias = arrays["log_bias"][indexes]
    samples = selected_sites["cell"].astype(str).to_numpy()
    tfs = selected_sites["tf"].astype(str).to_numpy()
    result = model.predict(counts, log_bias, samples, tfs)
    strengths = np.asarray([model.bias_strengths_[sample] for sample in samples])
    methods = score_methods(
        counts,
        arrays["parametric_expected"][indexes],
        arrays["baseline_expected"][indexes],
        arrays["pwm_expected"][indexes] if "pwm_expected" in arrays else None,
        log_bias,
        strengths,
        result,
        model.positions,
        model.total_dispersion_,
        residuals,
    )
    metrics = []
    for (cell, tf), _group in selected_sites.groupby(["cell", "tf"], sort=True):
        for method, (score, profiles, probability) in methods.items():
            row = metric_row(
                selected_sites,
                score,
                profiles,
                cell=str(cell),
                tf=str(tf),
                method=method,
                split=split_name,
                positions=model.positions,
                probability=probability,
            )
            if row is not None:
                metrics.append(row)
    metrics_frame = pd.DataFrame(metrics)
    bootstrap_rows = []
    baseline_score = methods["DWM"][0]
    for (cell, tf), group in selected_sites.groupby(["cell", "tf"], sort=True):
        task_indexes = group.index.to_numpy(dtype=int)
        for method, (score, _profiles, _probability) in methods.items():
            if method == "DWM":
                continue
            bootstrap_rows.append(
                {
                    "cell": cell,
                    "tf": tf,
                    "split": split_name,
                    "method": method,
                    **block_bootstrap_delta(
                        selected_sites.iloc[task_indexes],
                        score[task_indexes],
                        baseline_score[task_indexes],
                        iterations=bootstrap_iterations,
                        seed=seed,
                    ),
                }
            )
    return metrics_frame, pd.DataFrame(bootstrap_rows), methods


def _configuration_document(
    *,
    output: Path,
    selected_residual: str,
    model_path: Path,
    calibration_path: Path,
    study_path: Path,
    inputs: list[Path],
) -> dict:
    document = {
        "schema": "fp-tools-parametric-factorization-configuration-freeze-v1",
        "selected_residual": selected_residual,
        "factorization_model": {
            "path": str(model_path),
            "sha256": sha256_file(model_path),
        },
        "factorization_model_metadata": {
            "path": str(model_path.with_suffix(".json")),
            "sha256": sha256_file(model_path.with_suffix(".json")),
        },
        "bias_calibration": {
            "path": str(calibration_path),
            "sha256": sha256_file(calibration_path),
        },
        "bias_calibration_metadata": {
            "path": str(calibration_path.with_suffix(".json")),
            "sha256": sha256_file(calibration_path.with_suffix(".json")),
        },
        "study": {"path": str(study_path), "sha256": sha256_file(study_path)},
        "inputs": [{"path": str(path), "sha256": sha256_file(path)} for path in inputs],
        "test_labels_opened": False,
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["configuration_id"] = sha256(canonical.encode()).hexdigest()
    output.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return document


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument(
        "--candidate", type=parse_name_path, action="append", required=True
    )
    parser.add_argument(
        "--baseline", type=parse_name_path, action="append", required=True
    )
    parser.add_argument("--pwm-baseline", type=parse_name_path, action="append")
    parser.add_argument("--control", type=parse_name_path, action="append")
    parser.add_argument("--bias-model", type=parse_name_path, action="append")
    parser.add_argument("--mode", choices=("tune", "test"), required=True)
    parser.add_argument("--configuration-freeze", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--maximum-training-sites-per-tf", type=int, default=10000)
    parser.add_argument("--max-iter", type=int, default=25)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    candidates = dict(args.candidate)
    baselines = dict(args.baseline)
    pwm_baselines = dict(args.pwm_baseline or [])
    if set(candidates) != set(baselines):
        raise ValueError("candidate and baseline cells must match")
    if pwm_baselines and set(pwm_baselines) != set(candidates):
        raise ValueError("PWM baseline and candidate cells must match")
    arrays, sites, baseline_inputs = collect_datasets(
        candidates, baselines, pwm_baselines or None
    )
    valid = arrays["valid"] & (arrays["counts"].sum(axis=1) > 0)
    sites = sites.copy()
    positions = (
        np.arange(arrays["counts"].shape[1], dtype=float)
        - arrays["counts"].shape[1] // 2
    )
    args.outdir.mkdir(parents=True, exist_ok=True)

    input_paths = []
    for prefix in candidates.values():
        input_paths.extend(artifact_paths(prefix))
    input_paths.extend(baseline_inputs)

    if args.mode == "tune":
        if args.configuration_freeze is not None:
            raise ValueError(
                "tune mode creates the configuration freeze; it does not accept one"
            )
        controls = dict(args.control or [])
        bias_models = dict(args.bias_model or [])
        if set(controls) != set(candidates) or set(bias_models) != set(candidates):
            raise ValueError(
                "tune mode requires one control and bias model per candidate cell"
            )
        calibration = calibration_from_controls(
            sorted(candidates), controls, bias_models
        )
        calibration_path, _ = calibration.save(
            args.outdir / "frozen_bias_calibration",
            {"chromosome_split": "train", "region_class": "low-signal nonpeak"},
        )
        train_mask = (
            valid & sites["chromosome_split"].astype(str).eq("train").to_numpy()
        )
        train_sites = sites.loc[train_mask].copy().reset_index()
        train_sites["site_hash"] = (
            arrays["site_hash"][train_mask]
            if "site_hash" in arrays
            else sites.loc[train_mask, "site_hash"].to_numpy()
        )
        selected_local = select_label_free_sites(
            train_sites, args.maximum_training_sites_per_tf
        )
        selected = train_sites.loc[selected_local, "index"].to_numpy(dtype=int)
        model = FrozenParametricFactorization(positions, seed=args.seed).fit(
            arrays["counts"][selected],
            arrays["log_bias"][selected],
            sites.iloc[selected]["cell"].astype(str),
            sites.iloc[selected]["tf"].astype(str),
            sites.iloc[selected]["motif_family"].astype(str),
            calibration,
            max_iter=args.max_iter,
        )
        model_path, _ = model.save(
            args.outdir / "frozen_parametric_factorization",
            {
                "training_chromosomes": study["chromosome_split"]["train"],
                "labels_used": False,
            },
        )
        tune_mask = (
            valid & sites["chromosome_split"].astype(str).eq("validation").to_numpy()
        )
        tune_indexes = np.flatnonzero(tune_mask)
        metrics, bootstrap, _methods = evaluate_split(
            arrays,
            sites,
            tune_indexes,
            model,
            split_name="validation",
            residuals=RESIDUALS,
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
        )
        selected_residual, residual_summary = select_residual(metrics)
        metrics.to_csv(
            args.outdir / "factorization_tune_metrics.tsv", sep="\t", index=False
        )
        bootstrap.to_csv(
            args.outdir / "factorization_tune_bootstrap.tsv", sep="\t", index=False
        )
        residual_summary.to_csv(
            args.outdir / "factorization_residual_selection.tsv", sep="\t", index=False
        )
        _configuration_document(
            output=args.outdir / "factorization_configuration.freeze.json",
            selected_residual=selected_residual,
            model_path=model_path,
            calibration_path=calibration_path,
            study_path=args.study,
            inputs=input_paths + list(controls.values()) + list(bias_models.values()),
        )
    else:
        if args.configuration_freeze is None:
            raise ValueError("test mode requires --configuration-freeze")
        configuration = json.loads(
            args.configuration_freeze.read_text(encoding="utf-8")
        )
        if (
            configuration.get("schema")
            != "fp-tools-parametric-factorization-configuration-freeze-v1"
        ):
            raise ValueError("unsupported factorization configuration freeze")
        for record in (
            configuration["factorization_model"],
            configuration["factorization_model_metadata"],
            configuration["bias_calibration"],
            configuration["bias_calibration_metadata"],
            configuration["study"],
            *configuration["inputs"],
        ):
            if sha256_file(record["path"]) != record["sha256"]:
                raise ValueError(
                    f"frozen configuration input changed: {record['path']}"
                )
        model = FrozenParametricFactorization.load(
            configuration["factorization_model"]["path"]
        )
        test_mask = valid & sites["chromosome_split"].astype(str).eq("test").to_numpy()
        test_indexes = np.flatnonzero(test_mask)
        residual = str(configuration["selected_residual"])
        metrics, bootstrap, _methods = evaluate_split(
            arrays,
            sites,
            test_indexes,
            model,
            split_name="test",
            residuals=[residual],
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
        )
        metrics.to_csv(
            args.outdir / "factorization_test_metrics.tsv", sep="\t", index=False
        )
        bootstrap.to_csv(
            args.outdir / "factorization_test_bootstrap.tsv", sep="\t", index=False
        )
        test_record = {
            "schema": "fp-tools-parametric-factorization-test-v1",
            "configuration_id": configuration["configuration_id"],
            "configuration_freeze": str(args.configuration_freeze),
            "configuration_freeze_sha256": sha256_file(args.configuration_freeze),
            "metrics_sha256": sha256_file(
                args.outdir / "factorization_test_metrics.tsv"
            ),
            "bootstrap_sha256": sha256_file(
                args.outdir / "factorization_test_bootstrap.tsv"
            ),
            "refitted": False,
        }
        (args.outdir / "factorization_test_manifest.json").write_text(
            json.dumps(test_record, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
