#!/usr/bin/env python3
"""Evaluate label-free spline/FDA/GP footprint models on locked motif sites.

The label-free models are fitted on motif-site pools selected without ChIP
labels.  ChIP labels are read only for the matched evaluation tables and the
separate supervised information ceiling.  Test data cannot be scored unless
``--unlock-test`` is supplied after model choices have been frozen.
"""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from search_tf_footprint_models import extract_profiles  # noqa: E402
from fp_tools.tools.functional_footprints import (  # noqa: E402
    BiasAwareFunctionalMixture,
    FdaMixtureModel,
    FunctionalPCA,
    HybridFdaGpModel,
    deviance_profiles,
    normalize_functional_profiles,
    orient_profiles,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    ConditionalSequenceBiasModel,
    calibrated_residuals,
    center_flank_likelihood_score,
    combined_strand_log_bias,
    estimate_nb_dispersion,
)
from fp_tools.utils.fasta import open_fasta  # noqa: E402


LABEL_FREE_MODELS = ("spline", "gp", "fda", "hybrid")
RESIDUAL_MODES = ("difference", "pearson", "deviance", "log-ratio", "nb-likelihood")
REQUIRED_SITE_COLUMNS = {
    "cell",
    "tf",
    "motif_family",
    "TFBS_chr",
    "TFBS_start",
    "TFBS_end",
    "TFBS_strand",
    "motif_score",
    "accessibility",
    "chip_label",
    "chromosome_split",
}


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*values: object, seed: int = 2026) -> int:
    digest = blake2b(digest_size=8)
    digest.update(str(seed).encode())
    for value in values:
        digest.update(b"\0")
        digest.update(str(value).encode())
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def site_hashes(frame: pd.DataFrame) -> np.ndarray:
    values = []
    for row in frame.itertuples(index=False):
        digest = blake2b(digest_size=8)
        for value in (
            row.TFBS_chr,
            int(row.TFBS_start),
            int(row.TFBS_end),
            row.TFBS_strand,
            getattr(row, "tf", ""),
        ):
            digest.update(str(value).encode())
            digest.update(b"\0")
        values.append(int.from_bytes(digest.digest(), "little"))
    return np.asarray(values, dtype=np.uint64)


def validate_sites(frame: pd.DataFrame, path: str | Path) -> pd.DataFrame:
    frame = frame.copy()
    if "accessibility" not in frame:
        if "central_accessibility" in frame:
            frame["accessibility"] = pd.to_numeric(
                frame["central_accessibility"], errors="coerce"
            ).fillna(0.0)
        else:
            # The evaluator replaces this compatibility placeholder with the
            # summed raw profile before fitting any coverage-aware model.
            frame["accessibility"] = 0.0
    missing = REQUIRED_SITE_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing site columns: {', '.join(sorted(missing))}")
    output = frame.copy()
    output["chip_label"] = pd.to_numeric(output["chip_label"], errors="raise").astype(int)
    if not set(output["chip_label"].unique()).issubset({0, 1}):
        raise ValueError(f"{path} contains non-binary ChIP labels")
    return output


def chromosome_split(chromosome: str, study: dict) -> str:
    for split, chromosomes in study["chromosome_split"].items():
        if chromosome in chromosomes:
            return split
    return "excluded"


def parse_cell_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("unlabeled site sources must use CELL=PATH")
    cell, path = value.split("=", 1)
    if not cell or not path:
        raise argparse.ArgumentTypeError("unlabeled site sources must use CELL=PATH")
    return cell, Path(path)


def build_unlabeled_training_sites(
    source: str | Path,
    cell: str,
    tasks: pd.DataFrame,
    study: dict,
    *,
    maximum_per_tf: int,
    seed: int,
) -> pd.DataFrame:
    """Sample training motif sites without reading occupancy labels."""

    source_frame = pd.read_csv(source, sep="\t")
    required = {"motif", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"}
    missing = required.difference(source_frame.columns)
    if missing:
        raise ValueError(f"{source} is missing motif-site columns: {', '.join(sorted(missing))}")
    if any("label" in column.lower() or "chip" in column.lower() for column in source_frame.columns):
        raise ValueError("unlabeled training-site files must not contain ChIP/label columns")
    rows = []
    for task in tasks[tasks["cell"] == cell].itertuples(index=False):
        motif_id = str(task.motif_id)
        selected = source_frame[source_frame["motif"].astype(str).str.contains(motif_id, regex=False)].copy()
        selected["chromosome_split"] = selected["TFBS_chr"].map(lambda value: chromosome_split(str(value), study))
        selected = selected[selected["chromosome_split"] == "train"]
        if len(selected) > maximum_per_tf:
            selected = selected.sample(
                n=maximum_per_tf,
                random_state=stable_seed(cell, task.tf, seed=seed),
                replace=False,
            )
        selected["cell"] = cell
        selected["tf"] = str(task.tf)
        selected["motif_id"] = motif_id
        selected["motif_family"] = str(task.motif_family)
        if "TFBS_score" in selected:
            selected["motif_score"] = pd.to_numeric(selected["TFBS_score"], errors="coerce")
        elif "score" in selected:
            selected["motif_score"] = pd.to_numeric(selected["score"], errors="coerce")
        else:
            selected["motif_score"] = np.nan
        selected["accessibility"] = 0.0
        rows.append(
            selected[
                [
                    "cell",
                    "tf",
                    "motif_id",
                    "motif_family",
                    "TFBS_chr",
                    "TFBS_start",
                    "TFBS_end",
                    "TFBS_strand",
                    "motif_score",
                    "accessibility",
                    "chromosome_split",
                ]
            ]
        )
    if not rows:
        raise ValueError(f"no development tasks or motif sites were available for {cell}")
    return pd.concat(rows, ignore_index=True)


def load_or_extract_profiles(
    sites: pd.DataFrame,
    signal: str | Path,
    cache: str | Path,
    flank: int,
) -> np.ndarray:
    cache_path = Path(cache)
    expected_hashes = site_hashes(sites)
    signal_path = Path(signal).expanduser().resolve()
    signal_stat = signal_path.stat()
    signal_identity = f"{signal_path}:{signal_stat.st_size}:{signal_stat.st_mtime_ns}"
    if cache_path.is_file():
        with np.load(cache_path, allow_pickle=False) as arrays:
            profiles = np.asarray(arrays["profiles"], dtype=np.float32)
            hashes = np.asarray(arrays["site_hash"], dtype=np.uint64)
            cached_identity = str(arrays["signal_identity"].item()) if "signal_identity" in arrays else ""
        if (
            profiles.shape == (len(sites), flank * 2 + 1)
            and np.array_equal(hashes, expected_hashes)
            and cached_identity == signal_identity
        ):
            return profiles
    profiles, valid = extract_profiles(sites, signal_path, flank)
    if not np.all(valid):
        profiles = np.nan_to_num(profiles, nan=0.0)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        profiles=profiles,
        valid=valid,
        site_hash=expected_hashes,
        signal_identity=np.asarray(signal_identity),
    )
    return profiles


def derive_parametric_expected_profiles(
    sites: pd.DataFrame,
    observed: np.ndarray,
    model_path: str | Path,
    genome: str | Path,
    cache: str | Path,
    flank: int,
) -> np.ndarray:
    """Predict site-level expected cuts from a frozen parametric bias model."""

    cache_path = Path(cache)
    expected_hashes = site_hashes(sites)
    model_path = Path(model_path).expanduser().resolve()
    model_digest = file_sha256(model_path)
    raw_digest = sha256(np.ascontiguousarray(observed, dtype=np.float32).tobytes()).hexdigest()
    genome_path = Path(genome).expanduser().resolve()
    genome_stat = genome_path.stat()
    genome_identity = f"{genome_path}:{genome_stat.st_size}:{genome_stat.st_mtime_ns}"
    if cache_path.is_file():
        with np.load(cache_path, allow_pickle=False) as arrays:
            profiles = np.asarray(arrays["profiles"], dtype=np.float32)
            hashes = np.asarray(arrays["site_hash"], dtype=np.uint64)
            cached_model = str(arrays["model_sha256"].item())
            cached_raw = str(arrays["raw_sha256"].item())
            cached_genome = str(arrays["genome_identity"].item())
        if (
            profiles.shape == observed.shape
            and np.array_equal(hashes, expected_hashes)
            and cached_model == model_digest
            and cached_raw == raw_digest
            and cached_genome == genome_identity
        ):
            return profiles

    model = ConditionalSequenceBiasModel.load(model_path)
    width = flank * 2 + 1
    if observed.shape != (len(sites), width):
        raise ValueError("observed profiles do not match sites and flank")
    margin = max(41, model.feature_spec.context_length // 2 + 1)
    positions = margin + np.arange(width)
    expected = np.zeros_like(observed, dtype=np.float32)
    with open_fasta(genome_path) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for index, row in enumerate(sites.itertuples(index=False)):
            chromosome = str(row.TFBS_chr)
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            start = center - flank - margin
            end = center + flank + margin + 1
            if chromosome not in lengths or start < 0 or end > lengths[chromosome]:
                continue
            sequence = fasta.fetch(chromosome, start, end).upper()
            log_bias = combined_strand_log_bias(model, sequence, positions)
            finite = np.isfinite(log_bias)
            total = float(np.sum(observed[index]))
            if total <= 0 or not finite.any():
                continue
            centered = np.where(finite, log_bias - np.max(log_bias[finite]), -np.inf)
            propensity = np.exp(centered)
            propensity /= max(float(propensity.sum()), np.finfo(float).tiny)
            expected[index] = (total * propensity).astype(np.float32)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        profiles=expected,
        site_hash=expected_hashes,
        model_sha256=np.asarray(model_digest),
        raw_sha256=np.asarray(raw_digest),
        genome_identity=np.asarray(genome_identity),
    )
    return expected


def binary_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float | int]:
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    finite = np.isfinite(scores)
    labels = labels[finite]
    scores = scores[finite]
    output: dict[str, float | int] = {
        "n_sites": int(len(labels)),
        "positive_sites": int(np.sum(labels == 1)),
        "negative_sites": int(np.sum(labels == 0)),
        "prevalence": float(np.mean(labels)) if len(labels) else np.nan,
    }
    if len(labels) and len(np.unique(labels)) == 2:
        output.update(
            {
                "auroc": float(roc_auc_score(labels, scores)),
                "auprc": float(average_precision_score(labels, scores)),
                "brier": float(brier_score_loss(labels, np.clip(scores, 0.0, 1.0)))
                if np.all((scores >= 0) & (scores <= 1))
                else np.nan,
            }
        )
    else:
        output.update({"auroc": np.nan, "auprc": np.nan, "brier": np.nan})
    return output


def residual_score(
    observed: np.ndarray,
    expected: np.ndarray,
    mode: str,
    dispersion: float,
) -> tuple[np.ndarray, np.ndarray]:
    if mode == "nb-likelihood":
        score = center_flank_likelihood_score(
            observed,
            expected,
            center_width=15,
            flank_width=30,
            gap=5,
            dispersion=dispersion,
        )
        residual = deviance_profiles(observed, expected, dispersion)
        return residual, score
    residual = calibrated_residuals(observed, expected, mode, dispersion=dispersion)
    midpoint = residual.shape[1] // 2
    center = residual[:, midpoint - 7:midpoint + 8].mean(axis=1)
    flanks = np.concatenate(
        [residual[:, midpoint - 42:midpoint - 12], residual[:, midpoint + 13:midpoint + 43]],
        axis=1,
    ).mean(axis=1)
    return residual, flanks - center


def _score_to_probability(values: np.ndarray) -> np.ndarray:
    scores = np.asarray(values, dtype=float)
    finite = np.isfinite(scores)
    if not finite.any():
        return np.full_like(scores, 0.5)
    location = np.median(scores[finite])
    scale = 1.4826 * np.median(np.abs(scores[finite] - location))
    if scale <= 0:
        scale = np.std(scores[finite]) or 1.0
    return 1.0 / (1.0 + np.exp(-np.clip((scores - location) / scale, -30.0, 30.0)))


def fit_supervised_ceiling(
    train_residual: np.ndarray,
    validation_residual: np.ndarray,
    train_sites: pd.DataFrame,
    validation_sites: pd.DataFrame,
    *,
    seed: int,
) -> tuple[np.ndarray, dict[str, object], FunctionalPCA]:
    train_residual = normalize_functional_profiles(train_residual)
    validation_residual = normalize_functional_profiles(validation_residual)
    fpca = FunctionalPCA(variance_threshold=0.95, max_components=20, seed=seed)
    train_scores = fpca.fit_transform(train_residual)
    validation_scores = fpca.transform(validation_residual)
    train_covariates = np.column_stack(
        [
            pd.to_numeric(train_sites["motif_score"], errors="coerce"),
            np.log1p(pd.to_numeric(train_sites["accessibility"], errors="coerce").clip(lower=0)),
        ]
    )
    validation_covariates = np.column_stack(
        [
            pd.to_numeric(validation_sites["motif_score"], errors="coerce"),
            np.log1p(pd.to_numeric(validation_sites["accessibility"], errors="coerce").clip(lower=0)),
        ]
    )
    train_x = np.column_stack([train_covariates, train_scores])
    validation_x = np.column_stack([validation_covariates, validation_scores])
    labels = train_sites["chip_label"].to_numpy(dtype=int)
    validation_labels = validation_sites["chip_label"].to_numpy(dtype=int)
    candidates = []
    for c_value in (0.03, 0.1, 0.3, 1.0, 3.0):
        for l1_ratio in (0.0, 0.5, 1.0):
            model = Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    (
                        "classifier",
                        LogisticRegression(
                            solver="saga",
                            C=c_value,
                            l1_ratio=l1_ratio,
                            class_weight="balanced",
                            max_iter=10000,
                            tol=1e-3,
                            random_state=seed,
                        ),
                    ),
                ]
            )
            model.fit(train_x, labels)
            probabilities = model.predict_proba(validation_x)[:, 1]
            metrics = binary_metrics(validation_labels, probabilities)
            prevalence = float(metrics["prevalence"])
            adjusted_auprc = (float(metrics["auprc"]) - prevalence) / max(1.0 - prevalence, 1e-6)
            selection = float(metrics["auroc"]) + adjusted_auprc
            candidates.append((selection, c_value, l1_ratio, model, probabilities, metrics))
    selected = max(candidates, key=lambda item: (item[0], float(item[5]["auprc"])))
    return selected[4], {
        "C": selected[1],
        "l1_ratio": selected[2],
        "selection_score": selected[0],
        **selected[5],
    }, fpca


def supervised_baseline(
    train_sites: pd.DataFrame,
    validation_sites: pd.DataFrame,
    seed: int,
) -> np.ndarray:
    columns = ["motif_score", "accessibility"]
    model = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    C=1.0,
                    class_weight="balanced",
                    max_iter=2000,
                    random_state=seed,
                ),
            ),
        ]
    )
    model.fit(train_sites[columns], train_sites["chip_label"].to_numpy(dtype=int))
    return model.predict_proba(validation_sites[columns])[:, 1]


def _evaluation_profiles(
    sites: pd.DataFrame,
    tracks: pd.DataFrame,
    cell: str,
    correction: str,
    cache_dir: Path,
    cache_label: str,
    flank: int,
    genome: Path | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    selected = tracks[(tracks["cell"] == cell) & (tracks["model"] == correction)]
    paths = {str(row.track): Path(row.signal) for row in selected.itertuples(index=False)}
    if "raw" not in paths:
        raise ValueError(f"tracks are missing {cell}/{correction}: raw")
    raw = load_or_extract_profiles(
        sites,
        paths["raw"],
        cache_dir / f"{cache_label}.{cell}.{correction}.raw.flank{flank}.npz",
        flank,
    )
    if "expected" in paths:
        expected = load_or_extract_profiles(
            sites,
            paths["expected"],
            cache_dir / f"{cache_label}.{cell}.{correction}.expected.flank{flank}.npz",
            flank,
        )
    elif "parametric_model" in paths:
        if genome is None:
            raise ValueError("--genome is required when a parametric_model track is supplied")
        expected = derive_parametric_expected_profiles(
            sites,
            raw,
            paths["parametric_model"],
            genome,
            cache_dir / f"{cache_label}.{cell}.{correction}.parametric_expected.flank{flank}.npz",
            flank,
        )
    else:
        raise ValueError(f"tracks are missing {cell}/{correction}: expected or parametric_model")
    strands = sites["TFBS_strand"].astype(str).to_numpy()
    return orient_profiles(raw, strands), orient_profiles(expected, strands)


def _label_free_models(
    train_observed: np.ndarray,
    train_expected: np.ndarray,
    train_sites: pd.DataFrame,
    positions: np.ndarray,
    dispersion: float,
    models: tuple[str, ...],
    seed: int,
    prior_profiles: dict[str, np.ndarray] | None = None,
) -> dict[str, object]:
    fitted: dict[str, object] = {}
    residual = deviance_profiles(train_observed, train_expected, dispersion)
    for name in models:
        if name in {"spline", "gp"}:
            model = BiasAwareFunctionalMixture(
                positions,
                smoother=name,
                dispersion=dispersion,
                max_iter=75,
                tolerance=1e-5,
                shrinkage=50.0,
            )
            model.fit(
                train_observed,
                train_expected,
                motif_score=train_sites["motif_score"].to_numpy(dtype=float),
                accessibility=train_observed.sum(axis=1),
                prior_profile=(prior_profiles or {}).get(name),
            )
        elif name == "fda":
            model = FdaMixtureModel(max_components=20, seed=seed).fit(
                residual,
                positions=positions,
                sample_weight=np.sqrt(np.maximum(train_observed.sum(axis=1), 1.0)),
            )
        elif name == "hybrid":
            model = HybridFdaGpModel(positions, max_components=20, seed=seed).fit(
                residual,
                sample_weight=np.sqrt(np.maximum(train_observed.sum(axis=1), 1.0)),
            )
        else:
            raise ValueError(f"unknown label-free model: {name}")
        fitted[name] = model
    return fitted


def _predict_label_free(
    model: object,
    name: str,
    observed: np.ndarray,
    expected: np.ndarray,
    sites: pd.DataFrame,
    dispersion: float,
) -> np.ndarray:
    if name in {"spline", "gp"}:
        assert isinstance(model, BiasAwareFunctionalMixture)
        return model.predict(
            observed,
            expected,
            motif_score=sites["motif_score"].to_numpy(dtype=float),
            accessibility=observed.sum(axis=1),
        )
    residual = deviance_profiles(observed, expected, dispersion)
    return model.predict_proba(residual)


def classify_failures(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (split, cell, tf), group in metrics.groupby(["split", "cell", "tf"], sort=True):
        baseline_rows = group[group["method"] == "supervised_baseline"]
        ceiling_rows = group[group["method"] == "supervised_fpca"]
        label_free = group[group["method"].isin(LABEL_FREE_MODELS)]
        if baseline_rows.empty or ceiling_rows.empty or label_free.empty:
            continue
        baseline = baseline_rows.sort_values("auprc", ascending=False).iloc[0]
        ceiling = ceiling_rows.sort_values("auprc", ascending=False).iloc[0]
        winner = label_free.sort_values(["auprc", "auroc"], ascending=False).iloc[0]
        denominator = max(abs(float(baseline.auprc)), 1e-8)
        ceiling_gain = (float(ceiling.auprc) - float(baseline.auprc)) / denominator
        label_free_gain = (float(winner.auprc) - float(baseline.auprc)) / denominator
        prevalence = max(float(winner.prevalence), 1e-8)
        label_free_over_chance = (float(winner.auprc) - prevalence) / max(1.0 - prevalence, 1e-8)
        if int(winner.positive_sites) < 100:
            status = "power_limited"
        elif float(winner.auroc) >= 0.65 and label_free_over_chance >= 0.10:
            status = "detectable"
        elif ceiling_gain < 0.10 and float(ceiling.auroc) < 0.65:
            status = "assay_limited_or_motif_ambiguous"
        elif str(winner.correction) != "DWM":
            status = "bias_limited"
        else:
            residual_rows = group[group["method"].str.startswith("residual_")]
            best_residual = residual_rows.sort_values("auprc", ascending=False).iloc[0] if len(residual_rows) else None
            if best_residual is not None and float(best_residual.auprc) > float(winner.auprc):
                status = "correction_score_limited"
            else:
                status = "shape_model_limited"
        rows.append(
            {
                "split": split,
                "cell": cell,
                "tf": tf,
                "motif_family": winner.motif_family,
                "positive_sites": int(winner.positive_sites),
                "supervised_baseline_auprc": float(baseline.auprc),
                "supervised_ceiling_auprc": float(ceiling.auprc),
                "supervised_relative_gain": ceiling_gain,
                "best_label_free_method": winner.method,
                "best_label_free_correction": winner.correction,
                "best_label_free_auroc": float(winner.auroc),
                "best_label_free_auprc": float(winner.auprc),
                "label_free_relative_gain": label_free_gain,
                "label_free_chance_adjusted_auprc": label_free_over_chance,
                "classification": status,
            }
        )
    return pd.DataFrame(rows)


def evaluate(
    study: dict,
    development_sites: pd.DataFrame,
    test_sites: pd.DataFrame | None,
    unlabeled_sources: dict[str, Path],
    tracks: pd.DataFrame,
    outdir: Path,
    *,
    corrections: tuple[str, ...],
    models: tuple[str, ...],
    maximum_unlabeled_sites_per_tf: int,
    flank: int,
    minimum_evaluation_sites: int,
    seed: int,
    genome: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == "development"].copy()
    positions = np.arange(-flank, flank + 1, dtype=float)
    cache_dir = outdir / "profile_cache"
    model_dir = outdir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, object]] = []
    score_rows: list[pd.DataFrame] = []
    frozen_rows: list[dict[str, object]] = []

    evaluation_sets: list[tuple[str, pd.DataFrame]] = [("validation", development_sites)]
    if test_sites is not None:
        evaluation_sets.append(("test", test_sites))

    for cell in sorted(tasks["cell"].unique()):
        if cell not in unlabeled_sources:
            raise ValueError(f"no unlabeled motif-site source was supplied for {cell}")
        unlabeled = build_unlabeled_training_sites(
            unlabeled_sources[cell],
            cell,
            tasks,
            study,
            maximum_per_tf=maximum_unlabeled_sites_per_tf,
            seed=seed,
        )
        unlabeled.to_csv(outdir / f"{cell}.unlabeled_training_sites.tsv.gz", sep="\t", index=False)
        cell_development = development_sites[development_sites["cell"] == cell].reset_index(drop=True)
        cell_test = (
            test_sites[test_sites["cell"] == cell].reset_index(drop=True)
            if test_sites is not None
            else None
        )

        for correction in corrections:
            train_observed, train_expected = _evaluation_profiles(
                unlabeled,
                tracks,
                cell,
                correction,
                cache_dir,
                "unlabeled",
                flank,
                genome,
            )
            dispersion = estimate_nb_dispersion(train_observed, train_expected)
            development_observed, development_expected = _evaluation_profiles(
                cell_development,
                tracks,
                cell,
                correction,
                cache_dir,
                "development",
                flank,
                genome,
            )
            cell_development["accessibility"] = development_observed.sum(axis=1)
            test_profiles = (
                _evaluation_profiles(
                    cell_test,
                    tracks,
                    cell,
                    correction,
                    cache_dir,
                    "test",
                    flank,
                    genome,
                )
                if cell_test is not None and len(cell_test)
                else None
            )
            if cell_test is not None and test_profiles is not None:
                cell_test["accessibility"] = test_profiles[0].sum(axis=1)

            global_prior: dict[str, np.ndarray] = {}
            global_limit = min(len(unlabeled), maximum_unlabeled_sites_per_tf * 2)
            global_indexes = np.random.default_rng(stable_seed(cell, correction, "global", seed=seed)).choice(
                len(unlabeled), size=global_limit, replace=False
            )
            for smoother in set(models).intersection({"spline", "gp"}):
                global_model = BiasAwareFunctionalMixture(
                    positions,
                    smoother=smoother,
                    dispersion=dispersion,
                    max_iter=60,
                    shrinkage=0.0,
                )
                global_result = global_model.fit(
                    train_observed[global_indexes],
                    train_expected[global_indexes],
                    motif_score=unlabeled.iloc[global_indexes]["motif_score"].to_numpy(dtype=float),
                    accessibility=train_observed[global_indexes].sum(axis=1),
                )
                global_prior[smoother] = global_result.footprint_profile

            family_priors: dict[tuple[str, str], np.ndarray] = {}
            for family, family_sites in unlabeled.groupby("motif_family", sort=True):
                indexes = family_sites.index.to_numpy(dtype=int)
                for smoother in set(models).intersection({"spline", "gp"}):
                    family_model = BiasAwareFunctionalMixture(
                        positions,
                        smoother=smoother,
                        dispersion=dispersion,
                        max_iter=60,
                        shrinkage=50.0,
                    )
                    result = family_model.fit(
                        train_observed[indexes],
                        train_expected[indexes],
                        motif_score=unlabeled.iloc[indexes]["motif_score"].to_numpy(dtype=float),
                        accessibility=train_observed[indexes].sum(axis=1),
                        prior_profile=global_prior.get(smoother),
                    )
                    family_priors[(str(family), smoother)] = result.footprint_profile

            for task in tasks[tasks["cell"] == cell].itertuples(index=False):
                tf = str(task.tf)
                family = str(task.motif_family)
                train_indexes = np.flatnonzero(unlabeled["tf"].to_numpy() == tf)
                if len(train_indexes) < 100:
                    continue
                start_time = perf_counter()
                prior_profiles = {
                    name: family_priors.get((family, name), global_prior.get(name))
                    for name in set(models).intersection({"spline", "gp"})
                }
                fitted = _label_free_models(
                    train_observed[train_indexes],
                    train_expected[train_indexes],
                    unlabeled.iloc[train_indexes],
                    positions,
                    dispersion,
                    models,
                    stable_seed(cell, tf, correction, seed=seed),
                    prior_profiles=prior_profiles,
                )
                fit_seconds = perf_counter() - start_time
                for name, model in fitted.items():
                    if isinstance(model, BiasAwareFunctionalMixture):
                        model.save(
                            model_dir / f"{cell}.{tf}.{correction}.{name}.npz",
                            metadata={
                                "cell": cell,
                                "tf": tf,
                                "motif_family": family,
                                "correction": correction,
                                "training_labels_used": False,
                            },
                        )
                    frozen_rows.append(
                        {
                            "cell": cell,
                            "tf": tf,
                            "motif_family": family,
                            "correction": correction,
                            "method": name,
                            "training_sites": int(len(train_indexes)),
                            "dispersion": dispersion,
                            "fit_seconds_all_models": fit_seconds,
                            "training_labels_used": False,
                            "status": "frozen_before_test",
                        }
                    )

                for split, all_sites in evaluation_sets:
                    cell_sites = cell_development if split == "validation" else cell_test
                    if cell_sites is None:
                        continue
                    split_mask = cell_sites["chromosome_split"].to_numpy() == split
                    tf_mask = cell_sites["tf"].astype(str).to_numpy() == tf
                    indexes = np.flatnonzero(split_mask & tf_mask)
                    if len(indexes) < minimum_evaluation_sites:
                        continue
                    if split == "validation":
                        observed, expected = development_observed[indexes], development_expected[indexes]
                    else:
                        assert test_profiles is not None
                        observed, expected = test_profiles[0][indexes], test_profiles[1][indexes]
                    selected_sites = cell_sites.iloc[indexes].reset_index(drop=True)
                    labels = selected_sites["chip_label"].to_numpy(dtype=int)

                    for name, model in fitted.items():
                        probabilities = _predict_label_free(
                            model, name, observed, expected, selected_sites, dispersion
                        )
                        metrics = binary_metrics(labels, probabilities)
                        metric_rows.append(
                            {
                                "split": split,
                                "cell": cell,
                                "tf": tf,
                                "motif_family": family,
                                "correction": correction,
                                "method": name,
                                "training_labels_used": False,
                                "training_sites": int(len(train_indexes)),
                                "fit_seconds_all_models": fit_seconds,
                                "dispersion": dispersion,
                                **metrics,
                            }
                        )
                        score_rows.append(
                            pd.DataFrame(
                                {
                                    "split": split,
                                    "cell": cell,
                                    "tf": tf,
                                    "correction": correction,
                                    "method": name,
                                    "site_index": indexes,
                                    "chip_label": labels,
                                    "score": probabilities,
                                }
                            )
                        )

                    for residual_mode in RESIDUAL_MODES:
                        _residual, score = residual_score(
                            observed, expected, residual_mode, dispersion
                        )
                        probabilities = _score_to_probability(score)
                        metric_rows.append(
                            {
                                "split": split,
                                "cell": cell,
                                "tf": tf,
                                "motif_family": family,
                                "correction": correction,
                                "method": f"residual_{residual_mode}",
                                "training_labels_used": False,
                                "training_sites": int(len(train_indexes)),
                                "fit_seconds_all_models": 0.0,
                                "dispersion": dispersion,
                                **binary_metrics(labels, probabilities),
                            }
                        )

            # The supervised ceiling is fitted once per TF/correction on matched
            # development sites and never used to train a deployable model.
            for task in tasks[tasks["cell"] == cell].itertuples(index=False):
                tf = str(task.tf)
                tf_mask = cell_development["tf"].astype(str).to_numpy() == tf
                train_indexes = np.flatnonzero(
                    tf_mask & (cell_development["chromosome_split"].to_numpy() == "train")
                )
                validation_indexes = np.flatnonzero(
                    tf_mask & (cell_development["chromosome_split"].to_numpy() == "validation")
                )
                if min(len(train_indexes), len(validation_indexes)) < minimum_evaluation_sites:
                    continue
                train_residual = deviance_profiles(
                    development_observed[train_indexes], development_expected[train_indexes], dispersion
                )
                validation_residual = deviance_profiles(
                    development_observed[validation_indexes], development_expected[validation_indexes], dispersion
                )
                train_sites = cell_development.iloc[train_indexes].reset_index(drop=True)
                validation_sites = cell_development.iloc[validation_indexes].reset_index(drop=True)
                ceiling, selected, _fpca = fit_supervised_ceiling(
                    train_residual,
                    validation_residual,
                    train_sites,
                    validation_sites,
                    seed=stable_seed(cell, tf, correction, "ceiling", seed=seed),
                )
                baseline = supervised_baseline(
                    train_sites,
                    validation_sites,
                    stable_seed(cell, tf, correction, "baseline", seed=seed),
                )
                for name, probabilities in (
                    ("supervised_baseline", baseline),
                    ("supervised_fpca", ceiling),
                ):
                    metric_rows.append(
                        {
                            "split": "validation",
                            "cell": cell,
                            "tf": tf,
                            "motif_family": str(task.motif_family),
                            "correction": correction,
                            "method": name,
                            "training_labels_used": True,
                            "training_sites": int(len(train_indexes)),
                            "fit_seconds_all_models": np.nan,
                            "dispersion": dispersion,
                            **binary_metrics(
                                validation_sites["chip_label"].to_numpy(dtype=int), probabilities
                            ),
                        }
                    )

    metrics = pd.DataFrame(metric_rows)
    scores = pd.concat(score_rows, ignore_index=True) if score_rows else pd.DataFrame()
    frozen = pd.DataFrame(frozen_rows)
    return metrics, scores, frozen


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"),
    )
    parser.add_argument("--development-sites", type=Path, required=True)
    parser.add_argument("--test-sites", type=Path)
    parser.add_argument("--unlock-test", action="store_true")
    parser.add_argument(
        "--unlabeled-sites",
        action="append",
        type=parse_cell_path,
        required=True,
        metavar="CELL=PATH",
    )
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--genome", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--corrections", nargs="+", default=["DWM", "PWM"])
    parser.add_argument("--models", nargs="+", choices=LABEL_FREE_MODELS, default=list(LABEL_FREE_MODELS))
    parser.add_argument("--maximum-unlabeled-sites-per-tf", type=int, default=10000)
    parser.add_argument("--minimum-evaluation-sites", type=int, default=100)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    if args.test_sites is not None and not args.unlock_test:
        raise SystemExit("test evaluation requires --unlock-test after the model manifest is frozen")
    if args.maximum_unlabeled_sites_per_tf < 100:
        raise SystemExit("--maximum-unlabeled-sites-per-tf must be at least 100")

    study = json.loads(args.study.read_text(encoding="utf-8"))
    development_sites = validate_sites(pd.read_csv(args.development_sites, sep="\t"), args.development_sites)
    development_sites["chromosome_split"] = development_sites["TFBS_chr"].map(
        lambda value: chromosome_split(str(value), study)
    )
    test_sites = (
        validate_sites(pd.read_csv(args.test_sites, sep="\t"), args.test_sites)
        if args.test_sites is not None
        else None
    )
    if test_sites is not None:
        test_sites["chromosome_split"] = test_sites["TFBS_chr"].map(
            lambda value: chromosome_split(str(value), study)
        )
    tracks = pd.read_csv(args.tracks, sep="\t")
    required_tracks = {"cell", "model", "track", "signal"}
    if not required_tracks.issubset(tracks.columns):
        raise SystemExit("tracks TSV must contain cell, model, track, and signal")
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics, scores, frozen = evaluate(
        study,
        development_sites,
        test_sites,
        dict(args.unlabeled_sites),
        tracks,
        args.outdir,
        corrections=tuple(args.corrections),
        models=tuple(args.models),
        maximum_unlabeled_sites_per_tf=args.maximum_unlabeled_sites_per_tf,
        flank=args.flank,
        minimum_evaluation_sites=args.minimum_evaluation_sites,
        seed=args.seed,
        genome=args.genome,
    )
    metrics.to_csv(args.outdir / "functional_metrics.tsv", sep="\t", index=False)
    scores.to_csv(args.outdir / "functional_site_scores.tsv.gz", sep="\t", index=False)
    frozen.to_csv(args.outdir / "frozen_functional_models.tsv", sep="\t", index=False)
    classification = classify_failures(metrics)
    classification.to_csv(args.outdir / "functional_failure_classification.tsv", sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-functional-benchmark-v1",
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "development_sites": str(args.development_sites),
        "development_sites_sha256": file_sha256(args.development_sites),
        "test_sites": str(args.test_sites) if args.test_sites else None,
        "test_sites_sha256": file_sha256(args.test_sites) if args.test_sites else None,
        "test_unlocked": bool(args.test_sites),
        "unlabeled_sources": {
            cell: {"path": str(path), "sha256": file_sha256(path)}
            for cell, path in args.unlabeled_sites
        },
        "tracks": str(args.tracks),
        "tracks_sha256": file_sha256(args.tracks),
        "genome": str(args.genome) if args.genome else None,
        "genome_sha256": file_sha256(args.genome) if args.genome else None,
        "corrections": args.corrections,
        "models": args.models,
        "maximum_unlabeled_sites_per_tf": args.maximum_unlabeled_sites_per_tf,
        "minimum_evaluation_sites": args.minimum_evaluation_sites,
        "flank": args.flank,
        "seed": args.seed,
        "metrics_rows": int(len(metrics)),
        "score_rows": int(len(scores)),
    }
    (args.outdir / "functional_benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(metrics.to_string(index=False))
    print("\nFailure classification")
    print(classification.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
