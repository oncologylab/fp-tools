#!/usr/bin/env python3
"""Evaluate frozen TF functional settings across depths and randomizations.

Hyperparameters are read from a 10M development winner table and are never
reselected at higher depth. The script reports shape-only template transfer,
the matched scalar residual baseline, depth stability, power classifications,
and per-TF aggregate curves with uncertainty.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import (  # noqa: E402
    _evaluation_profiles,
    _score_to_probability,
    binary_metrics,
    chromosome_split,
    residual_score,
    stable_seed,
    validate_sites,
)
from evaluate_functional_template_transfer import (  # noqa: E402
    balanced_training_indexes,
    selection_score,
    training_indexes,
)
from run_footprint_ablation_plan import output_is_complete  # noqa: E402
from fp_tools.tools.functional_footprints import (  # noqa: E402
    FunctionalTemplateDetector,
    normalize_functional_profiles,
    profile_descriptors,
    site_accessibility_background,
    standardized_functional_separation,
)
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


DEPTH_ORDER = {"10000000": 0, "25000000": 1, "50000000": 2, "full": 3}
DEPTH_LABEL = {"10000000": "10M", "25000000": "25M", "50000000": "50M", "full": "Full"}


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover_signal_matrix(
    plan: pd.DataFrame,
    *,
    depths: tuple[str, ...] = tuple(DEPTH_ORDER),
    randomization_seeds: tuple[int, ...] | None = None,
    correction: str = "fp_tools_dwm",
    require_complete: bool = True,
    completed_jobs: set[str] | None = None,
) -> pd.DataFrame:
    selected = plan[
        (plan["stage"].astype(str) == "correction")
        & (plan["correction"].astype(str) == correction)
        & (plan["depth"].astype(str).isin(depths))
    ].copy()
    if randomization_seeds is not None:
        selected = selected[
            pd.to_numeric(selected["seed"], errors="raise").astype(int).isin(randomization_seeds)
        ].copy()
    rows = []
    for row in selected.itertuples(index=False):
        corrected = Path(str(row.expected_output))
        prefix = corrected.name.removesuffix("_corrected.bw")
        raw = corrected.parent / f"{prefix}_uncorrected.bw"
        expected = corrected.parent / f"{prefix}_expected.bw"
        files_complete = all(output_is_complete(path) for path in (corrected, raw, expected))
        runner_complete = completed_jobs is None or str(row.job_id) in completed_jobs
        complete = files_complete and runner_complete
        if require_complete and not complete:
            raise FileNotFoundError(
                f"incomplete correction output for {row.job_id}: {corrected.parent}"
            )
        if complete:
            rows.append(
                {
                    "sample": str(row.sample),
                    "cell": str(row.cell),
                    "depth": str(row.depth),
                    "seed": int(row.seed),
                    "correction": correction,
                    "raw": str(raw),
                    "expected": str(expected),
                    "corrected": str(corrected),
                    "job_id": str(row.job_id),
                }
            )
    output = pd.DataFrame(rows)
    if len(output):
        output["depth_order"] = output["depth"].map(DEPTH_ORDER)
        output = output.sort_values(["depth_order", "seed", "cell"]).reset_index(drop=True)
    return output


def completed_correction_jobs(status: pd.DataFrame) -> set[str]:
    """Return correction jobs whose runner reached a terminal successful state."""

    required = {"job_id", "stage", "state"}
    missing = required.difference(status.columns)
    if missing:
        raise ValueError(
            "runner status lacks columns: " + ", ".join(sorted(missing))
        )
    latest = status.drop_duplicates("job_id", keep="last")
    successful = latest[
        latest["stage"].astype(str).eq("correction")
        & latest["state"].astype(str).isin({"completed", "skipped_existing"})
    ]
    return set(successful["job_id"].astype(str))


def prepare_profile(
    observed: np.ndarray,
    expected: np.ndarray,
    positions: np.ndarray,
    *,
    residual_mode: str,
    background: str,
    dispersion: float,
) -> tuple[np.ndarray, np.ndarray]:
    adjusted = site_accessibility_background(
        observed,
        expected,
        positions,
        method=background,
    )
    return residual_score(observed, adjusted, residual_mode, dispersion)


def _global_to_local(sites: pd.DataFrame, cell: str, indexes: np.ndarray) -> np.ndarray:
    cell_rows = np.flatnonzero(sites["cell"].astype(str).to_numpy() == cell)
    lookup = pd.Series(np.arange(len(cell_rows)), index=cell_rows)
    return lookup.loc[indexes].to_numpy(dtype=int)


def _training_profile_matrix(
    sites: pd.DataFrame,
    indexes: np.ndarray,
    profiles_by_cell: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    parts = []
    labels = []
    selected_sites = sites.iloc[indexes]
    for cell in sorted(selected_sites["cell"].astype(str).unique()):
        global_indexes = indexes[selected_sites["cell"].astype(str).to_numpy() == cell]
        parts.append(profiles_by_cell[cell][_global_to_local(sites, cell, global_indexes)])
        labels.append(sites.iloc[global_indexes]["chip_label"].to_numpy(dtype=int))
    return np.vstack(parts), np.concatenate(labels)


def _training_score_vector(
    sites: pd.DataFrame,
    indexes: np.ndarray,
    scores_by_cell: dict[str, np.ndarray],
) -> np.ndarray:
    parts = []
    selected_sites = sites.iloc[indexes]
    for cell in sorted(selected_sites["cell"].astype(str).unique()):
        global_indexes = indexes[
            selected_sites["cell"].astype(str).to_numpy() == cell
        ]
        parts.append(
            scores_by_cell[cell][_global_to_local(sites, cell, global_indexes)]
        )
    return np.concatenate(parts)


def calibrated_scalar_probability(
    training_scores: np.ndarray,
    training_labels: np.ndarray,
    evaluation_scores: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Calibrate scalar direction from the same training split as the template."""

    values = np.asarray(training_scores, dtype=float)
    labels = np.asarray(training_labels, dtype=int)
    evaluation = np.asarray(evaluation_scores, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() < 4 or np.unique(labels[finite]).size != 2:
        return np.full_like(evaluation, 0.5), 0.0
    positive = values[finite & (labels == 1)]
    negative = values[finite & (labels == 0)]
    direction = float(np.sign(np.mean(positive) - np.mean(negative)))
    if direction == 0.0:
        direction = 1.0
    oriented = direction * values[finite]
    location = float(np.median(oriented))
    scale = float(1.4826 * np.median(np.abs(oriented - location)))
    if scale <= 0:
        scale = float(np.std(oriented)) or 1.0
    probability = 1.0 / (
        1.0
        + np.exp(
            -np.clip((direction * evaluation - location) / scale, -30.0, 30.0)
        )
    )
    return probability, direction


def profile_curve_rows(
    profiles: np.ndarray,
    labels: np.ndarray,
    positions: np.ndarray,
    metadata: dict,
) -> pd.DataFrame:
    normalized = normalize_functional_profiles(profiles, positions)
    data = {"position": positions.astype(int)}
    for label, name in ((0, "negative"), (1, "positive")):
        group = normalized[labels == label]
        data[f"{name}_mean"] = np.mean(group, axis=0)
        data[f"{name}_se"] = np.std(group, axis=0, ddof=1) / np.sqrt(max(len(group), 1))
    data["difference"] = data["positive_mean"] - data["negative_mean"]
    return pd.DataFrame({**metadata, **data})


def evaluate_signal_pair(
    *,
    study: dict,
    sites: pd.DataFrame,
    signal_rows: pd.DataFrame,
    winner_rows: pd.DataFrame,
    cache_dir: Path,
    genome: Path | None,
    flank: int,
    minimum_sites: int,
    maximum_train_per_tf_class: int,
    seed: int,
    skip_invalid: bool = False,
) -> tuple[list[dict], list[pd.DataFrame], list[dict]]:
    depth = str(signal_rows["depth"].iloc[0])
    randomization = int(signal_rows["seed"].iloc[0])
    positions = np.arange(-flank, flank + 1, dtype=float)
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == "development"]
    observed_by_cell = {}
    expected_by_cell = {}
    dispersion_by_cell = {}
    for signal in signal_rows.itertuples(index=False):
        cell_sites = sites[sites["cell"].astype(str) == str(signal.cell)].reset_index(drop=True)
        tracks = pd.DataFrame(
            [
                {"cell": signal.cell, "model": "DWM", "track": "raw", "signal": signal.raw},
                {"cell": signal.cell, "model": "DWM", "track": "expected", "signal": signal.expected},
            ]
        )
        cache_label = f"depth_{depth}.seed_{randomization}"
        observed, expected = _evaluation_profiles(
            cell_sites,
            tracks,
            str(signal.cell),
            "DWM",
            cache_dir,
            cache_label,
            flank,
            genome,
        )
        observed_by_cell[str(signal.cell)] = observed
        expected_by_cell[str(signal.cell)] = expected
        dispersion_by_cell[str(signal.cell)] = estimate_nb_dispersion(observed, expected)

    unique_configs = winner_rows[["residual", "background"]].drop_duplicates()
    profile_cache = {}
    score_cache = {}
    for config in unique_configs.itertuples(index=False):
        for cell in observed_by_cell:
            residual, scalar = prepare_profile(
                observed_by_cell[cell],
                expected_by_cell[cell],
                positions,
                residual_mode=str(config.residual),
                background=str(config.background),
                dispersion=dispersion_by_cell[cell],
            )
            key = (cell, str(config.residual), str(config.background))
            profile_cache[key] = residual
            score_cache[key] = scalar

    metric_rows = []
    curve_frames = []
    failure_rows = []
    for winner in winner_rows.itertuples(index=False):
        target_cell = str(winner.cell)
        target_tf = str(winner.tf)
        family = str(winner.motif_family)
        scope = str(winner.training_scope)
        validation = np.flatnonzero(
            (sites["cell"].astype(str).to_numpy() == target_cell)
            & (sites["tf"].astype(str).to_numpy() == target_tf)
            & (sites["chromosome_split"].astype(str).to_numpy() == "validation")
        )
        train = training_indexes(
            sites,
            target_cell=target_cell,
            target_tf=target_tf,
            target_family=family,
            scope=scope,
        )
        train = balanced_training_indexes(
            sites,
            train,
            maximum_per_tf_class=maximum_train_per_tf_class,
            seed=stable_seed(
                target_cell,
                target_tf,
                scope,
                depth,
                randomization,
                seed=seed,
            ),
        )
        train_labels = sites.iloc[train]["chip_label"].to_numpy(dtype=int)
        labels = sites.iloc[validation]["chip_label"].to_numpy(dtype=int)
        if (
            min(np.sum(train_labels == 0), np.sum(train_labels == 1)) < minimum_sites
            or min(np.sum(labels == 0), np.sum(labels == 1)) < minimum_sites
        ):
            continue
        profiles_by_cell = {
            cell: profile_cache[(cell, str(winner.residual), str(winner.background))]
            for cell in observed_by_cell
        }
        train_profiles, train_labels = _training_profile_matrix(
            sites,
            train,
            profiles_by_cell,
        )
        validation_profiles = profiles_by_cell[target_cell][
            _global_to_local(sites, target_cell, validation)
        ]
        started = perf_counter()
        try:
            detector = FunctionalTemplateDetector(
                positions,
                smoother=str(winner.smoother),
                window_limit=float(winner.window_limit),
            ).fit(train_profiles, train_labels)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
            if not skip_invalid:
                raise
            failure_rows.append(
                {
                    "cell": target_cell,
                    "tf": target_tf,
                    "motif_family": family,
                    "depth": depth,
                    "depth_label": DEPTH_LABEL[depth],
                    "seed": randomization,
                    "training_scope": scope,
                    "residual": str(winner.residual),
                    "background": str(winner.background),
                    "smoother": str(winner.smoother),
                    "window_limit": float(winner.window_limit),
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
            continue
        probabilities = detector.predict_proba(validation_profiles)
        metrics = binary_metrics(labels, probabilities)
        descriptor = asdict(profile_descriptors(detector.footprint_template_, positions))
        common = {
            "cell": target_cell,
            "tf": target_tf,
            "motif_family": family,
            "depth": depth,
            "depth_label": DEPTH_LABEL[depth],
            "seed": randomization,
            "training_scope": scope,
            "residual": str(winner.residual),
            "background": str(winner.background),
            "smoother": str(winner.smoother),
            "window_limit": float(winner.window_limit),
            "configuration_frozen_at_10m": True,
            "training_labels_used": True,
            "motif_or_accessibility_features_used": False,
            "train_sites": int(len(train)),
            "validation_sites": int(len(validation)),
            "fit_seconds": perf_counter() - started,
            "functional_separation": standardized_functional_separation(
                validation_profiles,
                labels,
                positions,
            ),
            **descriptor,
        }
        metric_rows.append(
            {
                **common,
                "method": "functional_template",
                "selection_score": selection_score(metrics),
                **metrics,
            }
        )
        scalar = score_cache[(target_cell, str(winner.residual), str(winner.background))][
            _global_to_local(sites, target_cell, validation)
        ]
        scalar_probabilities = _score_to_probability(scalar)
        scalar_metrics = binary_metrics(labels, scalar_probabilities)
        metric_rows.append(
            {
                **common,
                "method": "scalar_center_flank",
                "selection_score": selection_score(scalar_metrics),
                **scalar_metrics,
            }
        )
        scores_by_cell = {
            cell: score_cache[(cell, str(winner.residual), str(winner.background))]
            for cell in observed_by_cell
        }
        train_scalar = _training_score_vector(sites, train, scores_by_cell)
        calibrated_probabilities, scalar_direction = calibrated_scalar_probability(
            train_scalar,
            train_labels,
            scalar,
        )
        calibrated_metrics = binary_metrics(labels, calibrated_probabilities)
        metric_rows.append(
            {
                **common,
                "method": "tf_calibrated_scalar",
                "scalar_direction": scalar_direction,
                "selection_score": selection_score(calibrated_metrics),
                **calibrated_metrics,
            }
        )
        curve_frames.append(
            profile_curve_rows(
                validation_profiles,
                labels,
                positions,
                {
                    "cell": target_cell,
                    "tf": target_tf,
                    "motif_family": family,
                    "depth": depth,
                    "depth_label": DEPTH_LABEL[depth],
                    "seed": randomization,
                    "training_scope": scope,
                    "residual": str(winner.residual),
                    "background": str(winner.background),
                },
            )
        )
    return metric_rows, curve_frames, failure_rows


def summarize_depth(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return metrics.copy()
    keys = [
        "cell",
        "tf",
        "motif_family",
        "training_scope",
        "method",
        "depth",
        "depth_label",
    ]
    summary = (
        metrics.groupby(keys, sort=True)
        .agg(
            seeds=("seed", "nunique"),
            auroc_mean=("auroc", "mean"),
            auroc_sd=("auroc", "std"),
            auprc_mean=("auprc", "mean"),
            auprc_sd=("auprc", "std"),
            brier_mean=("brier", "mean"),
            functional_separation_mean=("functional_separation", "mean"),
            functional_separation_sd=("functional_separation", "std"),
            runtime_mean=("fit_seconds", "mean"),
        )
        .reset_index()
    )
    summary["depth_order"] = summary["depth"].map(DEPTH_ORDER)
    return summary.sort_values(keys[:-2] + ["depth_order"]).reset_index(drop=True)


def classify_depth_limits(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    selected = summary[
        (summary["training_scope"] == "same_cell_ceiling")
        & (summary["method"] == "functional_template")
    ]
    for (cell, tf, family), group in selected.groupby(["cell", "tf", "motif_family"], sort=True):
        values = group.set_index("depth")
        low = values.loc["10000000"] if "10000000" in values.index else None
        high_name = "full" if "full" in values.index else ("50000000" if "50000000" in values.index else None)
        if low is None or high_name is None:
            continue
        high = values.loc[high_name]
        gain = float(high.auroc_mean - low.auroc_mean)
        if float(high.auroc_mean) >= 0.65:
            status = "shape_detectable_at_high_depth"
        elif gain >= 0.03:
            status = "power_limited"
        else:
            status = "assay_limited_or_motif_ambiguous"
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "motif_family": family,
                "low_depth": "10000000",
                "high_depth": high_name,
                "low_auroc": float(low.auroc_mean),
                "high_auroc": float(high.auroc_mean),
                "auroc_gain": gain,
                "low_auprc": float(low.auprc_mean),
                "high_auprc": float(high.auprc_mean),
                "relative_auprc_gain": (
                    float(high.auprc_mean - low.auprc_mean) / max(float(low.auprc_mean), 1e-8)
                ),
                "classification": status,
            }
        )
    return pd.DataFrame(rows)


def render_depth_plots(summary: pd.DataFrame, profiles: pd.DataFrame, output: Path) -> None:
    if summary.empty or profiles.empty:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output) as pdf:
        for (cell, tf), task in summary.groupby(["cell", "tf"], sort=True):
            figure, axes = plt.subplots(1, 3, figsize=(15, 4.2))
            same = task[task["training_scope"] == "same_cell_ceiling"]
            for method, color in (
                ("functional_template", "#C23B33"),
                ("tf_calibrated_scalar", "#E07A2D"),
                ("scalar_center_flank", "#2A6FBB"),
            ):
                line = same[same["method"] == method].sort_values("depth_order")
                if line.empty:
                    continue
                x = line["depth_order"].to_numpy(dtype=float)
                axes[0].errorbar(x, line["auroc_mean"], yerr=line["auroc_sd"].fillna(0), marker="o", color=color, label=method)
                axes[1].errorbar(x, line["auprc_mean"], yerr=line["auprc_sd"].fillna(0), marker="o", color=color, label=method)
            for axis, label in zip(axes[:2], ("AUROC", "AUPRC")):
                axis.set_xticks(range(4), ["10M", "25M", "50M", "Full"])
                axis.set_ylabel(label)
                axis.set_xlabel("Fragments")
                axis.grid(alpha=0.2)
            curves = profiles[
                (profiles["cell"] == cell)
                & (profiles["tf"] == tf)
                & (profiles["training_scope"] == "same_cell_ceiling")
            ]
            for depth, color in (("10000000", "#888888"), ("50000000", "#E07A2D"), ("full", "#C23B33")):
                curve = curves[curves["depth"] == depth]
                if curve.empty:
                    continue
                means = curve.groupby("position")["difference"].mean()
                axes[2].plot(means.index, means.values, color=color, label=DEPTH_LABEL[depth])
            axes[2].axvline(0, color="#666666", linewidth=0.7, linestyle="--")
            axes[2].axhline(0, color="#999999", linewidth=0.6)
            axes[2].set_xlabel("Position from motif center (bp)")
            axes[2].set_ylabel("Positive − matched-negative residual")
            handles, labels = axes[0].get_legend_handles_labels()
            if handles:
                axes[0].legend(handles, labels, frameon=False, fontsize=8)
            axes[2].legend(frameon=False, fontsize=8)
            figure.suptitle(f"{cell} — {tf}: frozen shape detector across depth")
            figure.tight_layout()
            pdf.savefig(figure, bbox_inches="tight")
            plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-run", type=Path, required=True)
    parser.add_argument("--signal-plan", type=Path, required=True)
    parser.add_argument("--winner-table", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--depths", nargs="+", default=list(DEPTH_ORDER))
    parser.add_argument(
        "--randomization-seeds",
        nargs="+",
        type=int,
        help="Optional depth/downsampling seeds to evaluate; defaults to every seed in the plan",
    )
    parser.add_argument("--minimum-sites-per-class", type=int, default=100)
    parser.add_argument("--maximum-train-per-tf-class", type=int, default=2000)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument(
        "--runner-status",
        type=Path,
        help="Optional ablation-runner status TSV used to exclude files still being written",
    )
    parser.add_argument("--seed", type=int)
    args = parser.parse_args(argv)
    unknown_depths = set(args.depths).difference(DEPTH_ORDER)
    if unknown_depths:
        raise SystemExit("unknown depths: " + ", ".join(sorted(unknown_depths)))
    base_manifest_path = args.base_run / "functional_benchmark_manifest.json"
    base_manifest = json.loads(base_manifest_path.read_text(encoding="utf-8"))
    if base_manifest.get("test_unlocked"):
        raise SystemExit("depth evaluation refuses a base run that opened test labels")
    study_path = Path(base_manifest["study"])
    sites_path = Path(base_manifest["development_sites"])
    genome = Path(base_manifest["genome"]) if base_manifest.get("genome") else None
    study = json.loads(study_path.read_text(encoding="utf-8"))
    sites = validate_sites(pd.read_csv(sites_path, sep="\t"), sites_path)
    sites["chromosome_split"] = sites["TFBS_chr"].map(
        lambda chromosome: chromosome_split(str(chromosome), study)
    )
    sites = sites[sites["chromosome_split"].isin(["train", "validation"])].reset_index(drop=True)
    signal_plan = pd.read_csv(args.signal_plan, sep="\t", keep_default_na=False)
    completed_jobs = None
    if args.runner_status is not None:
        completed_jobs = completed_correction_jobs(
            pd.read_csv(args.runner_status, sep="\t", keep_default_na=False)
        )
    signals = discover_signal_matrix(
        signal_plan,
        depths=tuple(args.depths),
        randomization_seeds=(
            None
            if args.randomization_seeds is None
            else tuple(args.randomization_seeds)
        ),
        require_complete=not args.allow_incomplete,
        completed_jobs=completed_jobs,
    )
    if signals.empty:
        raise SystemExit("no complete depth signals were found")
    winners = pd.read_csv(args.winner_table, sep="\t")
    winners = winners[winners["correction"].astype(str) == "DWM"].copy()
    seed = int(args.seed if args.seed is not None else base_manifest["seed"])
    flank = int(base_manifest["flank"])
    args.outdir.mkdir(parents=True, exist_ok=True)
    metric_rows = []
    curve_frames = []
    failure_rows = []
    for (depth, randomization), signal_rows in signals.groupby(["depth", "seed"], sort=False):
        required_cells = set(pd.DataFrame(study["tasks"]).query("split == 'development'")["cell"].astype(str))
        if set(signal_rows["cell"].astype(str)) != required_cells:
            if args.allow_incomplete:
                continue
            raise SystemExit(f"depth {depth}/seed {randomization} is missing a development cell")
        rows, curves, failures = evaluate_signal_pair(
            study=study,
            sites=sites,
            signal_rows=signal_rows,
            winner_rows=winners,
            cache_dir=args.outdir / "profile_cache",
            genome=genome,
            flank=flank,
            minimum_sites=args.minimum_sites_per_class,
            maximum_train_per_tf_class=args.maximum_train_per_tf_class,
            seed=seed,
            skip_invalid=args.allow_incomplete,
        )
        metric_rows.extend(rows)
        curve_frames.extend(curves)
        failure_rows.extend(failures)
        print(f"completed depth={depth} seed={randomization}: {len(rows):,} metric rows", flush=True)
    metrics = pd.DataFrame(metric_rows)
    profiles = pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame()
    summary = summarize_depth(metrics)
    classification = classify_depth_limits(summary)
    metrics_path = args.outdir / "functional_depth_metrics.tsv.gz"
    profiles_path = args.outdir / "functional_depth_profiles.tsv.gz"
    summary_path = args.outdir / "functional_depth_summary.tsv"
    classification_path = args.outdir / "functional_depth_classification.tsv"
    failures_path = args.outdir / "functional_depth_failures.tsv"
    plot_path = args.outdir / "functional_depth_panels.pdf"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    profiles.to_csv(profiles_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    classification.to_csv(classification_path, sep="\t", index=False)
    pd.DataFrame(failure_rows).to_csv(failures_path, sep="\t", index=False)
    render_depth_plots(summary, profiles, plot_path)
    manifest = {
        "schema": "fp-tools-functional-depth-matrix-v1",
        "locked_test_labels_read": False,
        "configuration_frozen_at_10m": True,
        "training_labels_used": True,
        "motif_or_accessibility_features_used": False,
        "base_manifest": str(base_manifest_path),
        "base_manifest_sha256": file_sha256(base_manifest_path),
        "signal_plan": str(args.signal_plan),
        "signal_plan_sha256": file_sha256(args.signal_plan),
        "runner_status": (
            None
            if args.runner_status is None
            else {
                "path": str(args.runner_status),
                "sha256": file_sha256(args.runner_status),
            }
        ),
        "winner_table": str(args.winner_table),
        "winner_table_sha256": file_sha256(args.winner_table),
        "depths": args.depths,
        "requested_randomization_seeds": args.randomization_seeds,
        "seeds": sorted(int(value) for value in signals["seed"].unique()),
        "allow_incomplete": args.allow_incomplete,
        "minimum_sites_per_class": args.minimum_sites_per_class,
        "maximum_train_per_tf_class": args.maximum_train_per_tf_class,
        "seed": seed,
        "metrics_rows": int(len(metrics)),
        "profile_rows": int(len(profiles)),
        "summary_rows": int(len(summary)),
        "classification_rows": int(len(classification)),
        "failure_rows": int(len(failure_rows)),
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {"path": str(profiles_path), "sha256": file_sha256(profiles_path)},
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "classification": {"path": str(classification_path), "sha256": file_sha256(classification_path)},
            "failures": {"path": str(failures_path), "sha256": file_sha256(failures_path)},
            "panels": {"path": str(plot_path), "sha256": file_sha256(plot_path)} if plot_path.is_file() else None,
        },
    }
    (args.outdir / "functional_depth_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(classification.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
