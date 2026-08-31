#!/usr/bin/env python3
"""Render blinded per-TF aggregates for DWM and label-free functional routes."""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import norm

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_strand_functional_templates import load_artifact  # noqa: E402
from fp_tools.tools.functional_footprints import normalize_functional_profiles  # noqa: E402


TASK_KEYS = ("cell", "tf", "motif_family")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*values: object, seed: int = 2026) -> int:
    digest = blake2b(digest_size=8)
    for value in (seed, *values):
        digest.update(str(value).encode())
        digest.update(b"\0")
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def normal_mean_band(
    profiles: np.ndarray, *, confidence: float = 0.95
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(profiles, dtype=float)
    if values.ndim != 2 or len(values) == 0:
        raise ValueError("profiles must be a nonempty site-by-position matrix")
    mean = np.mean(values, axis=0)
    if len(values) == 1:
        return mean, np.full_like(mean, np.nan), np.full_like(mean, np.nan)
    critical = float(norm.ppf(0.5 + confidence / 2.0))
    standard_error = np.std(values, axis=0, ddof=1) / np.sqrt(len(values))
    return mean, mean - critical * standard_error, mean + critical * standard_error


def normal_difference_band(
    positive: np.ndarray,
    negative: np.ndarray,
    *,
    confidence: float = 0.95,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    first = np.asarray(positive, dtype=float)
    second = np.asarray(negative, dtype=float)
    if first.ndim != 2 or second.ndim != 2 or first.shape[1] != second.shape[1]:
        raise ValueError("positive and negative profiles must share positions")
    if len(first) == 0 or len(second) == 0:
        raise ValueError("both groups must contain sites")
    difference = np.mean(first, axis=0) - np.mean(second, axis=0)
    if min(len(first), len(second)) < 2:
        return difference, np.full_like(difference, np.nan), np.full_like(difference, np.nan)
    critical = float(norm.ppf(0.5 + confidence / 2.0))
    variance = np.var(first, axis=0, ddof=1) / len(first) + np.var(
        second, axis=0, ddof=1
    ) / len(second)
    standard_error = np.sqrt(np.maximum(variance, 0.0))
    return (
        difference,
        difference - critical * standard_error,
        difference + critical * standard_error,
    )


def conservative_aggregate_difference(
    positive: pd.DataFrame, negative: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Difference of stored means with a conservative independent-CI envelope."""

    first = positive.sort_values("position")
    second = negative.sort_values("position")
    if not np.array_equal(first["position"].to_numpy(), second["position"].to_numpy()):
        raise ValueError("stored aggregate positions do not agree")
    difference = first["mean"].to_numpy() - second["mean"].to_numpy()
    lower = first["lower_95"].to_numpy() - second["upper_95"].to_numpy()
    upper = first["upper_95"].to_numpy() - second["lower_95"].to_numpy()
    return difference, lower, upper


def curve_summary(
    positions: np.ndarray,
    difference: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
) -> dict[str, float]:
    x = np.asarray(positions, dtype=float)
    delta = np.asarray(difference, dtype=float)
    center = np.abs(x) <= 10
    shoulders = (np.abs(x) >= 20) & (np.abs(x) <= 50)
    core = np.abs(x) <= 50
    excludes_zero = ((np.asarray(lower) > 0) | (np.asarray(upper) < 0)) & np.isfinite(
        lower
    ) & np.isfinite(upper)
    return {
        "visual_rms_difference": float(np.sqrt(np.mean(np.square(delta[core])))),
        "center_mean_difference": float(np.mean(delta[center])),
        "shoulder_mean_difference": float(np.mean(delta[shoulders])),
        "protection_shape_contrast": float(
            np.mean(delta[shoulders]) - np.mean(delta[center])
        ),
        "band_exclusion_fraction_within_50bp": float(np.mean(excludes_zero[core])),
    }


def development_best(metrics: pd.DataFrame) -> pd.DataFrame:
    required = {
        *TASK_KEYS,
        "bias_configuration",
        "candidate_id",
        "family",
        "channel",
        "status",
        "converged",
        "selection_score",
        "auroc",
        "auprc",
    }
    missing = sorted(required.difference(metrics.columns))
    if missing:
        raise ValueError(f"evaluation metrics lack columns: {missing}")
    passing = metrics[metrics["status"].eq("ok") & metrics["converged"].eq(True)].copy()
    output = (
        passing.sort_values(
            list(TASK_KEYS)
            + ["selection_score", "auprc", "auroc", "bias_configuration", "candidate_id"],
            ascending=[True] * len(TASK_KEYS) + [False, False, False, True, True],
            kind="mergesort",
        )
        .groupby(list(TASK_KEYS), as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )
    output["display_method"] = "Development best (diagnostic)"
    output["selection_labels_used"] = True
    return output


def unlabeled_routes(rule_rows: pd.DataFrame, rule: str) -> pd.DataFrame:
    output = rule_rows[rule_rows["selection_rule"].eq(rule)].copy()
    if output.empty:
        raise ValueError(f"selection-rule table has no rows for {rule}")
    if output.duplicated(list(TASK_KEYS)).any():
        raise ValueError("unlabeled rule has duplicate task routes")
    output["family"] = output.get("family_evaluation", output.get("family_cv"))
    output["display_method"] = f"Unlabeled: {rule}"
    if output["selection_labels_used"].astype(bool).any():
        raise ValueError("unlabeled route unexpectedly records label use")
    return output


def fixed_dwm_routes(
    candidate_matrix: pd.DataFrame, fixed_candidate_id: str
) -> pd.DataFrame:
    output = candidate_matrix[
        candidate_matrix["bias_configuration"].eq("DWM")
        & candidate_matrix["candidate_id"].eq(fixed_candidate_id)
    ].copy()
    if output.duplicated(list(TASK_KEYS)).any():
        raise ValueError("DWM route has duplicate tasks")
    output["display_method"] = "Current DWM"
    output["selection_labels_used"] = False
    output["channel"] = "normalized_deviance"
    return output


def artifact_paths(manifest: dict) -> dict[tuple[str, str], Path]:
    paths: dict[tuple[str, str], Path] = {}
    for row in manifest["artifacts"]:
        if row.get("purpose") != "evaluation":
            continue
        key = (str(row["model"]), str(row["cell"]))
        if key in paths:
            raise ValueError(f"duplicate evaluation artifact {key}")
        paths[key] = Path(row["path"])
    return paths


def extract_route_profiles(
    sites: pd.DataFrame,
    profiles: dict[str, np.ndarray],
    *,
    tf: str,
    channel: str,
    positions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if channel not in profiles:
        raise ValueError(f"artifact lacks requested profile channel {channel}")
    mask = sites["tf"].astype(str).eq(tf) & sites["chromosome_split"].eq("validation")
    labels = sites.loc[mask, "chip_label"].to_numpy(dtype=int)
    values = normalize_functional_profiles(profiles[channel][mask.to_numpy()], positions)
    if len(values) != len(labels):
        raise AssertionError("profile and label counts differ")
    return values, labels


def _smooth(values: np.ndarray, sigma: float) -> np.ndarray:
    if sigma < 0:
        raise ValueError("plot smoothing must be non-negative")
    if sigma == 0:
        return np.asarray(values, dtype=float).copy()
    return gaussian_filter1d(np.asarray(values, dtype=float), sigma=sigma, axis=-1, mode="nearest")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label-free-run", type=Path, required=True)
    parser.add_argument("--selection-rule-rows", type=Path, required=True)
    parser.add_argument("--selection-rule", default="likelihood_plus_depletion")
    parser.add_argument("--candidate-matrix", type=Path, required=True)
    parser.add_argument("--policy-manifest", type=Path, required=True)
    parser.add_argument("--dwm-aggregate-profiles", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--key-out", type=Path)
    parser.add_argument("--summary-out", type=Path)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--plot-sigma-bp", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    if not 0.5 < args.confidence < 1.0:
        raise SystemExit("--confidence must be between 0.5 and 1")
    if args.plot_sigma_bp < 0:
        raise SystemExit("--plot-sigma-bp must be non-negative")

    manifest_path = args.label_free_run / "strand_label_free_manifest.json"
    label_free_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if label_free_manifest.get("locked_test_labels_read"):
        raise SystemExit("renderer refuses a run that opened locked test labels")
    study_path = Path(label_free_manifest["study"])
    study = json.loads(study_path.read_text(encoding="utf-8"))
    flank = int(study["profile_flank_bp"])
    positions = np.arange(-flank, flank + 1, dtype=float)
    metrics_path = args.label_free_run / "strand_label_free_metrics.tsv.gz"
    metrics = pd.read_csv(metrics_path, sep="\t")
    rule_rows = pd.read_csv(args.selection_rule_rows, sep="\t")
    candidate_matrix = pd.read_csv(args.candidate_matrix, sep="\t")
    policy = json.loads(args.policy_manifest.read_text(encoding="utf-8"))
    fixed_dwm_candidate = policy["global_dwm_fallback"]["candidate_id"]
    dwm_aggregates = pd.read_csv(args.dwm_aggregate_profiles, sep="\t")
    dwm_aggregates = dwm_aggregates[
        dwm_aggregates["correction"].eq("DWM")
        & dwm_aggregates["signal"].eq("normalized_deviance")
        & dwm_aggregates["split"].eq("validation")
    ].copy()

    routes = pd.concat(
        [
            fixed_dwm_routes(candidate_matrix, fixed_dwm_candidate),
            unlabeled_routes(rule_rows, args.selection_rule),
            development_best(metrics),
        ],
        ignore_index=True,
        sort=False,
    )
    paths = artifact_paths(label_free_manifest)
    artifact_cache: dict[tuple[str, str], tuple[pd.DataFrame, dict[str, np.ndarray]]] = {}

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    args.out.parent.mkdir(parents=True, exist_ok=True)
    key_out = args.key_out or args.out.with_name(args.out.stem + "_blinding_key.tsv")
    summary_out = args.summary_out or args.out.with_name(args.out.stem + "_summary.tsv")
    method_order = [
        "Current DWM",
        f"Unlabeled: {args.selection_rule}",
        "Development best (diagnostic)",
    ]
    tasks = routes[list(TASK_KEYS)].drop_duplicates().sort_values(list(TASK_KEYS))
    key_rows: list[dict] = []
    summary_rows: list[dict] = []
    with PdfPages(args.out) as pdf:
        for task in tasks.itertuples(index=False):
            cell, tf, motif_family = str(task.cell), str(task.tf), str(task.motif_family)
            task_routes = routes[
                routes["cell"].astype(str).eq(cell)
                & routes["tf"].astype(str).eq(tf)
                & routes["motif_family"].astype(str).eq(motif_family)
            ]
            figure, axes = plt.subplots(
                2,
                len(method_order),
                figsize=(15.6, 6.7),
                sharex="col",
                sharey="row",
                squeeze=False,
                gridspec_kw={"height_ratios": [2.2, 1.0]},
            )
            swap = bool(stable_seed(cell, tf, "blind", seed=args.seed) % 2)
            blinded = {
                0: "Group B" if swap else "Group A",
                1: "Group A" if swap else "Group B",
            }
            key_rows.extend(
                {
                    "cell": cell,
                    "tf": tf,
                    "blinded_group": blinded[label],
                    "actual_group": "chip_positive" if label else "matched_negative",
                }
                for label in (0, 1)
            )
            for column, method in enumerate(method_order):
                axis = axes[0, column]
                difference_axis = axes[1, column]
                route = task_routes[task_routes["display_method"].eq(method)]
                if route.empty:
                    axis.text(0.5, 0.5, "Insufficient sites", transform=axis.transAxes, ha="center")
                    difference_axis.set_visible(False)
                    continue
                winner = route.iloc[0]
                if method == "Current DWM":
                    stored = dwm_aggregates[
                        dwm_aggregates["cell"].astype(str).eq(cell)
                        & dwm_aggregates["tf"].astype(str).eq(tf)
                    ]
                    positive = stored[stored["group"].eq("chip_positive")].sort_values("position")
                    negative = stored[stored["group"].eq("matched_negative")].sort_values("position")
                    if positive.empty or negative.empty:
                        axis.set_visible(False)
                        difference_axis.set_visible(False)
                        continue
                    group_curves = {
                        1: (
                            positive["mean"].to_numpy(),
                            positive["lower_95"].to_numpy(),
                            positive["upper_95"].to_numpy(),
                            int(positive["sites"].iloc[0]),
                        ),
                        0: (
                            negative["mean"].to_numpy(),
                            negative["lower_95"].to_numpy(),
                            negative["upper_95"].to_numpy(),
                            int(negative["sites"].iloc[0]),
                        ),
                    }
                    difference, difference_lower, difference_upper = conservative_aggregate_difference(
                        positive, negative
                    )
                    detail = str(winner["candidate_id"])
                else:
                    model = str(winner["bias_configuration"])
                    key = (model, cell)
                    if key not in artifact_cache:
                        if key not in paths:
                            raise SystemExit(f"missing evaluation artifact for {model}/{cell}")
                        sites, profile_map, _document = load_artifact(paths[key], cell, study)
                        artifact_cache[key] = (sites, profile_map)
                    sites, profile_map = artifact_cache[key]
                    values, labels = extract_route_profiles(
                        sites,
                        profile_map,
                        tf=tf,
                        channel=str(winner["channel"]),
                        positions=positions,
                    )
                    if min(np.sum(labels == 0), np.sum(labels == 1)) < 2:
                        axis.text(0.5, 0.5, "Insufficient labeled sites", transform=axis.transAxes, ha="center")
                        difference_axis.set_visible(False)
                        continue
                    group_curves = {}
                    for label in (0, 1):
                        mean, lower, upper = normal_mean_band(
                            values[labels == label], confidence=args.confidence
                        )
                        group_curves[label] = (
                            mean,
                            lower,
                            upper,
                            int(np.sum(labels == label)),
                        )
                    difference, difference_lower, difference_upper = normal_difference_band(
                        values[labels == 1], values[labels == 0], confidence=args.confidence
                    )
                    detail = f"{winner['bias_configuration']}; {winner['candidate_id']}"

                colors = {0: "#2A6FBB", 1: "#C23B33"}
                for label in (0, 1):
                    mean, lower, upper, _sites = group_curves[label]
                    mean = _smooth(mean, args.plot_sigma_bp)
                    lower = _smooth(lower, args.plot_sigma_bp)
                    upper = _smooth(upper, args.plot_sigma_bp)
                    axis.plot(positions, mean, color=colors[label], linewidth=1.7, label=blinded[label])
                    axis.fill_between(positions, lower, upper, color=colors[label], alpha=0.18, linewidth=0)
                plotted_difference = _smooth(difference, args.plot_sigma_bp)
                plotted_lower = _smooth(difference_lower, args.plot_sigma_bp)
                plotted_upper = _smooth(difference_upper, args.plot_sigma_bp)
                if swap:
                    plotted_difference = -plotted_difference
                    plotted_lower, plotted_upper = -plotted_upper, -plotted_lower
                difference_axis.fill_between(
                    positions, plotted_lower, plotted_upper, color="#6A3D9A", alpha=0.18, linewidth=0
                )
                difference_axis.plot(positions, plotted_difference, color="#6A3D9A", linewidth=1.7)
                summary_rows.append(
                    {
                        "cell": cell,
                        "tf": tf,
                        "motif_family": motif_family,
                        "method": method,
                        "bias_configuration": winner.get("bias_configuration", "DWM"),
                        "candidate_id": winner["candidate_id"],
                        "channel": winner.get("channel", "normalized_deviance"),
                        "selection_labels_used": bool(winner["selection_labels_used"]),
                        "auroc": float(winner["auroc"]),
                        "auprc": float(winner["auprc"]),
                        "positive_sites": int(group_curves[1][3]),
                        "negative_sites": int(group_curves[0][3]),
                        **curve_summary(
                            positions, difference, difference_lower, difference_upper
                        ),
                    }
                )
                axis.axvline(0, color="#666666", linewidth=0.7, linestyle="--")
                axis.axhline(0, color="#999999", linewidth=0.6)
                difference_axis.axvspan(-15, 15, color="#777777", alpha=0.06, linewidth=0)
                difference_axis.axvline(0, color="#666666", linewidth=0.7, linestyle="--")
                difference_axis.axhline(0, color="#999999", linewidth=0.6)
                axis.set_title(
                    f"{method}\nAUROC {float(winner['auroc']):.3f}; AUPRC {float(winner['auprc']):.3f}\n{detail}",
                    fontsize=8,
                )
                difference_axis.set_xlabel("Position from motif center (bp)")
            axes[0, 0].set_ylabel("Normalized functional residual")
            axes[1, 0].set_ylabel("Group A − Group B")
            handles, labels_text = axes[0, 0].get_legend_handles_labels()
            if handles:
                figure.legend(
                    handles,
                    labels_text,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 0.923),
                    ncol=2,
                    frameon=False,
                )
            figure.suptitle(f"{cell} — {tf}: blinded aggregate audit", y=0.992)
            figure.subplots_adjust(
                top=0.79, bottom=0.09, left=0.06, right=0.99, hspace=0.08, wspace=0.04
            )
            pdf.savefig(figure)
            plt.close(figure)

    key_frame = pd.DataFrame(key_rows).drop_duplicates()
    summary = pd.DataFrame(summary_rows)
    key_frame.to_csv(key_out, sep="\t", index=False)
    summary.to_csv(summary_out, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-label-free-per-tf-aggregate-audit-v1",
        "locked_holdout_labels_read": False,
        "blinded_groups": True,
        "development_best_is_diagnostic_only": True,
        "unlabeled_selection_rule": args.selection_rule,
        "confidence_method_new_routes": "site-wise normal mean/difference interval",
        "confidence_method_dwm": "stored group bootstrap intervals; conservative difference envelope",
        "display_only_gaussian_sigma_bp": args.plot_sigma_bp,
        "inputs": {
            str(path): file_sha256(path)
            for path in (
                manifest_path,
                metrics_path,
                args.selection_rule_rows,
                args.candidate_matrix,
                args.policy_manifest,
                args.dwm_aggregate_profiles,
            )
        },
        "outputs": {
            "pdf": {"path": str(args.out), "sha256": file_sha256(args.out)},
            "blinding_key": {"path": str(key_out), "sha256": file_sha256(key_out)},
            "summary": {"path": str(summary_out), "sha256": file_sha256(summary_out)},
        },
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
