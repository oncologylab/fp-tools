#!/usr/bin/env python3
"""Evaluate a covariate-balanced sensitivity subset without outcome tuning.

The locked evaluator deliberately fails a task when propensity matching leaves
large standardized differences.  This diagnostic applies a fixed grid of
coarsened exact matches to the already-scored cohort.  It selects the largest
subset passing the balance and power thresholds using covariates only; model
scores are not consulted until after that subset is fixed.

This is a post-unblinding sensitivity analysis and never changes the original
locked result or its promotion status.
"""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path
import sys
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import binary_metrics  # noqa: E402
from evaluate_locked_holdout_policy import (  # noqa: E402
    aggregate_rows,
    bootstrap_delta,
    standardized_difference,
)
from fp_tools.tools.functional_footprints import normalize_functional_profiles  # noqa: E402
from pool_functional_profile_artifacts import load_artifact  # noqa: E402


MATCH_FEATURES = (
    "motif_score",
    "log_accessibility",
    "gc_fraction",
    "peak_position_signed",
    "peak_position_abs",
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _stable_key(row: Any, seed: int, namespace: str) -> int:
    digest = blake2b(digest_size=8)
    for value in (
        seed,
        namespace,
        row.TFBS_chr,
        int(row.TFBS_start),
        int(row.TFBS_end),
        row.TFBS_strand,
        int(row.label),
    ):
        digest.update(str(value).encode())
        digest.update(b"\0")
    return int.from_bytes(digest.digest(), "little")


def _quantile_bin(values: pd.Series, bins: int) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.nunique(dropna=True) <= 1:
        return pd.Series(np.zeros(len(values), dtype=int), index=values.index)
    return pd.qcut(numeric, bins, labels=False, duplicates="drop").astype("Int64")


def coarsened_exact_match(
    frame: pd.DataFrame,
    *,
    accessibility_bins: int,
    gc_bins: int,
    motif_bins: int = 3,
    peak_position_bins: int = 3,
    seed: int = 2026,
) -> pd.DataFrame:
    """Return a deterministic 1:1 match using pooled, outcome-blind bins."""
    required = {
        "label",
        "TFBS_chr",
        "TFBS_start",
        "TFBS_end",
        "TFBS_strand",
        *MATCH_FEATURES,
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError("score table is missing columns: " + ", ".join(missing))
    working = frame[frame["label"].isin([0, 1])].dropna(
        subset=list(MATCH_FEATURES)
    ).copy()
    working["_motif_bin"] = _quantile_bin(working["motif_score"], motif_bins)
    working["_access_bin"] = _quantile_bin(
        working["log_accessibility"], accessibility_bins
    )
    working["_gc_bin"] = _quantile_bin(working["gc_fraction"], gc_bins)
    working["_position_bin"] = _quantile_bin(
        working["peak_position_abs"], peak_position_bins
    )
    working["_position_sign"] = (
        working["peak_position_signed"].to_numpy(dtype=float) >= 0
    ).astype(int)
    strata = (
        "_motif_bin",
        "_access_bin",
        "_gc_bin",
        "_position_bin",
        "_position_sign",
    )
    selected: list[int] = []
    for key, group in working.groupby(list(strata), sort=True, observed=True):
        key_text = "|".join(str(value) for value in key)
        positive = group[group["label"].eq(1)].copy()
        negative = group[group["label"].eq(0)].copy()
        count = min(len(positive), len(negative))
        if not count:
            continue
        for label, subset in ((1, positive), (0, negative)):
            subset["_stable_key"] = [
                _stable_key(row, seed, f"{key_text}|{label}")
                for row in subset.itertuples(index=False)
            ]
            selected.extend(
                subset.sort_values("_stable_key", kind="mergesort").index[:count]
            )
    output = working.loc[selected].copy()
    return output.drop(columns=[*strata]).sort_values(
        ["TFBS_chr", "TFBS_start", "label"], kind="mergesort"
    ).reset_index(drop=True)


def matching_diagnostics(frame: pd.DataFrame) -> dict[str, float]:
    values = {
        feature: standardized_difference(frame, feature) for feature in MATCH_FEATURES
    }
    return {
        **{f"smd_{feature}": float(value) for feature, value in values.items()},
        "maximum_absolute_smd": float(max(abs(value) for value in values.values())),
    }


def evaluate_grid(
    frame: pd.DataFrame,
    *,
    accessibility_bins: Sequence[int],
    gc_bins: Sequence[int],
    motif_bins: int,
    peak_position_bins: int,
    minimum_pairs: int,
    maximum_smd: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    subsets: dict[tuple[int, int], pd.DataFrame] = {}
    for access_count in accessibility_bins:
        for gc_count in gc_bins:
            matched = coarsened_exact_match(
                frame,
                accessibility_bins=int(access_count),
                gc_bins=int(gc_count),
                motif_bins=motif_bins,
                peak_position_bins=peak_position_bins,
                seed=seed,
            )
            pairs = int(matched["label"].value_counts().min()) if len(matched) else 0
            diagnostics = matching_diagnostics(matched) if pairs else {
                **{f"smd_{feature}": np.nan for feature in MATCH_FEATURES},
                "maximum_absolute_smd": np.nan,
            }
            row: dict[str, Any] = {
                "accessibility_bins": int(access_count),
                "gc_bins": int(gc_count),
                "motif_bins": int(motif_bins),
                "peak_position_bins": int(peak_position_bins),
                "matched_per_class": pairs,
                **diagnostics,
            }
            row["balance_pass"] = bool(
                pairs >= minimum_pairs
                and np.isfinite(row["maximum_absolute_smd"])
                and row["maximum_absolute_smd"] <= maximum_smd
            )
            if pairs:
                labels = matched["label"].to_numpy(dtype=int)
                candidate = binary_metrics(
                    labels, matched["candidate_probability"].to_numpy(dtype=float)
                )
                reference = binary_metrics(
                    labels, matched["reference_probability"].to_numpy(dtype=float)
                )
                row.update(
                    {
                        "candidate_auroc": float(candidate["auroc"]),
                        "reference_auroc": float(reference["auroc"]),
                        "auroc_gain": float(candidate["auroc"] - reference["auroc"]),
                        "candidate_auprc": float(candidate["auprc"]),
                        "reference_auprc": float(reference["auprc"]),
                        "relative_auprc_gain": float(
                            (candidate["auprc"] - reference["auprc"])
                            / max(reference["auprc"], 1e-8)
                        ),
                    }
                )
            rows.append(row)
            subsets[(int(access_count), int(gc_count))] = matched
    grid = pd.DataFrame(rows)
    eligible = grid[grid["balance_pass"]].sort_values(
        ["matched_per_class", "accessibility_bins", "gc_bins"],
        ascending=[False, True, True],
        kind="mergesort",
    )
    if eligible.empty:
        raise ValueError("no common-support grid passes the balance and power thresholds")
    selected_row = eligible.iloc[0]
    selected = subsets[
        (int(selected_row["accessibility_bins"]), int(selected_row["gc_bins"]))
    ]
    summary = selected_row.to_dict()
    summary["selection_rule"] = (
        "largest matched subset with maximum absolute SMD at or below the threshold; "
        "ties prefer fewer accessibility and GC bins"
    )
    summary["selection_uses_model_scores"] = False
    return grid, selected, summary


def profile_rows(
    selected: pd.DataFrame,
    candidate_artifact: Path,
    reference_artifact: Path,
    *,
    profile_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    candidate_doc, candidate_sites, candidate_arrays = load_artifact(candidate_artifact)
    reference_doc, reference_sites, reference_arrays = load_artifact(reference_artifact)
    if len(candidate_sites) != len(reference_sites):
        raise ValueError("candidate and reference artifacts have different site counts")
    indexes = selected["artifact_index"].to_numpy(dtype=int)
    positions = np.arange(
        -(candidate_arrays["shared_strand_residual"].shape[1] // 2),
        candidate_arrays["shared_strand_residual"].shape[1] // 2 + 1,
    )
    labels = selected["label"].to_numpy(dtype=int)
    candidate = normalize_functional_profiles(
        candidate_arrays["shared_strand_residual"][indexes], positions
    )
    reference = normalize_functional_profiles(
        reference_arrays["combined_residual"][indexes], positions
    )
    rows = aggregate_rows(
        cell=str(selected["cell"].iloc[0]),
        tf=str(selected["tf"].iloc[0]),
        method="candidate",
        profiles=candidate,
        labels=labels,
        positions=positions,
        iterations=profile_bootstrap,
        seed=seed,
    )
    rows.extend(
        aggregate_rows(
            cell=str(selected["cell"].iloc[0]),
            tf=str(selected["tf"].iloc[0]),
            method="reference",
            profiles=reference,
            labels=labels,
            positions=positions,
            iterations=profile_bootstrap,
            seed=seed,
        )
    )
    return pd.DataFrame(rows)


def render_pdf(summary: dict[str, Any], profiles: pd.DataFrame, output: Path) -> None:
    colors = {"chip_positive": "#C23B33", "matched_negative": "#2A6FBB"}
    figure, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), constrained_layout=True)
    for axis, method, title in (
        (axes[0], "candidate", "LOG81 + TF-aware anchored FDA"),
        (axes[1], "reference", "Conventional DWM (TOBIAS-style)"),
    ):
        for group, label in (
            ("chip_positive", "ChIP-supported"),
            ("matched_negative", "Matched negative"),
        ):
            curve = profiles[
                profiles["method"].eq(method) & profiles["group"].eq(group)
            ].sort_values("position")
            axis.plot(curve["position"], curve["mean"], color=colors[group], label=label)
            axis.fill_between(
                curve["position"].to_numpy(dtype=float),
                curve["lower_95"].to_numpy(dtype=float),
                curve["upper_95"].to_numpy(dtype=float),
                color=colors[group],
                alpha=0.18,
            )
        axis.axvline(0, color="#777777", lw=0.8, ls="--")
        axis.set_title(title, fontsize=10)
        axis.set_xlabel("Position relative to MAX motif (bp)")
        axis.set_ylabel("Normalized signed-deviance residual")
        axis.legend(frameon=False, fontsize=8)
    axes[2].axis("off")
    text = (
        "WTC11 MAX external sensitivity\n\n"
        f"Balanced sites/class: {int(summary['matched_per_class']):,}\n"
        f"Maximum SMD: {summary['maximum_absolute_smd']:.3f}\n\n"
        f"AUROC: {summary['reference_auroc']:.3f} → {summary['candidate_auroc']:.3f}\n"
        f"AUROC gain: {summary['auroc_gain']:+.3f}\n"
        f"AUPRC: {summary['reference_auprc']:.3f} → {summary['candidate_auprc']:.3f}\n"
        f"Relative AUPRC gain: {summary['relative_auprc_gain']:+.1%}\n\n"
        f"Bootstrap P(AUROC gain > 0):\n"
        f"{summary['auroc_gain_bootstrap_probability_positive']:.3f}\n\n"
        "Post-unblinding common-support sensitivity;\n"
        "candidate and detector were frozen beforehand."
    )
    axes[2].text(0.01, 0.98, text, ha="left", va="top", family="monospace", fontsize=8.5)
    figure.suptitle("MAX footprint separation: before versus after", fontsize=12)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output)
    plt.close(figure)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--candidate-artifact", type=Path)
    parser.add_argument("--reference-artifact", type=Path)
    parser.add_argument("--accessibility-bins", nargs="+", type=int, default=[3, 4, 5, 6, 8, 10])
    parser.add_argument("--gc-bins", nargs="+", type=int, default=[2, 3, 4, 5])
    parser.add_argument("--motif-bins", type=int, default=3)
    parser.add_argument("--peak-position-bins", type=int, default=3)
    parser.add_argument("--minimum-pairs", type=int, default=200)
    parser.add_argument("--maximum-smd", type=float, default=0.1)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--profile-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if bool(args.candidate_artifact) != bool(args.reference_artifact):
        raise SystemExit("provide both --candidate-artifact and --reference-artifact")
    frame = pd.read_csv(args.scores, sep="\t")
    grid, selected, summary = evaluate_grid(
        frame,
        accessibility_bins=args.accessibility_bins,
        gc_bins=args.gc_bins,
        motif_bins=args.motif_bins,
        peak_position_bins=args.peak_position_bins,
        minimum_pairs=args.minimum_pairs,
        maximum_smd=args.maximum_smd,
        seed=args.seed,
    )
    bootstrap = bootstrap_delta(
        selected,
        selected["candidate_probability"].to_numpy(dtype=float),
        selected["reference_probability"].to_numpy(dtype=float),
        iterations=args.bootstrap,
        seed=args.seed,
    )
    summary.update(bootstrap)
    summary.update(
        {
            "schema": "fp-tools-common-support-sensitivity-v1",
            "post_unblinding_sensitivity": True,
            "scores": str(args.scores),
            "scores_sha256": file_sha256(args.scores),
            "minimum_pairs": int(args.minimum_pairs),
            "maximum_smd_threshold": float(args.maximum_smd),
            "seed": int(args.seed),
        }
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    grid.to_csv(args.outdir / "common_support_grid.tsv", sep="\t", index=False)
    selected.to_csv(args.outdir / "common_support_selected_sites.tsv.gz", sep="\t", index=False)
    (args.outdir / "common_support_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if args.candidate_artifact:
        profiles = profile_rows(
            selected,
            args.candidate_artifact,
            args.reference_artifact,
            profile_bootstrap=args.profile_bootstrap,
            seed=args.seed,
        )
        profiles.to_csv(
            args.outdir / "common_support_aggregate_profiles.tsv.gz", sep="\t", index=False
        )
        render_pdf(summary, profiles, args.outdir / "common_support_before_after.pdf")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
