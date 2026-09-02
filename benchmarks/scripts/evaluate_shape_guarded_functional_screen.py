#!/usr/bin/env python3
"""Compare validation-frozen functional policies on raw-matched test support.

This is a secondary development diagnostic.  It separates discrimination from
footprint morphology by requiring central protection before a numerical gain
can be called a shape-guarded footprint gain.  It cannot produce a promotion
decision because the development test labels have already been opened.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_parametric_factorization import block_bootstrap_delta  # noqa: E402
from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from freeze_label_free_functional_models import immutable_write_json  # noqa: E402
from summarize_frozen_parametric_tf_evidence import load_output  # noqa: E402


SCHEMA = "fp-tools-shape-guarded-functional-screen-v1"
TEST_SCHEMA = "fp-tools-frozen-functional-test-results-v1"
NAKED_SCHEMA = "fp-tools-frozen-functional-naked-dna-v1"
RAW_SCHEMA = "fp-tools-frozen-bias-shrinkage-test-v1"
SITE_KEYS = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand"]


def parse_named_path(value: str) -> tuple[str, Path]:
    fields = value.split("=", 1)
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError("value must use NAME=PATH")
    return fields[0], Path(fields[1])


def unique_named_paths(
    values: list[tuple[str, Path]], label: str
) -> dict[str, Path]:
    output: dict[str, Path] = {}
    for name, path in values:
        if name in output:
            raise ValueError(f"duplicate {label}: {name}")
        output[name] = path
    return output


def canonical_id(document: dict) -> str:
    encoded = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return sha256(encoded.encode()).hexdigest()


def task_shapes(profiles: pd.DataFrame) -> pd.DataFrame:
    selected = profiles[
        profiles["method"].astype(str).str.startswith("frozen_")
    ].copy()
    rows = []
    for (cell, tf), group in selected.groupby(["cell", "tf"], sort=True):
        methods = group["method"].astype(str).unique()
        if len(methods) != 1:
            raise ValueError(f"multiple frozen profile methods for {cell}/{tf}")
        descriptors = {}
        for name in (
            "depletion",
            "width",
            "shoulder_distance",
            "asymmetry",
            "periodicity",
        ):
            values = group[name].dropna().astype(float).unique()
            if len(values) > 1:
                raise ValueError(f"nonconstant {name} for {cell}/{tf}")
            descriptors[name] = values[0] if len(values) else np.nan
        center = group[group["position"].abs() <= 5]
        center_contrast = float(center["positive_minus_negative"].mean())
        significant_points = int((center["upper_95"] < 0).sum())
        protected = bool(
            descriptors["depletion"] > 0
            and center_contrast < 0
            and significant_points > 0
        )
        strong = bool(
            descriptors["depletion"] >= 0.02
            and center_contrast < 0
            and significant_points >= 2
        )
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "method": methods[0],
                **{f"shape_{key}": value for key, value in descriptors.items()},
                "shape_center_label_contrast": center_contrast,
                "shape_significant_center_points": significant_points,
                "shape_has_central_protection": protected,
                "shape_has_strong_canonical_protection": strong,
            }
        )
    return pd.DataFrame(rows)


def common_support(
    candidate: pd.DataFrame,
    raw: pd.DataFrame,
    *,
    minimum_coverage: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_candidate = set(SITE_KEYS + ["label", "candidate_probability", "dwm_score"])
    required_raw = set(SITE_KEYS + ["chip_label", "raw_score", "dwm_score"])
    if missing := required_candidate.difference(candidate.columns):
        raise ValueError("candidate scores lack: " + ", ".join(sorted(missing)))
    if missing := required_raw.difference(raw.columns):
        raise ValueError("raw scores lack: " + ", ".join(sorted(missing)))
    if candidate.duplicated(SITE_KEYS).any() or raw.duplicated(SITE_KEYS).any():
        raise ValueError("score tables contain duplicate motif sites")
    # Newer frozen-policy score tables retain their own raw diagnostic.  The
    # independently frozen shrinkage table is the raw guardrail for this audit,
    # so remove any candidate-side copy before merging to keep its column name
    # and provenance unambiguous.
    candidate_for_merge = candidate.drop(
        columns=[column for column in ("raw_score", "chip_label") if column in candidate]
    )
    merged = candidate_for_merge.merge(
        raw[SITE_KEYS + ["chip_label", "raw_score", "dwm_score"]],
        on=SITE_KEYS,
        how="inner",
        suffixes=("_candidate", "_raw"),
        validate="one_to_one",
    )
    if not np.array_equal(
        merged["label"].to_numpy(dtype=int),
        merged["chip_label"].to_numpy(dtype=int),
    ):
        raise ValueError("candidate and raw labels differ on common support")
    finite = np.isfinite(
        merged[
            ["candidate_probability", "raw_score", "dwm_score_candidate"]
        ].to_numpy(dtype=float)
    ).all(axis=1)
    merged = merged.loc[finite].copy()
    coverage_rows = []
    for (cell, tf), source in candidate.groupby(["cell", "tf"], sort=True):
        kept = merged[merged["cell"].eq(cell) & merged["tf"].eq(tf)]
        coverage = len(kept) / max(len(source), 1)
        if coverage < minimum_coverage:
            raise ValueError(
                f"raw common support below {minimum_coverage:.1%} for {cell}/{tf}: "
                f"{coverage:.1%}"
            )
        coverage_rows.append(
            {
                "cell": cell,
                "tf": tf,
                "candidate_sites": int(len(source)),
                "common_finite_sites": int(len(kept)),
                "common_support_fraction": float(coverage),
            }
        )
    return merged, pd.DataFrame(coverage_rows)


def _metric(labels: np.ndarray, scores: np.ndarray) -> tuple[float, float]:
    return (
        float(roc_auc_score(labels, scores)),
        float(average_precision_score(labels, scores)),
    )


def task_metrics(
    scores: pd.DataFrame,
    shapes: pd.DataFrame,
    safety: pd.DataFrame,
    *,
    family: str,
    minimum_sites_per_class: int,
    bootstrap_iterations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    shape_index = shapes.set_index(["cell", "tf"])
    safety_index = safety.set_index(["cell", "tf"])
    rows, bootstrap_rows = [], []
    for (cell, tf), group in scores.groupby(["cell", "tf"], sort=True):
        labels = group["chip_label"].to_numpy(dtype=int)
        candidate = group["candidate_probability"].to_numpy(dtype=float)
        raw = group["raw_score"].to_numpy(dtype=float)
        dwm = group["dwm_score_candidate"].to_numpy(dtype=float)
        positive = int((labels == 1).sum())
        negative = int((labels == 0).sum())
        candidate_auc, candidate_ap = _metric(labels, candidate)
        raw_auc, raw_ap = _metric(labels, raw)
        dwm_auc, dwm_ap = _metric(labels, dwm)
        shape = shape_index.loc[(cell, tf)].to_dict()
        method = str(shape.pop("method"))
        naked = (
            safety_index.loc[(cell, tf)].to_dict()
            if (cell, tf) in safety_index.index
            else {}
        )
        row = {
            "family": family,
            "cell": cell,
            "tf": tf,
            "motif_family": str(group["motif_family"].iloc[0]),
            "candidate_id": str(group["candidate_id"].iloc[0]),
            "method": method,
            "status": (
                "eligible"
                if min(positive, negative) >= minimum_sites_per_class
                else "underpowered"
            ),
            "n_sites": int(len(group)),
            "n_positive": positive,
            "n_negative": negative,
            "auroc": candidate_auc,
            "auprc": candidate_ap,
            "raw_auroc": raw_auc,
            "raw_auprc": raw_ap,
            "dwm_auroc": dwm_auc,
            "dwm_auprc": dwm_ap,
            "auroc_gain_over_raw": candidate_auc - raw_auc,
            "relative_auprc_gain_over_raw": (candidate_ap - raw_ap)
            / max(raw_ap, 1e-8),
            "auroc_gain_over_dwm": candidate_auc - dwm_auc,
            "relative_auprc_gain_over_dwm": (candidate_ap - dwm_ap)
            / max(dwm_ap, 1e-8),
            **shape,
            **{f"naked_{key}": value for key, value in naked.items()},
        }
        for baseline_name, baseline in (("raw", raw), ("dwm", dwm)):
            bootstrap = block_bootstrap_delta(
                group,
                candidate,
                baseline,
                iterations=bootstrap_iterations,
                seed=(
                    int.from_bytes(
                        sha256(f"{family}|{cell}|{tf}|{baseline_name}|{seed}".encode()).digest()[:4],
                        "little",
                    )
                ),
            )
            bootstrap_rows.append(
                {
                    "family": family,
                    "cell": cell,
                    "tf": tf,
                    "baseline": baseline_name,
                    **bootstrap,
                }
            )
            for key, value in bootstrap.items():
                row[f"{baseline_name}_{key}"] = value
        row["classification"] = classify(row)
        rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(bootstrap_rows)


def classify(row: dict | pd.Series) -> str:
    if str(row.get("status")) != "eligible":
        return "underpowered"
    if row.get("naked_passes_safety") is not True:
        return "naked_dna_safety_failed_or_missing"
    raw_gain = (
        float(row.get("auroc_gain_over_raw", -np.inf)) > 0
        and float(row.get("relative_auprc_gain_over_raw", -np.inf)) > 0
    )
    raw_significant = (
        float(row.get("raw_auroc_gain_lower_95", -np.inf)) > 0
        and float(row.get("raw_relative_auprc_gain_lower_95", -np.inf)) > 0
    )
    protected = bool(row.get("shape_has_central_protection", False))
    strong = bool(row.get("shape_has_strong_canonical_protection", False))
    if raw_gain and raw_significant and strong:
        return "strong_shape_guarded_gain"
    if raw_gain and raw_significant and protected:
        return "weak_shape_guarded_gain"
    if raw_gain and raw_significant:
        return "occupancy_gain_without_protection"
    if (
        float(row.get("auroc_gain_over_dwm", -np.inf)) > 0
        and float(row.get("relative_auprc_gain_over_dwm", -np.inf)) > 0
    ):
        return "dwm_recovery_or_uncertain_raw_gain"
    return "no_joint_discrimination_gain"


def paired_safety(rates: pd.DataFrame) -> pd.DataFrame:
    selected = rates[rates["method"].astype(str).eq("paired_safety")].copy()
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError("duplicate paired naked-DNA safety row")
    keep = [
        "cell",
        "tf",
        "candidate_false_positive_rate",
        "candidate_wilson_upper_95",
        "candidate_minus_dwm",
        "passes_safety",
    ]
    return selected[[column for column in keep if column in selected]]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-manifest", type=Path, required=True)
    parser.add_argument(
        "--candidate", type=parse_named_path, action="append", required=True
    )
    parser.add_argument("--naked", type=parse_named_path, action="append", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimum-sites-per-class", type=int, default=200)
    parser.add_argument("--minimum-common-support", type=float, default=0.95)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    candidates = unique_named_paths(args.candidate, "candidate family")
    naked_paths = unique_named_paths(args.naked, "naked-DNA family")
    if set(candidates) != set(naked_paths):
        raise ValueError("candidate and naked-DNA family names differ")
    raw_doc, raw, raw_path = load_output(
        args.raw_manifest, expected_schema=RAW_SCHEMA, output="site_scores"
    )
    all_metrics, all_bootstrap, all_coverage = [], [], []
    inputs = [args.raw_manifest, raw_path]
    policies = {}
    for family in sorted(candidates):
        test_doc, scores, scores_path = load_output(
            candidates[family], expected_schema=TEST_SCHEMA, output="site_scores"
        )
        _test_doc, profiles, profiles_path = load_output(
            candidates[family], expected_schema=TEST_SCHEMA, output="profiles"
        )
        naked_doc, rates, rates_path = load_output(
            naked_paths[family], expected_schema=NAKED_SCHEMA, output="rates"
        )
        if test_doc.get("models_refitted_on_test") is not False:
            raise ValueError(f"{family} does not certify frozen test models")
        if naked_doc.get("models_refitted_on_naked_dna") is not False:
            raise ValueError(f"{family} refit models on naked DNA")
        if test_doc.get("policy_id") != naked_doc.get("policy_id"):
            raise ValueError(f"{family} test and naked-DNA policies differ")
        policies[family] = str(test_doc["policy_id"])
        merged, coverage = common_support(
            scores, raw, minimum_coverage=args.minimum_common_support
        )
        shapes = task_shapes(profiles)
        metrics, bootstrap = task_metrics(
            merged,
            shapes,
            paired_safety(rates),
            family=family,
            minimum_sites_per_class=args.minimum_sites_per_class,
            bootstrap_iterations=args.bootstrap_iterations,
            seed=args.seed,
        )
        coverage.insert(0, "family", family)
        all_metrics.append(metrics)
        all_bootstrap.append(bootstrap)
        all_coverage.append(coverage)
        inputs.extend(
            [
                candidates[family],
                scores_path,
                profiles_path,
                naked_paths[family],
                rates_path,
            ]
        )

    metrics = pd.concat(all_metrics, ignore_index=True)
    bootstrap = pd.concat(all_bootstrap, ignore_index=True)
    coverage = pd.concat(all_coverage, ignore_index=True)
    summary = (
        metrics.groupby(["family", "classification"], sort=True)
        .size()
        .rename("tasks")
        .reset_index()
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "metrics": args.outdir / "shape_guarded_metrics.tsv",
        "bootstrap": args.outdir / "shape_guarded_bootstrap.tsv",
        "coverage": args.outdir / "shape_guarded_common_support.tsv",
        "summary": args.outdir / "shape_guarded_summary.tsv",
    }
    for frame, name in (
        (metrics, "metrics"),
        (bootstrap, "bootstrap"),
        (coverage, "coverage"),
        (summary, "summary"),
    ):
        frame.to_csv(paths[name], sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "development_test_labels_previously_opened": True,
        "secondary_diagnostic_only": True,
        "promotion_decision": False,
        "raw_signal_guardrail": True,
        "shape_guard": {
            "center": "absolute positions <=5 bp",
            "central_protection": (
                "depletion > 0, mean positive-minus-negative center contrast < 0, "
                "and at least one pointwise upper 95% band < 0"
            ),
            "strong_canonical_protection": (
                "depletion >= 0.02, mean center contrast < 0, and at least two "
                "pointwise upper 95% bands < 0"
            ),
        },
        "policies": policies,
        "raw_manifest_schema": raw_doc["schema"],
        "inputs": [
            {"path": str(path), "sha256": file_sha256(path)} for path in inputs
        ],
        "outputs": {
            name: {"path": str(path), "sha256": file_sha256(path)}
            for name, path in paths.items()
        },
    }
    manifest["screen_id"] = canonical_id(manifest)
    immutable_write_json(args.outdir / "shape_guarded_manifest.json", manifest)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
