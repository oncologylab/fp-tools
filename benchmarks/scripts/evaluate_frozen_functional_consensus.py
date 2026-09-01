#!/usr/bin/env python3
"""Freeze and evaluate a calibrated count--FDA footprint consensus.

The score transform and primary operator are selected on development validation
chromosomes only.  Test chromosomes and naked-DNA controls are separate
subcommands so neither can influence the frozen consensus.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import binary_metrics  # noqa: E402
from evaluate_functional_template_transfer import selection_score  # noqa: E402
from evaluate_naked_dna_functional_policy import wilson_interval  # noqa: E402
from evaluate_parametric_factorization import block_bootstrap_delta  # noqa: E402
from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from freeze_functional_call_thresholds import (  # noqa: E402
    SCHEMA as THRESHOLD_SCHEMA,
    upper_tail_threshold,
)
from freeze_label_free_functional_models import immutable_write_json  # noqa: E402


FREEZE_SCHEMA = "fp-tools-frozen-functional-consensus-v1"
TEST_SCHEMA = "fp-tools-frozen-functional-consensus-test-v1"
NAKED_SCHEMA = "fp-tools-frozen-functional-consensus-naked-dna-v1"
FUNCTIONAL_TEST_SCHEMA = "fp-tools-frozen-functional-test-results-v1"
FUNCTIONAL_NAKED_SCHEMA = "fp-tools-frozen-functional-naked-dna-v1"

TASK_COLUMNS = ["cell", "tf"]
OPERATORS = (
    {"operator": "count_only", "kind": "weighted_mean", "count_weight": 1.0},
    {"operator": "rank_min", "kind": "minimum"},
    {
        "operator": "rank_mean_count_0p75",
        "kind": "weighted_mean",
        "count_weight": 0.75,
    },
    {"operator": "rank_geometric_mean", "kind": "geometric_mean"},
    {
        "operator": "rank_mean_count_0p50",
        "kind": "weighted_mean",
        "count_weight": 0.50,
    },
    {
        "operator": "rank_mean_count_0p25",
        "kind": "weighted_mean",
        "count_weight": 0.25,
    },
    {"operator": "fda_only", "kind": "weighted_mean", "count_weight": 0.0},
)
SIMPLICITY_ORDER = tuple(record["operator"] for record in OPERATORS)


def canonical_identifier(document: dict, field: str) -> str:
    content = dict(document)
    observed = str(content.pop(field))
    canonical = json.dumps(content, sort_keys=True, separators=(",", ":"))
    expected = sha256(canonical.encode()).hexdigest()
    if observed != expected:
        raise ValueError(f"{field} does not match the frozen document")
    return observed


def require_checksum(record: dict, purpose: str) -> Path:
    path = Path(record["path"])
    if file_sha256(path) != str(record["sha256"]):
        raise ValueError(f"{purpose} checksum mismatch: {path}")
    return path


def empirical_negative_cdf(
    scores: np.ndarray,
    sorted_negative_scores: np.ndarray,
) -> np.ndarray:
    """Map scores to a smoothed empirical matched-negative percentile."""

    values = np.asarray(scores, dtype=float)
    negative = np.asarray(sorted_negative_scores, dtype=float)
    if negative.ndim != 1 or not len(negative) or not np.isfinite(negative).all():
        raise ValueError("negative calibration scores must be a nonempty finite vector")
    if np.any(negative[1:] < negative[:-1]):
        raise ValueError("negative calibration scores must be sorted")
    output = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    ranks = np.searchsorted(negative, values[finite], side="right")
    output[finite] = (ranks + 0.5) / (len(negative) + 1.0)
    return output


def consensus_score(
    count_rank: np.ndarray,
    fda_rank: np.ndarray,
    operator: dict,
) -> np.ndarray:
    count = np.asarray(count_rank, dtype=float)
    fda = np.asarray(fda_rank, dtype=float)
    if count.shape != fda.shape:
        raise ValueError("count and FDA ranks have different shapes")
    kind = str(operator["kind"])
    if kind == "minimum":
        return np.minimum(count, fda)
    if kind == "geometric_mean":
        return np.sqrt(np.clip(count, 0.0, 1.0) * np.clip(fda, 0.0, 1.0))
    if kind == "weighted_mean":
        weight = float(operator["count_weight"])
        if not 0.0 <= weight <= 1.0:
            raise ValueError("count weight must be between zero and one")
        return weight * count + (1.0 - weight) * fda
    raise ValueError(f"unsupported consensus operator: {kind}")


def score_standardized_separation(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    finite = np.isfinite(scores)
    positive = scores[finite & (labels == 1)]
    negative = scores[finite & (labels == 0)]
    if not len(positive) or not len(negative):
        return np.nan
    pooled = np.sqrt(0.5 * (np.var(positive) + np.var(negative)))
    return float((np.mean(positive) - np.mean(negative)) / max(pooled, 1e-12))


def immutable_savez(path: Path, arrays: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    normalized = {key: np.asarray(value) for key, value in arrays.items()}
    if path.exists():
        with np.load(path, allow_pickle=False) as observed:
            if set(observed.files) != set(normalized):
                raise ValueError(f"immutable NPZ members differ: {path}")
            for key, expected in normalized.items():
                if not np.array_equal(observed[key], expected, equal_nan=True):
                    raise ValueError(f"immutable NPZ array differs: {path}::{key}")
        return
    np.savez_compressed(path, **normalized)


def write_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        frame.to_csv(
            path,
            sep="\t",
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
    else:
        frame.to_csv(path, sep="\t", index=False)


def load_threshold_scores(path: Path) -> tuple[dict, pd.DataFrame]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != THRESHOLD_SCHEMA:
        raise ValueError(f"unsupported threshold freeze: {path}")
    canonical_identifier(document, "threshold_id")
    if document.get("naked_dna_read") is not False:
        raise ValueError("threshold freeze was not created before naked-DNA access")
    require_checksum(document["policy"], "functional policy")
    require_checksum(document["reference_configuration"], "reference configuration")
    require_checksum(document["thresholds"], "threshold table")
    scores_record = document.get("validation_site_scores")
    if scores_record is None:
        raise ValueError("threshold freeze lacks validation site scores")
    scores_path = require_checksum(scores_record, "validation site scores")
    scores = pd.read_csv(scores_path, sep="\t")
    required = {
        "cell",
        "tf",
        "motif_family",
        "bias_configuration",
        "candidate_id",
        "artifact_index",
        "site_hash",
        "TFBS_chr",
        "label",
        "candidate_score",
        "dwm_score",
    }
    missing = required.difference(scores.columns)
    if missing:
        raise ValueError(
            "validation scores lack columns: " + ", ".join(sorted(missing))
        )
    return document, scores


def align_validation_scores(count: pd.DataFrame, fda: pd.DataFrame) -> pd.DataFrame:
    keys = ["cell", "tf", "artifact_index"]
    for name, frame in (("count", count), ("FDA", fda)):
        if frame.duplicated(keys).any():
            raise ValueError(f"{name} validation scores contain duplicate sites")
    merged = count.merge(
        fda,
        on=keys,
        how="outer",
        suffixes=("_count", "_fda"),
        indicator=True,
        validate="one_to_one",
    )
    if not merged["_merge"].eq("both").all():
        raise ValueError("count and FDA validation sites do not align")
    exact = [
        "motif_family",
        "bias_configuration",
        "site_hash",
        "TFBS_chr",
        "label",
    ]
    for column in exact:
        if not merged[f"{column}_count"].equals(merged[f"{column}_fda"]):
            raise ValueError(f"count and FDA validation {column} values differ")
    if not np.allclose(
        merged["dwm_score_count"],
        merged["dwm_score_fda"],
        rtol=0.0,
        atol=1e-12,
        equal_nan=True,
    ):
        raise ValueError("count and FDA validation DWM scores differ")
    return pd.DataFrame(
        {
            "cell": merged["cell"].astype(str),
            "tf": merged["tf"].astype(str),
            "motif_family": merged["motif_family_count"].astype(str),
            "bias_configuration": merged["bias_configuration_count"].astype(str),
            "artifact_index": merged["artifact_index"].astype(int),
            "site_hash": merged["site_hash_count"].astype("uint64"),
            "TFBS_chr": merged["TFBS_chr_count"].astype(str),
            "label": merged["label_count"].astype(int),
            "count_score": merged["candidate_score_count"].astype(float),
            "fda_score": merged["candidate_score_fda"].astype(float),
            "dwm_score": merged["dwm_score_count"].astype(float),
        }
    ).sort_values(keys, kind="mergesort").reset_index(drop=True)


def task_metric_record(
    frame: pd.DataFrame,
    score: np.ndarray,
    *,
    operator: str,
    minimum_sites_per_class: int,
) -> dict:
    labels = frame["label"].to_numpy(dtype=int)
    candidate = binary_metrics(labels, score)
    baseline = binary_metrics(labels, frame["dwm_score"].to_numpy(dtype=float))
    positive = int(candidate["positive_sites"])
    negative = int(candidate["negative_sites"])
    candidate_selection = selection_score(candidate)
    baseline_selection = selection_score(baseline)
    return {
        "cell": str(frame["cell"].iloc[0]),
        "tf": str(frame["tf"].iloc[0]),
        "motif_family": str(frame["motif_family"].iloc[0]),
        "operator": operator,
        "status": (
            "eligible"
            if min(positive, negative) >= minimum_sites_per_class
            else "underpowered"
        ),
        "n_sites": int(candidate["n_sites"]),
        "n_positive": positive,
        "n_negative": negative,
        "auroc": float(candidate["auroc"]),
        "auprc": float(candidate["auprc"]),
        "prevalence": float(candidate["prevalence"]),
        "dwm_auroc": float(baseline["auroc"]),
        "dwm_auprc": float(baseline["auprc"]),
        "auroc_gain_over_dwm": float(candidate["auroc"] - baseline["auroc"]),
        "relative_auprc_gain_over_dwm": float(
            (candidate["auprc"] - baseline["auprc"])
            / max(float(baseline["auprc"]), 1e-8)
        ),
        "selection_score": candidate_selection,
        "dwm_selection_score": baseline_selection,
        "selection_score_gain_over_dwm": candidate_selection - baseline_selection,
        "score_standardized_separation": score_standardized_separation(
            labels, score
        ),
    }


def summarize_operators(metrics: pd.DataFrame) -> pd.DataFrame:
    eligible = metrics[metrics["status"].astype(str).eq("eligible")]
    rows = []
    for operator in SIMPLICITY_ORDER:
        selected = eligible[eligible["operator"].astype(str).eq(operator)]
        ctcf = selected[selected["tf"].astype(str).eq("CTCF")]
        values = selected["selection_score_gain_over_dwm"].to_numpy(dtype=float)
        standard_error = (
            float(np.std(values, ddof=1) / np.sqrt(len(values)))
            if len(values) > 1
            else 0.0
        )
        rows.append(
            {
                "operator": operator,
                "simplicity_rank": SIMPLICITY_ORDER.index(operator),
                "eligible_tasks": len(selected),
                "eligible_families": selected["motif_family"].nunique(),
                "eligible_contexts": selected["cell"].nunique(),
                "mean_selection_score_gain_over_dwm": float(np.mean(values)),
                "selection_score_gain_standard_error": standard_error,
                "mean_auroc_gain_over_dwm": float(
                    selected["auroc_gain_over_dwm"].mean()
                ),
                "mean_relative_auprc_gain_over_dwm": float(
                    selected["relative_auprc_gain_over_dwm"].mean()
                ),
                "tasks_with_positive_auroc_gain": int(
                    selected["auroc_gain_over_dwm"].gt(0).sum()
                ),
                "tasks_with_positive_auprc_gain": int(
                    selected["relative_auprc_gain_over_dwm"].gt(0).sum()
                ),
                "minimum_ctcf_auroc_gain": (
                    float(ctcf["auroc_gain_over_dwm"].min()) if len(ctcf) else np.nan
                ),
                "passes_ctcf_nonregression": bool(
                    len(ctcf) and ctcf["auroc_gain_over_dwm"].ge(-0.02).all()
                ),
            }
        )
    return pd.DataFrame(rows)


def select_primary_operator(summary: pd.DataFrame) -> tuple[str, float, str]:
    candidates = summary[summary["passes_ctcf_nonregression"].eq(True)].copy()
    if candidates.empty:
        raise ValueError("no consensus operator passes CTCF non-regression")
    best = candidates.sort_values(
        ["mean_selection_score_gain_over_dwm", "simplicity_rank"],
        ascending=[False, True],
        kind="mergesort",
    ).iloc[0]
    cutoff = float(
        best["mean_selection_score_gain_over_dwm"]
        - best["selection_score_gain_standard_error"]
    )
    within_one_se = candidates[
        candidates["mean_selection_score_gain_over_dwm"].ge(cutoff)
    ].sort_values("simplicity_rank", kind="mergesort")
    selected = str(within_one_se.iloc[0]["operator"])
    return selected, cutoff, str(best["operator"])


def freeze_consensus(args: argparse.Namespace) -> int:
    count_document, count = load_threshold_scores(args.count_threshold_freeze)
    fda_document, fda = load_threshold_scores(args.fda_threshold_freeze)
    aligned = align_validation_scores(count, fda)
    arrays: dict[str, np.ndarray] = {}
    task_records = []
    metric_rows = []
    validation_frames = []
    for task_index, ((_cell, _tf), task) in enumerate(
        aligned.groupby(TASK_COLUMNS, sort=True)
    ):
        task = task.reset_index(drop=True)
        labels = task["label"].to_numpy(dtype=int)
        negative = labels == 0
        if not np.any(negative) or np.unique(labels).size != 2:
            raise ValueError(f"validation task lacks both classes: {_cell} {_tf}")
        count_key = f"count_negative_{task_index:03d}"
        fda_key = f"fda_negative_{task_index:03d}"
        arrays[count_key] = np.sort(task.loc[negative, "count_score"].to_numpy(float))
        arrays[fda_key] = np.sort(task.loc[negative, "fda_score"].to_numpy(float))
        count_rank = empirical_negative_cdf(task["count_score"], arrays[count_key])
        fda_rank = empirical_negative_cdf(task["fda_score"], arrays[fda_key])
        transformed = task.copy()
        transformed["count_negative_percentile"] = count_rank
        transformed["fda_negative_percentile"] = fda_rank
        thresholds = {}
        for operator in OPERATORS:
            name = str(operator["operator"])
            score = consensus_score(count_rank, fda_rank, operator)
            transformed[name] = score
            threshold, calls = upper_tail_threshold(
                score[negative], args.target_negative_call_rate
            )
            thresholds[name] = {
                "threshold": threshold,
                "validation_negative_calls": calls,
                "validation_negative_sites": int(np.sum(negative)),
            }
            metric_rows.append(
                task_metric_record(
                    task,
                    score,
                    operator=name,
                    minimum_sites_per_class=args.minimum_sites_per_class,
                )
            )
        dwm_threshold, dwm_calls = upper_tail_threshold(
            task.loc[negative, "dwm_score"].to_numpy(float),
            args.target_negative_call_rate,
        )
        task_records.append(
            {
                "cell": str(_cell),
                "tf": str(_tf),
                "motif_family": str(task["motif_family"].iloc[0]),
                "count_negative_array": count_key,
                "fda_negative_array": fda_key,
                "thresholds": thresholds,
                "dwm_threshold": dwm_threshold,
                "dwm_validation_negative_calls": dwm_calls,
                "dwm_validation_negative_sites": int(np.sum(negative)),
            }
        )
        validation_frames.append(transformed)

    metrics = pd.DataFrame(metric_rows)
    summary = summarize_operators(metrics)
    selected, cutoff, unconstrained_best = select_primary_operator(summary)
    summary["selected_for_primary"] = summary["operator"].astype(str).eq(selected)
    args.outdir.mkdir(parents=True, exist_ok=True)
    calibration_path = args.outdir / "functional_consensus_calibration.npz"
    metrics_path = args.outdir / "functional_consensus_validation_metrics.tsv"
    summary_path = args.outdir / "functional_consensus_validation_summary.tsv"
    scores_path = args.outdir / "functional_consensus_validation_scores.tsv.gz"
    immutable_savez(calibration_path, arrays)
    write_frame(metrics_path, metrics)
    write_frame(summary_path, summary)
    write_frame(scores_path, pd.concat(validation_frames, ignore_index=True))
    document = {
        "schema": FREEZE_SCHEMA,
        "count_threshold_freeze": {
            "path": str(args.count_threshold_freeze),
            "sha256": file_sha256(args.count_threshold_freeze),
            "threshold_id": count_document["threshold_id"],
            "policy_id": count_document["policy_id"],
        },
        "fda_threshold_freeze": {
            "path": str(args.fda_threshold_freeze),
            "sha256": file_sha256(args.fda_threshold_freeze),
            "threshold_id": fda_document["threshold_id"],
            "policy_id": fda_document["policy_id"],
        },
        "operators": list(OPERATORS),
        "selection_rule": {
            "metric": "mean label-adjusted AP plus AUROC gain over DWM",
            "one_standard_error": True,
            "simplicity_order": list(SIMPLICITY_ORDER),
            "ctcf_max_auroc_loss": 0.02,
            "unconstrained_best": unconstrained_best,
            "one_se_cutoff": cutoff,
        },
        "selected_operator": selected,
        "target_negative_call_rate": args.target_negative_call_rate,
        "minimum_sites_per_class": args.minimum_sites_per_class,
        "tasks": task_records,
        "calibration": {
            "path": str(calibration_path),
            "sha256": file_sha256(calibration_path),
            "kind": "smoothed empirical matched-negative CDF",
        },
        "validation_metrics": {
            "path": str(metrics_path),
            "sha256": file_sha256(metrics_path),
        },
        "validation_summary": {
            "path": str(summary_path),
            "sha256": file_sha256(summary_path),
        },
        "validation_scores": {
            "path": str(scores_path),
            "sha256": file_sha256(scores_path),
        },
        "training_labels_used_for_score_models": False,
        "validation_labels_used_for_operator_selection": True,
        "test_scores_read": False,
        "naked_dna_scores_read": False,
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["consensus_id"] = sha256(canonical.encode()).hexdigest()
    immutable_write_json(args.outdir / "functional_consensus.freeze.json", document)
    columns = [
        "operator",
        "eligible_tasks",
        "mean_auroc_gain_over_dwm",
        "mean_relative_auprc_gain_over_dwm",
        "minimum_ctcf_auroc_gain",
        "selected_for_primary",
    ]
    print(summary[columns].to_string(index=False))
    return 0


def load_consensus(path: Path) -> tuple[dict, dict[str, np.ndarray]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != FREEZE_SCHEMA:
        raise ValueError("unsupported functional consensus freeze")
    canonical_identifier(document, "consensus_id")
    for key in (
        "count_threshold_freeze",
        "fda_threshold_freeze",
        "calibration",
        "validation_metrics",
        "validation_summary",
        "validation_scores",
    ):
        require_checksum(document[key], key.replace("_", " "))
    if document.get("test_scores_read") is not False:
        raise ValueError("consensus was not frozen before test access")
    if document.get("naked_dna_scores_read") is not False:
        raise ValueError("consensus was not frozen before naked-DNA access")
    with np.load(document["calibration"]["path"], allow_pickle=False) as source:
        arrays = {key: np.asarray(source[key]) for key in source.files}
    for task in document["tasks"]:
        for key in (task["count_negative_array"], task["fda_negative_array"]):
            if key not in arrays:
                raise ValueError(f"consensus calibration lacks array: {key}")
    return document, arrays


def load_result_scores(
    path: Path,
    *,
    schema: str,
    output_name: str,
) -> tuple[dict, pd.DataFrame]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != schema:
        raise ValueError(f"unsupported result manifest: {path}")
    for name, record in document["outputs"].items():
        require_checksum(record, f"{name} result")
    score_path = Path(document["outputs"][output_name]["path"])
    return document, pd.read_csv(score_path, sep="\t")


def compare_columns(
    merged: pd.DataFrame,
    columns: list[str],
    *,
    left: str,
    right: str,
) -> None:
    for column in columns:
        first = merged[f"{column}_{left}"]
        second = merged[f"{column}_{right}"]
        if pd.api.types.is_numeric_dtype(first) and pd.api.types.is_numeric_dtype(second):
            equal = np.allclose(first, second, rtol=0.0, atol=1e-12, equal_nan=True)
        else:
            equal = first.fillna("<NA>").astype(str).equals(
                second.fillna("<NA>").astype(str)
            )
        if not equal:
            raise ValueError(f"count and FDA {column} values differ")


def align_test_scores(count: pd.DataFrame, fda: pd.DataFrame) -> pd.DataFrame:
    keys = ["cell", "tf", "artifact_index"]
    for name, frame in (("count", count), ("FDA", fda)):
        if frame.duplicated(keys).any():
            raise ValueError(f"{name} test scores contain duplicate sites")
    merged = count.merge(
        fda,
        on=keys,
        how="outer",
        suffixes=("_count", "_fda"),
        indicator=True,
        validate="one_to_one",
    )
    if not merged["_merge"].eq("both").all():
        raise ValueError("count and FDA test sites do not align")
    compare_columns(
        merged,
        [
            "motif_family",
            "bias_configuration",
            "TFBS_chr",
            "TFBS_start",
            "TFBS_end",
            "TFBS_strand",
            "motif_score",
            "accessibility",
            "label",
            "dwm_score",
            "direct_score",
            "log_accessibility",
            "motif_id",
        ],
        left="count",
        right="fda",
    )
    output = pd.DataFrame({column: merged[column] for column in keys})
    for column in (
        "motif_family",
        "bias_configuration",
        "TFBS_chr",
        "TFBS_start",
        "TFBS_end",
        "TFBS_strand",
        "motif_score",
        "accessibility",
        "label",
        "dwm_score",
        "direct_score",
        "log_accessibility",
        "motif_id",
    ):
        output[column] = merged[f"{column}_count"]
    output["count_score"] = merged["candidate_probability_count"].astype(float)
    output["fda_score"] = merged["candidate_probability_fda"].astype(float)
    return output.sort_values(keys, kind="mergesort").reset_index(drop=True)


def task_lookup(document: dict) -> dict[tuple[str, str], dict]:
    tasks = {(str(row["cell"]), str(row["tf"])): row for row in document["tasks"]}
    if len(tasks) != len(document["tasks"]):
        raise ValueError("consensus freeze contains duplicate tasks")
    return tasks


def apply_consensus_scores(
    frame: pd.DataFrame,
    document: dict,
    arrays: dict[str, np.ndarray],
) -> pd.DataFrame:
    tasks = task_lookup(document)
    frames = []
    for key, selected in frame.groupby(TASK_COLUMNS, sort=True):
        normalized_key = (str(key[0]), str(key[1]))
        if normalized_key not in tasks:
            raise ValueError(f"scores contain task absent from consensus: {normalized_key}")
        task = tasks[normalized_key]
        output = selected.copy()
        count_rank = empirical_negative_cdf(
            output["count_score"], arrays[task["count_negative_array"]]
        )
        fda_rank = empirical_negative_cdf(
            output["fda_score"], arrays[task["fda_negative_array"]]
        )
        output["count_negative_percentile"] = count_rank
        output["fda_negative_percentile"] = fda_rank
        for operator in document["operators"]:
            output[str(operator["operator"])] = consensus_score(
                count_rank, fda_rank, operator
            )
        frames.append(output)
    observed = set(zip(frame["cell"].astype(str), frame["tf"].astype(str), strict=True))
    if observed != set(tasks):
        missing = sorted(set(tasks).difference(observed))
        raise ValueError(f"consensus tasks absent from score input: {missing}")
    return pd.concat(frames, ignore_index=True)


def evaluate_test(args: argparse.Namespace) -> int:
    consensus, arrays = load_consensus(args.consensus_freeze)
    count_document, count = load_result_scores(
        args.count_manifest,
        schema=FUNCTIONAL_TEST_SCHEMA,
        output_name="site_scores",
    )
    fda_document, fda = load_result_scores(
        args.fda_manifest,
        schema=FUNCTIONAL_TEST_SCHEMA,
        output_name="site_scores",
    )
    if count_document["policy_id"] != consensus["count_threshold_freeze"]["policy_id"]:
        raise ValueError("count test policy does not match consensus")
    if fda_document["policy_id"] != consensus["fda_threshold_freeze"]["policy_id"]:
        raise ValueError("FDA test policy does not match consensus")
    scored = apply_consensus_scores(align_test_scores(count, fda), consensus, arrays)
    metrics_rows = []
    bootstrap_rows = []
    for (_cell, _tf), task in scored.groupby(TASK_COLUMNS, sort=True):
        task = task.reset_index(drop=True)
        labels = task["label"].to_numpy(dtype=int)
        positive = int(np.sum(labels == 1))
        negative = int(np.sum(labels == 0))
        for operator in consensus["operators"]:
            name = str(operator["operator"])
            metrics_rows.append(
                task_metric_record(
                    task,
                    task[name].to_numpy(dtype=float),
                    operator=name,
                    minimum_sites_per_class=args.minimum_sites_per_class,
                )
            )
            if min(positive, negative) >= args.minimum_sites_per_class:
                sites = task.rename(columns={"label": "chip_label"})
                bootstrap_rows.append(
                    {
                        "cell": str(_cell),
                        "tf": str(_tf),
                        "motif_family": str(task["motif_family"].iloc[0]),
                        "operator": name,
                        **block_bootstrap_delta(
                            sites,
                            task[name].to_numpy(dtype=float),
                            task["dwm_score"].to_numpy(dtype=float),
                            iterations=args.bootstrap_iterations,
                            seed=args.seed,
                        ),
                    }
                )
    metrics = pd.DataFrame(metrics_rows)
    summary = summarize_operators(metrics)
    summary["selected_for_primary"] = summary["operator"].astype(str).eq(
        str(consensus["selected_operator"])
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.outdir / "functional_consensus_test_metrics.tsv"
    bootstrap_path = args.outdir / "functional_consensus_test_bootstrap.tsv"
    summary_path = args.outdir / "functional_consensus_test_summary.tsv"
    scores_path = args.outdir / "functional_consensus_test_site_scores.tsv.gz"
    write_frame(metrics_path, metrics)
    write_frame(bootstrap_path, pd.DataFrame(bootstrap_rows))
    write_frame(summary_path, summary)
    write_frame(scores_path, scored)
    manifest = {
        "schema": TEST_SCHEMA,
        "consensus_id": consensus["consensus_id"],
        "selected_operator": consensus["selected_operator"],
        "consensus_freeze": {
            "path": str(args.consensus_freeze),
            "sha256": file_sha256(args.consensus_freeze),
        },
        "count_test_manifest": {
            "path": str(args.count_manifest),
            "sha256": file_sha256(args.count_manifest),
        },
        "fda_test_manifest": {
            "path": str(args.fda_manifest),
            "sha256": file_sha256(args.fda_manifest),
        },
        "test_labels_used_for_operator_selection": False,
        "models_refitted_on_test": False,
        "minimum_sites_per_class": args.minimum_sites_per_class,
        "bootstrap_iterations": args.bootstrap_iterations,
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "bootstrap": {
                "path": str(bootstrap_path),
                "sha256": file_sha256(bootstrap_path),
            },
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "scores": {"path": str(scores_path), "sha256": file_sha256(scores_path)},
        },
    }
    immutable_write_json(args.outdir / "functional_consensus_test_manifest.json", manifest)
    print(summary.to_string(index=False))
    return 0


def align_naked_scores(count: pd.DataFrame, fda: pd.DataFrame) -> pd.DataFrame:
    keys = ["cell", "tf", "replicate", "site_hash"]
    frames = {}
    for family, frame in (("count", count), ("fda", fda)):
        for method in ("candidate", "DWM"):
            selected = frame[frame["method"].astype(str).eq(method)].copy()
            if selected.duplicated(keys).any():
                raise ValueError(f"{family} {method} naked-DNA scores contain duplicates")
            frames[(family, method)] = selected
    candidate = frames[("count", "candidate")].merge(
        frames[("fda", "candidate")],
        on=keys,
        how="outer",
        suffixes=("_count", "_fda"),
        indicator=True,
        validate="one_to_one",
    )
    if not candidate["_merge"].eq("both").all():
        raise ValueError("count and FDA naked-DNA candidate sites do not align")
    compare_columns(
        candidate,
        ["valid", "informative"],
        left="count",
        right="fda",
    )
    dwm = frames[("count", "DWM")].merge(
        frames[("fda", "DWM")],
        on=keys,
        how="outer",
        suffixes=("_count", "_fda"),
        indicator=True,
        validate="one_to_one",
    )
    if not dwm["_merge"].eq("both").all():
        raise ValueError("count and FDA naked-DNA DWM sites do not align")
    compare_columns(
        dwm,
        ["score", "valid", "informative"],
        left="count",
        right="fda",
    )
    output = candidate[keys].copy()
    output["count_score"] = candidate["score_count"].astype(float)
    output["fda_score"] = candidate["score_fda"].astype(float)
    output["candidate_valid"] = candidate["valid_count"].astype(bool)
    output["candidate_informative"] = candidate["informative_count"].astype(bool)
    dwm_values = dwm[keys + ["score_count", "valid_count", "informative_count"]].rename(
        columns={
            "score_count": "dwm_score",
            "valid_count": "dwm_valid",
            "informative_count": "dwm_informative",
        }
    )
    output = output.merge(dwm_values, on=keys, validate="one_to_one")
    return output.sort_values(keys, kind="mergesort").reset_index(drop=True)


def safety_record(
    candidate_score: np.ndarray,
    candidate_valid: np.ndarray,
    candidate_informative: np.ndarray,
    candidate_threshold: float,
    dwm_score: np.ndarray,
    dwm_valid: np.ndarray,
    dwm_informative: np.ndarray,
    dwm_threshold: float,
) -> tuple[dict, np.ndarray]:
    candidate_finite = np.asarray(candidate_valid, bool) & np.isfinite(candidate_score)
    dwm_finite = np.asarray(dwm_valid, bool) & np.isfinite(dwm_score)
    candidate_calls = (
        candidate_finite
        & np.asarray(candidate_informative, bool)
        & (candidate_score >= candidate_threshold)
    )
    dwm_calls = (
        dwm_finite
        & np.asarray(dwm_informative, bool)
        & (dwm_score >= dwm_threshold)
    )
    candidate_total = int(np.sum(candidate_finite))
    candidate_count = int(np.sum(candidate_calls))
    _low, candidate_upper = wilson_interval(candidate_count, candidate_total)
    point_rate = candidate_count / candidate_total if candidate_total else np.nan
    paired = candidate_finite & dwm_finite
    paired_total = int(np.sum(paired))
    paired_candidate_count = int(np.sum(candidate_calls & paired))
    paired_dwm_count = int(np.sum(dwm_calls & paired))
    paired_candidate_rate = (
        paired_candidate_count / paired_total if paired_total else np.nan
    )
    paired_dwm_rate = paired_dwm_count / paired_total if paired_total else np.nan
    increase = paired_candidate_rate - paired_dwm_rate
    return (
        {
            "valid_sites": candidate_total,
            "informative_sites": int(
                np.sum(candidate_finite & np.asarray(candidate_informative, bool))
            ),
            "calls": candidate_count,
            "false_positive_rate": point_rate,
            "wilson_upper_95": candidate_upper,
            "paired_sites": paired_total,
            "paired_candidate_calls": paired_candidate_count,
            "paired_dwm_calls": paired_dwm_count,
            "paired_candidate_false_positive_rate": paired_candidate_rate,
            "paired_dwm_false_positive_rate": paired_dwm_rate,
            "candidate_minus_dwm": increase,
            "passes_point_rate": bool(point_rate <= 0.05),
            "passes_wilson": bool(candidate_upper <= 0.05),
            "passes_increase": bool(increase <= 0.01),
            "passes_safety": bool(
                point_rate <= 0.05
                and candidate_upper <= 0.05
                and increase <= 0.01
            ),
        },
        candidate_calls,
    )


def evaluate_naked(args: argparse.Namespace) -> int:
    consensus, arrays = load_consensus(args.consensus_freeze)
    count_document, count = load_result_scores(
        args.count_manifest,
        schema=FUNCTIONAL_NAKED_SCHEMA,
        output_name="scores",
    )
    fda_document, fda = load_result_scores(
        args.fda_manifest,
        schema=FUNCTIONAL_NAKED_SCHEMA,
        output_name="scores",
    )
    if count_document["policy_id"] != consensus["count_threshold_freeze"]["policy_id"]:
        raise ValueError("count naked-DNA policy does not match consensus")
    if fda_document["policy_id"] != consensus["fda_threshold_freeze"]["policy_id"]:
        raise ValueError("FDA naked-DNA policy does not match consensus")
    scored = apply_consensus_scores(align_naked_scores(count, fda), consensus, arrays)
    tasks = task_lookup(consensus)
    rate_rows = []
    call_columns = []
    for key, task_frame in scored.groupby(TASK_COLUMNS, sort=True):
        normalized_key = (str(key[0]), str(key[1]))
        frozen_task = tasks[normalized_key]
        indexes = task_frame.index.to_numpy(dtype=int)
        for operator in consensus["operators"]:
            name = str(operator["operator"])
            threshold = float(frozen_task["thresholds"][name]["threshold"])
            rate, calls = safety_record(
                task_frame[name].to_numpy(dtype=float),
                task_frame["candidate_valid"].to_numpy(dtype=bool),
                task_frame["candidate_informative"].to_numpy(dtype=bool),
                threshold,
                task_frame["dwm_score"].to_numpy(dtype=float),
                task_frame["dwm_valid"].to_numpy(dtype=bool),
                task_frame["dwm_informative"].to_numpy(dtype=bool),
                float(frozen_task["dwm_threshold"]),
            )
            column = f"called_{name}"
            if column not in scored:
                scored[column] = False
                call_columns.append(column)
            scored.loc[indexes, column] = calls
            rate_rows.append(
                {
                    "cell": normalized_key[0],
                    "tf": normalized_key[1],
                    "motif_family": frozen_task["motif_family"],
                    "replicate": str(task_frame["replicate"].iloc[0]),
                    "operator": name,
                    "threshold": threshold,
                    "selected_for_primary": name == consensus["selected_operator"],
                    **rate,
                }
            )
    rates = pd.DataFrame(rate_rows)
    summary = (
        rates.groupby(["operator", "selected_for_primary"], as_index=False)
        .agg(
            tasks=("tf", "size"),
            tasks_passing_safety=("passes_safety", "sum"),
            mean_false_positive_rate=("false_positive_rate", "mean"),
            maximum_false_positive_rate=("false_positive_rate", "max"),
            maximum_candidate_minus_dwm=("candidate_minus_dwm", "max"),
        )
        .sort_values("operator", kind="mergesort")
    )
    summary["all_tasks_pass_safety"] = (
        summary["tasks"] == summary["tasks_passing_safety"]
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    rates_path = args.outdir / "functional_consensus_naked_dna_rates.tsv"
    summary_path = args.outdir / "functional_consensus_naked_dna_summary.tsv"
    scores_path = args.outdir / "functional_consensus_naked_dna_site_scores.tsv.gz"
    write_frame(rates_path, rates)
    write_frame(summary_path, summary)
    write_frame(scores_path, scored)
    primary = rates[rates["selected_for_primary"].eq(True)]
    manifest = {
        "schema": NAKED_SCHEMA,
        "consensus_id": consensus["consensus_id"],
        "selected_operator": consensus["selected_operator"],
        "consensus_freeze": {
            "path": str(args.consensus_freeze),
            "sha256": file_sha256(args.consensus_freeze),
        },
        "count_naked_manifest": {
            "path": str(args.count_manifest),
            "sha256": file_sha256(args.count_manifest),
        },
        "fda_naked_manifest": {
            "path": str(args.fda_manifest),
            "sha256": file_sha256(args.fda_manifest),
        },
        "models_refitted_on_naked_dna": False,
        "thresholds_changed_on_naked_dna": False,
        "operators_selected_on_naked_dna": False,
        "primary_tasks": len(primary),
        "primary_tasks_passing_safety": int(primary["passes_safety"].sum()),
        "primary_all_tasks_pass_safety": bool(primary["passes_safety"].all()),
        "outputs": {
            "rates": {"path": str(rates_path), "sha256": file_sha256(rates_path)},
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "scores": {"path": str(scores_path), "sha256": file_sha256(scores_path)},
        },
    }
    immutable_write_json(args.outdir / "functional_consensus_naked_dna_manifest.json", manifest)
    print(summary.to_string(index=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    freeze = subparsers.add_parser("freeze", help="freeze validation calibration")
    freeze.add_argument("--count-threshold-freeze", type=Path, required=True)
    freeze.add_argument("--fda-threshold-freeze", type=Path, required=True)
    freeze.add_argument("--outdir", type=Path, required=True)
    freeze.add_argument("--target-negative-call-rate", type=float, default=0.025)
    freeze.add_argument("--minimum-sites-per-class", type=int, default=200)
    freeze.set_defaults(handler=freeze_consensus)

    test = subparsers.add_parser("test", help="evaluate frozen consensus on test")
    test.add_argument("--consensus-freeze", type=Path, required=True)
    test.add_argument("--count-manifest", type=Path, required=True)
    test.add_argument("--fda-manifest", type=Path, required=True)
    test.add_argument("--outdir", type=Path, required=True)
    test.add_argument("--minimum-sites-per-class", type=int, default=200)
    test.add_argument("--bootstrap-iterations", type=int, default=1000)
    test.add_argument("--seed", type=int, default=2026)
    test.set_defaults(handler=evaluate_test)

    naked = subparsers.add_parser(
        "naked",
        help="evaluate frozen consensus on independent naked DNA",
    )
    naked.add_argument("--consensus-freeze", type=Path, required=True)
    naked.add_argument("--count-manifest", type=Path, required=True)
    naked.add_argument("--fda-manifest", type=Path, required=True)
    naked.add_argument("--outdir", type=Path, required=True)
    naked.set_defaults(handler=evaluate_naked)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
