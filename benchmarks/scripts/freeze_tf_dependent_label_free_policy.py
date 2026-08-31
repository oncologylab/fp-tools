#!/usr/bin/env python3
"""Freeze fail-closed TF-family detector routes before holdout scoring."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd


TASK_KEYS = ("cell", "tf", "motif_family")
CONFIG_KEYS = ("bias_configuration", "candidate_id")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def choose_configuration(frame: pd.DataFrame, *, require_configuration: str | None = None) -> pd.Series:
    values = frame.copy()
    if require_configuration is not None:
        values = values[values["bias_configuration"].eq(require_configuration)]
    task_count = values[list(TASK_KEYS)].drop_duplicates().shape[0]
    summary = (
        values.groupby(list(CONFIG_KEYS), sort=True)
        .agg(
            contexts=("tf", "size"),
            mean_selection_score=("selection_score", "mean"),
            minimum_selection_score=("selection_score", "min"),
            selection_sd=("selection_score", "std"),
            mean_auroc=("auroc", "mean"),
            mean_auprc=("auprc", "mean"),
        )
        .reset_index()
    )
    summary = summary[summary["contexts"].eq(task_count)].copy()
    if summary.empty:
        raise ValueError("no candidate configuration covers every requested context")
    summary["robust_selection_score"] = summary["mean_selection_score"] - (
        0.5
        * summary["selection_sd"].fillna(0.0)
        / np.sqrt(summary["contexts"].clip(lower=1))
    )
    return summary.sort_values(
        ["robust_selection_score", "mean_selection_score", "minimum_selection_score", *CONFIG_KEYS],
        ascending=[False, False, False, True, True],
        kind="mergesort",
    ).iloc[0]


def configuration_rows(frame: pd.DataFrame, selected: pd.Series) -> pd.DataFrame:
    return frame[
        frame["bias_configuration"].eq(selected["bias_configuration"])
        & frame["candidate_id"].eq(selected["candidate_id"])
    ].copy()


def freeze_policy(
    matrix: pd.DataFrame,
    *,
    minimum_contexts: int,
    minimum_mean_auroc_gain: float,
    minimum_relative_auprc_gain: float,
    maximum_context_auroc_loss: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    required = {*TASK_KEYS, *CONFIG_KEYS, "selection_score", "auroc", "auprc"}
    missing = sorted(required.difference(matrix.columns))
    if missing:
        raise ValueError(f"candidate matrix lacks columns: {missing}")
    global_dwm = choose_configuration(matrix, require_configuration="DWM")
    family_rows = []
    context_rows = []
    for motif_family, family in matrix.groupby("motif_family", sort=True):
        new_rows = family[~family["bias_configuration"].eq("DWM")]
        if new_rows.empty:
            reference = choose_configuration(family, require_configuration="DWM")
            family_rows.append(
                {
                    "motif_family": motif_family,
                    "development_contexts": 0,
                    "candidate_bias_configuration": None,
                    "candidate_id": None,
                    "reference_bias_configuration": reference["bias_configuration"],
                    "reference_candidate_id": reference["candidate_id"],
                    "mean_auroc_gain": np.nan,
                    "relative_auprc_gain": np.nan,
                    "maximum_context_auroc_loss": np.nan,
                    "contexts_with_positive_auroc_gain": 0,
                    "new_route_passes_development_gates": False,
                    "failure_reasons": "no_eligible_new_candidate_context",
                    "recommended_bias_configuration": reference["bias_configuration"],
                    "recommended_candidate_id": reference["candidate_id"],
                }
            )
            continue
        candidate_tasks = new_rows[list(TASK_KEYS)].drop_duplicates()
        comparable = family.merge(candidate_tasks, on=list(TASK_KEYS), how="inner")
        contexts = comparable[list(TASK_KEYS)].drop_duplicates()
        reference = choose_configuration(comparable, require_configuration="DWM")
        candidate = choose_configuration(
            comparable[~comparable["bias_configuration"].eq("DWM")]
        )
        reference_rows = configuration_rows(comparable, reference)
        candidate_rows = configuration_rows(comparable, candidate)
        paired = reference_rows[list(TASK_KEYS) + ["auroc", "auprc"]].merge(
            candidate_rows[list(TASK_KEYS) + ["auroc", "auprc"]],
            on=list(TASK_KEYS),
            suffixes=("_dwm", "_candidate"),
            validate="one_to_one",
        )
        paired["auroc_gain"] = paired["auroc_candidate"] - paired["auroc_dwm"]
        paired["auprc_gain"] = paired["auprc_candidate"] - paired["auprc_dwm"]
        paired["motif_family"] = motif_family
        mean_auroc_gain = float(paired["auroc_gain"].mean())
        mean_dwm_ap = float(paired["auprc_dwm"].mean())
        relative_ap_gain = float(paired["auprc_candidate"].mean() / max(mean_dwm_ap, 1e-8) - 1.0)
        maximum_loss = float(np.maximum(-paired["auroc_gain"], 0.0).max())
        eligible = len(contexts) >= minimum_contexts
        passes = bool(
            eligible
            and mean_auroc_gain >= minimum_mean_auroc_gain
            and relative_ap_gain >= minimum_relative_auprc_gain
            and maximum_loss <= maximum_context_auroc_loss
        )
        reasons = []
        if not eligible:
            reasons.append("insufficient_independent_contexts")
        if mean_auroc_gain < minimum_mean_auroc_gain:
            reasons.append("mean_auroc_gain")
        if relative_ap_gain < minimum_relative_auprc_gain:
            reasons.append("relative_auprc_gain")
        if maximum_loss > maximum_context_auroc_loss:
            reasons.append("context_nonregression")
        recommended = candidate if passes else reference
        family_rows.append(
            {
                "motif_family": motif_family,
                "development_contexts": int(len(contexts)),
                "candidate_bias_configuration": candidate["bias_configuration"],
                "candidate_id": candidate["candidate_id"],
                "reference_bias_configuration": reference["bias_configuration"],
                "reference_candidate_id": reference["candidate_id"],
                "mean_auroc_gain": mean_auroc_gain,
                "relative_auprc_gain": relative_ap_gain,
                "maximum_context_auroc_loss": maximum_loss,
                "contexts_with_positive_auroc_gain": int(np.sum(paired["auroc_gain"] > 0)),
                "new_route_passes_development_gates": passes,
                "failure_reasons": ",".join(reasons),
                "recommended_bias_configuration": recommended["bias_configuration"],
                "recommended_candidate_id": recommended["candidate_id"],
            }
        )
        context_rows.append(paired)
    return (
        pd.DataFrame(family_rows),
        pd.concat(context_rows, ignore_index=True) if context_rows else pd.DataFrame(),
        global_dwm,
    )


def map_holdout_routes(study: dict, policy: pd.DataFrame, global_dwm: pd.Series) -> pd.DataFrame:
    routes = policy.set_index("motif_family")
    rows = []
    for task in study["tasks"]:
        if task["split"] != "locked_holdout":
            continue
        family = str(task["motif_family"])
        if family in routes.index:
            route = routes.loc[family]
            source = "development_family"
            bias_configuration = route["recommended_bias_configuration"]
            candidate_id = route["recommended_candidate_id"]
        else:
            source = "unseen_family_dwm_fallback"
            bias_configuration = global_dwm["bias_configuration"]
            candidate_id = global_dwm["candidate_id"]
        rows.append(
            {
                "cell": task["cell"],
                "tf": task["tf"],
                "motif_id": task["motif_id"],
                "motif_family": family,
                "role": task["role"],
                "route_source": source,
                "bias_configuration": bias_configuration,
                "candidate_id": candidate_id,
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"),
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimum-contexts", type=int, default=2)
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    gates = study["promotion_gates"]
    matrix = pd.read_csv(args.matrix, sep="\t")
    policy, context, global_dwm = freeze_policy(
        matrix,
        minimum_contexts=args.minimum_contexts,
        minimum_mean_auroc_gain=float(gates["minimum_mean_auroc_gain"]),
        minimum_relative_auprc_gain=float(gates["minimum_relative_auprc_gain"]),
        maximum_context_auroc_loss=float(gates["maximum_positive_control_auroc_loss"]),
    )
    holdout = map_holdout_routes(study, policy, global_dwm)
    args.outdir.mkdir(parents=True, exist_ok=True)
    policy_path = args.outdir / "tf_family_label_free_policy.tsv"
    context_path = args.outdir / "tf_family_development_contexts.tsv"
    holdout_path = args.outdir / "locked_holdout_routes.tsv"
    policy.to_csv(policy_path, sep="\t", index=False)
    context.to_csv(context_path, sep="\t", index=False)
    holdout.to_csv(holdout_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-tf-dependent-label-free-policy-v1",
        "locked_holdout_labels_read": False,
        "policy_frozen_before_holdout": True,
        "unseen_family_behavior": "DWM fallback",
        "matrix": str(args.matrix),
        "matrix_sha256": file_sha256(args.matrix),
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "gates": {
            "minimum_contexts": args.minimum_contexts,
            "minimum_mean_auroc_gain": gates["minimum_mean_auroc_gain"],
            "minimum_relative_auprc_gain": gates["minimum_relative_auprc_gain"],
            "maximum_context_auroc_loss": gates["maximum_positive_control_auroc_loss"],
        },
        "global_dwm_fallback": {
            "bias_configuration": global_dwm["bias_configuration"],
            "candidate_id": global_dwm["candidate_id"],
        },
        "outputs": {
            path.name: {"path": str(path), "sha256": file_sha256(path)}
            for path in (policy_path, context_path, holdout_path)
        },
    }
    (args.outdir / "tf_family_label_free_policy_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(policy.to_string(index=False))
    print("\nLocked holdout routes")
    print(holdout.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
