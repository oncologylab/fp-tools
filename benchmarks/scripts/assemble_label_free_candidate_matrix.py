#!/usr/bin/env python3
"""Assemble prior-free DWM and strand-aware label-free detector metrics."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd


IDENTITY = ("cell", "tf", "motif_family")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selection_score(auroc: pd.Series, auprc: pd.Series, prevalence: pd.Series) -> pd.Series:
    denominator = np.maximum(1.0 - prevalence.to_numpy(dtype=float), 1e-8)
    return auroc.to_numpy(dtype=float) + (
        auprc.to_numpy(dtype=float) - prevalence.to_numpy(dtype=float)
    ) / denominator


def prepare_dwm(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        *IDENTITY,
        "correction",
        "candidate_id",
        "family",
        "prior_constraint",
        "training_labels_used",
        "status",
        "shape_auroc",
        "shape_auprc",
        "shape_brier",
        "prevalence",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"DWM metrics lack columns: {missing}")
    output = frame[
        frame["status"].eq("ok")
        & frame["correction"].eq("DWM")
        & frame["training_labels_used"].eq(False)
        & frame["prior_constraint"].eq("none")
        & ~frame["family"].eq("anchored-fda")
    ].copy()
    output["bias_configuration"] = "DWM"
    output["auroc"] = output["shape_auroc"]
    output["auprc"] = output["shape_auprc"]
    output["brier"] = output["shape_brier"]
    output["selection_score"] = selection_score(
        output["auroc"], output["auprc"], output["prevalence"]
    )
    output["motif_or_accessibility_features_used"] = False
    output["source"] = "dwm_prior_free_shape"
    return output


def prepare_strand(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        *IDENTITY,
        "bias_configuration",
        "candidate_id",
        "family",
        "training_labels_used",
        "motif_or_accessibility_features_used",
        "status",
        "auroc",
        "auprc",
        "brier",
        "prevalence",
        "selection_score",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"strand metrics lack columns: {missing}")
    output = frame[frame["status"].eq("ok")].copy()
    if output["training_labels_used"].astype(bool).any():
        raise ValueError("strand candidate matrix refuses models trained with labels")
    if output["motif_or_accessibility_features_used"].astype(bool).any():
        raise ValueError("strand candidate matrix refuses motif/accessibility features")
    output["source"] = "strand_label_free"
    return output


def assemble(dwm: pd.DataFrame, strand: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    left = prepare_dwm(dwm)
    right = prepare_strand(strand)
    common_columns = sorted(set(left.columns).intersection(right.columns))
    preferred = [
        *IDENTITY,
        "bias_configuration",
        "candidate_id",
        "family",
        "source",
        "training_labels_used",
        "motif_or_accessibility_features_used",
        "n_sites",
        "positive_sites",
        "negative_sites",
        "prevalence",
        "auroc",
        "auprc",
        "brier",
        "selection_score",
    ]
    columns = [name for name in preferred if name in common_columns]
    columns += [name for name in common_columns if name not in columns]
    matrix = pd.concat([left[columns], right[columns]], ignore_index=True)
    duplicate = list(IDENTITY) + ["bias_configuration", "candidate_id"]
    if matrix.duplicated(duplicate).any():
        raise ValueError("candidate matrix has duplicate task/configuration rows")
    audit_rows = []
    for keys, group in matrix.groupby(list(IDENTITY), sort=True):
        audit_rows.append(
            {
                **dict(zip(IDENTITY, keys)),
                "configurations": int(group[["bias_configuration", "candidate_id"]].drop_duplicates().shape[0]),
                "bias_configurations": ",".join(sorted(group["bias_configuration"].unique())),
                "site_count_min": int(group["n_sites"].min()) if "n_sites" in group else 0,
                "site_count_max": int(group["n_sites"].max()) if "n_sites" in group else 0,
                "prevalence_min": float(group["prevalence"].min()),
                "prevalence_max": float(group["prevalence"].max()),
            }
        )
    audit = pd.DataFrame(audit_rows)
    if np.any(audit["prevalence_max"] - audit["prevalence_min"] > 1e-10):
        raise ValueError("DWM and strand evaluations do not use identical task prevalence")
    if "n_sites" in matrix and np.any(audit["site_count_min"] != audit["site_count_max"]):
        raise ValueError("DWM and strand evaluations do not use identical site counts")
    return matrix.sort_values(list(IDENTITY) + ["bias_configuration", "candidate_id"]), audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dwm-metrics", type=Path, required=True)
    parser.add_argument("--strand-metrics", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    matrix, audit = assemble(
        pd.read_csv(args.dwm_metrics, sep="\t"),
        pd.read_csv(args.strand_metrics, sep="\t"),
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix_path = args.outdir / "label_free_candidate_matrix.tsv.gz"
    audit_path = args.outdir / "label_free_candidate_audit.tsv"
    matrix.to_csv(matrix_path, sep="\t", index=False)
    audit.to_csv(audit_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-label-free-candidate-matrix-v1",
        "locked_test_labels_read": False,
        "training_labels_used": False,
        "motif_or_accessibility_features_used": False,
        "dwm_metrics": str(args.dwm_metrics),
        "dwm_metrics_sha256": file_sha256(args.dwm_metrics),
        "strand_metrics": str(args.strand_metrics),
        "strand_metrics_sha256": file_sha256(args.strand_metrics),
        "rows": int(len(matrix)),
        "tasks": int(matrix[list(IDENTITY)].drop_duplicates().shape[0]),
        "outputs": {
            path.name: {"path": str(path), "sha256": file_sha256(path)}
            for path in (matrix_path, audit_path)
        },
    }
    (args.outdir / "label_free_candidate_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(audit.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
