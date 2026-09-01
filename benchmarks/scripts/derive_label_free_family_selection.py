#!/usr/bin/env python3
"""Derive a checksum-locked detector-family selection from one completed grid."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_strand_label_free_models import (  # noqa: E402
    file_sha256,
    select_winners,
)


SCHEMA = "fp-tools-strand-label-free-evaluation-v1"


def load_grid(path: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    manifest_path = path / "strand_label_free_manifest.json"
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    if document.get("schema") != SCHEMA:
        raise ValueError("unsupported label-free selection manifest")
    if document.get("locked_test_labels_read") is not False:
        raise ValueError("source grid opened locked test labels")
    if document.get("training_labels_used") is not False:
        raise ValueError("source grid used labels for model fitting")
    outputs = document.get("outputs", {})
    for name in ("metrics", "profiles", "winners"):
        record = outputs.get(name, {})
        if file_sha256(record.get("path", "")) != record.get("sha256"):
            raise ValueError(f"source grid {name} checksum mismatch")
    metrics = pd.read_csv(outputs["metrics"]["path"], sep="\t")
    profiles = pd.read_csv(outputs["profiles"]["path"], sep="\t")
    return document, metrics, profiles


def select_families(
    metrics: pd.DataFrame,
    profiles: pd.DataFrame,
    families: set[str],
    *,
    minimum_sites_per_class: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = metrics[metrics["family"].astype(str).isin(families)].copy()
    if selected.empty:
        raise ValueError("requested detector families have no candidates")
    selected["evaluation_status"] = "underpowered"
    eligible = (
        selected["validation_positive_sites"].fillna(0).astype(int)
        >= minimum_sites_per_class
    ) & (
        selected["validation_negative_sites"].fillna(0).astype(int)
        >= minimum_sites_per_class
    )
    selected.loc[eligible, "evaluation_status"] = "eligible"
    winners = select_winners(selected[eligible].copy())
    expected = selected.loc[eligible, ["cell", "tf", "bias_configuration"]].drop_duplicates()
    observed = winners[["cell", "tf", "bias_configuration"]].drop_duplicates()
    coverage = expected.merge(observed, how="outer", indicator=True)
    if not coverage["_merge"].eq("both").all():
        raise ValueError("a requested family lacks a usable winner for an eligible TF")
    candidate_ids = set(selected["candidate_id"].astype(str))
    selected_profiles = profiles[
        profiles["candidate_id"].astype(str).isin(candidate_ids)
    ].copy()
    return selected, selected_profiles, winners


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--family", action="append", required=True)
    parser.add_argument("--minimum-sites-per-class", type=int, default=200)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    families = {str(value) for value in args.family}
    source, metrics, profiles = load_grid(args.run)
    selected, selected_profiles, winners = select_families(
        metrics,
        profiles,
        families,
        minimum_sites_per_class=args.minimum_sites_per_class,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.outdir / "strand_label_free_metrics.tsv.gz"
    profiles_path = args.outdir / "strand_label_free_profiles.tsv.gz"
    winners_path = args.outdir / "strand_label_free_winners.tsv"
    selected.to_csv(metrics_path, sep="\t", index=False)
    selected_profiles.to_csv(profiles_path, sep="\t", index=False)
    winners.to_csv(winners_path, sep="\t", index=False)
    source_manifest = args.run / "strand_label_free_manifest.json"
    document = {
        "schema": SCHEMA,
        "derived_family_selection": True,
        "locked_test_labels_read": False,
        "training_labels_used": False,
        "study": source["study"],
        "study_sha256": source["study_sha256"],
        "artifacts": source["artifacts"],
        "families": sorted(families),
        "candidate_count": int(selected["candidate_id"].nunique()),
        "candidates": [
            row.dropna().to_dict()
            for _, row in selected.drop_duplicates("candidate_id")[
                [
                    "candidate_id",
                    "family",
                    "smoother",
                    "background",
                    "window",
                    "channel",
                    "training_pool",
                    "anchor_strength",
                    "covariate_ridge",
                ]
            ].iterrows()
        ],
        "minimum_evaluation_sites": args.minimum_sites_per_class,
        "seed": int(source["seed"]),
        "metrics_rows": int(len(selected)),
        "winner_rows": int(len(winners)),
        "source_run": {
            "path": str(source_manifest),
            "sha256": file_sha256(source_manifest),
        },
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {
                "path": str(profiles_path),
                "sha256": file_sha256(profiles_path),
            },
            "winners": {"path": str(winners_path), "sha256": file_sha256(winners_path)},
        },
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["selection_id"] = sha256(canonical.encode()).hexdigest()
    manifest_path = args.outdir / "strand_label_free_manifest.json"
    manifest_path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(winners[["cell", "tf", "candidate_id", "auroc", "auprc"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
