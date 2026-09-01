#!/usr/bin/env python3
"""Combine disjoint validation grids and select eligible per-TF winners."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from evaluate_strand_label_free_models import file_sha256, select_winners


SCHEMA = "fp-tools-strand-label-free-evaluation-v1"


def validate_run(path: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    manifest_path = path / "strand_label_free_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA:
        raise ValueError(f"unsupported functional-grid manifest: {manifest_path}")
    if manifest.get("locked_test_labels_read") is not False:
        raise ValueError(f"functional grid opened test labels: {manifest_path}")
    if manifest.get("training_labels_used") is not False:
        raise ValueError(f"functional grid used fitting labels: {manifest_path}")
    outputs = manifest["outputs"]
    for record in outputs.values():
        if file_sha256(record["path"]) != record["sha256"]:
            raise ValueError(f"functional-grid output changed: {record['path']}")
    metrics = pd.read_csv(outputs["metrics"]["path"], sep="\t")
    profiles = pd.read_csv(outputs["profiles"]["path"], sep="\t")
    return manifest, metrics, profiles


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimum-sites-per-class", type=int, default=200)
    args = parser.parse_args(argv)

    manifests = []
    metric_frames = []
    profile_frames = []
    source_records = []
    for path in args.run:
        manifest, metrics, profiles = validate_run(path)
        manifests.append(manifest)
        metric_frames.append(metrics)
        profile_frames.append(profiles)
        manifest_path = path / "strand_label_free_manifest.json"
        source_records.append(
            {
                "path": str(manifest_path),
                "sha256": file_sha256(manifest_path),
            }
        )
    study_hashes = {manifest["study_sha256"] for manifest in manifests}
    seeds = {int(manifest["seed"]) for manifest in manifests}
    artifact_signatures = {
        json.dumps(manifest["artifacts"], sort_keys=True, separators=(",", ":"))
        for manifest in manifests
    }
    if len(study_hashes) != 1 or len(seeds) != 1 or len(artifact_signatures) != 1:
        raise ValueError("functional grids do not share study, seed and artifacts")

    metrics = pd.concat(metric_frames, ignore_index=True)
    duplicate = metrics.duplicated(
        ["cell", "tf", "bias_configuration", "candidate_id"]
    )
    if duplicate.any():
        repeated = metrics.loc[duplicate, "candidate_id"].astype(str).unique()
        raise ValueError("functional grids overlap: " + ", ".join(sorted(repeated)))
    metrics["evaluation_status"] = "underpowered"
    eligible = (
        metrics["validation_positive_sites"].fillna(0).astype(int)
        >= args.minimum_sites_per_class
    ) & (
        metrics["validation_negative_sites"].fillna(0).astype(int)
        >= args.minimum_sites_per_class
    )
    metrics.loc[eligible, "evaluation_status"] = "eligible"
    winners = select_winners(metrics[eligible].copy())
    profiles = pd.concat(profile_frames, ignore_index=True)

    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.outdir / "strand_label_free_metrics.tsv.gz"
    profiles_path = args.outdir / "strand_label_free_profiles.tsv.gz"
    winners_path = args.outdir / "strand_label_free_winners.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    profiles.to_csv(profiles_path, sep="\t", index=False)
    winners.to_csv(winners_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "combined_selection": True,
        "locked_test_labels_read": False,
        "training_labels_used": False,
        "study": manifests[0]["study"],
        "study_sha256": manifests[0]["study_sha256"],
        "artifacts": manifests[0]["artifacts"],
        "candidate_count": int(metrics["candidate_id"].nunique()),
        "candidates": [
            row.dropna().to_dict()
            for _, row in metrics.drop_duplicates("candidate_id")[
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
        "seed": next(iter(seeds)),
        "metrics_rows": int(len(metrics)),
        "winner_rows": int(len(winners)),
        "source_runs": source_records,
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {
                "path": str(profiles_path),
                "sha256": file_sha256(profiles_path),
            },
            "winners": {"path": str(winners_path), "sha256": file_sha256(winners_path)},
        },
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["selection_id"] = sha256(canonical.encode()).hexdigest()
    manifest_path = args.outdir / "strand_label_free_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        winners[
            ["cell", "tf", "candidate_id", "auroc", "auprc"]
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
