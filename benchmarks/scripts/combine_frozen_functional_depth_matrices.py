#!/usr/bin/env python3
"""Combine integrity-checked frozen depth runs without rebuilding profiles."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_frozen_functional_depth_matrix import (  # noqa: E402
    SCHEMA as DEPTH_SCHEMA,
    classify_depth,
    summarize_metrics,
    summarize_replicates,
)
from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from freeze_label_free_functional_models import immutable_write_json  # noqa: E402


SCHEMA = "fp-tools-combined-frozen-functional-depth-matrix-v2"


def load_manifest(path: Path) -> tuple[dict, dict[str, pd.DataFrame]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != DEPTH_SCHEMA:
        raise ValueError(f"unsupported frozen depth manifest: {path}")
    if document.get("raw_signal_guardrail") is not True:
        raise ValueError(f"frozen depth manifest lacks the raw-signal guardrail: {path}")
    frames = {}
    for name in ("metrics", "profiles", "artifacts"):
        record = document["outputs"][name]
        if file_sha256(record["path"]) != record["sha256"]:
            raise ValueError(f"frozen depth {name} checksum mismatch: {path}")
        frames[name] = pd.read_csv(record["path"], sep="\t")
    return document, frames


def reject_duplicates(frame: pd.DataFrame, columns: list[str], name: str) -> None:
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f"{name} lacks columns: {', '.join(missing)}")
    duplicated = frame.duplicated(columns, keep=False)
    if duplicated.any():
        example = frame.loc[duplicated, columns].iloc[0].to_dict()
        raise ValueError(f"duplicate {name} row: {example}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", action="append", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if len(set(args.manifest)) != len(args.manifest):
        raise SystemExit("duplicate depth manifest")

    documents = []
    grouped: dict[str, list[pd.DataFrame]] = {
        "metrics": [],
        "profiles": [],
        "artifacts": [],
    }
    for path in args.manifest:
        document, frames = load_manifest(path)
        documents.append((path, document))
        for name, frame in frames.items():
            frame = frame.copy()
            frame["source_depth_manifest"] = str(path)
            grouped[name].append(frame)
    policy_ids = {str(document["policy_id"]) for _path, document in documents}
    if len(policy_ids) != 1:
        raise ValueError("depth manifests use different frozen policies")
    dispersions = {float(document["dispersion"]) for _path, document in documents}
    if len(dispersions) != 1:
        raise ValueError("depth manifests use different dispersions")

    metrics = pd.concat(grouped["metrics"], ignore_index=True)
    profiles = pd.concat(grouped["profiles"], ignore_index=True)
    artifacts = pd.concat(grouped["artifacts"], ignore_index=True)
    reject_duplicates(
        metrics,
        ["cell", "sample", "tf", "candidate_id", "method", "depth", "seed"],
        "metric",
    )
    reject_duplicates(
        profiles,
        [
            "cell",
            "sample",
            "tf",
            "candidate_id",
            "method",
            "depth",
            "seed",
            "position",
        ],
        "profile",
    )
    reject_duplicates(artifacts, ["cell", "sample", "depth", "seed"], "artifact")
    summary = summarize_metrics(metrics)
    replicate_summary = summarize_replicates(metrics)
    classification = classify_depth(replicate_summary)

    args.outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "metrics": args.outdir / "frozen_functional_depth_metrics.tsv.gz",
        "profiles": args.outdir / "frozen_functional_depth_profiles.tsv.gz",
        "artifacts": args.outdir / "frozen_functional_depth_artifacts.tsv",
        "summary": args.outdir / "frozen_functional_depth_summary.tsv",
        "replicate_summary": (
            args.outdir / "frozen_functional_depth_replicate_summary.tsv"
        ),
        "classification": args.outdir / "frozen_functional_depth_classification.tsv",
    }
    metrics.to_csv(paths["metrics"], sep="\t", index=False)
    profiles.to_csv(paths["profiles"], sep="\t", index=False)
    artifacts.to_csv(paths["artifacts"], sep="\t", index=False)
    summary.to_csv(paths["summary"], sep="\t", index=False)
    replicate_summary.to_csv(paths["replicate_summary"], sep="\t", index=False)
    classification.to_csv(paths["classification"], sep="\t", index=False)
    document = {
        "schema": SCHEMA,
        "policy_id": next(iter(policy_ids)),
        "dispersion": next(iter(dispersions)),
        "source_manifests": [
            {"path": str(path), "sha256": file_sha256(path)}
            for path, _document in documents
        ],
        "depths": sorted(metrics["depth"].astype(str).unique()),
        "seeds": sorted(int(value) for value in metrics["seed"].unique()),
        "models_refitted_by_depth": False,
        "raw_signal_guardrail": True,
        "outputs": {
            name: {"path": str(path), "sha256": file_sha256(path)}
            for name, path in paths.items()
        },
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["combined_depth_id"] = sha256(canonical.encode()).hexdigest()
    immutable_write_json(args.outdir / "combined_depth_manifest.json", document)
    print(classification.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
