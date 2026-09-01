#!/usr/bin/env python3
"""Screen frozen bias predictions for unexplained naked-DNA motif structure.

The input profiles must be generated from an independent enzyme-only library.
This command compares observed and predicted cuts at motif sites, rejects any
site table containing ChIP labels, and applies the frozen absolute 0.25
center/flank residual rule with a bootstrap confidence interval.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from fp_tools.tools.frozen_bias_evaluation import motif_residual_effect  # noqa: E402


SCHEMA = "fp-tools-frozen-bias-motif-residual-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("profiles must use LABEL=PREFIX")
    label, path = value.split("=", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError("profiles must use LABEL=PREFIX")
    return label, Path(path)


def artifact_paths(prefix: Path) -> tuple[Path, Path, Path]:
    if prefix.suffix in {".npz", ".json"}:
        prefix = prefix.with_suffix("")
    return (
        Path(str(prefix) + ".npz"),
        Path(str(prefix) + ".json"),
        Path(str(prefix) + ".sites.tsv.gz"),
    )


def load_profile_artifact(prefix: Path) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    npz_path, json_path, sites_path = artifact_paths(prefix)
    document = json.loads(json_path.read_text(encoding="utf-8"))
    if document.get("profiles_sha256") != file_sha256(npz_path):
        raise ValueError("profile checksum mismatch")
    if document.get("sites_sha256") != file_sha256(sites_path):
        raise ValueError("site-table checksum mismatch")
    sites = pd.read_csv(sites_path, sep="\t")
    forbidden = [
        column
        for column in sites.columns
        if "chip" in column.lower() or "label" in column.lower()
    ]
    if forbidden:
        raise ValueError(
            "naked-DNA motif residual inputs cannot contain ChIP labels: "
            + ", ".join(forbidden)
        )
    with np.load(npz_path, allow_pickle=False) as source:
        arrays = {name: np.asarray(source[name]) for name in source.files}
    required = {
        "plus_observed",
        "minus_observed",
        "plus_expected",
        "minus_expected",
        "valid",
    }
    missing = required.difference(arrays)
    if missing:
        raise ValueError("profile artifact is missing: " + ", ".join(sorted(missing)))
    if len(sites) != len(arrays["valid"]):
        raise ValueError("profile arrays and site table have different lengths")
    return arrays, sites


def summarize_artifact(
    label: str,
    arrays: dict[str, np.ndarray],
    sites: pd.DataFrame,
    *,
    threshold: float,
    bootstraps: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed = np.asarray(
        arrays["plus_observed"] + arrays["minus_observed"], dtype=float
    )
    expected = np.asarray(
        arrays["plus_expected"] + arrays["minus_expected"], dtype=float
    )
    valid = np.asarray(arrays["valid"], dtype=bool)
    positions = np.arange(observed.shape[1], dtype=float) - observed.shape[1] // 2
    if "tf" not in sites.columns:
        raise ValueError("site table must contain tf")
    summary_rows: list[dict[str, object]] = []
    curve_frames: list[pd.DataFrame] = []
    for tf, group in sites.groupby("tf", sort=True):
        indexes = group.index.to_numpy(dtype=int)
        indexes = indexes[valid[indexes] & (observed[indexes].sum(axis=1) > 0)]
        if len(indexes) < 2:
            continue
        result = motif_residual_effect(
            observed[indexes],
            expected[indexes],
            positions,
            bootstraps=bootstraps,
            seed=seed,
            threshold=threshold,
        )
        first = sites.loc[indexes[0]]
        summary_rows.append(
            {
                "candidate_id": label,
                "tf": str(tf),
                "motif_family": str(first.get("motif_family", "")),
                "motif_id": str(first.get("motif", first.get("motif_id", ""))),
                **result,
            }
        )
        mean_observed = observed[indexes].mean(axis=0)
        mean_expected = expected[indexes].mean(axis=0)
        curve_frames.append(
            pd.DataFrame(
                {
                    "candidate_id": label,
                    "tf": str(tf),
                    "position": positions.astype(int),
                    "mean_observed": mean_observed,
                    "mean_expected": mean_expected,
                    "mean_log_ratio_residual": np.log(
                        (mean_observed + 0.5) / (mean_expected + 0.5)
                    ),
                }
            )
        )
    return pd.DataFrame(summary_rows), (
        pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame()
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument(
        "--profile", type=parse_named_path, action="append", required=True
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    if study.get("status") != "development_locked_holdout_unscored":
        raise ValueError("motif residual screening requires the locked, unscored study")
    threshold = float(
        study["promotion_gates"]["motif_residual_absolute_center_flank_limit"]
    )
    summary_frames, curve_frames = [], []
    for label, prefix in args.profile:
        arrays, sites = load_profile_artifact(prefix)
        summary, curves = summarize_artifact(
            label,
            arrays,
            sites,
            threshold=threshold,
            bootstraps=args.bootstraps,
            seed=args.seed,
        )
        summary_frames.append(summary)
        curve_frames.append(curves)
    summaries = pd.concat(summary_frames, ignore_index=True)
    curves = pd.concat(curve_frames, ignore_index=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "motif_residual_summary.tsv"
    curves_path = args.outdir / "motif_residual_profiles.tsv.gz"
    summaries.to_csv(summary_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "profiles": [
            {
                "candidate_id": label,
                "prefix": str(prefix),
                "npz_sha256": file_sha256(artifact_paths(prefix)[0]),
                "sites_sha256": file_sha256(artifact_paths(prefix)[2]),
            }
            for label, prefix in args.profile
        ],
        "threshold": threshold,
        "chipped_labels_used": False,
        "passed": bool(not summaries["motif_residual_flag"].any()),
        "outputs": {
            summary_path.name: file_sha256(summary_path),
            curves_path.name: file_sha256(curves_path),
        },
    }
    (args.outdir / "motif_residual_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(summaries.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
