#!/usr/bin/env python3
"""Test condition-specific changes in complete TF footprint functions."""

from __future__ import annotations

import argparse
from dataclasses import asdict
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

from fp_tools.tools.functional_footprints import (  # noqa: E402
    functional_differential_test,
    normalize_functional_profiles,
)


REQUIRED_MANIFEST_COLUMNS = {"sample", "condition", "replicate", "profiles_npz", "sites_tsv"}
CHANNELS = {
    "combined": "combined_residual",
    "shared": "shared_strand_residual",
    "antisymmetric": "antisymmetric_strand_residual",
}


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_contrast(value: str) -> tuple[str, str]:
    fields = value.split(",")
    if len(fields) != 2 or not all(fields) or fields[0] == fields[1]:
        raise argparse.ArgumentTypeError("contrasts must use CONDITION_A,CONDITION_B")
    return fields[0], fields[1]


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    output = np.full_like(values, np.nan)
    finite = np.flatnonzero(np.isfinite(values))
    if len(finite) == 0:
        return output
    order = finite[np.argsort(values[finite])]
    ranked = values[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    output[order] = np.minimum(adjusted, 1.0)
    return output


def load_manifest_profiles(
    manifest: pd.DataFrame,
    channel: str,
    *,
    positions: np.ndarray,
) -> dict[str, list[tuple[np.ndarray, str, str, str]]]:
    """Load valid per-site curves grouped by TF without unsafe pickle."""

    array_name = CHANNELS[channel]
    grouped: dict[str, list[tuple[np.ndarray, str, str, str]]] = {}
    for row in manifest.itertuples(index=False):
        sites = pd.read_csv(row.sites_tsv, sep="\t")
        if "tf" not in sites:
            raise ValueError(f"{row.sites_tsv} must contain a tf column")
        with np.load(row.profiles_npz, allow_pickle=False) as arrays:
            if array_name not in arrays or "valid" not in arrays:
                raise ValueError(f"{row.profiles_npz} is missing {array_name} or valid")
            profiles = np.asarray(arrays[array_name], dtype=float)
            valid = np.asarray(arrays["valid"], dtype=bool)
        if profiles.shape != (len(sites), len(positions)) or valid.shape != (len(sites),):
            raise ValueError(f"{row.profiles_npz} does not match its sites table or profile flank")
        profiles = normalize_functional_profiles(profiles, positions)
        for tf, indexes in sites.groupby("tf", sort=True).groups.items():
            selected = np.asarray(list(indexes), dtype=int)
            selected = selected[valid[selected]]
            if len(selected):
                grouped.setdefault(str(tf), []).append(
                    (profiles[selected], str(row.condition), str(row.replicate), str(row.sample))
                )
    return grouped


def replicate_consistency(entries: list[tuple[np.ndarray, str, str, str]], contrast: tuple[str, str]) -> float:
    curves: dict[str, list[np.ndarray]] = {contrast[0]: [], contrast[1]: []}
    for profiles, condition, _replicate, _sample in entries:
        if condition in curves:
            curves[condition].append(np.mean(profiles, axis=0))
    correlations = []
    for condition_curves in curves.values():
        for first in range(len(condition_curves)):
            for second in range(first + 1, len(condition_curves)):
                correlations.append(np.corrcoef(condition_curves[first], condition_curves[second])[0, 1])
    return float(np.nanmean(correlations)) if correlations else np.nan


def evaluate_differential_profiles(
    manifest: pd.DataFrame,
    contrasts: Sequence[tuple[str, str]],
    *,
    channel: str,
    flank: int,
    bootstraps: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    positions = np.arange(-flank, flank + 1, dtype=float)
    by_tf = load_manifest_profiles(manifest, channel, positions=positions)
    summary_rows: list[dict] = []
    curve_rows: list[dict] = []
    for contrast in contrasts:
        for tf, entries in sorted(by_tf.items()):
            profiles = []
            conditions = []
            replicates = []
            sites_by_condition = {contrast[0]: 0, contrast[1]: 0}
            replicate_names = {contrast[0]: set(), contrast[1]: set()}
            for values, condition, replicate, _sample in entries:
                if condition not in contrast:
                    continue
                profiles.append(values)
                conditions.extend([condition] * len(values))
                replicates.extend([replicate] * len(values))
                sites_by_condition[condition] += len(values)
                replicate_names[condition].add(replicate)
            if any(len(replicate_names[condition]) < 2 for condition in contrast):
                continue
            result = functional_differential_test(
                np.concatenate(profiles),
                conditions,
                replicates,
                contrast,
                positions=positions,
                n_bootstrap=bootstraps,
                seed=seed,
            )
            descriptors = asdict(result.descriptor_change)
            summary_rows.append(
                {
                    "tf": tf,
                    "condition_a": contrast[0],
                    "condition_b": contrast[1],
                    "channel": channel,
                    "replicates_a": len(replicate_names[contrast[0]]),
                    "replicates_b": len(replicate_names[contrast[1]]),
                    "sites_a": sites_by_condition[contrast[0]],
                    "sites_b": sites_by_condition[contrast[1]],
                    "global_pvalue": result.global_pvalue,
                    "replicate_profile_correlation": replicate_consistency(entries, contrast),
                    **{f"change_{name}": value for name, value in descriptors.items()},
                }
            )
            for index, position in enumerate(result.positions):
                curve_rows.append(
                    {
                        "tf": tf,
                        "condition_a": contrast[0],
                        "condition_b": contrast[1],
                        "channel": channel,
                        "position": int(position),
                        "difference": float(result.difference[index]),
                        "pointwise_lower_95": float(result.pointwise_lower[index]),
                        "pointwise_upper_95": float(result.pointwise_upper[index]),
                        "simultaneous_lower_95": float(result.simultaneous_lower[index]),
                        "simultaneous_upper_95": float(result.simultaneous_upper[index]),
                    }
                )
    summary = pd.DataFrame(summary_rows)
    if len(summary):
        summary["global_qvalue"] = benjamini_hochberg(summary["global_pvalue"].to_numpy())
    return summary, pd.DataFrame(curve_rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--contrast", type=parse_contrast, action="append", required=True)
    parser.add_argument("--channel", choices=tuple(CHANNELS), default="combined")
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = pd.read_csv(args.manifest, sep="\t")
    missing = REQUIRED_MANIFEST_COLUMNS.difference(manifest.columns)
    if missing:
        raise SystemExit("manifest is missing columns: " + ", ".join(sorted(missing)))
    summary, curves = evaluate_differential_profiles(
        manifest,
        args.contrast,
        channel=args.channel,
        flank=args.flank,
        bootstraps=args.bootstraps,
        seed=args.seed,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "differential_functional_summary.tsv"
    curves_path = args.outdir / "differential_functional_curves.tsv.gz"
    summary.to_csv(summary_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    document = {
        "schema": "fp-tools-differential-functional-benchmark-v1",
        "manifest": str(args.manifest),
        "manifest_sha256": file_sha256(args.manifest),
        "contrasts": [list(value) for value in args.contrast],
        "channel": args.channel,
        "flank": args.flank,
        "bootstraps": args.bootstraps,
        "seed": args.seed,
        "outputs": {
            "summary": {"path": str(summary_path), "sha256": file_sha256(summary_path)},
            "curves": {"path": str(curves_path), "sha256": file_sha256(curves_path)},
        },
    }
    (args.outdir / "differential_functional_manifest.json").write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
