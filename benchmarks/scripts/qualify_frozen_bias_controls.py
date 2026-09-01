#!/usr/bin/env python3
"""Qualify frozen bias ensembles on independent enzyme-control windows.

Candidate selection uses only conditional cut likelihood, calibration,
cross-library stability, paired block-bootstrap intervals, runtime, and model
size.  No ChIP labels or motif identities are accepted by this command.
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
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_parametric_bias import ControlWindowDataset  # noqa: E402
from fp_tools.tools.frozen_bias_evaluation import (  # noqa: E402
    TOBIAS_DWM_SCHEMA,
    TobiasDwmReferenceModel,
    conditional_control_scores,
    paired_block_bootstrap_gain,
    retain_control_candidates,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    ConditionalSequenceBiasModel,
    encode_sequence,
    reverse_complement_contexts,
)


SCHEMA = "fp-tools-frozen-bias-control-qualification-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("values must use NAME=PATH")
    name, raw_path = value.split("=", 1)
    if not name or not raw_path:
        raise argparse.ArgumentTypeError("values must use NAME=PATH")
    return name, Path(raw_path)


def load_control_model(path: Path):
    document = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    if document.get("schema") == TOBIAS_DWM_SCHEMA:
        return TobiasDwmReferenceModel.load(path)
    return ConditionalSequenceBiasModel.load(path)


def model_arrays_with_keys(
    dataset: ControlWindowDataset,
    model,
    library: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce ``model_arrays`` while retaining stable window identities."""

    full_sequence_valid = np.asarray(
        [set(str(value).upper()).issubset({"A", "C", "G", "T"}) for value in dataset.sequences],
        dtype=bool,
    )
    if not np.any(full_sequence_valid):
        raise ValueError(f"control library {library!r} has no fully resolved windows")
    original_indexes = np.flatnonzero(full_sequence_valid)
    dataset = dataset.subset(original_indexes, split=dataset.split)
    context_length = int(model.feature_spec.context_length)
    left = context_length // 2
    first = int(dataset.margin) - left
    if first < 0:
        raise ValueError("dataset margin is too small for a candidate model")
    encoded = np.stack([encode_sequence(str(value)) for value in dataset.sequences])
    windows = np.lib.stride_tricks.sliding_window_view(encoded, context_length, axis=1)
    forward = np.asarray(
        windows[:, first : first + dataset.window_size], dtype=np.uint8
    )
    reverse = reverse_complement_contexts(forward)
    contexts = np.concatenate([forward, reverse], axis=0)
    counts = np.concatenate(
        [dataset.forward_counts, dataset.reverse_counts], axis=0
    ).astype(float)
    chromosomes = np.tile(np.asarray(dataset.chromosomes, dtype=str), 2)
    strands = np.repeat(["forward", "reverse"], len(dataset.starts))
    indexes = np.tile(original_indexes, 2)
    keys = np.asarray(
        [f"{library}|{strand}|{index}" for strand, index in zip(strands, indexes)]
    )
    blocks = np.asarray(
        [f"{library}|{chromosome}" for chromosome in chromosomes], dtype=str
    )
    valid = np.all(contexts < 4, axis=2)
    counts[~valid] = 0.0
    keep = (counts.sum(axis=1) > 0) & valid.any(axis=1)
    return contexts[keep], counts[keep], keys[keep], blocks[keep]


def score_candidate(
    candidate_id: str,
    model_path: Path,
    datasets: dict[str, ControlWindowDataset],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    model = load_control_model(model_path)
    expected_shift = tuple(int(value) for value in model.metadata.get("read_shift", ()))
    window_frames: list[pd.DataFrame] = []
    library_rows: list[dict[str, object]] = []
    for library, dataset in sorted(datasets.items()):
        if expected_shift and tuple(dataset.shift) != expected_shift:
            continue
        contexts, counts, keys, blocks = model_arrays_with_keys(dataset, model, library)
        scores = conditional_control_scores(model, contexts, counts)
        window_frames.append(
            pd.DataFrame(
                {
                    "candidate_id": candidate_id,
                    "library": library,
                    "window_key": keys,
                    "block": blocks,
                    "log_likelihood": scores.log_likelihood,
                    "null_log_likelihood": scores.null_log_likelihood,
                    "deviance": scores.deviance,
                    "cuts": scores.totals,
                }
            )
        )
        library_rows.append(
            {
                "candidate_id": candidate_id,
                "library": library,
                "windows": len(scores.totals),
                "cuts": int(np.sum(scores.totals)),
                "conditional_nll": scores.conditional_nll,
                "nll_gain": scores.nll_gain,
                "multinomial_deviance_per_cut": scores.deviance_per_cut,
                "calibration_error": scores.calibration_error,
                "aggregate_jsd": scores.aggregate_jsd,
            }
        )
    if not window_frames:
        raise ValueError(f"candidate {candidate_id} matched no control datasets")
    metadata = {
        "candidate_id": candidate_id,
        "model_npz": str(model_path),
        "model_sha256": file_sha256(model_path),
        "model_json": str(model_path.with_suffix(".json")),
        "context_length": int(model.feature_spec.context_length),
        "feature_name": model.feature_spec.name,
        "shift_forward": int(expected_shift[0]) if expected_shift else np.nan,
        "shift_reverse": int(expected_shift[1]) if expected_shift else np.nan,
        "model_size_mb": model_path.stat().st_size / (1024 * 1024),
    }
    return (
        pd.concat(window_frames, ignore_index=True),
        pd.DataFrame(library_rows),
        metadata,
    )


def summarize_candidates(
    windows: pd.DataFrame,
    libraries: pd.DataFrame,
    model_metadata: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for candidate_id, group in windows.groupby("candidate_id", sort=True):
        block = (
            group.groupby("block", as_index=False)
            .agg(log_likelihood=("log_likelihood", "sum"), cuts=("cuts", "sum"))
            .assign(block_nll=lambda value: -value["log_likelihood"] / value["cuts"])
        )
        library = libraries[libraries["candidate_id"].eq(candidate_id)]
        metadata = model_metadata[model_metadata["candidate_id"].eq(candidate_id)].iloc[
            0
        ]
        rows.append(
            {
                **metadata.to_dict(),
                "evaluated_libraries": int(library["library"].nunique()),
                "windows": int(len(group)),
                "cuts": int(group["cuts"].sum()),
                "mean_conditional_nll": float(
                    -group["log_likelihood"].sum() / group["cuts"].sum()
                ),
                "standard_error_conditional_nll": float(
                    block["block_nll"].std(ddof=1) / np.sqrt(len(block))
                    if len(block) > 1
                    else np.nan
                ),
                "mean_nll_gain": float(library["nll_gain"].mean()),
                "minimum_library_nll_gain": float(library["nll_gain"].min()),
                "mean_deviance": float(
                    np.sum(group["deviance"]) / np.sum(group["cuts"])
                ),
                "mean_calibration_error": float(library["calibration_error"].mean()),
                "mean_aggregate_jsd": float(library["aggregate_jsd"].mean()),
                "cross_library_nll_sd": float(library["conditional_nll"].std(ddof=1)),
            }
        )
    return pd.DataFrame(rows)


def paired_candidate_gain(
    windows: pd.DataFrame,
    candidate_id: str,
    reference_id: str,
    *,
    bootstraps: int,
    seed: int,
) -> dict[str, object]:
    candidate = windows[windows["candidate_id"].eq(candidate_id)].drop(
        columns="candidate_id"
    )
    reference = windows[windows["candidate_id"].eq(reference_id)].drop(
        columns="candidate_id"
    )
    paired = candidate.merge(
        reference,
        on=["library", "window_key", "block"],
        suffixes=("_candidate", "_reference"),
        how="inner",
        validate="one_to_one",
    )
    paired_support = np.isclose(
        paired["cuts_candidate"], paired["cuts_reference"], rtol=0.0, atol=1e-8
    )
    paired = paired.loc[paired_support].copy()
    if paired.empty:
        raise ValueError(
            f"candidate {candidate_id} and reference {reference_id} lack paired windows"
        )
    result = paired_block_bootstrap_gain(
        paired["log_likelihood_candidate"].to_numpy(),
        paired["log_likelihood_reference"].to_numpy(),
        paired["cuts_candidate"].to_numpy(),
        paired["block"].to_numpy(),
        bootstraps=bootstraps,
        seed=seed,
    )
    library_rows = []
    for library, group in paired.groupby("library", sort=True):
        library_rows.append(
            (
                library,
                float(
                    np.sum(
                        group["log_likelihood_candidate"]
                        - group["log_likelihood_reference"]
                    )
                    / np.sum(group["cuts_candidate"])
                ),
            )
        )
    return {
        "candidate_id": candidate_id,
        "reference_id": reference_id,
        "candidate_finite_windows": int(len(candidate)),
        "reference_finite_windows": int(len(reference)),
        "paired_finite_windows": int(len(paired)),
        "paired_support_fraction": float(
            len(paired) / max(len(candidate), len(reference), 1)
        ),
        "paired_cuts": float(paired["cuts_candidate"].sum()),
        "paired_libraries": len(library_rows),
        "minimum_library_gain": min(value for _name, value in library_rows),
        **result,
    }


def qualify(
    candidates: Sequence[tuple[str, Path]],
    datasets: Sequence[tuple[str, Path]],
    *,
    reference_id: str | None,
    bootstraps: int,
    seed: int,
    maximum_model_size_mb: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    loaded_datasets = {name: ControlWindowDataset.load(path) for name, path in datasets}
    if len(loaded_datasets) != len(datasets):
        raise ValueError("dataset names must be unique")
    window_frames = []
    library_frames = []
    metadata_rows = []
    for candidate_id, path in candidates:
        candidate_windows, candidate_libraries, metadata = score_candidate(
            candidate_id, path, loaded_datasets
        )
        window_frames.append(candidate_windows)
        library_frames.append(candidate_libraries)
        metadata_rows.append(metadata)
    windows = pd.concat(window_frames, ignore_index=True)
    libraries = pd.concat(library_frames, ignore_index=True)
    metadata = pd.DataFrame(metadata_rows)
    summary = summarize_candidates(windows, libraries, metadata)

    paired_rows: list[dict[str, object]] = []
    if reference_id is not None:
        if reference_id not in set(metadata["candidate_id"]):
            raise ValueError(f"unknown reference candidate: {reference_id}")
        for candidate_id in metadata["candidate_id"]:
            if candidate_id == reference_id:
                continue
            paired_rows.append(
                paired_candidate_gain(
                    windows,
                    candidate_id,
                    reference_id,
                    bootstraps=bootstraps,
                    seed=seed,
                )
            )
    paired = pd.DataFrame(paired_rows)

    eligible = summary[
        (summary["minimum_library_nll_gain"] > 0)
        & (summary["model_size_mb"] <= maximum_model_size_mb)
    ]
    if eligible.empty:
        selection = retain_control_candidates(
            summary, maximum_model_size_mb=maximum_model_size_mb
        )
        return windows, libraries, paired, selection
    best = eligible.sort_values(
        ["mean_conditional_nll", "context_length", "candidate_id"],
        kind="mergesort",
    ).iloc[0]
    threshold = float(best["mean_conditional_nll"]) + float(
        best["standard_error_conditional_nll"]
    )
    smallest = (
        eligible[eligible["mean_conditional_nll"] <= threshold]
        .sort_values(
            ["context_length", "model_size_mb", "mean_conditional_nll", "candidate_id"],
            kind="mergesort",
        )
        .iloc[0]
    )
    gains_over_smallest = []
    for candidate_id in summary["candidate_id"]:
        if candidate_id == smallest["candidate_id"]:
            gains_over_smallest.append(np.nan)
        else:
            gain = paired_candidate_gain(
                windows,
                str(candidate_id),
                str(smallest["candidate_id"]),
                bootstraps=bootstraps,
                seed=seed + 1,
            )
            gains_over_smallest.append(float(gain["paired_gain_lower_95"]))
    summary["gain_over_smallest_lower_95"] = gains_over_smallest
    selection = retain_control_candidates(
        summary, maximum_model_size_mb=maximum_model_size_mb
    )
    if reference_id is not None and not paired.empty:
        reference_metrics = paired.set_index("candidate_id")
        selection["passed_reference_gate"] = selection["candidate_id"].map(
            lambda value: bool(
                value == reference_id
                or (
                    value in reference_metrics.index
                    and reference_metrics.at[value, "minimum_library_gain"] > 0
                    and reference_metrics.at[value, "paired_gain_lower_95"] > 0
                )
            )
        )
        selection.loc[
            selection["retained"] & ~selection["passed_reference_gate"],
            ["retained", "retention_reason"],
        ] = [False, "failed paired reference likelihood gate"]
    return windows, libraries, paired, selection


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument(
        "--candidate", type=parse_named_path, action="append", required=True
    )
    parser.add_argument(
        "--dataset", type=parse_named_path, action="append", required=True
    )
    parser.add_argument("--reference-id")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bootstraps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    if study.get("status") != "development_locked_holdout_unscored":
        raise ValueError("control qualification requires the locked, unscored study")
    args.outdir.mkdir(parents=True, exist_ok=True)
    windows, libraries, paired, selection = qualify(
        args.candidate,
        args.dataset,
        reference_id=args.reference_id,
        bootstraps=args.bootstraps,
        seed=args.seed,
        maximum_model_size_mb=float(study["promotion_gates"]["maximum_model_size_mb"]),
    )
    paths = {
        "control_window_scores.tsv.gz": windows,
        "control_library_metrics.tsv": libraries,
        "paired_control_likelihood.tsv": paired,
        "control_candidate_selection.tsv": selection,
    }
    for name, frame in paths.items():
        frame.to_csv(args.outdir / name, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "candidates": [
            {"candidate_id": name, "path": str(path), "sha256": file_sha256(path)}
            for name, path in args.candidate
        ],
        "datasets": [
            {"library": name, "path": str(path), "sha256": file_sha256(path)}
            for name, path in args.dataset
        ],
        "reference_id": args.reference_id,
        "chipped_labels_used": False,
        "outputs": {
            name: {
                "path": str(args.outdir / name),
                "sha256": file_sha256(args.outdir / name),
            }
            for name in paths
        },
    }
    (args.outdir / "control_qualification_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(selection.sort_values("mean_conditional_nll").to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
