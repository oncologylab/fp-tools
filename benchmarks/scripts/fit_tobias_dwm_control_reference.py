#!/usr/bin/env python3
"""Fit a safe TOBIAS-style DWM reference on frozen control windows.

The output contains numerical arrays and provenance in NPZ plus JSON; it never
serializes a Python pickle.  This is a research baseline for paired likelihood
qualification, not a replacement for the current production DWM code.
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
from fit_combined_control_bias import load_datasets, parse_named_path  # noqa: E402
from fp_tools.tools.frozen_bias_evaluation import (  # noqa: E402
    TobiasDwmReferenceModel,
    conditional_control_scores,
)


SCHEMA = "fp-tools-tobias-dwm-control-fit-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fit_references(
    training: dict[tuple[tuple[int, int], str], ControlWindowDataset],
    validation: dict[tuple[tuple[int, int], str], ControlWindowDataset],
    outdir: Path,
    *,
    context_length: int = 25,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)
    artifact_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    shifts = sorted({key[0] for key in training})
    if set(shifts) != {key[0] for key in validation}:
        raise ValueError("training and validation must contain the same shifts")
    for shift in shifts:
        model = TobiasDwmReferenceModel(context_length=context_length)
        train_names = sorted(
            name for candidate_shift, name in training if candidate_shift == shift
        )
        arrays = [
            training[(shift, name)].model_arrays(model.feature_spec)
            for name in train_names
        ]
        contexts = np.concatenate([value[0] for value in arrays], axis=0)
        counts = np.concatenate([value[1] for value in arrays], axis=0)
        model.fit(contexts, counts)
        stem = outdir / f"tobias_dwm.shift_{shift[0]}_{shift[1]}"
        npz_path, json_path = model.save(
            stem,
            metadata={
                "training_source": "combined_naked_dna_mitochondrial",
                "configuration": "conventional_tobias_style_dwm",
                "read_shift": list(shift),
                "training_datasets": train_names,
                "background": "uniform valid candidate contexts",
                "k_flank": (context_length - 1) // 2,
            },
        )
        artifact_rows.append(
            {
                "candidate_id": f"DWM_{shift[0]}_{shift[1]}",
                "shift_forward": shift[0],
                "shift_reverse": shift[1],
                "model_npz": str(npz_path),
                "model_json": str(json_path),
                "model_size_mb": npz_path.stat().st_size / (1024 * 1024),
                "training_contexts": int(len(contexts) * contexts.shape[1]),
                "training_cuts": int(np.sum(counts)),
            }
        )
        for candidate_shift, name in sorted(validation):
            if candidate_shift != shift:
                continue
            candidate_contexts, candidate_counts = validation[
                (candidate_shift, name)
            ].model_arrays(model.feature_spec)
            scores = conditional_control_scores(
                model, candidate_contexts, candidate_counts
            )
            metric_rows.append(
                {
                    "candidate_id": f"DWM_{shift[0]}_{shift[1]}",
                    "library": name,
                    "shift_forward": shift[0],
                    "shift_reverse": shift[1],
                    "windows": len(scores.totals),
                    "cuts": int(np.sum(scores.totals)),
                    "conditional_nll": scores.conditional_nll,
                    "nll_gain": scores.nll_gain,
                    "multinomial_deviance_per_cut": scores.deviance_per_cut,
                    "calibration_error": scores.calibration_error,
                    "aggregate_jsd": scores.aggregate_jsd,
                }
            )
    return pd.DataFrame(artifact_rows), pd.DataFrame(metric_rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument(
        "--training-dataset", type=parse_named_path, action="append", required=True
    )
    parser.add_argument(
        "--validation-dataset", type=parse_named_path, action="append", required=True
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--context-length",
        type=int,
        default=25,
        help="odd DWM context length; 25 matches the production ±12-bp default",
    )
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    if study.get("status") != "development_locked_holdout_unscored":
        raise ValueError("DWM qualification requires the locked, unscored study")
    training = load_datasets(args.training_dataset, required_split="train")
    validation = load_datasets(args.validation_dataset, required_split="validation")
    artifacts, metrics = fit_references(
        training,
        validation,
        args.outdir,
        context_length=args.context_length,
    )
    artifacts_path = args.outdir / "dwm_reference_models.tsv"
    metrics_path = args.outdir / "dwm_reference_metrics.tsv"
    artifacts.to_csv(artifacts_path, sep="\t", index=False)
    metrics.to_csv(metrics_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "training_datasets": [
            {"name": name, "path": str(path), "sha256": file_sha256(path)}
            for name, path in args.training_dataset
        ],
        "validation_datasets": [
            {"name": name, "path": str(path), "sha256": file_sha256(path)}
            for name, path in args.validation_dataset
        ],
        "outputs": {
            path.name: {"path": str(path), "sha256": file_sha256(path)}
            for path in (artifacts_path, metrics_path)
        },
        "chipped_labels_used": False,
        "context_length": args.context_length,
    }
    (args.outdir / "dwm_reference_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(metrics.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
