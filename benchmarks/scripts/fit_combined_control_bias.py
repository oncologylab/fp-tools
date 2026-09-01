#!/usr/bin/env python3
"""Fit frozen sequence-bias models jointly across enzyme-control sources.

This research-only command consumes checksum-protected control-window artifacts
created by ``evaluate_parametric_bias.py``.  It concatenates the requested
training libraries, fits the fixed parametric grid without ChIP labels, and
geometrically ensembles matching seed fits.  The fitted coefficient artifact
remains one compact NPZ plus JSON pair.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter
from typing import Sequence

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_parametric_bias import (  # noqa: E402
    ControlWindowDataset,
    _fit_one,
    conditional_metrics,
    feature_spec,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    ConditionalSequenceBiasModel,
    ensemble_sequence_bias_models,
)


SCHEMA = "fp-tools-combined-control-bias-v1"
MODEL_NAMES = ("selma10", "loglinear21", "loglinear41", "loglinear81")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("datasets must use NAME=CONTROL.npz")
    name, raw_path = value.split("=", 1)
    if not name or not raw_path:
        raise argparse.ArgumentTypeError("datasets must use NAME=CONTROL.npz")
    return name, Path(raw_path)


def safe_token(value: object) -> str:
    return str(value).replace("-", "m").replace(".", "p").replace("/", "_")


def load_datasets(
    values: Sequence[tuple[str, Path]],
    *,
    required_split: str,
) -> dict[tuple[tuple[int, int], str], ControlWindowDataset]:
    datasets: dict[tuple[tuple[int, int], str], ControlWindowDataset] = {}
    for name, path in values:
        dataset = ControlWindowDataset.load(path)
        if dataset.split != required_split:
            raise ValueError(
                f"{name} has split {dataset.split!r}; expected {required_split!r}"
            )
        key = (tuple(dataset.shift), name)
        if key in datasets:
            raise ValueError(
                f"duplicate dataset identity: {name}, shift={dataset.shift}"
            )
        datasets[key] = dataset
    return datasets


def fit_combined_grid(
    training: dict[tuple[tuple[int, int], str], ControlWindowDataset],
    validation: dict[tuple[tuple[int, int], str], ControlWindowDataset],
    outdir: Path,
    *,
    models: Sequence[str],
    l2_values: Sequence[float],
    seeds: Sequence[int],
    epochs: int,
    batch_windows: int,
    jobs: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Fit seed models and compact coefficient-mean ensembles."""

    if jobs < 1:
        raise ValueError("jobs must be positive")

    shifts = sorted({key[0] for key in training})
    if set(shifts) != {key[0] for key in validation}:
        raise ValueError("training and validation must contain the same read shifts")
    model_dir = outdir / "models"
    ensemble_dir = outdir / "ensembles"
    model_dir.mkdir(parents=True, exist_ok=True)
    ensemble_dir.mkdir(parents=True, exist_ok=True)
    artifact_rows: list[dict[str, object]] = []
    ensemble_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []

    for shift in shifts:
        training_names = sorted(
            name for candidate_shift, name in training if candidate_shift == shift
        )
        validation_names = sorted(
            name for candidate_shift, name in validation if candidate_shift == shift
        )
        if len(training_names) < 2:
            raise ValueError(
                "combined-source fitting requires at least two training datasets"
            )
        for model_name in models:
            spec = feature_spec(model_name)
            arrays = [
                training[(shift, name)].model_arrays(spec) for name in training_names
            ]
            contexts = np.concatenate([value[0] for value in arrays], axis=0)
            counts = np.concatenate([value[1] for value in arrays], axis=0)
            validation_arrays = {
                name: validation[(shift, name)].model_arrays(spec)
                for name in validation_names
            }
            for l2 in l2_values:
                members: list[ConditionalSequenceBiasModel] = []
                member_paths: list[Path] = []
                total_runtime = 0.0
                maximum_memory = 0.0

                def fit_seed(seed: int):
                    return _fit_one(
                        spec,
                        contexts,
                        counts,
                        l2=float(l2),
                        epochs=epochs,
                        batch_windows=batch_windows,
                        seed=int(seed),
                    )

                if jobs == 1 or len(seeds) == 1:
                    fitted_seeds = [fit_seed(int(seed)) for seed in seeds]
                else:
                    with ThreadPoolExecutor(
                        max_workers=min(int(jobs), len(seeds)),
                        thread_name_prefix="frozen-bias-seed",
                    ) as executor:
                        fitted_seeds = list(executor.map(fit_seed, seeds))

                for seed, (model, runtime, memory) in zip(seeds, fitted_seeds):
                    stem = model_dir / (
                        f"combined.{model_name}.shift_{shift[0]}_{shift[1]}."
                        f"l2_{safe_token(l2)}.seed_{seed}"
                    )
                    npz_path, json_path = model.save(
                        stem,
                        metadata={
                            "training_source": "combined_naked_dna_mitochondrial",
                            "configuration": "frozen_cross_control_pooled",
                            "read_shift": list(shift),
                            "training_datasets": training_names,
                            "seed": int(seed),
                        },
                    )
                    members.append(model)
                    member_paths.append(npz_path)
                    total_runtime += runtime
                    maximum_memory = max(maximum_memory, memory)
                    artifact_rows.append(
                        {
                            "source": "combined_naked_dna_mitochondrial",
                            "shift_forward": shift[0],
                            "shift_reverse": shift[1],
                            "model": model_name,
                            "context_length": spec.context_length,
                            "l2": float(l2),
                            "seed": int(seed),
                            "model_npz": str(npz_path),
                            "model_json": str(json_path),
                            "runtime_seconds": runtime,
                            "peak_memory_increment_mb": memory,
                            "model_size_mb": npz_path.stat().st_size / (1024 * 1024),
                        }
                    )
                started = perf_counter()
                ensemble = ensemble_sequence_bias_models(members)
                ensemble_stem = ensemble_dir / (
                    f"combined.{model_name}.shift_{shift[0]}_{shift[1]}."
                    f"l2_{safe_token(l2)}.seed_ensemble"
                )
                ensemble_npz, ensemble_json = ensemble.save(
                    ensemble_stem,
                    metadata={
                        "training_source": "combined_naked_dna_mitochondrial",
                        "configuration": "frozen_cross_control_pooled",
                        "read_shift": list(shift),
                        "l2": float(l2),
                        "member_seeds": [int(seed) for seed in seeds],
                        "member_models": [
                            {"path": str(path), "sha256": file_sha256(path)}
                            for path in member_paths
                        ],
                    },
                )
                ensemble_runtime = perf_counter() - started
                ensemble_rows.append(
                    {
                        "source": "combined_naked_dna_mitochondrial",
                        "shift_forward": shift[0],
                        "shift_reverse": shift[1],
                        "model": model_name,
                        "context_length": spec.context_length,
                        "l2": float(l2),
                        "member_count": len(members),
                        "member_seeds": ",".join(str(seed) for seed in seeds),
                        "model_npz": str(ensemble_npz),
                        "model_json": str(ensemble_json),
                        "training_runtime_seconds": total_runtime,
                        "ensemble_runtime_seconds": ensemble_runtime,
                        "peak_memory_increment_mb": maximum_memory,
                        "model_size_mb": ensemble_npz.stat().st_size / (1024 * 1024),
                    }
                )
                for name, (
                    candidate_contexts,
                    candidate_counts,
                ) in validation_arrays.items():
                    metric_rows.append(
                        {
                            "source": "combined_naked_dna_mitochondrial",
                            "library": name,
                            "shift_forward": shift[0],
                            "shift_reverse": shift[1],
                            "model": model_name,
                            "context_length": spec.context_length,
                            "l2": float(l2),
                            "model_npz": str(ensemble_npz),
                            **conditional_metrics(
                                ensemble, candidate_contexts, candidate_counts
                            ),
                        }
                    )
    return (
        pd.DataFrame(artifact_rows),
        pd.DataFrame(ensemble_rows),
        pd.DataFrame(metric_rows),
    )


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
    parser.add_argument("--model", choices=MODEL_NAMES, action="append", dest="models")
    parser.add_argument("--l2", type=float, action="append", dest="l2_values")
    parser.add_argument("--fit-seed", type=int, action="append", dest="seeds")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-windows", type=int, default=64)
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="fit deterministic seed members concurrently (research training only)",
    )
    args = parser.parse_args(argv)

    study = json.loads(args.study.read_text(encoding="utf-8"))
    if study.get("status") != "development_locked_holdout_unscored":
        raise ValueError("combined fitting requires the locked, unscored study")
    models = args.models or list(MODEL_NAMES)
    l2_values = args.l2_values or [1e-4, 1e-3, 1e-2]
    seeds = args.seeds or [2026, 2027, 2028, 2029, 2030]
    if seeds != list(study["random_seeds"]):
        raise ValueError("fit seeds must exactly match the frozen study")
    training = load_datasets(args.training_dataset, required_split="train")
    validation = load_datasets(args.validation_dataset, required_split="validation")
    args.outdir.mkdir(parents=True, exist_ok=True)
    artifacts, ensembles, metrics = fit_combined_grid(
        training,
        validation,
        args.outdir,
        models=models,
        l2_values=l2_values,
        seeds=seeds,
        epochs=args.epochs,
        batch_windows=args.batch_windows,
        jobs=args.jobs,
    )
    artifacts_path = args.outdir / "combined_bias_model_artifacts.tsv"
    ensembles_path = args.outdir / "combined_bias_model_ensembles.tsv"
    metrics_path = args.outdir / "combined_bias_model_metrics.tsv"
    artifacts.to_csv(artifacts_path, sep="\t", index=False)
    ensembles.to_csv(ensembles_path, sep="\t", index=False)
    metrics.to_csv(metrics_path, sep="\t", index=False)
    inputs = [path for _name, path in args.training_dataset + args.validation_dataset]
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
        "models": models,
        "l2_values": l2_values,
        "seeds": seeds,
        "jobs": args.jobs,
        "chipped_labels_used": False,
        "outputs": {
            path.name: {"path": str(path), "sha256": file_sha256(path)}
            for path in (artifacts_path, ensembles_path, metrics_path)
        },
        "input_count": len(inputs),
    }
    (args.outdir / "combined_bias_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(metrics.sort_values("conditional_nll").head(12).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
