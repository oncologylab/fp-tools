#!/usr/bin/env python3
"""Evaluate shared and antisymmetric strand footprint shapes per TF.

Artifacts are constructed without labels by ``build_strand_functional_profiles``.
This development-only evaluator then measures same-cell ceilings and cross-cell
or leave-TF-out transfer using only functional cut-shape channels.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import site_hashes  # noqa: E402
from evaluate_functional_footprints import (  # noqa: E402
    binary_metrics,
    chromosome_split,
    stable_seed,
)
from evaluate_functional_template_transfer import (  # noqa: E402
    TRAINING_SCOPES,
    balanced_training_indexes,
    selection_score,
    training_indexes,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    FunctionalTemplateDetector,
    MultichannelFunctionalTemplateDetector,
)


CHANNEL_SETS = {
    "combined": ("combined_residual",),
    "shared": ("shared_strand_residual",),
    "antisymmetric": ("antisymmetric_strand_residual",),
    "shared_antisymmetric": ("shared_strand_residual", "antisymmetric_strand_residual"),
    "all": (
        "combined_residual",
        "shared_strand_residual",
        "antisymmetric_strand_residual",
    ),
}


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_artifact(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("artifact must use MODEL,CELL,JSON")
    return fields[0], fields[1], Path(fields[2])


def load_artifact(path: Path, expected_cell: str, study: dict) -> tuple[pd.DataFrame, dict[str, np.ndarray], dict]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != "fp-tools-strand-functional-profiles-v1":
        raise ValueError(f"unsupported strand artifact schema: {path}")
    if document.get("metadata", {}).get("labels_used") is not False:
        raise ValueError(f"strand artifact does not certify label-free construction: {path}")
    profiles_path = Path(document["profiles_npz"])
    sites_path = Path(document["sites"])
    if file_sha256(profiles_path) != document["profiles_sha256"]:
        raise ValueError(f"profile checksum mismatch: {profiles_path}")
    if file_sha256(sites_path) != document["sites_sha256"]:
        raise ValueError(f"site checksum mismatch: {sites_path}")
    sites = pd.read_csv(sites_path, sep="\t").reset_index(drop=True)
    if "cell" not in sites or set(sites["cell"].astype(str)) != {expected_cell}:
        raise ValueError(f"artifact sites do not exclusively contain {expected_cell}: {path}")
    sites["chromosome_split"] = sites["TFBS_chr"].map(
        lambda chromosome: chromosome_split(str(chromosome), study)
    )
    with np.load(profiles_path, allow_pickle=False) as arrays:
        if not np.array_equal(arrays["site_hash"], site_hashes(sites)):
            raise ValueError(f"site order hash mismatch: {path}")
        valid = np.asarray(arrays["valid"], dtype=bool)
        profiles = {
            name: np.asarray(arrays[name], dtype=np.float64)
            for names in CHANNEL_SETS.values()
            for name in names
        }
    eligible = valid & sites["chromosome_split"].isin(["train", "validation"]).to_numpy()
    return sites.loc[eligible].reset_index(drop=True), {
        name: values[eligible] for name, values in profiles.items()
    }, document


def stack_channels(profiles: dict[str, np.ndarray], channel_set: str) -> np.ndarray:
    names = CHANNEL_SETS[channel_set]
    return np.stack([profiles[name] for name in names], axis=1)


def _evaluate_one(
    *,
    sites: pd.DataFrame,
    cell_profiles: dict[str, np.ndarray],
    target_cell: str,
    target_tf: str,
    target_family: str,
    bias_configuration: str,
    channel_set: str,
    scope: str,
    smoother: str,
    window: float,
    positions: np.ndarray,
    minimum_sites: int,
    maximum_train_per_tf_class: int,
    seed: int,
) -> tuple[dict, pd.DataFrame | None]:
    started = perf_counter()
    validation = np.flatnonzero(
        (sites["cell"].astype(str).to_numpy() == target_cell)
        & (sites["tf"].astype(str).to_numpy() == target_tf)
        & (sites["chromosome_split"].astype(str).to_numpy() == "validation")
    )
    train = training_indexes(
        sites,
        target_cell=target_cell,
        target_tf=target_tf,
        target_family=target_family,
        scope=scope,
    )
    train = balanced_training_indexes(
        sites,
        train,
        maximum_per_tf_class=maximum_train_per_tf_class,
        seed=stable_seed(
            bias_configuration,
            channel_set,
            target_cell,
            target_tf,
            scope,
            smoother,
            window,
            seed=seed,
        ),
    )
    train_labels = sites.iloc[train]["chip_label"].to_numpy(dtype=int)
    validation_labels = sites.iloc[validation]["chip_label"].to_numpy(dtype=int)
    base = {
        "cell": target_cell,
        "tf": target_tf,
        "motif_family": target_family,
        "bias_configuration": bias_configuration,
        "channel_set": channel_set,
        "training_scope": scope,
        "smoother": smoother,
        "window_limit": float(window),
        "training_labels_used": True,
        "motif_or_accessibility_features_used": False,
        "train_sites": int(len(train)),
        "validation_sites": int(len(validation)),
        "train_positive_sites": int(np.sum(train_labels == 1)),
        "train_negative_sites": int(np.sum(train_labels == 0)),
        "validation_positive_sites": int(np.sum(validation_labels == 1)),
        "validation_negative_sites": int(np.sum(validation_labels == 0)),
    }
    if (
        min(np.sum(train_labels == 0), np.sum(train_labels == 1)) < minimum_sites
        or min(np.sum(validation_labels == 0), np.sum(validation_labels == 1)) < minimum_sites
    ):
        return {**base, "status": "insufficient_sites"}, None

    train_parts = []
    train_label_parts = []
    for cell in sorted(sites.iloc[train]["cell"].astype(str).unique()):
        global_cell_rows = np.flatnonzero(sites["cell"].astype(str).to_numpy() == cell)
        lookup = pd.Series(np.arange(len(global_cell_rows)), index=global_cell_rows)
        selected_global = train[sites.iloc[train]["cell"].astype(str).to_numpy() == cell]
        train_parts.append(cell_profiles[cell][lookup.loc[selected_global].to_numpy(dtype=int)])
        train_label_parts.append(sites.iloc[selected_global]["chip_label"].to_numpy(dtype=int))
    train_profiles = np.vstack(train_parts)
    train_labels = np.concatenate(train_label_parts)
    target_global_rows = np.flatnonzero(sites["cell"].astype(str).to_numpy() == target_cell)
    target_lookup = pd.Series(np.arange(len(target_global_rows)), index=target_global_rows)
    validation_profiles = cell_profiles[target_cell][target_lookup.loc[validation].to_numpy(dtype=int)]

    if train_profiles.shape[1] == 1:
        detector = FunctionalTemplateDetector(
            positions,
            smoother=smoother,
            window_limit=window,
        ).fit(train_profiles[:, 0, :], train_labels)
        probabilities = detector.predict_proba(validation_profiles[:, 0, :])
        templates = np.asarray([detector.footprint_template_])
        standard_errors = np.asarray([detector.template_standard_error_])
        channel_weights = np.ones(1)
    else:
        detector = MultichannelFunctionalTemplateDetector(
            positions,
            smoother=smoother,
            window_limit=window,
        ).fit(train_profiles, train_labels)
        probabilities = detector.predict_proba(validation_profiles)
        templates = np.asarray(
            [model.footprint_template_ for model in detector.channel_models_]
        )
        standard_errors = np.asarray(
            [model.template_standard_error_ for model in detector.channel_models_]
        )
        channel_weights = detector.channel_discriminant_
    metrics = binary_metrics(validation_labels, probabilities)
    row = {
        **base,
        "status": "ok",
        "fit_seconds": perf_counter() - started,
        "selection_score": selection_score(metrics),
        "channel_weights": ",".join(f"{value:.8g}" for value in channel_weights),
        **metrics,
    }
    names = CHANNEL_SETS[channel_set]
    curves = []
    for index, name in enumerate(names):
        curves.append(
            pd.DataFrame(
                {
                    "cell": target_cell,
                    "tf": target_tf,
                    "motif_family": target_family,
                    "bias_configuration": bias_configuration,
                    "channel_set": channel_set,
                    "channel": name,
                    "training_scope": scope,
                    "smoother": smoother,
                    "window_limit": float(window),
                    "position": positions.astype(int),
                    "template": templates[index],
                    "template_standard_error": standard_errors[index],
                    "channel_weight": float(channel_weights[index]),
                }
            )
        )
    return row, pd.concat(curves, ignore_index=True)


def select_winners(metrics: pd.DataFrame) -> pd.DataFrame:
    passing = metrics[metrics["status"] == "ok"].copy()
    keys = ["cell", "tf", "bias_configuration", "training_scope"]
    if passing.empty:
        return passing
    return (
        passing.sort_values(
            keys + ["selection_score", "auprc", "auroc"],
            ascending=[True] * len(keys) + [False, False, False],
            kind="mergesort",
        )
        .groupby(keys, sort=True, as_index=False)
        .head(1)
        .reset_index(drop=True)
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact",
        action="append",
        type=parse_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"),
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--channel-sets", nargs="+", choices=tuple(CHANNEL_SETS), default=list(CHANNEL_SETS))
    parser.add_argument("--training-scopes", nargs="+", choices=TRAINING_SCOPES, default=list(TRAINING_SCOPES))
    parser.add_argument("--smoothers", nargs="+", choices=("spline", "gp"), default=["spline", "gp"])
    parser.add_argument("--windows", nargs="+", type=float, default=[30.0, 50.0, 80.0])
    parser.add_argument("--minimum-sites-per-class", type=int, default=100)
    parser.add_argument("--maximum-train-per-tf-class", type=int, default=2000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == "development"]
    loaded: dict[str, dict[str, tuple[pd.DataFrame, dict[str, np.ndarray], dict]]] = {}
    input_manifest = []
    for model, cell, path in args.artifact:
        if cell in loaded.setdefault(model, {}):
            raise SystemExit(f"duplicate artifact for {model}/{cell}")
        loaded[model][cell] = load_artifact(path, cell, study)
        input_manifest.append(
            {"model": model, "cell": cell, "path": str(path), "sha256": file_sha256(path)}
        )
    required_cells = set(tasks["cell"].astype(str))
    incomplete = {model: required_cells.difference(cells) for model, cells in loaded.items()}
    incomplete = {model: missing for model, missing in incomplete.items() if missing}
    if incomplete:
        raise SystemExit("artifact models are missing development cells: " + str(incomplete))

    rows = []
    curves = []
    futures = {}
    positions = np.arange(-int(study["profile_flank_bp"]), int(study["profile_flank_bp"]) + 1, dtype=float)
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=args.workers) as executor:
        for model, cells in loaded.items():
            global_sites = pd.concat(
                [cells[cell][0] for cell in sorted(cells)],
                ignore_index=True,
            )
            for channel_set in args.channel_sets:
                cell_profiles = {
                    cell: stack_channels(cells[cell][1], channel_set)
                    for cell in cells
                }
                for task in tasks.itertuples(index=False):
                    for scope in args.training_scopes:
                        for smoother in args.smoothers:
                            for window in args.windows:
                                future = executor.submit(
                                    _evaluate_one,
                                    sites=global_sites,
                                    cell_profiles=cell_profiles,
                                    target_cell=str(task.cell),
                                    target_tf=str(task.tf),
                                    target_family=str(task.motif_family),
                                    bias_configuration=model,
                                    channel_set=channel_set,
                                    scope=scope,
                                    smoother=smoother,
                                    window=window,
                                    positions=positions,
                                    minimum_sites=args.minimum_sites_per_class,
                                    maximum_train_per_tf_class=args.maximum_train_per_tf_class,
                                    seed=args.seed,
                                )
                                futures[future] = (model, str(task.cell), str(task.tf), channel_set, scope)
        for future in as_completed(futures):
            identity = futures[future]
            try:
                row, curve = future.result()
            except Exception as error:
                model, cell, tf, channel_set, scope = identity
                row = {
                    "bias_configuration": model,
                    "cell": cell,
                    "tf": tf,
                    "channel_set": channel_set,
                    "training_scope": scope,
                    "status": "error",
                    "error": f"{type(error).__name__}: {error}",
                }
                curve = None
            rows.append(row)
            if curve is not None:
                curves.append(curve)
    metrics = pd.DataFrame(rows)
    profiles = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    winners = select_winners(metrics)
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.outdir / "strand_functional_metrics.tsv.gz"
    profiles_path = args.outdir / "strand_functional_profiles.tsv.gz"
    winners_path = args.outdir / "strand_functional_winners.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    profiles.to_csv(profiles_path, sep="\t", index=False)
    winners.to_csv(winners_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-strand-functional-template-evaluation-v1",
        "locked_test_labels_read": False,
        "selection_split": "validation",
        "profile_construction_labels_used": False,
        "template_training_labels_used": True,
        "motif_or_accessibility_features_used": False,
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "artifacts": input_manifest,
        "channel_sets": args.channel_sets,
        "training_scopes": args.training_scopes,
        "smoothers": args.smoothers,
        "windows": args.windows,
        "minimum_sites_per_class": args.minimum_sites_per_class,
        "maximum_train_per_tf_class": args.maximum_train_per_tf_class,
        "workers": args.workers,
        "seed": args.seed,
        "metrics_rows": int(len(metrics)),
        "winner_rows": int(len(winners)),
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {"path": str(profiles_path), "sha256": file_sha256(profiles_path)},
            "winners": {"path": str(winners_path), "sha256": file_sha256(winners_path)},
        },
    }
    (args.outdir / "strand_functional_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    columns = [
        "cell",
        "tf",
        "bias_configuration",
        "training_scope",
        "channel_set",
        "smoother",
        "window_limit",
        "auroc",
        "auprc",
    ]
    print(winners[columns].to_string(index=False) if len(winners) else "no eligible winners")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
