#!/usr/bin/env python3
"""Render blinded per-TF aggregate curves for DWM and strand candidates."""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_depth_matrix import prepare_profile  # noqa: E402
from evaluate_functional_footprints import (  # noqa: E402
    _evaluation_profiles,
    chromosome_split,
    validate_sites,
)
from evaluate_strand_functional_templates import (  # noqa: E402
    CHANNEL_SETS,
    load_artifact,
    stack_channels,
)
from fp_tools.tools.functional_footprints import normalize_functional_profiles  # noqa: E402
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*values: object, seed: int = 2026) -> int:
    digest = blake2b(digest_size=8)
    for value in (seed, *values):
        digest.update(str(value).encode())
        digest.update(b"\0")
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def mean_band(profiles: np.ndarray, *, bootstraps: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(profiles, dtype=float)
    mean = np.mean(values, axis=0)
    if len(values) < 2 or bootstraps < 2:
        return mean, mean.copy(), mean.copy()
    rng = np.random.default_rng(seed)
    sampled = np.empty((bootstraps, values.shape[1]), dtype=np.float32)
    for index in range(bootstraps):
        sampled[index] = np.mean(values[rng.integers(0, len(values), len(values))], axis=0)
    lower, upper = np.quantile(sampled, [0.025, 0.975], axis=0)
    return mean, lower, upper


def combined_strand_shape(
    profiles: dict[str, np.ndarray],
    winner: pd.Series,
    positions: np.ndarray,
) -> np.ndarray:
    values = stack_channels(profiles, str(winner.channel_set))
    normalized = np.stack(
        [normalize_functional_profiles(values[:, channel, :], positions) for channel in range(values.shape[1])],
        axis=1,
    )
    weights = np.fromstring(str(winner.channel_weights), sep=",")
    if len(weights) != values.shape[1]:
        raise ValueError("winner channel weights do not match channel set")
    weights = weights / max(float(np.sum(np.abs(weights))), np.finfo(float).eps)
    return np.einsum("ncp,c->np", normalized, weights)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-run", type=Path, required=True)
    parser.add_argument("--dwm-winners", type=Path, required=True)
    parser.add_argument("--strand-run", type=Path, required=True)
    parser.add_argument("--models", nargs="+", default=["MT_LOG81_4m5", "MT_SELMA10_4m4"])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--key-out", type=Path)
    parser.add_argument("--summary-out", type=Path)
    parser.add_argument("--bootstraps", type=int, default=200)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    if args.bootstraps < 2:
        raise SystemExit("--bootstraps must be at least two")
    base_manifest = json.loads((args.base_run / "functional_benchmark_manifest.json").read_text(encoding="utf-8"))
    if base_manifest.get("test_unlocked"):
        raise SystemExit("aggregate renderer refuses a base run that opened test labels")
    study_path = Path(base_manifest["study"])
    study = json.loads(study_path.read_text(encoding="utf-8"))
    sites_path = Path(base_manifest["development_sites"])
    sites = validate_sites(pd.read_csv(sites_path, sep="\t"), sites_path)
    sites["chromosome_split"] = sites["TFBS_chr"].map(
        lambda chromosome: chromosome_split(str(chromosome), study)
    )
    sites = sites[sites["chromosome_split"].isin(["train", "validation"])].reset_index(drop=True)
    tracks_path = Path(base_manifest["tracks"])
    tracks = pd.read_csv(tracks_path, sep="\t")
    genome = Path(base_manifest["genome"]) if base_manifest.get("genome") else None
    flank = int(base_manifest["flank"])
    positions = np.arange(-flank, flank + 1, dtype=float)
    dwm_winners = pd.read_csv(args.dwm_winners, sep="\t")
    dwm_winners = dwm_winners[
        (dwm_winners["correction"] == "DWM")
        & (dwm_winners["training_scope"] == "same_cell_ceiling")
    ].copy()
    strand_manifest_path = args.strand_run / "strand_functional_manifest.json"
    strand_manifest = json.loads(strand_manifest_path.read_text(encoding="utf-8"))
    strand_winners = pd.read_csv(args.strand_run / "strand_functional_winners.tsv", sep="\t")
    strand_winners = strand_winners[
        strand_winners["training_scope"].eq("same_cell_ceiling")
        & strand_winners["bias_configuration"].isin(args.models)
    ].copy()
    artifact_paths = {
        (str(row["model"]), str(row["cell"])): Path(row["path"])
        for row in strand_manifest["artifacts"]
    }
    artifacts = {}
    for model in args.models:
        for cell in sorted(strand_winners[strand_winners["bias_configuration"] == model]["cell"].unique()):
            key = (model, str(cell))
            if key not in artifact_paths:
                raise SystemExit(f"strand run lacks artifact for {model}/{cell}")
            artifacts[key] = load_artifact(artifact_paths[key], str(cell), study)

    dwm_profiles = {}
    for cell in sorted(dwm_winners["cell"].unique()):
        cell_sites = sites[sites["cell"].astype(str) == cell].reset_index(drop=True)
        observed, expected = _evaluation_profiles(
            cell_sites,
            tracks,
            cell,
            "DWM",
            args.base_run / "profile_cache",
            "development",
            flank,
            genome,
        )
        dwm_profiles[cell] = (cell_sites, observed, expected, estimate_nb_dispersion(observed, expected))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    args.out.parent.mkdir(parents=True, exist_ok=True)
    key_out = args.key_out or args.out.with_name(args.out.stem + "_blinding_key.tsv")
    summary_out = args.summary_out or args.out.with_name(args.out.stem + "_summary.tsv")
    key_rows = []
    summary_rows = []
    methods = ["DWM", *args.models]
    tasks = sorted(set(zip(dwm_winners["cell"].astype(str), dwm_winners["tf"].astype(str))))
    with PdfPages(args.out) as pdf:
        for cell, tf in tasks:
            figure, axes = plt.subplots(1, len(methods), figsize=(5.0 * len(methods), 4.2), squeeze=False, sharey=True)
            swap = bool(stable_seed(cell, tf, "blind", seed=args.seed) % 2)
            group_names = {0: "Group B" if swap else "Group A", 1: "Group A" if swap else "Group B"}
            key_rows.extend(
                [
                    {"cell": cell, "tf": tf, "blinded_group": group_names[label], "actual_group": "chip_positive" if label else "matched_negative"}
                    for label in (0, 1)
                ]
            )
            for axis, method in zip(axes[0], methods):
                if method == "DWM":
                    candidate = dwm_winners[(dwm_winners["cell"] == cell) & (dwm_winners["tf"] == tf)]
                    if candidate.empty:
                        axis.set_visible(False)
                        continue
                    winner = candidate.iloc[0]
                    cell_sites, observed, expected, dispersion = dwm_profiles[cell]
                    residual, _ = prepare_profile(
                        observed,
                        expected,
                        positions,
                        residual_mode=str(winner.residual),
                        background=str(winner.background),
                        dispersion=dispersion,
                    )
                    mask = (cell_sites["tf"].astype(str) == tf) & (cell_sites["chromosome_split"] == "validation")
                    profiles = normalize_functional_profiles(residual[mask.to_numpy()], positions)
                    labels = cell_sites.loc[mask, "chip_label"].to_numpy(dtype=int)
                    auc, ap = float(winner.auroc), float(winner.auprc)
                    detail = f"{winner.residual}, {winner.background}"
                else:
                    candidate = strand_winners[
                        (strand_winners["cell"] == cell)
                        & (strand_winners["tf"] == tf)
                        & (strand_winners["bias_configuration"] == method)
                    ]
                    if candidate.empty:
                        axis.set_visible(False)
                        continue
                    winner = candidate.iloc[0]
                    artifact_sites, artifact_profiles, _document = artifacts[(method, cell)]
                    mask = (artifact_sites["tf"].astype(str) == tf) & (artifact_sites["chromosome_split"] == "validation")
                    profiles = combined_strand_shape(
                        {name: values[mask.to_numpy()] for name, values in artifact_profiles.items()},
                        winner,
                        positions,
                    )
                    labels = artifact_sites.loc[mask, "chip_label"].to_numpy(dtype=int)
                    auc, ap = float(winner.auroc), float(winner.auprc)
                    detail = str(winner.channel_set)
                colors = {0: "#2A6FBB", 1: "#C23B33"}
                for label in (0, 1):
                    mean, lower, upper = mean_band(
                        profiles[labels == label],
                        bootstraps=args.bootstraps,
                        seed=stable_seed(cell, tf, method, label, seed=args.seed),
                    )
                    axis.plot(positions, mean, color=colors[label], linewidth=1.7, label=group_names[label])
                    axis.fill_between(positions, lower, upper, color=colors[label], alpha=0.18, linewidth=0)
                difference = np.mean(profiles[labels == 1], axis=0) - np.mean(profiles[labels == 0], axis=0)
                separation = float(np.sqrt(np.mean(np.square(difference[np.abs(positions) <= 50]))))
                summary_rows.append(
                    {
                        "cell": cell,
                        "tf": tf,
                        "method": method,
                        "auroc": auc,
                        "auprc": ap,
                        "visual_rms_difference": separation,
                        "detail": detail,
                    }
                )
                axis.axvline(0, color="#666666", linewidth=0.7, linestyle="--")
                axis.axhline(0, color="#999999", linewidth=0.6)
                axis.set_title(f"{method}\nAUROC {auc:.3f}; AUPRC {ap:.3f}\n{detail}", fontsize=9)
                axis.set_xlabel("Position from motif center (bp)")
            axes[0, 0].set_ylabel("Normalized functional residual")
            handles, labels_text = axes[0, 0].get_legend_handles_labels()
            if handles:
                figure.legend(handles, labels_text, loc="upper center", ncol=2, frameon=False)
            figure.suptitle(f"{cell} — {tf}: blinded matched aggregate profiles", y=1.03)
            figure.tight_layout()
            pdf.savefig(figure, bbox_inches="tight")
            plt.close(figure)
    pd.DataFrame(key_rows).drop_duplicates().to_csv(key_out, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(summary_out, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-functional-aggregate-comparison-v1",
        "locked_test_labels_read": False,
        "blinded": True,
        "base_run": str(args.base_run),
        "dwm_winners": str(args.dwm_winners),
        "dwm_winners_sha256": file_sha256(args.dwm_winners),
        "strand_manifest": str(strand_manifest_path),
        "strand_manifest_sha256": file_sha256(strand_manifest_path),
        "models": args.models,
        "bootstraps": args.bootstraps,
        "seed": args.seed,
        "outputs": {
            "pdf": {"path": str(args.out), "sha256": file_sha256(args.out)},
            "blinding_key": {"path": str(key_out), "sha256": file_sha256(key_out)},
            "summary": {"path": str(summary_out), "sha256": file_sha256(summary_out)},
        },
    }
    args.out.with_suffix(".json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(pd.DataFrame(summary_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
