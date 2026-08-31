#!/usr/bin/env python3
"""Search TF-specific footprint geometries without touching locked test data.

The search is deliberately staged.  A broad, inexpensive geometry screen is
performed on training chromosomes, the best geometries are expanded across
normalization and symmetry choices, and validation chromosomes are used only
to select among the resulting shortlist.  Test chromosomes are not read.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from fp_tools.utils import bigwig as pyBigWig


SITE_COLUMNS = {
    "cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end", "chip_label",
}
SIGNAL_COLUMNS = {"cell", "correction", "signal"}


@dataclass(frozen=True)
class Candidate:
    correction: str
    center_width: int
    flank_width: int
    gap: int
    shoulder: str = "mean"
    center: str = "mean"
    normalization: str = "none"
    asymmetry_penalty: float = 0.0

    @property
    def identifier(self) -> str:
        return (
            f"{self.correction}.c{self.center_width}.f{self.flank_width}.g{self.gap}"
            f".{self.shoulder}.{self.center}.{self.normalization}"
            f".a{self.asymmetry_penalty:g}"
        )


def read_sites(paths: Iterable[Path]) -> pd.DataFrame:
    frame = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    missing = sorted(SITE_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError("site table is missing columns: " + ", ".join(missing))
    frame = frame.copy()
    frame["chip_label"] = pd.to_numeric(frame["chip_label"], errors="raise").astype(int)
    if not set(frame["chip_label"].unique()).issubset({0, 1}):
        raise ValueError("chip_label must contain only 0 and 1")
    key = ["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end"]
    return frame.drop_duplicates(key, keep="first").reset_index(drop=True)


def read_signals(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    if "correction" not in frame and "method" in frame:
        frame = frame.rename(columns={"method": "correction"})
    missing = sorted(SIGNAL_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError("signal manifest is missing columns: " + ", ".join(missing))
    if frame.duplicated(["cell", "correction"]).any():
        raise ValueError("signal manifest has duplicate cell/correction rows")
    missing_files = [str(value) for value in frame["signal"] if not Path(value).is_file()]
    if missing_files:
        raise FileNotFoundError("signal files do not exist: " + ", ".join(missing_files[:3]))
    return frame


def attach_chromosome_split(sites: pd.DataFrame, study: dict[str, object]) -> pd.DataFrame:
    by_chromosome = {
        chromosome: split
        for split, chromosomes in study["chromosome_split"].items()
        for chromosome in chromosomes
    }
    output = sites.copy()
    output["chromosome_split"] = output["TFBS_chr"].map(by_chromosome).fillna("excluded")
    return output


def deterministic_class_sample(
    sites: pd.DataFrame,
    maximum_per_class: int,
    seed: int,
    negative_pool_multiplier: int = 1,
) -> pd.DataFrame:
    """Cap each cell/TF/split/label stratum reproducibly."""

    if maximum_per_class <= 0:
        return sites.reset_index(drop=True)
    key = ["cell", "tf", "chromosome_split", "chip_label"]
    ranked = sites.copy()
    hashes = pd.util.hash_pandas_object(
        ranked[["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end"]],
        index=False,
    ).to_numpy(dtype=np.uint64)
    ranked["_sample_order"] = hashes ^ np.uint64(seed)
    ranked = ranked.sort_values(key + ["_sample_order"], kind="mergesort")
    ranked["_maximum"] = np.where(
        ranked["chip_label"].to_numpy(dtype=int) == 0,
        maximum_per_class * max(1, int(negative_pool_multiplier)),
        maximum_per_class,
    )
    ranked["_rank"] = ranked.groupby(key, sort=False).cumcount()
    ranked = ranked[ranked["_rank"] < ranked["_maximum"]]
    return ranked.drop(columns=["_sample_order", "_maximum", "_rank"]).sort_index().reset_index(drop=True)


def write_merged_regions(sites: pd.DataFrame, path: Path, padding: int) -> int:
    """Write merged windows that cover all sampled motif sites."""

    intervals = []
    for row in sites.itertuples(index=False):
        center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
        intervals.append((str(row.TFBS_chr), max(0, center - padding), center + padding + 1))
    intervals.sort(key=lambda item: (item[0], item[1], item[2]))
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            previous = merged[-1]
            merged[-1] = (chrom, previous[1], max(previous[2], end))
        else:
            merged.append((chrom, start, end))
    with path.open("w", encoding="utf-8") as handle:
        for index, (chrom, start, end) in enumerate(merged, start=1):
            handle.write(f"{chrom}\t{start}\t{end}\tevaluation_region_{index}\n")
    return len(merged)


def extract_profiles(sites: pd.DataFrame, signal: Path, flank: int) -> tuple[np.ndarray, np.ndarray]:
    """Extract fixed-width motif-centered profiles and a finite-row mask."""

    width = flank * 2 + 1
    profiles = np.full((len(sites), width), np.nan, dtype=np.float32)
    handle = pyBigWig.open(str(signal))
    try:
        chromosomes = handle.chroms()
        for index, row in enumerate(sites.itertuples(index=False)):
            chrom = str(row.TFBS_chr)
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            start, end = center - flank, center + flank + 1
            if chrom not in chromosomes or start < 0 or end > int(chromosomes[chrom]):
                continue
            values = np.asarray(handle.values(chrom, start, end, numpy=True), dtype=float)
            if len(values) != width:
                continue
            # A missing bigWig interval means no recorded cut, not a missing site.
            profiles[index] = np.nan_to_num(values, nan=0.0).astype(np.float32)
    finally:
        handle.close()
    valid = np.isfinite(profiles).all(axis=1)
    return profiles, valid


def _segment(profile: np.ndarray, middle: int, left: int, right: int) -> np.ndarray:
    return profile[:, middle + left:middle + right]


def _summary(values: np.ndarray, statistic: str) -> np.ndarray:
    if statistic == "mean":
        return np.mean(values, axis=1)
    if statistic == "minimum":
        return np.min(values, axis=1)
    if statistic == "q25":
        return np.quantile(values, 0.25, axis=1)
    raise ValueError(f"unknown summary statistic: {statistic}")


def score_candidate(profiles: np.ndarray, candidate: Candidate) -> np.ndarray:
    """Calculate one TF-shape hypothesis for every profile."""

    middle = profiles.shape[1] // 2
    left_center = -(candidate.center_width // 2)
    right_center = left_center + candidate.center_width
    left = _segment(
        profiles,
        middle,
        left_center - candidate.gap - candidate.flank_width,
        left_center - candidate.gap,
    ).mean(axis=1)
    right = _segment(
        profiles,
        middle,
        right_center + candidate.gap,
        right_center + candidate.gap + candidate.flank_width,
    ).mean(axis=1)
    center_values = _segment(profiles, middle, left_center, right_center)
    center = _summary(center_values, candidate.center)
    if candidate.shoulder == "mean":
        shoulder = (left + right) / 2.0
    elif candidate.shoulder == "minimum":
        shoulder = np.minimum(left, right)
    elif candidate.shoulder == "maximum":
        shoulder = np.maximum(left, right)
    else:
        raise ValueError(f"unknown shoulder statistic: {candidate.shoulder}")
    score = shoulder - center - candidate.asymmetry_penalty * np.abs(left - right)

    span_left = left_center - candidate.gap - candidate.flank_width
    span_right = right_center + candidate.gap + candidate.flank_width
    local = _segment(profiles, middle, span_left, span_right)
    epsilon = max(float(np.nanmedian(np.abs(local))) * 1e-6, 1e-8)
    if candidate.normalization == "none":
        scale = 1.0
    elif candidate.normalization == "sqrt_abs":
        scale = np.sqrt(np.mean(np.abs(local), axis=1) + epsilon)
    elif candidate.normalization == "rms":
        scale = np.sqrt(np.mean(np.square(local), axis=1) + epsilon)
    elif candidate.normalization == "flank_sd":
        flanks = np.concatenate(
            [
                _segment(profiles, middle, span_left, left_center - candidate.gap),
                _segment(profiles, middle, right_center + candidate.gap, span_right),
            ],
            axis=1,
        )
        scale = np.std(flanks, axis=1) + epsilon
    elif candidate.normalization == "flank_mad":
        flanks = np.concatenate(
            [
                _segment(profiles, middle, span_left, left_center - candidate.gap),
                _segment(profiles, middle, right_center + candidate.gap, span_right),
            ],
            axis=1,
        )
        median = np.median(flanks, axis=1, keepdims=True)
        scale = 1.4826 * np.median(np.abs(flanks - median), axis=1) + epsilon
    else:
        raise ValueError(f"unknown normalization: {candidate.normalization}")
    return np.asarray(score / scale, dtype=float)


def binary_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float | int]:
    finite = np.isfinite(scores)
    labels, scores = labels[finite], scores[finite]
    if len(labels) == 0 or len(np.unique(labels)) != 2:
        return {"n_sites": int(len(labels)), "positive_sites": int(labels.sum()), "auroc": np.nan, "auprc": np.nan}
    return {
        "n_sites": int(len(labels)),
        "positive_sites": int(labels.sum()),
        "auroc": float(roc_auc_score(labels, scores)),
        "auprc": float(average_precision_score(labels, scores)),
    }


def geometry_grid(correction: str, flank: int, quick: bool = False) -> list[Candidate]:
    centers = [7, 11, 15, 21, 27, 33, 41] if quick else [5, 7, 9, 11, 13, 15, 17, 21, 25, 29, 33, 41, 49]
    shoulders = [8, 16, 32, 48] if quick else [4, 8, 12, 16, 24, 32, 48, 64]
    gaps = [0, 4, 8] if quick else [0, 2, 4, 8, 12]
    candidates = []
    for center_width in centers:
        for flank_width in shoulders:
            for gap in gaps:
                if center_width // 2 + gap + flank_width > flank:
                    continue
                for shoulder in ("mean", "minimum"):
                    candidates.append(Candidate(correction, center_width, flank_width, gap, shoulder))
    return candidates


def expanded_grid(candidate: Candidate) -> list[Candidate]:
    return [
        replace(
            candidate,
            center=center,
            normalization=normalization,
            asymmetry_penalty=penalty,
        )
        for center in ("mean", "minimum", "q25")
        for normalization in ("none", "sqrt_abs", "rms", "flank_sd", "flank_mad")
        for penalty in (0.0, 0.25, 0.5, 1.0)
    ]


def evaluate_candidates(
    profiles: np.ndarray,
    labels: np.ndarray,
    candidates: Iterable[Candidate],
    *,
    cell: str,
    tf: str,
    split: str,
    stage: str,
) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        metrics = binary_metrics(labels, score_candidate(profiles, candidate))
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "chromosome_split": split,
                "stage": stage,
                "candidate": candidate.identifier,
                **asdict(candidate),
                **metrics,
            }
        )
    return pd.DataFrame(rows)


def shortlist(metrics: pd.DataFrame, count: int) -> pd.DataFrame:
    ranked = metrics.copy()
    prevalence = ranked["positive_sites"] / ranked["n_sites"].clip(lower=1)
    # AUROC and chance-adjusted AUPRC have equal, bounded influence.
    adjusted_auprc = (ranked["auprc"] - prevalence) / (1.0 - prevalence).clip(lower=1e-6)
    ranked["selection_score"] = ranked["auroc"] + adjusted_auprc
    return ranked.sort_values(
        ["selection_score", "auprc", "auroc"], ascending=False, kind="mergesort"
    ).head(count)


def candidate_from_row(row) -> Candidate:
    return Candidate(
        correction=str(row.correction),
        center_width=int(row.center_width),
        flank_width=int(row.flank_width),
        gap=int(row.gap),
        shoulder=str(row.shoulder),
        center=str(row.center),
        normalization=str(row.normalization),
        asymmetry_penalty=float(row.asymmetry_penalty),
    )


def run_search(
    sites: pd.DataFrame,
    signals: pd.DataFrame,
    outdir: Path,
    flank: int,
    stage1_keep: int,
    validation_keep: int,
    quick: bool,
    profiles_only: bool = False,
    profile_cache_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_metrics = []
    winners = []
    cache_dir = profile_cache_dir or (outdir / "profile_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        cell_signals = signals[signals["cell"] == cell]
        if cell_signals.empty:
            continue
        profiles_by_correction = {}
        valid_by_correction = {}
        for signal_row in cell_signals.itertuples(index=False):
            cache = cache_dir / f"{cell}.{signal_row.correction}.flank{flank}.npz"
            site_hash = pd.util.hash_pandas_object(
                cell_sites[["cell", "tf", "TFBS_chr", "TFBS_start", "TFBS_end"]],
                index=False,
            ).to_numpy(dtype=np.uint64)
            if cache.is_file():
                payload = np.load(cache)
                profiles = payload["profiles"]
                valid = payload["valid"].astype(bool)
                cached_hash = payload["site_hash"] if "site_hash" in payload else np.array([])
                if len(profiles) != len(cell_sites) or (
                    len(cached_hash) and not np.array_equal(cached_hash, site_hash)
                ):
                    raise ValueError(f"stale profile cache does not match sampled sites: {cache}")
            else:
                profiles, valid = extract_profiles(cell_sites, Path(signal_row.signal), flank)
                np.savez_compressed(cache, profiles=profiles, valid=valid, site_hash=site_hash)
            profiles_by_correction[str(signal_row.correction)] = profiles
            valid_by_correction[str(signal_row.correction)] = valid

        if profiles_only:
            continue

        common = np.logical_and.reduce(list(valid_by_correction.values()))
        for tf, tf_sites in cell_sites.groupby("tf", sort=True):
            positions = tf_sites.index.to_numpy()
            positions = positions[common[positions]]
            selected_sites = cell_sites.iloc[positions]
            train_mask = selected_sites["chromosome_split"].to_numpy() == "train"
            validation_mask = selected_sites["chromosome_split"].to_numpy() == "validation"
            train_labels = selected_sites.loc[train_mask, "chip_label"].to_numpy(dtype=int)
            validation_labels = selected_sites.loc[validation_mask, "chip_label"].to_numpy(dtype=int)
            if len(np.unique(train_labels)) != 2 or len(np.unique(validation_labels)) != 2:
                continue

            if "footprint_score" in selected_sites:
                for split, mask, labels in (
                    ("train", train_mask, train_labels),
                    ("validation", validation_mask, validation_labels),
                ):
                    baseline = binary_metrics(
                        labels,
                        pd.to_numeric(
                            selected_sites.loc[mask, "footprint_score"], errors="coerce"
                        ).to_numpy(dtype=float),
                    )
                    all_metrics.append(
                        pd.DataFrame(
                            [
                                {
                                    "cell": str(cell),
                                    "tf": str(tf),
                                    "chromosome_split": split,
                                    "stage": "source_table_baseline",
                                    "candidate": "source_table_footprint_score",
                                    "correction": "source_table",
                                    **baseline,
                                }
                            ]
                        )
                    )

            stage1_parts = []
            expanded: list[Candidate] = []
            for correction, cell_profiles in profiles_by_correction.items():
                tf_profiles = cell_profiles[positions]
                first = evaluate_candidates(
                    tf_profiles[train_mask], train_labels, geometry_grid(correction, flank, quick),
                    cell=str(cell), tf=str(tf), split="train", stage="geometry",
                )
                stage1_parts.append(first)
                expanded.extend(
                    candidate
                    for row in shortlist(first, stage1_keep).itertuples(index=False)
                    for candidate in expanded_grid(candidate_from_row(row))
                )
            stage1 = pd.concat(stage1_parts, ignore_index=True)
            all_metrics.append(stage1)

            stage2_parts = []
            for correction, correction_candidates in pd.Series(
                expanded, dtype=object
            ).groupby(lambda index: expanded[index].correction):
                candidates = list(correction_candidates)
                tf_profiles = profiles_by_correction[str(correction)][positions]
                stage2_parts.append(
                    evaluate_candidates(
                        tf_profiles[train_mask], train_labels, candidates,
                        cell=str(cell), tf=str(tf), split="train", stage="expanded",
                    )
                )
            stage2 = pd.concat(stage2_parts, ignore_index=True).drop_duplicates("candidate")
            all_metrics.append(stage2)
            frozen_candidates = [
                candidate_from_row(row)
                for row in shortlist(stage2, validation_keep).itertuples(index=False)
            ]

            validation_parts = []
            for correction in sorted({candidate.correction for candidate in frozen_candidates}):
                candidates = [candidate for candidate in frozen_candidates if candidate.correction == correction]
                tf_profiles = profiles_by_correction[correction][positions]
                validation_parts.append(
                    evaluate_candidates(
                        tf_profiles[validation_mask], validation_labels, candidates,
                        cell=str(cell), tf=str(tf), split="validation", stage="frozen_validation",
                    )
                )
            validation = pd.concat(validation_parts, ignore_index=True)
            all_metrics.append(validation)
            winner = shortlist(validation, 1).iloc[0].to_dict()
            winner["selection_status"] = "frozen_without_test_evaluation"
            winners.append(winner)

    metrics = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    return metrics, pd.DataFrame(winners)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", nargs="+", type=Path, required=True)
    parser.add_argument("--signals", type=Path, required=True)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--max-sites-per-class", type=int, default=5000)
    parser.add_argument("--negative-pool-multiplier", type=int, default=1)
    parser.add_argument("--stage1-keep", type=int, default=20)
    parser.add_argument("--validation-keep", type=int, default=50)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--profiles-only", action="store_true")
    parser.add_argument("--profile-cache-dir", type=Path)
    parser.add_argument(
        "--chromosome-splits", nargs="+", default=["train", "validation"],
        choices=["train", "validation", "test"],
    )
    args = parser.parse_args(argv)

    if args.flank < 32:
        parser.error("--flank must be at least 32")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    sites = attach_chromosome_split(read_sites(args.sites), study)
    sites = sites[sites["chromosome_split"].isin(args.chromosome_splits)].copy()
    if not args.profiles_only and not {"train", "validation"}.issubset(args.chromosome_splits):
        parser.error("model search requires both train and validation chromosome splits")
    sites = deterministic_class_sample(
        sites, args.max_sites_per_class, args.seed, args.negative_pool_multiplier
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    sites.to_csv(args.outdir / "sampled_sites.tsv.gz", sep="\t", index=False, compression="gzip")
    write_merged_regions(sites, args.outdir / "sampled_regions.bed", padding=args.flank + 100)
    metrics, winners = run_search(
        sites,
        read_signals(args.signals),
        args.outdir,
        args.flank,
        args.stage1_keep,
        args.validation_keep,
        args.quick,
        args.profiles_only,
        args.profile_cache_dir,
    )
    metrics.to_csv(args.outdir / "per_tf_candidate_metrics.tsv.gz", sep="\t", index=False, compression="gzip")
    winners.to_csv(args.outdir / "per_tf_frozen_candidates.tsv", sep="\t", index=False)
    print(winners.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
