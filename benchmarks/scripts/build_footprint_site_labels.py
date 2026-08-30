#!/usr/bin/env python3
"""Build conservative motif-site labels from reproducible TF ChIP peaks.

The natural output contains positive and high-confidence negative motif sites.
Sites near a ChIP peak but lacking summit support are written as indeterminate
when requested and are never silently used as negatives.  An optional matched
output uses a deterministic propensity score to select negative controls with
similar measured covariates.
"""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
import gzip
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class Peak:
    start: int
    end: int
    summit: int | None


@dataclass(frozen=True)
class PeakIndex:
    """Compact chromosome index for logarithmic peak and summit queries."""

    peaks: tuple[Peak, ...]
    starts: tuple[int, ...]
    prefix_max_ends: tuple[int, ...]
    summit_positions: tuple[int, ...]
    summit_peaks: tuple[Peak, ...]


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_peaks(path: Path) -> dict[str, list[Peak]]:
    peaks: dict[str, list[Peak]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number} has fewer than three columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number} has invalid coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number} has invalid interval {start}-{end}")
            summit: int | None = None
            if len(fields) >= 10:
                try:
                    offset = int(fields[9])
                except ValueError:
                    offset = -1
                if offset >= 0:
                    summit = start + offset
            peaks.setdefault(fields[0], []).append(Peak(start, end, summit))
    for chrom in peaks:
        peaks[chrom].sort(key=lambda peak: (peak.start, peak.end, peak.summit or -1))
    return peaks


def read_sites(
    path: Path,
    site_id_column: int = 4,
    motif_score_column: int = 5,
    strand_column: int = 6,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number} has fewer than three columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number} has invalid coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number} has invalid interval {start}-{end}")

            def field(column: int, default: str) -> str:
                index = column - 1
                return fields[index] if column > 0 and index < len(fields) else default

            raw_score = field(motif_score_column, "nan")
            try:
                motif_score = float(raw_score)
            except ValueError:
                motif_score = np.nan
            site_id = field(site_id_column, f"{fields[0]}:{start}-{end}")
            strand = field(strand_column, ".")
            rows.append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "strand": strand if strand in {"+", "-"} else ".",
                    "site_id": site_id,
                    "motif_score": motif_score,
                }
            )
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError(f"{path} contains no motif sites")
    duplicated = frame.duplicated(["chrom", "start", "end", "strand"])
    if duplicated.any():
        row = frame.loc[duplicated].iloc[0]
        raise ValueError(
            f"{path} contains duplicate motif site {row.chrom}:{row.start}-{row.end}:{row.strand}"
        )
    return frame


def interval_distance(position: int, peak: Peak) -> int:
    if peak.start <= position < peak.end:
        return 0
    if position < peak.start:
        return peak.start - position
    return position - peak.end + 1


def index_peaks(peaks: dict[str, list[Peak]]) -> dict[str, PeakIndex]:
    indexes: dict[str, PeakIndex] = {}
    for chrom, chrom_peaks in peaks.items():
        ordered = tuple(sorted(chrom_peaks, key=lambda peak: (peak.start, peak.end)))
        prefix_max_ends: list[int] = []
        maximum_end = -1
        for peak in ordered:
            maximum_end = max(maximum_end, peak.end)
            prefix_max_ends.append(maximum_end)
        summit_pairs = sorted(
            (peak.summit, peak)
            for peak in ordered
            if peak.summit is not None
        )
        indexes[chrom] = PeakIndex(
            peaks=ordered,
            starts=tuple(peak.start for peak in ordered),
            prefix_max_ends=tuple(prefix_max_ends),
            summit_positions=tuple(int(summit) for summit, _ in summit_pairs),
            summit_peaks=tuple(peak for _, peak in summit_pairs),
        )
    return indexes


def nearest_peak_distance(position: int, index: PeakIndex) -> int:
    insertion = bisect_right(index.starts, position)
    if insertion and index.prefix_max_ends[insertion - 1] > position:
        return 0
    distances: list[int] = []
    if insertion < len(index.starts):
        distances.append(index.starts[insertion] - position)
    if insertion:
        distances.append(position - index.prefix_max_ends[insertion - 1] + 1)
    return min(distances)


def nearest_summit_distance(position: int, index: PeakIndex) -> float:
    if not index.summit_positions:
        return np.inf
    insertion = bisect_left(index.summit_positions, position)
    distances: list[int] = []
    if insertion < len(index.summit_positions):
        distances.append(index.summit_positions[insertion] - position)
    if insertion:
        distances.append(position - index.summit_positions[insertion - 1])
    return float(min(distances))


def has_supported_summit(
    position: int,
    index: PeakIndex,
    maximum_distance: int,
) -> bool:
    left = bisect_left(index.summit_positions, position - maximum_distance)
    right = bisect_right(index.summit_positions, position + maximum_distance)
    return any(
        peak.start <= position < peak.end for peak in index.summit_peaks[left:right]
    )


def label_site(
    chrom: str,
    start: int,
    end: int,
    peaks: dict[str, PeakIndex],
    positive_summit_distance: int,
    negative_peak_distance: int,
) -> tuple[int, str, float, float]:
    center = (start + end) // 2
    peak_index = peaks.get(chrom)
    if peak_index is None:
        return 0, "negative_no_chromosome_peaks", np.inf, np.inf
    peak_distance = nearest_peak_distance(center, peak_index)
    summit_distance = nearest_summit_distance(center, peak_index)
    if peak_distance == 0:
        if has_supported_summit(center, peak_index, positive_summit_distance):
            return 1, "positive_summit_supported", 0.0, summit_distance
        return -1, "indeterminate_peak_without_nearby_summit", 0.0, summit_distance
    if peak_distance > negative_peak_distance:
        return (
            0,
            "negative_far_from_peak",
            float(peak_distance),
            summit_distance,
        )
    return (
        -1,
        "indeterminate_near_peak",
        float(peak_distance),
        summit_distance,
    )


def label_sites(
    sites: pd.DataFrame,
    peaks: dict[str, list[Peak]],
    positive_summit_distance: int = 100,
    negative_peak_distance: int = 500,
) -> pd.DataFrame:
    if positive_summit_distance < 0 or negative_peak_distance < 0:
        raise ValueError("label distances must be non-negative")
    peak_indexes = index_peaks(peaks)
    labels = [
        label_site(
            row.chrom,
            int(row.start),
            int(row.end),
            peak_indexes,
            positive_summit_distance,
            negative_peak_distance,
        )
        for row in sites.itertuples(index=False)
    ]
    output = sites.copy()
    output[[
        "label",
        "label_reason",
        "nearest_peak_distance",
        "nearest_summit_distance",
    ]] = pd.DataFrame(labels, index=output.index)
    return output


def merge_features(sites: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    keys = ["chrom", "start", "end"]
    missing = [column for column in keys if column not in features]
    if missing:
        raise ValueError(f"feature table is missing columns: {', '.join(missing)}")
    if features.duplicated(keys).any():
        raise ValueError("feature table contains duplicate chrom/start/end rows")
    overlapping = [
        column for column in features.columns if column in sites.columns and column not in keys
    ]
    if overlapping:
        raise ValueError(
            "feature columns collide with site columns: " + ", ".join(overlapping)
        )
    return sites.merge(features, on=keys, how="left", validate="many_to_one")


def _pop_closest(available: list[tuple[float, int]], target: float) -> tuple[float, int]:
    position = bisect_left(available, (target, -1))
    candidates = []
    if position < len(available):
        candidates.append((abs(available[position][0] - target), position))
    if position:
        candidates.append((abs(available[position - 1][0] - target), position - 1))
    _, selected = min(candidates)
    return available.pop(selected)


def propensity_match(
    frame: pd.DataFrame,
    feature_columns: Iterable[str],
    negative_ratio: int = 1,
    seed: int = 2026,
) -> pd.DataFrame:
    feature_columns = list(feature_columns)
    if negative_ratio < 1:
        raise ValueError("negative_ratio must be at least one")
    missing = [column for column in feature_columns if column not in frame]
    if missing:
        raise ValueError(f"matching features are missing: {', '.join(missing)}")
    usable = frame[frame["label"].isin([0, 1])].copy()
    usable[feature_columns] = usable[feature_columns].apply(pd.to_numeric, errors="coerce")
    usable = usable.dropna(subset=feature_columns)
    positives = usable[usable["label"] == 1]
    negatives = usable[usable["label"] == 0]
    required_negatives = len(positives) * negative_ratio
    if positives.empty or len(negatives) < required_negatives:
        raise ValueError(
            f"matching requires {required_negatives} negatives for {len(positives)} positives; "
            f"found {len(negatives)}"
        )
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=seed),
    )
    model.fit(usable[feature_columns], usable["label"])
    usable["propensity_score"] = model.predict_proba(usable[feature_columns])[:, 1]
    positives = usable[usable["label"] == 1].sort_values(
        ["propensity_score", "chrom", "start", "end"], kind="mergesort"
    )
    available = sorted(
        (float(row.propensity_score), int(index))
        for index, row in usable[usable["label"] == 0].iterrows()
    )
    selected_negative_indices: list[int] = []
    for row in positives.itertuples():
        for _ in range(negative_ratio):
            _, selected_index = _pop_closest(available, float(row.propensity_score))
            selected_negative_indices.append(selected_index)
    matched = pd.concat([positives, usable.loc[selected_negative_indices]], ignore_index=True)
    return matched.sort_values(
        ["label", "chrom", "start", "end"], ascending=[False, True, True, True]
    ).reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--chip-peaks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--indeterminate-out", type=Path)
    parser.add_argument("--features", type=Path)
    parser.add_argument("--matched-out", type=Path)
    parser.add_argument("--match-columns", nargs="+", default=["motif_score"])
    parser.add_argument("--negative-ratio", type=int, default=1)
    parser.add_argument("--positive-summit-distance", type=int, default=100)
    parser.add_argument("--negative-peak-distance", type=int, default=500)
    parser.add_argument("--site-id-column", type=int, default=4)
    parser.add_argument("--motif-score-column", type=int, default=5)
    parser.add_argument("--strand-column", type=int, default=6)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    sites = read_sites(
        args.sites,
        site_id_column=args.site_id_column,
        motif_score_column=args.motif_score_column,
        strand_column=args.strand_column,
    )
    labelled = label_sites(
        sites,
        read_peaks(args.chip_peaks),
        positive_summit_distance=args.positive_summit_distance,
        negative_peak_distance=args.negative_peak_distance,
    )
    if args.features:
        labelled = merge_features(labelled, pd.read_csv(args.features, sep="\t"))
    natural = labelled[labelled["label"].isin([0, 1])].copy()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    natural.to_csv(args.out, sep="\t", index=False)
    if args.indeterminate_out:
        args.indeterminate_out.parent.mkdir(parents=True, exist_ok=True)
        labelled[labelled["label"] == -1].to_csv(
            args.indeterminate_out, sep="\t", index=False
        )
    if args.matched_out:
        matched = propensity_match(
            natural,
            args.match_columns,
            negative_ratio=args.negative_ratio,
            seed=args.seed,
        )
        args.matched_out.parent.mkdir(parents=True, exist_ok=True)
        matched.to_csv(args.matched_out, sep="\t", index=False)
    counts = labelled["label_reason"].value_counts().sort_index()
    print(counts.to_string())
    print(f"wrote {len(natural):,} labelled sites to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
