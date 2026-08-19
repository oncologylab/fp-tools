"""Small, dependency-free BED interval overlap helpers."""

from __future__ import annotations

import bisect
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator


def _iter_bed(path: str | Path) -> Iterator[tuple[str, int, int, str]]:
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                yield fields[0], int(fields[1]), int(fields[2]), line
            except ValueError:
                continue


class IntervalIndex:
    """Chromosome-partitioned interval index for overlap membership queries."""

    def __init__(self, intervals: Iterable[tuple[str, int, int]] = ()):
        grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in intervals:
            if end > start:
                grouped[str(chrom)].append((int(start), int(end)))
        self._intervals: dict[str, tuple[list[int], list[int]]] = {}
        for chrom, values in grouped.items():
            values.sort()
            merged: list[list[int]] = []
            for start, end in values:
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            self._intervals[chrom] = (
                [item[0] for item in merged],
                [item[1] for item in merged],
            )

    @classmethod
    def from_bed(cls, path: str | Path) -> "IntervalIndex":
        return cls((chrom, start, end) for chrom, start, end, _line in _iter_bed(path))

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        starts, ends = self._intervals.get(str(chrom), ([], []))
        if not starts or end <= start:
            return False
        index = bisect.bisect_left(starts, int(end)) - 1
        return index >= 0 and ends[index] > int(start)


def intersect_bed(
    query: str | Path,
    regions: str | Path,
    output: str | Path,
    *,
    invert: bool = False,
) -> Path:
    """Write query BED records that overlap (or do not overlap) regions."""

    index = IntervalIndex.from_bed(regions)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for chrom, start, end, line in _iter_bed(query):
            keep = index.overlaps(chrom, start, end)
            if keep != invert:
                handle.write(line if line.endswith("\n") else line + "\n")
    return output


def filter_regions(region_list, regions: str | Path, *, invert: bool = False):
    """Return a RegionList-like object filtered by overlap membership."""

    index = IntervalIndex.from_bed(regions)
    return region_list.__class__(
        region
        for region in region_list
        if index.overlaps(region.chrom, region.start, region.end) != invert
    )
