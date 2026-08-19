"""Cross-platform read-only BAM access used by core footprinting commands."""

from __future__ import annotations

import os
import gzip
from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path

try:
    import pysam as _pysam
except ImportError:  # pragma: no cover - Windows path
    _pysam = None

try:
    import bamnostic as _bamnostic
except ImportError:  # pragma: no cover - dependency metadata installs it
    _bamnostic = None


def backend_name() -> str:
    return "pysam" if _pysam is not None else "bamnostic"


def open_alignment(path: str | os.PathLike[str], mode: str = "rb"):
    """Open a BAM for reading with normalized record attributes."""

    path_text = str(path)
    if path_text.endswith((".tsv", ".tsv.gz", ".fragments", ".fragments.gz")):
        return FragmentAlignment(path_text)
    if mode not in {"r", "rb"}:
        if _pysam is None:
            raise RuntimeError("Writing BAM files requires pysam; use fragment input on Windows")
        return _pysam.AlignmentFile(str(path), mode)
    if _pysam is not None:
        return _pysam.AlignmentFile(str(path), mode)
    if _bamnostic is None:
        raise ImportError("BAM support requires pysam or bamnostic")
    return _BamnosticAlignment(str(path))


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


@dataclass
class _FragmentRecord:
    reference_name: str
    reference_start: int
    reference_end: int
    is_reverse: bool
    query_name: str
    template_length: int
    flag: int
    reference_id: int
    query_alignment_start: int = 0
    query_alignment_end: int = 1
    query_length: int = 1
    cigartuples: tuple[tuple[int, int], ...] = ((0, 1),)
    is_unmapped: bool = False
    is_duplicate: bool = False

    def get_tags(self):
        return []

    def infer_query_length(self):
        return self.query_length


class FragmentAlignment:
    """Read 10x-style fragments as paired one-base cut-site anchors."""

    def __init__(self, path: str):
        self.path = str(Path(path).expanduser().resolve())
        self.filename = self.path.encode("utf-8")
        raw: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
        lengths: dict[str, int] = {}
        with _open_text(self.path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                if end - start < 2:
                    continue
                copies = int(fields[4]) if len(fields) > 4 and fields[4].isdigit() else 1
                raw[chrom].append((start, end, copies))
                lengths[chrom] = max(lengths.get(chrom, 0), end)
        self.references = tuple(raw)
        self.lengths = tuple(lengths[chrom] for chrom in self.references)
        self._reference_ids = {chrom: index for index, chrom in enumerate(self.references)}
        self._raw = raw

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def close(self):
        return None

    def has_index(self):
        return True

    def check_index(self):
        return True

    def get_reference_name(self, reference_id: int) -> str:
        return self.references[int(reference_id)]

    def _records(self, chrom: str, start: int, end: int):
        ref_id = self._reference_ids[chrom]
        for row_index, (fragment_start, fragment_end, copies) in enumerate(self._raw.get(chrom, [])):
            if fragment_start >= end or fragment_end <= start:
                continue
            reverse_start = fragment_end - 2
            template_length = fragment_end - fragment_start
            for copy_index in range(copies):
                name = f"{chrom}:{fragment_start}-{fragment_end}:{row_index}:{copy_index}"
                yield _FragmentRecord(
                    chrom, fragment_start, fragment_start + 1, False, name,
                    template_length, 99, ref_id,
                )
                yield _FragmentRecord(
                    chrom, reverse_start, reverse_start + 1, True, name,
                    -template_length, 147, ref_id,
                )

    def fetch(self, reference=None, start=None, end=None, until_eof=False):
        if until_eof:
            for chrom, chrom_length in zip(self.references, self.lengths):
                yield from self._records(chrom, 0, chrom_length)
            return
        if reference not in self._raw:
            return
        yield from self._records(reference, int(start or 0), int(end or self.lengths[self._reference_ids[reference]]))


def index_alignment(path: str | os.PathLike[str]) -> bool:
    """Create a BAM index when pysam is available; otherwise use scan fallback."""

    if _pysam is None:
        return False
    _pysam.index(str(path))
    return True


class _BamnosticRecord:
    def __init__(self, record):
        self._record = record

    def __getattr__(self, name):
        if name == "template_length":
            return self._record.tlen
        if name == "query_alignment_start":
            total = 0
            for operation, size in (self._record.cigartuples or []):
                if operation == 4:
                    total += int(size)
                elif operation != 5:
                    break
            return total
        if name == "query_alignment_end":
            start = self.query_alignment_start
            length = sum(
                int(size)
                for operation, size in (self._record.cigartuples or [])
                if operation in {0, 1, 7, 8}
            )
            return start + length
        if name == "query_alignment_length":
            return self.query_alignment_end - self.query_alignment_start
        if name == "infer_query_length":
            return lambda: int(getattr(self._record, "query_length", 0) or self.query_alignment_end)
        return getattr(self._record, name)


class _BamnosticAlignment:
    """bamnostic reader with an indexed-query or one-pass in-memory fallback."""

    def __init__(self, path: str):
        self.path = str(Path(path).expanduser().resolve())
        self.filename = self.path.encode("utf-8")
        self._handle = _bamnostic.AlignmentFile(self.path, "rb", require_index=False)
        self.references = tuple(self._handle.references)
        self.lengths = tuple(self._handle.lengths)
        self._cache = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()
        return False

    def __iter__(self):
        for record in self._handle:
            yield _BamnosticRecord(record)

    def close(self):
        self._handle.close()

    def has_index(self) -> bool:
        try:
            return bool(self._handle.check_index())
        except (AttributeError, OSError, ValueError):
            return False

    def check_index(self):
        if not self.has_index():
            raise ValueError("BAM index is not available")
        return True

    def get_reference_name(self, reference_id: int) -> str:
        return self.references[int(reference_id)]

    def fetch(self, reference=None, start=None, end=None, until_eof=False):
        if until_eof:
            fresh = _bamnostic.AlignmentFile(self.path, "rb", require_index=False)
            try:
                for record in fresh:
                    yield _BamnosticRecord(record)
            finally:
                fresh.close()
            return
        if reference is None:
            raise ValueError("reference is required unless until_eof=True")
        try:
            for record in self._handle.fetch(reference, int(start or 0), int(end or self.lengths[self.references.index(reference)])):
                yield _BamnosticRecord(record)
            return
        except (AttributeError, OSError, ValueError):
            pass
        if self._cache is None:
            grouped = defaultdict(list)
            fresh = _bamnostic.AlignmentFile(self.path, "rb", require_index=False)
            try:
                for record in fresh:
                    if not record.is_unmapped:
                        grouped[record.reference_name].append(_BamnosticRecord(record))
            finally:
                fresh.close()
            self._cache = grouped
        lo = int(start or 0)
        hi = int(end if end is not None else self.lengths[self.references.index(reference)])
        for record in self._cache.get(reference, []):
            if record.reference_start < hi and record.reference_end > lo:
                yield record
