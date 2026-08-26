"""Cross-platform bigWig I/O with a small pyBigWig-compatible surface.

``pyBigWig`` is used when it is available.  Other platforms use
``pybigtools``, whose wheels cover Windows, macOS, and Linux.  The adapter is
deliberately limited to the methods used by fp-tools so scientific code does
not need platform checks.
"""

from __future__ import annotations

import math
import ntpath
import os
import queue
import threading
from pathlib import Path
from typing import Iterable

import numpy as np

# pyBigWig exposes this feature flag and legacy code uses it to request arrays.
numpy = 1

try:  # Fast native backend on platforms where a wheel is available.
    import pyBigWig as _pybigwig
except ImportError:  # pragma: no cover - exercised on Windows/macOS CI
    _pybigwig = None

try:
    import pybigtools as _pybigtools
except ImportError:  # pragma: no cover - dependency metadata installs it
    _pybigtools = None


def backend_name() -> str:
    """Return the selected bigWig backend name."""

    return "pyBigWig" if _pybigwig is not None else "pybigtools"


def _pybigtools_write_path(path: str | os.PathLike[str]) -> str:
    """Return a local path that pybigtools will not parse as a URL.

    pybigtools checks strings with Rust's URL parser before opening an output.
    A normal absolute Windows path (``C:\\...``) is consequently mistaken for
    a URL with scheme ``c``.  The equivalent Windows extended-length form is
    still a valid local filesystem path and is unambiguously not a URL.
    """

    if os.name != "nt":
        return os.path.abspath(os.fspath(path))
    resolved = ntpath.abspath(os.fspath(path))
    if resolved.startswith("\\\\?\\"):
        return resolved
    if resolved.startswith("\\\\"):
        return "\\\\?\\UNC\\" + resolved.lstrip("\\")
    return "\\\\?\\" + resolved


def open(path: str | os.PathLike[str], mode: str = "r"):  # noqa: A001 - API compatibility
    """Open a bigWig file using the best available backend."""

    path = str(path)
    normalized_mode = mode.replace("b", "")
    if _pybigwig is not None:
        return _pybigwig.open(path, normalized_mode)
    if _pybigtools is None:
        raise ImportError("bigWig support requires pybigtools or pyBigWig")
    if normalized_mode.startswith("w"):
        return _StreamingBigWigWriter(path)
    return _BigWigReader(path)


class _BigWigReader:
    def __init__(self, path: str):
        self.path = path
        self._handle = _pybigtools.open(path)

    def __bool__(self) -> bool:
        return self._handle is not None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()
        return False

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def chroms(self, chrom: str | None = None):
        result = dict(self._handle.chroms())
        return result.get(chrom) if chrom is not None else result

    def header(self) -> dict[str, float | int]:
        """Return the pyBigWig-compatible summary header."""

        info = self._handle.info()
        summary = info.get("summary") or {}
        covered = int(summary.get("basesCovered", 0) or 0)
        mean = float(summary.get("mean", 0.0) or 0.0)
        std = float(summary.get("std", 0.0) or 0.0)
        return {
            "version": int(info.get("version", 0) or 0),
            "nLevels": int(info.get("zoomLevels", 0) or 0),
            "nBasesCovered": covered,
            "minVal": float(summary.get("min", 0.0) or 0.0),
            "maxVal": float(summary.get("max", 0.0) or 0.0),
            "sumData": mean * covered,
            "sumSquared": (std * std + mean * mean) * covered,
        }

    def values(self, chrom: str, start: int, end: int, numpy: bool = False):
        values = list(self._handle.values(chrom, int(start), int(end), fillna=None))
        return np.asarray(values, dtype=float) if numpy else values

    def intervals(self, chrom: str, start: int | None = None, end: int | None = None):
        chrom_size = self.chroms(chrom)
        if chrom_size is None:
            return None
        lo = 0 if start is None else max(0, int(start))
        hi = int(chrom_size) if end is None else min(int(chrom_size), int(end))
        records = list(self._handle.records(chrom, lo, hi))
        return records or None

    def stats(
        self,
        chrom: str,
        start: int,
        end: int,
        type: str = "mean",  # noqa: A002 - pyBigWig API compatibility
        nBins: int = 1,
        **_kwargs,
    ):
        if type not in {"mean", "min", "max", "sum", "coverage", "std"}:
            raise ValueError(f"Unsupported bigWig summary type: {type}")
        edges = np.linspace(int(start), int(end), max(1, int(nBins)) + 1, dtype=int)
        summaries = []
        for left, right in zip(edges[:-1], edges[1:]):
            values = np.asarray(self.values(chrom, int(left), int(right), numpy=True), dtype=float)
            finite = values[np.isfinite(values)]
            if finite.size == 0:
                summaries.append(None)
            elif type == "mean":
                summaries.append(float(np.mean(finite)))
            elif type == "min":
                summaries.append(float(np.min(finite)))
            elif type == "max":
                summaries.append(float(np.max(finite)))
            elif type == "sum":
                summaries.append(float(np.sum(finite)))
            elif type == "std":
                summaries.append(float(np.std(finite)))
            else:
                summaries.append(float(finite.size / max(1, right - left)))
        return summaries


class _StreamingBigWigWriter:
    """Feed incremental pyBigWig-style entries to a pybigtools writer."""

    def __init__(self, path: str):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._header: dict[str, int] | None = None
        self._queue = queue.Queue(maxsize=100_000)
        self._sentinel = object()
        self._thread = None
        self._error = None
        self._closed = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        if exc_type is None:
            self.close()
        else:
            self._finish()
        return False

    def addHeader(self, header: Iterable[tuple[str, int]], **_kwargs) -> None:
        self._header = {str(chrom): int(size) for chrom, size in header}
        self._thread = threading.Thread(target=self._write, daemon=True)
        self._thread.start()

    def _write(self) -> None:
        try:
            writer = _pybigtools.open(_pybigtools_write_path(self.path), "w")
            writer.write(self._header, self._records())
        except Exception as exc:  # surfaced in close/addEntries
            self._error = exc

    def addEntries(self, chroms, starts, ends=None, values=None, span=1, **_kwargs) -> None:
        if values is None:
            raise ValueError("values are required when writing bigWig entries")
        values = list(values)
        if isinstance(chroms, str):
            if isinstance(starts, (int, np.integer)):
                first_start = int(starts)
                step = int(_kwargs.get("step", 1))
                starts = [first_start + index * step for index in range(len(values))]
            chrom_values = [chroms] * len(starts)
        else:
            chrom_values = list(chroms)
            starts = list(starts)
        if ends is None:
            end_values = [int(start) + int(span) for start in starts]
        else:
            end_values = list(ends)
        if not (len(chrom_values) == len(starts) == len(end_values) == len(values)):
            raise ValueError("chroms, starts, ends, and values must have equal lengths")
        for chrom, start, end, value in zip(chrom_values, starts, end_values, values):
            numeric = float(value)
            if math.isfinite(numeric):
                record = (str(chrom), int(start), int(end), numeric)
                while True:
                    if self._error is not None:
                        raise RuntimeError("pybigtools writer failed") from self._error
                    try:
                        self._queue.put(record, timeout=0.2)
                        break
                    except queue.Full:
                        continue

    def _records(self):
        while True:
            record = self._queue.get()
            if record is self._sentinel:
                return
            yield record

    def close(self) -> None:
        if self._closed:
            return
        if self._header is None:
            raise ValueError("addHeader must be called before closing a bigWig writer")
        self._finish()
        if self._error is not None:
            raise RuntimeError("pybigtools writer failed") from self._error
        self._closed = True

    def _finish(self) -> None:
        if self._thread is None:
            return
        while self._thread.is_alive():
            try:
                self._queue.put(self._sentinel, timeout=0.2)
                break
            except queue.Full:
                if self._error is not None:
                    break
        self._thread.join()
        self._thread = None
