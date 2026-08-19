"""Indexed FASTA access that keeps generated indexes outside source data."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pyfastx

try:
    import pysam
except ImportError:  # pragma: no cover - Windows path
    pysam = None


def _cache_root() -> Path:
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    path = base / "fp-tools" / "fasta"
    path.mkdir(parents=True, exist_ok=True)
    return path


class FastaFile:
    """Subset of :class:`pysam.FastaFile` used by fp-tools."""

    def __init__(self, path: str | os.PathLike[str]):
        self.filename = str(Path(path).expanduser().resolve())
        if pysam is not None:
            self.index_path = None
            self._fasta = pysam.FastaFile(self.filename)
            self._pysam = True
            return
        source = Path(self.filename)
        identity = f"{self.filename}:{source.stat().st_size}:{source.stat().st_mtime_ns}"
        index_name = hashlib.sha256(identity.encode("utf-8")).hexdigest() + ".fxi"
        self.index_path = _cache_root() / index_name
        self._fasta = pyfastx.Fasta(self.filename, index_file=str(self.index_path), full_name=False)
        self._pysam = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()
        return False

    @property
    def references(self) -> tuple[str, ...]:
        return tuple(self._fasta.references) if self._pysam else tuple(self._fasta.keys())

    @property
    def lengths(self) -> tuple[int, ...]:
        return tuple(self._fasta.lengths) if self._pysam else tuple(len(self._fasta[name]) for name in self.references)

    def fetch(self, reference: str, start: int | None = None, end: int | None = None) -> str:
        if self._pysam:
            return str(self._fasta.fetch(str(reference), start, end))
        record = self._fasta[str(reference)]
        lo = 0 if start is None else int(start)
        hi = len(record) if end is None else int(end)
        return str(record[lo:hi])

    def close(self) -> None:
        if self._pysam and self._fasta is not None:
            self._fasta.close()
        self._fasta = None


def open_fasta(path: str | os.PathLike[str]) -> FastaFile:
    return FastaFile(path)
