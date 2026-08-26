"""Indexed FASTA access that keeps generated indexes outside source data."""

from __future__ import annotations

import atexit
import hashlib
import errno
import os
from pathlib import Path
import time
import uuid

try:
    import pysam
except ImportError:  # pragma: no cover - Windows path
    pysam = None

if pysam is None:  # pragma: no cover - exercised by Windows CI
    import pyfastx
else:
    pyfastx = None


_PROCESS_INDEX_TOKEN = uuid.uuid4().hex
_PROCESS_INDEXES: set[Path] = set()


def _cleanup_process_indexes() -> None:
    for index_path in tuple(_PROCESS_INDEXES):
        try:
            index_path.unlink(missing_ok=True)
        except OSError:
            pass
        else:
            _PROCESS_INDEXES.discard(index_path)


atexit.register(_cleanup_process_indexes)


def _cache_root() -> Path:
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    path = base / "fp-tools" / "fasta"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _open_pyfastx(filename: str, index_path: Path):
    """Open pyfastx while serializing creation of its shared SQLite index."""

    lock_path = index_path.with_suffix(index_path.suffix + ".lock")
    deadline = time.monotonic() + 120
    lock_fd = None
    while lock_fd is None:
        try:
            lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except OSError as exc:
            if exc.errno not in {errno.EACCES, errno.EEXIST}:
                raise
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Timed out waiting for FASTA index lock: {lock_path}")
            time.sleep(0.05)
    try:
        return pyfastx.Fasta(filename, index_file=str(index_path), full_name=False)
    finally:
        os.close(lock_fd)
        lock_path.unlink(missing_ok=True)


def _use_process_local_index() -> bool:
    return os.name == "nt"


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
        self._process_local_index = _use_process_local_index()
        if self._process_local_index:
            index_name = f"{index_name}.{_PROCESS_INDEX_TOKEN}"
        self.index_path = _cache_root() / index_name
        self._fasta = _open_pyfastx(self.filename, self.index_path)
        if self._process_local_index:
            _PROCESS_INDEXES.add(self.index_path)
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

    def get_reference_length(self, reference: str) -> int:
        """Return the length of one reference sequence."""

        if self._pysam:
            return int(self._fasta.get_reference_length(str(reference)))
        return int(len(self._fasta[str(reference)]))

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
        elif self._fasta is not None:
            fasta = self._fasta
            close = getattr(fasta, "close", None)
            if close is not None:
                close()
            del fasta
        self._fasta = None


def open_fasta(path: str | os.PathLike[str]) -> FastaFile:
    return FastaFile(path)
