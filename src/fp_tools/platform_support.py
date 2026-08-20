"""Platform support boundaries shared by fp-tools command interfaces."""

from __future__ import annotations

import platform


RAW_READ_LINUX_ONLY_MESSAGE = (
    "FASTQ-to-BAM preparation is supported only by the Linux command line and "
    "the Linux fp-tools container. Native macOS and Windows installations, "
    "desktop applications, and every fp-tools GUI start from coordinate-sorted "
    "BAM/BAI files plus matching peak BED files."
)


def supports_raw_read_preparation() -> bool:
    """Return whether this process is running on a supported raw-read host."""

    return platform.system() == "Linux"


def require_raw_read_preparation_support() -> None:
    """Raise a concise error when raw-read preparation is unavailable."""

    if not supports_raw_read_preparation():
        raise ValueError(RAW_READ_LINUX_ONLY_MESSAGE)
