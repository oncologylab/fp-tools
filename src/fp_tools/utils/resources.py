"""Processor budgets shared by commands and their workflow wrappers."""

from __future__ import annotations

import os
from collections.abc import Callable


def available_cores() -> int:
    """Return process-visible logical processors, including affinity limits."""
    process_count = getattr(os, "process_cpu_count", None)
    if process_count is not None:
        count = process_count()
        if count:
            return max(1, int(count))
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError, NotImplementedError):
        return max(1, os.cpu_count() or 1)


def resolve_cores(given: int | None, *, warn: Callable[[str], None] | None = None) -> int:
    """Resolve automatic budgets, preserving the established core-limit policy."""
    available = available_cores()
    if given is None:
        return available
    if given < 1:
        message = f"Invalid '--cores' value {given}; using all available cores ({available})."
    elif given > available:
        message = f"Requested {given} cores, but only {available} are available; using {available}."
    else:
        return int(given)
    if warn is not None:
        warn(message)
    return available
