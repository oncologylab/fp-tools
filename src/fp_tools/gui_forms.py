"""Small Streamlit form helpers shared across GUI pages."""

from __future__ import annotations

from collections.abc import Iterable


def newline_values(text: str) -> list[str]:
    return [line.strip() for line in text.splitlines() if line.strip()]


def joined_lines(values: Iterable[str]) -> str:
    return "\n".join(str(value) for value in values if str(value).strip())
