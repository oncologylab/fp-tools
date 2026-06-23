#!/usr/bin/env python
"""Compatibility wrapper for the PBMC5k manuscript signature plots."""

from __future__ import annotations

import sys

from fp_tools.tools.find_signature_fp import main


if __name__ == "__main__":
    raise SystemExit(main(["--legacy-pbmc5k-names", *sys.argv[1:]]))
