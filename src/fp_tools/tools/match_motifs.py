#!/usr/bin/env python
"""Command entry point for motif matching against footprint score tracks."""

from fp_tools.tools.diff_footprints import match_motifs_cli


def main():
    return match_motifs_cli()


if __name__ == "__main__":
    main()
