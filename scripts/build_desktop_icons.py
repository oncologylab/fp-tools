#!/usr/bin/env python3
"""Build native desktop icon files from the canonical fp-tools logo."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image


def build_icons(source: Path, output_dir: Path) -> tuple[Path, Path]:
    """Write multi-resolution Windows ICO and macOS ICNS assets."""

    output_dir.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as opened:
        logo = opened.convert("RGBA")
        if logo.width != logo.height:
            raise ValueError("Desktop icon source must be square")
        ico_path = output_dir / "fp-tools.ico"
        icns_path = output_dir / "fp-tools.icns"
        logo.save(
            ico_path,
            format="ICO",
            sizes=[
                (16, 16),
                (24, 24),
                (32, 32),
                (48, 48),
                (64, 64),
                (128, 128),
                (256, 256),
            ],
        )
        logo.save(icns_path, format="ICNS")
    return ico_path, icns_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("docs/assets/fp_tools_logo_icon_1024.png"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("build/desktop-icons"))
    args = parser.parse_args()
    ico_path, icns_path = build_icons(args.source, args.output_dir)
    print(ico_path)
    print(icns_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
