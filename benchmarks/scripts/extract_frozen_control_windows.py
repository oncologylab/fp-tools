#!/usr/bin/env python3
"""Extract checksum-protected control windows without fitting a bias model.

This separates independent naked-DNA candidate/final controls from training:
the resulting windows can be scored by a frozen model, but this command never
estimates coefficients or reads TF labels.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_parametric_bias import (  # noqa: E402
    build_or_load_datasets,
    parse_name_path,
    parse_shift,
)


SCHEMA = "fp-tools-frozen-control-window-extraction-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_manifest(
    *,
    study: Path,
    source: str,
    samples: Sequence[tuple[str, Path]],
    window_manifest,
    output_table: Path,
) -> dict[str, object]:
    artifacts = []
    for path in sorted({Path(value) for value in window_manifest["cache_npz"]}):
        sidecar = path.with_suffix(".json")
        artifacts.append(
            {
                "path": str(path),
                "sha256": file_sha256(path),
                "metadata_path": str(sidecar),
                "metadata_sha256": file_sha256(sidecar),
            }
        )
    return {
        "schema": SCHEMA,
        "study": str(study),
        "study_sha256": file_sha256(study),
        "source": source,
        "samples": [
            {"name": name, "path": str(path), "sha256": file_sha256(path)}
            for name, path in samples
        ],
        "models_fitted": False,
        "chipped_labels_used": False,
        "window_artifacts": artifacts,
        "window_table": {
            "path": str(output_table),
            "sha256": file_sha256(output_table),
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument(
        "--sample", type=parse_name_path, action="append", required=True
    )
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument(
        "--source", choices=("naked_dna", "mitochondrial"), required=True
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument(
        "--read-shift", type=parse_shift, action="append", dest="read_shifts"
    )
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--window-size", type=int, default=200)
    parser.add_argument("--margin", type=int, default=41)
    parser.add_argument("--train-windows", type=int, default=3000)
    parser.add_argument("--validation-windows", type=int, default=1000)
    parser.add_argument("--candidate-factor", type=int, default=4)
    parser.add_argument("--minimum-mapq", type=int, default=30)
    parser.add_argument("--keep-duplicates", action="store_true")
    parser.add_argument(
        "--mitochondrial-chromosome",
        action="append",
        dest="mitochondrial_chromosomes",
    )
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    if study.get("status") != "development_locked_holdout_unscored":
        raise ValueError("control extraction requires the locked, unscored study")
    args.read_shifts = args.read_shifts or [(4, -5), (4, -4)]
    args.mitochondrial_chromosomes = args.mitochondrial_chromosomes or ["chrM", "MT"]
    args.peaks = None
    args.blacklist = None
    args.low_signal_quantile = 0.75
    args.outdir.mkdir(parents=True, exist_ok=True)
    _datasets, window_manifest = build_or_load_datasets(args, study)
    table_path = args.outdir / "control_windows.tsv"
    window_manifest.to_csv(table_path, sep="\t", index=False)
    manifest = build_manifest(
        study=args.study,
        source=args.source,
        samples=args.sample,
        window_manifest=window_manifest,
        output_table=table_path,
    )
    (args.outdir / "control_window_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(window_manifest.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
