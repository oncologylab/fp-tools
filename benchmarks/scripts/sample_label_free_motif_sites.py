#!/usr/bin/env python3
"""Sample larger label-free TF motif pools on frozen training chromosomes."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import (  # noqa: E402
    build_unlabeled_training_sites,
)


SCHEMA = "fp-tools-label-free-motif-pools-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_cell_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("sources must use CELL=PATH")
    cell, path = value.split("=", 1)
    if not cell or not path:
        raise argparse.ArgumentTypeError("sources must use CELL=PATH")
    return cell, Path(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument(
        "--source",
        action="append",
        type=parse_cell_path,
        required=True,
        metavar="CELL=TSV",
    )
    parser.add_argument("--maximum-per-tf", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.maximum_per_tf < 100:
        raise SystemExit("--maximum-per-tf must be at least 100")

    experiment_study = json.loads(args.study.read_text(encoding="utf-8"))
    task_study_path = args.study
    task_study = experiment_study
    if "tasks" not in task_study:
        task_study_path = Path(str(experiment_study["base_study"]))
        task_study = json.loads(task_study_path.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(task_study["tasks"])
    tasks = tasks[tasks["split"].astype(str).eq("development")]
    sources = dict(args.source)
    if len(sources) != len(args.source):
        raise SystemExit("duplicate source cell")
    expected_cells = set(tasks["cell"].astype(str))
    if set(sources) != expected_cells:
        raise SystemExit("sources must exactly match development cells")

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.outdir / "label_free_motif_pools.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing.get("schema") != SCHEMA:
            raise ValueError("unsupported existing motif-pool manifest")
        if existing.get("maximum_per_tf") != args.maximum_per_tf:
            raise ValueError("existing motif-pool size differs from requested size")
        if existing.get("seed") != args.seed:
            raise ValueError("existing motif-pool seed differs from requested seed")
        expected_sources = {
            cell: (str(path), file_sha256(path)) for cell, path in sources.items()
        }
        observed_sources = {
            str(record["cell"]): (
                str(record["source"]),
                str(record["source_sha256"]),
            )
            for record in existing["outputs"]
        }
        if observed_sources != expected_sources:
            raise ValueError("existing motif-pool sources differ from requested sources")
        for record in existing["outputs"]:
            if file_sha256(record["output"]) != record["output_sha256"]:
                raise ValueError("existing motif-pool output checksum mismatch")
        if file_sha256(existing["counts"]["path"]) != existing["counts"]["sha256"]:
            raise ValueError("existing motif-pool count checksum mismatch")
        print(pd.read_csv(existing["counts"]["path"], sep="\t").to_string(index=False))
        return 0

    records = []
    count_frames = []
    for cell, source in sorted(sources.items()):
        sites = build_unlabeled_training_sites(
            source,
            cell,
            tasks,
            experiment_study,
            maximum_per_tf=args.maximum_per_tf,
            seed=args.seed,
        )
        forbidden = [
            column
            for column in sites.columns
            if "label" in column.lower() or "chip" in column.lower()
        ]
        if forbidden:
            raise ValueError(
                "sampled label-free sites contain forbidden columns: "
                + ", ".join(forbidden)
            )
        if set(sites["chromosome_split"].astype(str)) != {"train"}:
            raise ValueError("sampled sites contain non-training chromosomes")
        output = args.outdir / f"{cell}.unlabeled_training_sites.tsv.gz"
        sites.to_csv(
            output,
            sep="\t",
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
        counts = (
            sites.groupby(["cell", "tf", "motif_id", "motif_family"], sort=True)
            .size()
            .rename("sites")
            .reset_index()
        )
        count_frames.append(counts)
        records.append(
            {
                "cell": cell,
                "source": str(source),
                "source_sha256": file_sha256(source),
                "output": str(output),
                "output_sha256": file_sha256(output),
                "sites": int(len(sites)),
                "tfs": int(sites["tf"].nunique()),
            }
        )

    counts = pd.concat(count_frames, ignore_index=True)
    counts_path = args.outdir / "label_free_motif_pool_counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "study": {"path": str(args.study), "sha256": file_sha256(args.study)},
        "task_study": {
            "path": str(task_study_path),
            "sha256": file_sha256(task_study_path),
        },
        "labels_used": False,
        "selection_split": "train",
        "maximum_per_tf": int(args.maximum_per_tf),
        "seed": int(args.seed),
        "outputs": records,
        "counts": {
            "path": str(counts_path),
            "sha256": file_sha256(counts_path),
        },
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["pool_id"] = sha256(canonical.encode()).hexdigest()
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(counts.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
