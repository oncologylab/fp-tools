#!/usr/bin/env python3
"""Discover and optionally download ENCODE ChIP peak files for study tasks."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import requests


ENCODE = "https://www.encodeproject.org"


def request_json(session: requests.Session, url: str, params: dict[str, str] | None = None) -> dict[str, Any]:
    response = session.get(url, params=params, headers={"accept": "application/json"}, timeout=60)
    response.raise_for_status()
    return response.json()


def candidate_rank(row: dict[str, object]) -> tuple[object, ...]:
    output_rank = {
        "optimal IDR thresholded peaks": 3,
        "IDR thresholded peaks": 2,
        "conservative IDR thresholded peaks": 1,
    }.get(str(row["output_type"]), 0)
    return (
        not bool(row["perturbed"]),
        output_rank,
        bool(row["preferred_default"]),
        int(row["biological_replicate_count"]),
        str(row["date_created"]),
        str(row["file_accession"]),
    )


def discover_task(
    session: requests.Session,
    cell: str,
    tf: str,
    assembly: str,
) -> list[dict[str, object]]:
    graph = request_json(
        session,
        f"{ENCODE}/search/",
        {
            "type": "Experiment",
            "assay_title": "TF ChIP-seq",
            "biosample_ontology.term_name": cell,
            "target.label": tf,
            "status": "released",
            "format": "json",
            "limit": "all",
        },
    ).get("@graph", [])
    rows = []
    for summary in graph:
        experiment = request_json(session, f"{ENCODE}{summary['@id']}?format=json")
        perturbed = bool(experiment.get("perturbed", False))
        for record in experiment.get("files", []):
            if not isinstance(record, dict):
                continue
            output_type = str(record.get("output_type", ""))
            if (
                record.get("status") != "released"
                or record.get("file_format") != "bed"
                or record.get("assembly") != assembly
                or "IDR thresholded peaks" not in output_type
            ):
                continue
            href = str(record.get("href") or f"/files/{record['accession']}/@@download/{record['accession']}.bed.gz")
            rows.append(
                {
                    "cell": cell,
                    "tf": tf,
                    "experiment_accession": str(experiment["accession"]),
                    "file_accession": str(record["accession"]),
                    "output_type": output_type,
                    "assembly": assembly,
                    "biological_replicates": ",".join(map(str, record.get("biological_replicates", []))),
                    "biological_replicate_count": len(record.get("biological_replicates", [])),
                    "preferred_default": bool(record.get("preferred_default", False)),
                    "perturbed": perturbed,
                    "date_created": str(record.get("date_created", "")),
                    "md5sum": str(record.get("md5sum", "")),
                    "url": href if href.startswith("http") else f"{ENCODE}{href}",
                }
            )
    return rows


def select_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    selected = []
    for _, group in candidates.groupby(["cell", "tf"], sort=True):
        records = group.to_dict("records")
        chosen = max(records, key=candidate_rank)
        chosen["selection_rank"] = "unperturbed,optimal-IDR,preferred,replicates,date"
        selected.append(chosen)
    return pd.DataFrame(selected)


def download_selected(session: requests.Session, selected: pd.DataFrame, directory: Path) -> pd.DataFrame:
    directory.mkdir(parents=True, exist_ok=True)
    output = selected.copy()
    paths = []
    for row in output.itertuples(index=False):
        path = directory / f"{row.cell}.{row.tf}.{row.file_accession}.bed.gz"
        if not path.is_file():
            response = session.get(str(row.url), stream=True, timeout=120)
            response.raise_for_status()
            with path.open("wb") as handle:
                for block in response.iter_content(8 * 1024 * 1024):
                    if block:
                        handle.write(block)
        digest = hashlib.md5(path.read_bytes()).hexdigest()  # noqa: S324 - verifies ENCODE metadata
        if row.md5sum and digest != row.md5sum:
            raise ValueError(f"MD5 mismatch for {row.file_accession}: {digest} != {row.md5sum}")
        paths.append(str(path))
    output["local_path"] = paths
    return output


def extract_task_motifs(study: dict[str, object], database: Path, output: Path, split: str) -> int:
    wanted = {str(task["motif_id"]) for task in study["tasks"] if task["split"] == split}
    blocks = []
    current: list[str] = []
    with database.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">") and current:
                blocks.append(current)
                current = []
            current.append(line)
    if current:
        blocks.append(current)
    selected = [block for block in blocks if block and block[0][1:].split()[0] in wanted]
    found = {block[0][1:].split()[0] for block in selected}
    missing = sorted(wanted.difference(found))
    if missing:
        raise ValueError("motif database is missing task motifs: " + ", ".join(missing))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for block in selected:
            handle.writelines(block)
    return len(selected)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--split", default="development")
    parser.add_argument("--motif-database", type=Path)
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args(argv)

    study = json.loads(args.study.read_text(encoding="utf-8"))
    tasks = [task for task in study["tasks"] if task["split"] == args.split]
    session = requests.Session()
    rows = []
    for task in tasks:
        rows.extend(discover_task(session, str(task["cell"]), str(task["tf"]), str(study["assembly"])))
    candidates = pd.DataFrame(rows)
    if candidates.empty:
        raise RuntimeError("ENCODE returned no matching peak files")
    missing_tasks = sorted(
        {(str(task["cell"]), str(task["tf"])) for task in tasks}
        .difference(set(zip(candidates["cell"], candidates["tf"])))
    )
    if missing_tasks:
        raise RuntimeError("ENCODE peak files are absent for: " + ", ".join(f"{cell}/{tf}" for cell, tf in missing_tasks))
    selected = select_candidates(candidates)
    args.outdir.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(args.outdir / "encode_chip_peak_candidates.tsv", sep="\t", index=False)
    if args.download:
        selected = download_selected(session, selected, args.outdir / "peaks")
    selected.to_csv(args.outdir / "encode_chip_peak_selected.tsv", sep="\t", index=False)
    if args.motif_database:
        count = extract_task_motifs(
            study, args.motif_database, args.outdir / f"{args.split}_task_motifs.jaspar", args.split
        )
        print(f"wrote {count} unique task motifs")
    print(selected[["cell", "tf", "experiment_accession", "file_accession", "output_type"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
