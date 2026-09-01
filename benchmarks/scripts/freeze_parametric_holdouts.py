#!/usr/bin/env python3
"""Freeze metadata-only ENCODE inputs for the parametric promotion holdout.

This selector deliberately does not download or inspect ChIP peak contents.
It locks file accessions, sizes, MD5 values, and deterministic selection
decisions before labels can be opened.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from typing import Any, Iterable

import requests


ENCODE = "https://www.encodeproject.org"
FULL_COLUMNS = [
    "source",
    "benchmark_tier",
    "cell_type",
    "donor",
    "tf",
    "assay",
    "experiment_accession",
    "file_accession",
    "assembly",
    "output_type",
    "file_format",
    "url",
    "checksum",
    "status",
    "local_path",
    "split",
    "notes",
]


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def request_json(
    session: requests.Session,
    url: str,
    params: dict[str, str] | None = None,
) -> dict[str, Any]:
    response = session.get(
        url,
        params=params,
        headers={"accept": "application/json"},
        timeout=60,
    )
    response.raise_for_status()
    return response.json()


def _analysis_step(record: dict[str, Any]) -> tuple[str, tuple[int, ...], str]:
    version = record.get("analysis_step_version") or {}
    if not isinstance(version, dict):
        return "", (), ""
    step = version.get("analysis_step") or {}
    if not isinstance(step, dict):
        step = {}
    label = str(step.get("step_label") or "")
    major = int(step.get("major_version") or 0)
    minor = int(version.get("minor_version") or 0)
    source_urls = []
    for pipeline in step.get("pipelines", []):
        if isinstance(pipeline, dict):
            source_urls.append(str(pipeline.get("source_url") or ""))
    numeric = tuple(
        int(value)
        for value in re.findall(r"\d+", " ".join(source_urls))[-4:]
    )
    return label, (major, minor, *numeric), str(version.get("name") or label)


def is_replicate_aware_idr(record: dict[str, Any]) -> bool:
    output_type = str(record.get("output_type") or "")
    replicates = record.get("biological_replicates") or []
    label, _version, _display = _analysis_step(record)
    if "IDR thresholded peaks" not in output_type or len(set(replicates)) < 2:
        return False
    if "pseudoreplicated" in label:
        return False
    return "replicated-idr" in label or label == "tf-idr-step"


def candidate_rank(record: dict[str, Any]) -> tuple[object, ...]:
    """Locked tie break: preferred, pipeline version, then accession."""

    _label, version, _display = _analysis_step(record)
    return (
        bool(record.get("preferred_default", False)),
        version,
        str(record.get("accession") or ""),
    )


def _eligible_file(record: Any, *, assembly: str, file_format: str) -> bool:
    return bool(
        isinstance(record, dict)
        and record.get("status") == "released"
        and record.get("assembly") == assembly
        and record.get("file_format") == file_format
    )


def _download_url(record: dict[str, Any]) -> str:
    accession = str(record["accession"])
    href = str(record.get("href") or "")
    if href:
        return href if href.startswith("http") else f"{ENCODE}{href}"
    suffix = "bam" if record.get("file_format") == "bam" else "bed.gz"
    return f"{ENCODE}/files/{accession}/@@download/{accession}.{suffix}"


def _manifest_row(
    record: dict[str, Any],
    *,
    cell: str,
    assay: str,
    experiment: str,
    tf: str = "",
    note: str,
) -> dict[str, str]:
    accession = str(record["accession"])
    suffix = "bam" if record.get("file_format") == "bam" else "bed.gz"
    target = (
        f"{cell}.{tf}.{accession}.{suffix}"
        if tf
        else f"{cell}.{accession}.{suffix}"
    )
    return {
        "source": "ENCODE",
        "benchmark_tier": "external_locked",
        "cell_type": cell,
        "donor": "",
        "tf": tf,
        "assay": assay,
        "experiment_accession": experiment,
        "file_accession": accession,
        "assembly": str(record.get("assembly") or "GRCh38"),
        "output_type": str(record.get("output_type") or ""),
        "file_format": suffix,
        "url": _download_url(record),
        "checksum": str(record.get("md5sum") or ""),
        "status": str(record.get("status") or "released"),
        "local_path": f"data/public/raw/frozen_parametric_holdout/{target}",
        "split": "external_locked",
        "notes": f"{note}; {int(record.get('file_size') or 0)} bytes",
    }


def select_atac_inputs(
    session: requests.Session,
    cell: str,
    experiment_accession: str,
    assembly: str,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    experiment = request_json(
        session,
        f"{ENCODE}/experiments/{experiment_accession}/?format=json",
    )
    if experiment.get("status") != "released" or experiment.get("assay_title") != "ATAC-seq":
        raise ValueError(f"{experiment_accession} is not a released ATAC-seq experiment")
    records = [record for record in experiment.get("files", []) if isinstance(record, dict)]
    alignments = [
        record
        for record in records
        if _eligible_file(record, assembly=assembly, file_format="bam")
        and record.get("output_type") == "alignments"
        and len(record.get("biological_replicates") or []) == 1
    ]
    by_replicate: dict[int, list[dict[str, Any]]] = {}
    for record in alignments:
        replicate = int((record.get("biological_replicates") or [0])[0])
        by_replicate.setdefault(replicate, []).append(record)
    if len(by_replicate) < 2:
        raise ValueError(f"{experiment_accession} has fewer than two released GRCh38 BAM replicates")
    selected_bams = [max(by_replicate[key], key=candidate_rank) for key in sorted(by_replicate)]
    peaks = [
        record
        for record in records
        if _eligible_file(record, assembly=assembly, file_format="bed")
        and record.get("output_type") == "conservative IDR thresholded peaks"
        and len(set(record.get("biological_replicates") or [])) >= 2
    ]
    if not peaks:
        raise ValueError(f"{experiment_accession} has no replicate-aware conservative IDR peak file")
    selected_peak = max(peaks, key=candidate_rank)
    rows = [
        _manifest_row(
            record,
            cell=cell,
            assay="ATAC-seq",
            experiment=experiment_accession,
            note=f"coordinate-sorted biological replicate {(record.get('biological_replicates') or [''])[0]}",
        )
        for record in selected_bams
    ]
    rows.append(
        _manifest_row(
            selected_peak,
            cell=cell,
            assay="ATAC-seq",
            experiment=experiment_accession,
            note="label-free replicate-aware conservative IDR peak universe",
        )
    )
    return rows, {
        "cell": cell,
        "experiment_accession": experiment_accession,
        "selected_bams": [str(record["accession"]) for record in selected_bams],
        "selected_peak": str(selected_peak["accession"]),
        "labels_read": False,
    }


def discover_chip_records(
    session: requests.Session,
    cell: str,
    tf: str,
) -> list[tuple[str, dict[str, Any]]]:
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
    candidates: list[tuple[str, dict[str, Any]]] = []
    for summary in graph:
        experiment = request_json(session, f"{ENCODE}{summary['@id']}?format=json")
        if experiment.get("status") != "released" or bool(experiment.get("perturbed", False)):
            continue
        if experiment.get("treatments"):
            continue
        for record in experiment.get("files", []):
            if _eligible_file(record, assembly="GRCh38", file_format="bed") and is_replicate_aware_idr(record):
                candidates.append((str(experiment["accession"]), record))
    return candidates


def select_chip_input(
    session: requests.Session,
    task: dict[str, Any],
) -> tuple[dict[str, str] | None, dict[str, Any]]:
    cell = str(task["cell"])
    tf = str(task["tf"])
    motif_id = task.get("motif_id")
    if not motif_id:
        return None, {
            "cell": cell,
            "tf": tf,
            "status": "ineligible_no_jaspar_motif",
            "labels_read": False,
        }
    candidates = discover_chip_records(session, cell, tf)
    if not candidates:
        return None, {
            "cell": cell,
            "tf": tf,
            "motif_id": motif_id,
            "status": "ineligible_no_untreated_replicate_aware_idr",
            "labels_read": False,
        }
    experiment, selected = max(candidates, key=lambda item: candidate_rank(item[1]))
    label, version, display = _analysis_step(selected)
    row = _manifest_row(
        selected,
        cell=cell,
        assay="TF ChIP-seq",
        experiment=experiment,
        tf=tf,
        note=(
            f"evaluation label only; motif {motif_id}; family {task['motif_family']}; "
            f"replicate-aware step {label}; selection preferred,pipeline,accession"
        ),
    )
    return row, {
        "cell": cell,
        "tf": tf,
        "motif_id": motif_id,
        "motif_family": str(task["motif_family"]),
        "role": str(task["role"]),
        "status": "selected_pending_power_check",
        "experiment_accession": experiment,
        "file_accession": str(selected["accession"]),
        "output_type": str(selected.get("output_type") or ""),
        "biological_replicates": list(selected.get("biological_replicates") or []),
        "preferred_default": bool(selected.get("preferred_default", False)),
        "pipeline_version": list(version),
        "pipeline_label": display,
        "candidate_count": len(candidates),
        "labels_read": False,
    }


def write_manifest(rows: Iterable[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FULL_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def freeze(
    study_path: Path,
    manifest_path: Path,
    freeze_path: Path,
    *,
    session: requests.Session | None = None,
) -> dict[str, Any]:
    study = json.loads(study_path.read_text(encoding="utf-8"))
    if study.get("status") != "development_locked_holdout_unscored":
        raise ValueError("holdout metadata can only be frozen for an unscored study")
    active_session = session or requests.Session()
    rows: list[dict[str, str]] = []
    atac_decisions = []
    for item in study["new_holdout_atac"]:
        selected_rows, decision = select_atac_inputs(
            active_session,
            str(item["cell"]),
            str(item["experiment_accession"]),
            str(study["assembly"]),
        )
        rows.extend(selected_rows)
        atac_decisions.append(decision)
    chip_decisions = []
    for task in study["new_holdout_tasks"]:
        row, decision = select_chip_input(active_session, task)
        if row is not None:
            rows.append(row)
        chip_decisions.append(decision)
    rows.sort(key=lambda row: (row["cell_type"], row["assay"], row["tf"], row["file_accession"]))
    write_manifest(rows, manifest_path)
    selector_path = Path(__file__).resolve()
    document = {
        "schema": "fp-tools-parametric-holdout-freeze-v1",
        "study": {"path": str(study_path), "sha256": sha256_file(study_path)},
        "selector": {"path": str(selector_path), "sha256": sha256_file(selector_path)},
        "manifest": {
            "path": str(manifest_path),
            "sha256": sha256_file(manifest_path),
            "rows": len(rows),
        },
        "selection_rule": study["holdout_selection"],
        "atac_decisions": atac_decisions,
        "chip_decisions": chip_decisions,
        "chipped_peak_contents_read": False,
        "holdout_labels_scored": False,
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["freeze_id"] = sha256(canonical.encode()).hexdigest()
    freeze_path.parent.mkdir(parents=True, exist_ok=True)
    freeze_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return document


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--freeze", dest="freeze_path", type=Path, required=True)
    args = parser.parse_args(argv)
    document = freeze(args.study, args.manifest, args.freeze_path)
    print(json.dumps({"freeze_id": document["freeze_id"], "rows": document["manifest"]["rows"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
