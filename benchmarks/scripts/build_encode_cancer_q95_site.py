#!/usr/bin/env python3
"""Publish compact, canonical differential-report payloads to the static site."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd

from fp_tools.utils.motif_databases import motif_db_path
from fp_tools.utils.motifs import MotifList


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "data/public/processed/encode_cancer_pairwise_q95_20260814"
SITE = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
MANIFEST = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.tsv"
COMPARISONS = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814_comparisons.tsv"
EXPECTED_MOTIFS = 1019
EXPECTED_PAIRS = 21
RELEASE_DATE = "2026-08-14"
REFERENCE_SCIENTIFIC_SHA256 = "72e545e4a1324edc5b172b3206105e60f8d5b77c7fca5032addcf81b9466a6ff"


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def scientific_digest(payload: dict) -> str:
    scientific = {
        key: payload[key]
        for key in ("title", "report_label", "conditions", "groups", "colors", "points", "aggregate", "change_label")
    }
    raw = json.dumps(scientific, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def read_payload(path: Path) -> dict:
    payload = json.loads(gzip.decompress(path.read_bytes()))
    if len(payload.get("points", [])) != EXPECTED_MOTIFS:
        raise ValueError(f"{path} does not contain {EXPECTED_MOTIFS} motifs")
    if len(payload.get("conditions", [])) != 2 or not payload.get("aggregate", {}).get("motifs"):
        raise ValueError(f"{path} is not a complete differential-report payload")
    return payload


def available_records(project: Path) -> list[dict]:
    comparisons = pd.read_csv(COMPARISONS, sep="\t", dtype=str)
    records = []
    for row in comparisons.itertuples(index=False):
        source = project / "pairs" / row.comparison / "results/report_payload.json.gz"
        marker_path = project / "pairs" / row.comparison / "complete.json"
        if not source.is_file() or not marker_path.is_file():
            continue
        payload = read_payload(source)
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
        if marker.get("payload_sha256") != sha256(source):
            raise ValueError(f"Checksum mismatch for {row.comparison}")
        records.append({
            "comparison": row.comparison,
            "condition1": payload["conditions"][0],
            "condition2": payload["conditions"][1],
            "source": source,
            "payload": payload,
            "marker": marker,
        })
    return records


def build(*, project: Path, site: Path, allow_partial: bool) -> None:
    records = available_records(project)
    if not allow_partial and len(records) != EXPECTED_PAIRS:
        raise ValueError(f"Expected {EXPECTED_PAIRS} completed pairs, found {len(records)}")
    reference = next((record for record in records if record["comparison"] == "HepG2_vs_K562"), None)
    if reference is None or scientific_digest(reference["payload"]) != REFERENCE_SCIENTIFIC_SHA256:
        raise ValueError("The browser reference comparison does not exactly match the preserved report")
    samples = pd.read_csv(MANIFEST, sep="\t", dtype=str, keep_default_na=False)
    data = site / "data"
    staging = site / "data.q95-staging"
    if staging.exists():
        shutil.rmtree(staging)
    reports = staging / "reports"
    reports.mkdir(parents=True)
    metadata_records = []
    long_rows = []
    for record in records:
        destination = reports / f"{record['comparison']}.json.gz"
        shutil.copy2(record["source"], destination)
        payload = record["payload"]
        aggregate = {item["prefix"]: item for item in payload["aggregate"]["motifs"]}
        for point in payload["points"]:
            long_rows.append({
                "comparison": record["comparison"],
                "condition1": payload["conditions"][0],
                "condition2": payload["conditions"][1],
                "prefix": point["prefix"],
                "name": point["name"],
                "motif_id": point.get("motif_id", ""),
                "group": point["group"],
                "n_sites": aggregate.get(point["prefix"], {}).get("n_sites", ""),
                "effect": point["change"],
                "pvalue": point["pvalue"],
                "qvalue": point["fdr"],
            })
        metadata_records.append({
            "comparison": record["comparison"],
            "condition1": payload["conditions"][0],
            "condition2": payload["conditions"][1],
            "file": f"data/reports/{record['comparison']}.json.gz",
            "payload_sha256": sha256(destination),
            "motifs": len(payload["points"]),
            "aggregate_motifs": len(payload["aggregate"]["motifs"]),
        })
    frame = pd.DataFrame(long_rows)
    if not frame.empty:
        frame.to_csv(staging / "all_pairwise_results.tsv.gz", sep="\t", index=False, compression="gzip")
    conditions = [
        {
            "name": condition,
            "samples": group.sort_values("display_order")["sample"].tolist(),
            "experiment": group["experiment"].iloc[0],
        }
        for condition, group in samples.groupby("condition", sort=True)
    ]
    metadata = {
        "schema": "fp-tools.encode-cancer-static-browser.v2",
        "release_date": RELEASE_DATE,
        "method": "Pair-specific released IDR-peak union; corrected cut-site q95 scaling; fp-tools footprint scoring",
        "conditions": conditions,
        "comparisons": metadata_records,
        "downloads": {"all_results": "data/all_pairwise_results.tsv.gz"},
        "reference": {
            "comparison": "K562 vs HepG2",
            "uncompressed_scientific_sha256": REFERENCE_SCIENTIFIC_SHA256,
            "points": 1019,
            "aggregate_motifs": 1009,
        },
    }
    matrices = {}
    for motif in MotifList().from_file(str(motif_db_path("jaspar2026_vertebrates"))):
        motif.set_prefix("name_id")
        matrices[motif.prefix] = [
            [round(float(value), 4) for value in row]
            for row in motif.counts
        ]
    if len(matrices) != EXPECTED_MOTIFS or (not frame.empty and set(matrices) != set(frame["prefix"])):
        raise ValueError("JASPAR motif matrices do not match the report payloads")
    (staging / "motif_matrices.json").write_text(
        json.dumps({"schema": "fp-tools.motif-matrices.v1", "motifs": matrices}, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    (staging / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if data.exists():
        shutil.rmtree(data)
    staging.replace(data)


def verify(site: Path, *, allow_partial: bool) -> None:
    metadata = json.loads((site / "data/metadata.json").read_text(encoding="utf-8"))
    expected = len(metadata["comparisons"])
    if not allow_partial and expected != EXPECTED_PAIRS:
        raise ValueError(f"Static site exposes {expected}, not {EXPECTED_PAIRS}, comparisons")
    if len(metadata["conditions"]) != 7 or sum(len(item["samples"]) for item in metadata["conditions"]) != 17:
        raise ValueError("Static site does not expose the locked seven-line, 17-sample design")
    matrices = json.loads((site / "data/motif_matrices.json").read_text(encoding="utf-8"))["motifs"]
    if len(matrices) != EXPECTED_MOTIFS:
        raise ValueError("Static site does not contain all 1,019 motif matrices")
    for record in metadata["comparisons"]:
        path = site / record["file"]
        if sha256(path) != record["payload_sha256"]:
            raise ValueError(f"Static payload checksum mismatch: {record['comparison']}")
        read_payload(path)
    reference = read_payload(site / "data/reports/HepG2_vs_K562.json.gz")
    if scientific_digest(reference) != REFERENCE_SCIENTIFIC_SHA256:
        raise ValueError("Static K562/HepG2 payload differs from the preserved report")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("build", "verify"))
    parser.add_argument("--project", type=Path, default=PROJECT)
    parser.add_argument("--site", type=Path, default=SITE)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "build":
        build(project=args.project, site=args.site, allow_partial=args.allow_partial)
    verify(args.site, allow_partial=args.allow_partial)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
