#!/usr/bin/env python3
"""Publish compact, canonical differential-report payloads to the static site."""

from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "data/public/processed/encode_cancer_pairwise_q95_20260814"
SITE = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
MANIFEST = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.tsv"
COMPARISONS = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814_comparisons.tsv"
REFERENCE_REPORT = ROOT / "docs/demos/reports/diff_footprints_K562_HepG2.html"
EXPECTED_MOTIFS = 1019
EXPECTED_PAIRS = 21
PROFILE_SHARDS = 16
RELEASE_DATE = "2026-08-14"
REFERENCE_SCIENTIFIC_SHA256 = "ae9c8abca29096a0f5b10bbb0952e1dd41c1c55a2344075ba0e17822b401812a"


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


def reference_logo_pngs() -> dict[str, bytes]:
    match = re.search(
        r'const reportPayloadB64="([^"]+)"',
        REFERENCE_REPORT.read_text(encoding="utf-8"),
    )
    if not match:
        raise ValueError("The preserved report does not contain its payload")
    payload = json.loads(gzip.decompress(base64.b64decode(match.group(1))))
    logos = {}
    for prefix, record in payload.get("logos", {}).items():
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", prefix):
            raise ValueError(f"Unsafe motif logo prefix: {prefix}")
        uri = record.get("png", "")
        if not uri.startswith("data:image/png;base64,"):
            raise ValueError(f"The preserved logo is not PNG: {prefix}")
        image = base64.b64decode(uri.split(",", 1)[1])
        if not image.startswith(b"\x89PNG\r\n\x1a\n"):
            raise ValueError(f"Invalid preserved PNG logo: {prefix}")
        logos[prefix] = image
    if len(logos) != EXPECTED_MOTIFS:
        raise ValueError(f"Expected {EXPECTED_MOTIFS} preserved logos, found {len(logos)}")
    return logos


def write_gzip_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    path.write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))


def profile_shard(prefix: str) -> int:
    return int(hashlib.sha256(prefix.encode("utf-8")).hexdigest()[:2], 16) % PROFILE_SHARDS


def split_browser_payload(payload: dict) -> tuple[dict, list[list[dict]]]:
    """Separate the fast comparison payload from lazily loaded profiles."""
    core = {key: value for key, value in payload.items() if key not in {"aggregate", "motif_matrices"}}
    summaries = []
    shards: list[list[dict]] = [[] for _ in range(PROFILE_SHARDS)]
    for motif in payload["aggregate"]["motifs"]:
        shard = profile_shard(motif["prefix"])
        conditions = []
        for condition in motif["conditions"]:
            conditions.append({
                **{key: value for key, value in condition.items() if key != "profile"},
                "samples": [
                    {key: value for key, value in sample.items() if key != "profile"}
                    for sample in condition["samples"]
                ],
            })
        summaries.append({
            **{key: value for key, value in motif.items() if key != "conditions"},
            "conditions": conditions,
            "profile_shard": shard,
        })
        shards[shard].append(motif)
    core["aggregate"] = {
        **{key: value for key, value in payload["aggregate"].items() if key != "motifs"},
        "motifs": summaries,
    }
    return core, shards


def read_browser_payload(path: Path) -> dict:
    payload = json.loads(gzip.decompress(path.read_bytes()))
    if len(payload.get("points", [])) != EXPECTED_MOTIFS:
        raise ValueError(f"{path} does not contain {EXPECTED_MOTIFS} motifs")
    motifs = payload.get("aggregate", {}).get("motifs", [])
    if not motifs or any("profile_shard" not in motif for motif in motifs):
        raise ValueError(f"{path} is not a complete compact browser payload")
    return payload


def reconstruct_browser_payload(site: Path, record: dict) -> dict:
    payload = read_browser_payload(site / record["file"])
    profiles = {}
    for shard in record["profile_shards"]:
        path = site / shard["file"]
        if sha256(path) != shard["sha256"]:
            raise ValueError(f"Static profile checksum mismatch: {record['comparison']} shard {shard['id']}")
        shard_payload = json.loads(gzip.decompress(path.read_bytes()))
        if len(shard_payload.get("motifs", [])) != shard["motifs"]:
            raise ValueError(f"Static profile count mismatch: {record['comparison']} shard {shard['id']}")
        profiles.update({motif["prefix"]: motif for motif in shard_payload["motifs"]})
    prefixes = [motif["prefix"] for motif in payload["aggregate"]["motifs"]]
    if set(prefixes) != set(profiles):
        raise ValueError(f"Static profiles are incomplete for {record['comparison']}")
    payload["aggregate"]["motifs"] = [profiles[prefix] for prefix in prefixes]
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
    profile_root = staging / "profiles"
    logo_root = staging / "logos"
    logo_root.mkdir(parents=True)
    metadata_records = []
    long_rows = []
    for record in records:
        payload = record["payload"]
        compact, profile_shards = split_browser_payload(payload)
        destination = reports / f"{record['comparison']}.json.gz"
        write_gzip_json(compact, destination)
        shard_records = []
        for shard_id, motifs in enumerate(profile_shards):
            shard_path = profile_root / record["comparison"] / f"{shard_id:02x}.json.gz"
            write_gzip_json({"motifs": motifs}, shard_path)
            shard_records.append({
                "id": shard_id,
                "file": str(shard_path.relative_to(staging)),
                "sha256": sha256(shard_path),
                "motifs": len(motifs),
            })
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
            "source_payload_sha256": record["marker"]["payload_sha256"],
            "profile_shards": [
                {**shard, "file": f"data/{shard['file']}"}
                for shard in shard_records
            ],
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
    logos = reference_logo_pngs()
    for prefix, image in logos.items():
        (logo_root / f"{prefix}.png").write_bytes(image)
    metadata = {
        "schema": "fp-tools.encode-cancer-static-browser.v2",
        "release_date": RELEASE_DATE,
        "method": "Pair-specific released IDR-peak union; corrected cut-site q95 scaling; fp-tools footprint scoring",
        "conditions": conditions,
        "comparisons": metadata_records,
        "downloads": {"all_results": "data/all_pairwise_results.tsv.gz"},
        "logos": {
            "base": "data/logos",
            "format": "png",
            "count": len(logos),
            "source": "preserved K562 versus HepG2 report",
        },
        "reference": {
            "comparison": "K562 vs HepG2",
            "uncompressed_scientific_sha256": REFERENCE_SCIENTIFIC_SHA256,
            "points": 1019,
            "aggregate_motifs": 1019,
        },
    }
    if not frame.empty and set(logos) != set(frame["prefix"]):
        raise ValueError("Preserved motif logos do not match the report payloads")
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
    expected_logos = reference_logo_pngs()
    if metadata.get("logos", {}).get("count") != EXPECTED_MOTIFS:
        raise ValueError("Static site does not declare all 1,019 motif logos")
    for prefix, expected in expected_logos.items():
        if (site / f"data/logos/{prefix}.png").read_bytes() != expected:
            raise ValueError(f"Static motif logo differs from the preserved report: {prefix}")
    for record in metadata["comparisons"]:
        path = site / record["file"]
        if sha256(path) != record["payload_sha256"]:
            raise ValueError(f"Static payload checksum mismatch: {record['comparison']}")
        read_browser_payload(path)
        reconstruct_browser_payload(site, record)
    reference_record = next(record for record in metadata["comparisons"] if record["comparison"] == "HepG2_vs_K562")
    reference = reconstruct_browser_payload(site, reference_record)
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
