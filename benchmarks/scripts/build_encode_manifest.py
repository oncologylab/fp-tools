#!/usr/bin/env python3
"""Build a conservative ENCODE manifest for fp-tools benchmarks.

The script queries released human GRCh38 ENCODE experiments and writes a TSV
manifest. It does not download data.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ENCODE = "https://www.encodeproject.org"
HEADERS = {"accept": "application/json"}
DEFAULT_CELL_TYPES = ["GM12878", "K562", "HepG2", "A549", "H1", "HCT116"]
DEFAULT_TFS = ["CTCF", "SPI1", "GATA1", "JUNB", "REST", "NRF1", "MAX"]
FIELDNAMES = [
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


def fetch_json(path_or_url: str) -> dict:
    if path_or_url.startswith("http"):
        url = path_or_url
    else:
        url = ENCODE + path_or_url
    separator = "&" if "?" in url else "?"
    request = Request(f"{url}{separator}format=json", headers=HEADERS)
    with urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def resolve_files(file_refs: list[dict], cache: dict[str, dict] | None = None) -> list[dict]:
    cache = cache if cache is not None else {}
    files = []
    for ref in file_refs:
        file_id = ref.get("@id") or ref.get("href")
        if not file_id:
            continue
        if file_id not in cache:
            cache[file_id] = fetch_json(file_id)
        files.append(cache[file_id])
    return files


def search(params: dict[str, str]) -> dict:
    query = urlencode({**params, "format": "json", "limit": "all"})
    request = Request(f"{ENCODE}/search/?{query}", headers=HEADERS)
    with urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def iter_experiments(cell_type: str, assay_title: str, tf: str | None = None):
    params = {
        "type": "Experiment",
        "assay_title": assay_title,
        "status": "released",
        "assembly": "GRCh38",
        "biosample_ontology.term_name": cell_type,
    }
    if tf:
        params["target.label"] = tf
    yield from search(params).get("@graph", [])


def choose_file(files: list[dict], output_type: str, file_format: str, assembly: str = "GRCh38") -> dict | None:
    candidates = []
    for item in files:
        if item.get("status") != "released":
            continue
        if item.get("output_type") != output_type:
            continue
        if item.get("file_format") != file_format:
            continue
        item_assembly = item.get("assembly")
        if item_assembly not in (assembly, None):
            continue
        candidates.append(item)
    candidates.sort(
        key=lambda item: (item.get("assembly") == assembly, item.get("date_created", ""), item.get("accession", "")),
        reverse=True,
    )
    return candidates[0] if candidates else None


def make_row(cell_type: str, tf: str, assay: str, experiment: dict, file_item: dict, output_type: str, split: str) -> dict[str, str]:
    href = file_item.get("href", "")
    accession = file_item.get("accession", "")
    file_format = file_item.get("file_format", "")
    download_name = href.rstrip("/").rsplit("/", 1)[-1] if href else f"{accession}.{file_format}"
    local_name = f"{cell_type}.{tf}.{assay}.{download_name}".replace(" ", "_")
    return {
        "source": "ENCODE",
        "benchmark_tier": "bulk",
        "cell_type": cell_type,
        "donor": "",
        "tf": tf,
        "assay": assay,
        "experiment_accession": experiment.get("accession", ""),
        "file_accession": accession,
        "assembly": "GRCh38",
        "output_type": output_type,
        "file_format": file_format,
        "url": ENCODE + href if href else "",
        "checksum": file_item.get("md5sum", ""),
        "status": file_item.get("status", ""),
        "local_path": f"data/public/raw/encode/{local_name}",
        "split": split,
        "notes": "auto_selected_latest_released_file",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-types", nargs="*", default=DEFAULT_CELL_TYPES)
    parser.add_argument("--tfs", nargs="*", default=DEFAULT_TFS)
    parser.add_argument("--out", default="benchmarks/manifests/encode_bulk_manifest.tsv")
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    file_cache: dict[str, dict] = {}
    for cell_type in args.cell_types:
        atac_experiments = list(iter_experiments(cell_type, "ATAC-seq"))
        for tf in args.tfs:
            chip_experiments = list(iter_experiments(cell_type, "TF ChIP-seq", tf=tf))
            if not atac_experiments or not chip_experiments:
                continue
            atac = atac_experiments[0]
            chip = chip_experiments[0]
            atac_files = resolve_files(atac.get("files", []), file_cache)
            chip_files = resolve_files(chip.get("files", []), file_cache)
            selected = [
                ("ATAC-seq", atac, choose_file(atac_files, "alignments", "bam"), "alignments"),
                ("ATAC-seq", atac, choose_file(atac_files, "IDR thresholded peaks", "bed"), "IDR thresholded peaks"),
                ("TF ChIP-seq", chip, choose_file(chip_files, "IDR thresholded peaks", "bed"), "IDR thresholded peaks"),
            ]
            if any(file_item is None for _, _, file_item, _ in selected):
                continue
            for assay, experiment, file_item, output_type in selected:
                rows.append(make_row(cell_type, tf, assay, experiment, file_item, output_type, "not_assigned"))
            time.sleep(0.2)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="	")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
