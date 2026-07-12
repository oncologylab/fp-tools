#!/usr/bin/env python3
"""Resolve the curated LCMV CD8 GSM collection to current ENA FASTQs."""

from __future__ import annotations

import argparse
import csv
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


EXPECTED = {
    "libraries": 58,
    "primary_atac": 27,
    "supplemental_atac": 6,
    "rna": 25,
    "condition_pairs": 9,
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {field: row.get(field, "") for field in fields} for row in rows
        )


def validate_selection(rows: list[dict[str, str]]) -> None:
    if len(rows) != EXPECTED["libraries"]:
        raise ValueError(f"Expected 58 curated GSM libraries, found {len(rows)}")
    gsms = [row["gsm_accession"] for row in rows]
    if len(set(gsms)) != len(gsms):
        raise ValueError("The curated selection contains duplicate GSM accessions")
    primary_atac = [
        row for row in rows if row["assay"] == "ATAC" and row["collection"] == "primary"
    ]
    supplemental_atac = [
        row
        for row in rows
        if row["assay"] == "ATAC" and row["collection"] == "supplemental"
    ]
    rna = [row for row in rows if row["assay"] == "RNA"]
    observed = {
        "primary_atac": len(primary_atac),
        "supplemental_atac": len(supplemental_atac),
        "rna": len(rna),
    }
    for key, expected in EXPECTED.items():
        if key in observed and observed[key] != expected:
            raise ValueError(
                f"Expected {expected} {key} libraries, found {observed[key]}"
            )
    groups: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["collection"] == "primary":
            groups[row["condition_pair_id"]].add(row["assay"])
        elif row["include_in_primary_paired_analysis"] != "false":
            raise ValueError(
                "Supplemental libraries must be excluded from paired analysis"
            )
    if len(groups) != EXPECTED["condition_pairs"]:
        raise ValueError(f"Expected 9 primary condition pairs, found {len(groups)}")
    missing = [group for group, assays in groups.items() if assays != {"ATAC", "RNA"}]
    if missing:
        raise ValueError(
            f"Condition pairs without both assays: {', '.join(sorted(missing))}"
        )


def resolve_sra(gsms: list[str]) -> dict[str, list[dict[str, str]]]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    term = " OR ".join(f"{gsm}[All Fields]" for gsm in gsms)
    search_url = (
        base
        + "esearch.fcgi?"
        + urllib.parse.urlencode(
            {"db": "sra", "term": term, "retmode": "json", "retmax": 200}
        )
    )
    with urllib.request.urlopen(search_url, timeout=120) as response:
        result = json.load(response)["esearchresult"]
    if int(result["count"]) != len(gsms):
        raise RuntimeError(
            f"NCBI returned {result['count']} SRA records for {len(gsms)} curated GSMs"
        )
    fetch_url = (
        base
        + "efetch.fcgi?"
        + urllib.parse.urlencode(
            {"db": "sra", "id": ",".join(result["idlist"]), "retmode": "xml"}
        )
    )
    with urllib.request.urlopen(fetch_url, timeout=180) as response:
        root = ET.fromstring(response.read())
    resolved: dict[str, list[dict[str, str]]] = defaultdict(list)
    for package in root.findall(".//EXPERIMENT_PACKAGE"):
        experiment = package.find("EXPERIMENT")
        if experiment is None:
            continue
        gsm = experiment.attrib.get("alias", "")
        layout_node = experiment.find(".//LIBRARY_LAYOUT/*")
        layout = layout_node.tag if layout_node is not None else ""
        strategy = experiment.findtext(".//LIBRARY_STRATEGY") or ""
        title = experiment.findtext("TITLE") or ""
        for run in package.findall("./RUN_SET/RUN"):
            resolved[gsm].append(
                {
                    "run_accession": run.attrib["accession"],
                    "library_layout": layout,
                    "library_strategy": strategy,
                    "library_title": title,
                    "sra_archive_bytes": run.attrib.get("size", ""),
                    "total_bases": run.attrib.get("total_bases", ""),
                }
            )
    missing = sorted(set(gsms) - set(resolved))
    if missing:
        raise RuntimeError(f"GSMs without SRA runs: {', '.join(missing)}")
    return resolved


def resolve_ena(runs: list[str]) -> dict[str, dict[str, str]]:
    query = " OR ".join(f'run_accession="{run}"' for run in runs)
    fields = "run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout"
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(
        {
            "result": "read_run",
            "query": query,
            "fields": fields,
            "format": "tsv",
            "limit": max(100, len(runs)),
        }
    )
    with urllib.request.urlopen(url, timeout=180) as response:
        rows = list(
            csv.DictReader(response.read().decode("utf-8").splitlines(), delimiter="\t")
        )
    result = {row["run_accession"]: row for row in rows}
    missing = sorted(set(runs) - set(result))
    if missing:
        raise RuntimeError(f"Runs without ENA FASTQs: {', '.join(missing)}")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--selection",
        type=Path,
        default=Path("benchmarks/manifests/compact/lcmv_cd8_libraries.tsv"),
    )
    parser.add_argument(
        "--root", type=Path, default=Path("data/public/raw/lcmv_cd8_bulk")
    )
    args = parser.parse_args(argv)
    selection = read_tsv(args.selection)
    validate_selection(selection)
    sra = resolve_sra([row["gsm_accession"] for row in selection])
    run_ids = [item["run_accession"] for gsm in sra.values() for item in gsm]
    ena = resolve_ena(run_ids)

    resolved_rows: list[dict[str, str]] = []
    download_rows: list[dict[str, str]] = []
    for library in selection:
        gsm = library["gsm_accession"]
        for run_index, run in enumerate(sra[gsm], start=1):
            ena_row = ena[run["run_accession"]]
            urls = [
                value if "://" in value else f"https://{value}"
                for value in ena_row["fastq_ftp"].split(";")
            ]
            md5s = ena_row["fastq_md5"].split(";")
            sizes = ena_row["fastq_bytes"].split(";")
            if (
                len(urls) not in {1, 2}
                or len(urls) != len(md5s)
                or len(urls) != len(sizes)
            ):
                raise RuntimeError(
                    f"Incomplete ENA FASTQ metadata for {run['run_accession']}"
                )
            if ena_row["library_layout"] != run["library_layout"]:
                raise RuntimeError(
                    f"NCBI/ENA layout disagreement for {run['run_accession']}"
                )
            base = (
                args.root
                / "fastq"
                / library["assay"].lower()
                / gsm
                / run["run_accession"]
            )
            local_paths = [
                str((base / Path(urllib.parse.urlparse(url).path).name).resolve())
                for url in urls
            ]
            resolved = {
                **library,
                **run,
                "run_index_within_gsm": str(run_index),
                "fastq_1": local_paths[0],
                "fastq_2": local_paths[1] if len(local_paths) == 2 else "",
                "fastq_1_md5": md5s[0],
                "fastq_2_md5": md5s[1] if len(md5s) == 2 else "",
                "fastq_1_bytes": sizes[0],
                "fastq_2_bytes": sizes[1] if len(sizes) == 2 else "",
            }
            resolved_rows.append(resolved)
            for index, (url, md5, size, local_path) in enumerate(
                zip(urls, md5s, sizes, local_paths, strict=True), start=1
            ):
                download_rows.append(
                    {
                        "file_accession": f"{run['run_accession']}_{index}",
                        "gsm_accession": gsm,
                        "run_accession": run["run_accession"],
                        "assay": library["assay"],
                        "url": url,
                        "checksum": md5,
                        "expected_bytes": size,
                        "local_path": local_path,
                    }
                )

    if len(resolved_rows) != 60 or len(download_rows) != 105:
        raise RuntimeError(
            f"Expected 60 runs/105 FASTQs, found {len(resolved_rows)}/{len(download_rows)}"
        )
    metadata = args.root / "metadata"
    resolved_fields = list(selection[0]) + [
        "run_accession",
        "run_index_within_gsm",
        "library_layout",
        "library_strategy",
        "library_title",
        "sra_archive_bytes",
        "total_bases",
        "fastq_1",
        "fastq_2",
        "fastq_1_md5",
        "fastq_2_md5",
        "fastq_1_bytes",
        "fastq_2_bytes",
    ]
    write_tsv(metadata / "resolved_runs.tsv", resolved_rows, resolved_fields)
    write_tsv(
        metadata / "download_manifest.tsv",
        download_rows,
        [
            "file_accession",
            "gsm_accession",
            "run_accession",
            "assay",
            "url",
            "checksum",
            "expected_bytes",
            "local_path",
        ],
    )

    preprocess_fields = [
        "run_accession",
        "sample",
        "condition",
        "replicate",
        "fastq_1",
        "fastq_2",
        "fastq_1_md5",
        "fastq_2_md5",
        *list(selection[0]),
    ]
    for collection in ("primary", "supplemental"):
        rows = []
        for row in resolved_rows:
            if row["assay"] != "ATAC" or row["collection"] != collection:
                continue
            condition = (
                row["original_condition"]
                if row["author"] == "Beltra"
                else row["harmonized_condition"]
            )
            rows.append(
                {
                    **row,
                    "sample": row["run_accession"],
                    "condition": condition,
                }
            )
        write_tsv(metadata / f"{collection}_atac_samples.tsv", rows, preprocess_fields)

    pair_rows = []
    for pair_id in sorted(
        {
            row["condition_pair_id"]
            for row in selection
            if row["collection"] == "primary"
        }
    ):
        libraries = [row for row in selection if row["condition_pair_id"] == pair_id]
        runs = [row for row in resolved_rows if row["condition_pair_id"] == pair_id]
        atac = [row for row in libraries if row["assay"] == "ATAC"]
        rna = [row for row in libraries if row["assay"] == "RNA"]
        pair_rows.append(
            {
                "condition_pair_id": pair_id,
                "author": libraries[0]["author"],
                "infection": libraries[0]["infection"],
                "day": libraries[0]["day"],
                "tissue": libraries[0]["tissue"],
                "harmonized_condition": libraries[0]["harmonized_condition"],
                "broader_state": libraries[0]["broader_state"],
                "atac_gsm_count": str(len(atac)),
                "rna_gsm_count": str(len(rna)),
                "atac_gsms": ",".join(row["gsm_accession"] for row in atac),
                "rna_gsms": ",".join(row["gsm_accession"] for row in rna),
                "atac_srrs": ",".join(
                    row["run_accession"] for row in runs if row["assay"] == "ATAC"
                ),
                "rna_srrs": ",".join(
                    row["run_accession"] for row in runs if row["assay"] == "RNA"
                ),
                "pairing_unit": "condition",
            }
        )
    pair_fields = list(pair_rows[0])
    write_tsv(metadata / "paired_conditions.tsv", pair_rows, pair_fields)
    print(
        f"Resolved {len(selection)} GSM libraries to {len(resolved_rows)} runs and "
        f"{len(download_rows)} FASTQ files under {metadata}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
