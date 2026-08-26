#!/usr/bin/env python3
"""Build a manifest-pinned, replicate-complete ENCODE ATAC project.

The runner is intentionally resumable and storage-conscious. Released GRCh38
alignment BAMs are handled one at a time, verified before use, and downloaded
BAMs are removed only after their corrected bigWig has been validated and
recorded in a content-hash marker. Existing local source BAMs are never removed.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from itertools import combinations
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from fp_tools.utils import bigwig as pyBigWig

try:
    import pysam
except ImportError:  # The public-data runner is Linux-only, but its helpers are portable.
    pysam = None


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / "benchmarks/manifests/encode_atac_8cell_20260813.tsv"
DEFAULT_SPEC = ROOT / "benchmarks/manifests/encode_atac_8cell_20260813.spec.json"
DEFAULT_PROJECT = ROOT / "data/public/processed/encode_atac_8cell_20260813"
DEFAULT_SCRATCH = ROOT / "data/public/scratch/encode_atac_8cell_20260813"
DEFAULT_GENOME = ROOT / "data/public/raw/genome/hg38.fa"
DEFAULT_BLACKLIST = ROOT / "data/public/raw/encode_hct116_reviewer_revision/ENCFF356LFX.bed.gz"
DEFAULT_HCT116_SOURCE = ROOT / "data/public/raw/encode_hct116_reviewer_revision"
ENCODE = "https://www.encodeproject.org"
CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
EXPECTED_SELECTED_BIOSAMPLES = {
    "ENCBS110BUW", "ENCBS978MLY", "ENCBS831RKD", "ENCBS851BDM",
    "ENCBS193JTP", "ENCBS041EUO", "ENCBS597AKG", "ENCBS611WWO",
    "ENCBS327ZDR", "ENCBS708GCD", "ENCBS820EFY", "ENCBS986JFA",
    "ENCBS533AYC", "ENCBS388XJL", "ENCBS378LMT", "ENCBS254ODH",
}
EXPECTED_CONDITIONS = {"GM12878", "HCT116", "HepG2", "IMR-90", "K562", "MCF-7", "PC-3", "Panc1"}
MANIFEST_COLUMNS = {
    "condition", "sample", "biological_replicate", "selected_biosample",
    "biosample", "library",
    "experiment", "bam_accession", "bam_size", "bam_md5", "peak_accession",
    "peak_size", "peak_md5", "note",
}
EXPECTED_MOTIFS = 1019
EXPECTED_MERGED_REGIONS = 229_924
EXPECTED_MERGED_BP = 183_445_261
EXPECTED_MERGED_MD5 = "83a7b7f59a50b899b0471ae5bfc20f8f"
ENCODE_S3_OVERRIDES = {
    # Official public S3 objects used when the ENCODE portal intermittently
    # returns HTTP 500. Each download is still checked against ENCODE MD5/size.
    "ENCFF646NWY.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/baf5cad2-f9e6-4adc-84d8-5ec034b49977/ENCFF646NWY.bam",
    "ENCFF987XOV.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/57ce8c5b-7229-44b1-b5dc-e62d82a5b7d2/ENCFF987XOV.bam",
    "ENCFF607OSL.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/53686683-9fd3-4615-99ad-e46222fdc38b/ENCFF607OSL.bam",
    "ENCFF772EFK.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/de34a7f6-2b63-4e29-bc3f-37df9075ed36/ENCFF772EFK.bam",
    "ENCFF024FNF.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/7bbe68ba-e540-424a-ac28-ad34eb84f994/ENCFF024FNF.bam",
    "ENCFF516GDK.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/e448f14b-5f2f-4a09-b60d-12a1815d5cf1/ENCFF516GDK.bam",
    "ENCFF822BKT.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/847ea9c3-7600-4e31-bf92-e1f0e4b4fae8/ENCFF822BKT.bam",
    "ENCFF836WDC.bam": "https://encode-public.s3.amazonaws.com/2021/02/24/04721c3d-c404-4890-9d63-b40cfb2ef407/ENCFF836WDC.bam",
}


class ProjectSpec(NamedTuple):
    project_id: str
    conditions: dict[str, tuple[str, ...]]
    selected_biosamples: frozenset[str]
    motifs_per_sample: int
    peak_universe: dict[str, object] | None

    @property
    def samples(self) -> frozenset[str]:
        return frozenset(sample for values in self.conditions.values() for sample in values)

    @property
    def replicate_counts(self) -> dict[str, int]:
        return {condition: len(samples) for condition, samples in self.conditions.items()}


def load_project_spec(path: Path) -> ProjectSpec:
    payload = json.loads(path.read_text(encoding="utf-8"))
    conditions = {
        str(condition): tuple(str(sample) for sample in samples)
        for condition, samples in payload["conditions"].items()
    }
    if not conditions or any(len(samples) < 2 for samples in conditions.values()):
        raise ValueError("Every project condition must declare at least two biological replicates")
    all_samples = [sample for samples in conditions.values() for sample in samples]
    if len(all_samples) != len(set(all_samples)):
        raise ValueError("Project specification sample IDs must be unique")
    selected_biosamples = frozenset(str(value) for value in payload["selected_biosamples"])
    if len(selected_biosamples) != len(all_samples):
        raise ValueError("Project specification must declare one unique selected biosample per sample")
    motifs = int(payload.get("motifs_per_sample", EXPECTED_MOTIFS))
    if motifs <= 0:
        raise ValueError("motifs_per_sample must be positive")
    peak_universe = payload.get("peak_universe")
    if peak_universe is not None:
        required = {"regions", "covered_bp", "md5"}
        if set(peak_universe) != required:
            raise ValueError(f"peak_universe must contain exactly {sorted(required)}")
        peak_universe = {
            "regions": int(peak_universe["regions"]),
            "covered_bp": int(peak_universe["covered_bp"]),
            "md5": str(peak_universe["md5"]),
        }
    return ProjectSpec(
        project_id=str(payload["project_id"]),
        conditions=conditions,
        selected_biosamples=selected_biosamples,
        motifs_per_sample=motifs,
        peak_universe=peak_universe,
    )


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def digest(path: Path, algorithm: str = "sha256") -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def resolve(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else ROOT / path


def read_manifest(path: Path, spec: ProjectSpec | None = None) -> pd.DataFrame:
    spec = spec or load_project_spec(DEFAULT_SPEC)
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = MANIFEST_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
    for column in ("biological_replicate", "bam_size", "peak_size"):
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype("int64")
    expected_rows = len(spec.samples)
    if len(frame) != expected_rows:
        raise ValueError(f"Manifest must contain exactly {expected_rows} rows, found {len(frame)}")
    if frame["sample"].duplicated().any() or frame["bam_accession"].duplicated().any():
        raise ValueError("Every sample and BAM accession must be unique")
    if set(frame["sample"]) != set(spec.samples):
        raise ValueError("Manifest sample set does not exactly match the project specification")
    if set(frame["selected_biosample"]) != set(spec.selected_biosamples):
        raise ValueError("Manifest selected-biosample set does not exactly match the project specification")
    if frame["biosample"].duplicated().any() or frame["library"].duplicated().any():
        raise ValueError("Every assayed biosample and library accession must be unique")
    if not frame["biosample"].str.fullmatch(r"ENCBS[A-Z0-9]+").all():
        raise ValueError("Assayed biosample accessions must use the ENCBS prefix")
    if not frame["library"].str.fullmatch(r"ENCLB[A-Z0-9]+").all():
        raise ValueError("Library accessions must use the ENCLB prefix")
    if set(frame["condition"]) != set(spec.conditions):
        raise ValueError("Manifest condition set does not exactly match the project specification")
    counts = frame.groupby("condition")["sample"].count().to_dict()
    if counts != spec.replicate_counts:
        raise ValueError(f"Condition replicate counts do not match the project specification: {counts}")
    for condition, subset in frame.groupby("condition"):
        expected_replicates = set(range(1, spec.replicate_counts[condition] + 1))
        if set(subset["biological_replicate"]) != expected_replicates:
            raise ValueError(f"{condition} must contain biological replicates {sorted(expected_replicates)}")
        expected_samples = set(spec.conditions[condition])
        if set(subset["sample"]) != expected_samples:
            raise ValueError(f"{condition} samples do not match the project specification")
        if subset["experiment"].nunique() != 1 or subset["peak_accession"].nunique() != 1:
            raise ValueError(f"{condition} replicates must use one experiment and conservative-IDR peak set")
    forbidden = {"auroc", "auprc", "precision", "recall", "f1", "mcc"}
    if forbidden.intersection(column.lower() for column in frame.columns):
        raise ValueError("Input membership must not depend on performance metrics")
    return frame.sort_values(["condition", "biological_replicate"]).reset_index(drop=True)


def encode_metadata(accession: str) -> dict[str, object]:
    request = urllib.request.Request(
        f"{ENCODE}/files/{accession}/?format=json",
        headers={"Accept": "application/json", "User-Agent": "fp-tools-encode-project/1"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def encode_experiment_metadata(accession: str) -> dict[str, object]:
    request = urllib.request.Request(
        f"{ENCODE}/experiments/{accession}/?format=json",
        headers={"Accept": "application/json", "User-Agent": "fp-tools-encode-project/1"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def _accession(value: object) -> str:
    if isinstance(value, dict):
        return str(value.get("accession") or value.get("@id") or "").strip("/").split("/")[-1]
    return str(value or "").strip("/").split("/")[-1]


def audit_encode_replicates(frame: pd.DataFrame, output: Path) -> pd.DataFrame:
    """Verify assayed biosample and library accessions against ENCODE records."""

    rows: list[dict[str, object]] = []
    for experiment, subset in frame.groupby("experiment", sort=True):
        metadata = encode_experiment_metadata(str(experiment))
        replicate_map = {
            int(item["biological_replicate_number"]): item
            for item in metadata.get("replicates", [])
            if item.get("biological_replicate_number") is not None
            and int(item.get("technical_replicate_number", 1)) == 1
        }
        for expected in subset.itertuples(index=False):
            observed = replicate_map.get(int(expected.biological_replicate))
            if observed is None:
                raise ValueError(
                    f"{experiment} lacks biological replicate {expected.biological_replicate}"
                )
            library = observed.get("library") or {}
            biosample = library.get("biosample") or {}
            observed_library = _accession(library)
            observed_biosample = _accession(biosample)
            part_of = _accession(biosample.get("part_of") if isinstance(biosample, dict) else None)
            if observed_library != expected.library or observed_biosample != expected.biosample:
                raise ValueError(
                    f"ENCODE library/biosample metadata differs for {expected.sample}: "
                    f"{observed_library}/{observed_biosample}"
                )
            if expected.selected_biosample not in {observed_biosample, part_of}:
                raise ValueError(
                    f"Selected biosample {expected.selected_biosample} is not the assayed "
                    f"biosample or its source for {expected.sample}"
                )
            rows.append(
                {
                    "sample": expected.sample,
                    "experiment": experiment,
                    "biological_replicate": int(expected.biological_replicate),
                    "selected_biosample": expected.selected_biosample,
                    "biosample": observed_biosample,
                    "library": observed_library,
                    "source_relation": "assayed biosample"
                    if expected.selected_biosample == observed_biosample
                    else "source biosample",
                    "audited_at": utc_now(),
                }
            )
    audit = pd.DataFrame(rows).sort_values(["experiment", "biological_replicate"])
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output, sep="\t", index=False)
    return audit


def audit_encode_files(frame: pd.DataFrame, output: Path) -> pd.DataFrame:
    cached = None
    if output.is_file():
        try:
            cached = pd.read_csv(output, sep="\t", dtype={"file_size": "int64"})
        except Exception:
            cached = None
    rows: list[dict[str, object]] = []
    for accession, expected_type in [
        *((accession, "alignments") for accession in frame["bam_accession"]),
        *((accession, "conservative IDR thresholded peaks") for accession in frame["peak_accession"].drop_duplicates()),
    ]:
        expected = frame.loc[
            frame["bam_accession"].eq(accession) if expected_type == "alignments" else frame["peak_accession"].eq(accession)
        ].iloc[0]
        expected_size = int(expected["bam_size"] if expected_type == "alignments" else expected["peak_size"])
        expected_md5 = str(expected["bam_md5"] if expected_type == "alignments" else expected["peak_md5"])
        try:
            metadata = encode_metadata(accession)
        except Exception:
            cached_row = None if cached is None else cached.loc[cached["accession"].eq(accession)]
            if cached_row is None or len(cached_row) != 1:
                raise
            prior = cached_row.iloc[0].to_dict()
            if prior.get("status") != "released" or prior.get("assembly") != "GRCh38" or prior.get("output_type") != expected_type:
                raise ValueError(f"Cached ENCODE audit is invalid for {accession}: {prior}")
            if int(prior.get("file_size", -1)) != expected_size or prior.get("md5sum") != expected_md5:
                raise ValueError(f"Cached ENCODE size/MD5 differs from the manifest for {accession}: {prior}")
            if expected_type == "alignments":
                observed_bio = {
                    int(float(value))
                    for value in str(prior.get("biological_replicates", "")).split(",")
                    if value and value.lower() != "nan"
                }
                if observed_bio != {int(expected["biological_replicate"])}:
                    raise ValueError(f"Cached ENCODE biological replicate differs for {accession}")
            rows.append(prior)
            continue
        biological_replicates = {
            int(value) for value in metadata.get("biological_replicates", [])
        }
        technical_replicates = {
            str(value) for value in metadata.get("technical_replicates", [])
        }
        observed = {
            "accession": accession,
            "status": metadata.get("status"),
            "assembly": metadata.get("assembly"),
            "output_type": metadata.get("output_type"),
            "file_size": metadata.get("file_size"),
            "md5sum": metadata.get("md5sum"),
            "dataset": metadata.get("dataset"),
            "biological_replicates": ",".join(str(value) for value in sorted(biological_replicates)),
            "technical_replicates": ",".join(sorted(technical_replicates)),
            "audited_at": utc_now(),
        }
        if observed["status"] != "released" or observed["assembly"] != "GRCh38":
            raise ValueError(f"{accession} is not a released GRCh38 file: {observed}")
        if observed["output_type"] != expected_type:
            raise ValueError(f"Unexpected output type for {accession}: {observed['output_type']}")
        if int(observed["file_size"] or -1) != expected_size or observed["md5sum"] != expected_md5:
            raise ValueError(f"ENCODE size/MD5 changed for {accession}: {observed}")
        if expected_type == "alignments":
            expected_bio = int(expected["biological_replicate"])
            if biological_replicates != {expected_bio}:
                raise ValueError(f"{accession} does not represent biological replicate {expected_bio}")
            if not any(value.startswith(f"{expected_bio}_") for value in technical_replicates):
                raise ValueError(f"{accession} lacks the expected technical-replicate annotation")
        rows.append(observed)
    audit = pd.DataFrame(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output, sep="\t", index=False)
    return audit


def run_logged(command: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("\n$ " + " ".join(command) + "\n")
        handle.flush()
        result = subprocess.run(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise RuntimeError(f"Command exited {result.returncode}; see {log}")


def encode_s3_url(filename: str, date_added: str) -> str:
    """Resolve one ENCODE file in the public requester-free S3 bucket."""
    prefix = f"{date_added[:10].replace('-', '/')}/"
    continuation: str | None = None
    for _page in range(100):
        query = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if continuation:
            query["continuation-token"] = continuation
        with urllib.request.urlopen(
            "https://encode-public.s3.amazonaws.com/?" + urllib.parse.urlencode(query), timeout=60
        ) as response:
            root = ET.fromstring(response.read())
        namespace = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
        for node in root.findall("s3:Contents/s3:Key", namespace):
            if node.text and node.text.endswith("/" + filename):
                return "https://encode-public.s3.amazonaws.com/" + urllib.parse.quote(node.text, safe="/")
        token = root.find("s3:NextContinuationToken", namespace)
        if token is None or not token.text:
            break
        continuation = token.text
    raise ValueError(f"Could not locate {filename} under ENCODE public S3 prefix {prefix}")


def download(accession: str, size: int, md5: str, destination: Path, log: Path, curl: str) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(str(destination) + ".part")
    if destination.is_file():
        if destination.stat().st_size == size and digest(destination, "md5") == md5:
            return destination
        raise ValueError(f"Existing file does not match manifest: {destination}")
    portal_url = f"{ENCODE}/files/{accession}/@@download/{destination.name}"
    metadata = None
    try:
        metadata = encode_metadata(accession)
    except Exception:
        pass
    urls: list[str] = []
    if destination.name in ENCODE_S3_OVERRIDES:
        urls.append(ENCODE_S3_OVERRIDES[destination.name])
    if metadata and metadata.get("date_created"):
        resolved_s3 = encode_s3_url(destination.name, str(metadata["date_created"]))
        if resolved_s3 not in urls:
            urls.append(resolved_s3)
    # Prefer ENCODE's official public object store. The portal download route
    # remains a verified fallback and resolves to the same released object.
    urls.append(portal_url)
    error = None
    for url in urls:
        try:
            run_logged(
                [
                    curl, "--fail", "--location", "--connect-timeout", "20", "--speed-limit", "1024",
                    "--speed-time", "120", "--retry", "8", "--retry-all-errors", "--retry-delay", "15",
                    "--retry-max-time", "1800", "--continue-at", "-", "--output", str(partial), url,
                ],
                log,
            )
            error = None
            break
        except RuntimeError as exc:
            error = exc
    if error is not None:
        raise error
    if partial.stat().st_size != size or digest(partial, "md5") != md5:
        raise ValueError(f"Downloaded file failed size/MD5 validation: {accession}")
    partial.replace(destination)
    return destination


def prepare_peaks(frame: pd.DataFrame, project: Path, curl: str, spec: ProjectSpec) -> tuple[Path, Path]:
    input_dir = project / "input" / "encode_conservative_idr_peaks"
    log = project / "logs" / "download_peaks.log"
    intervals: list[tuple[str, int, int]] = []
    for row in frame.drop_duplicates("peak_accession").sort_values("peak_accession").itertuples(index=False):
        gz = download(row.peak_accession, int(row.peak_size), row.peak_md5, input_dir / f"{row.peak_accession}.bed.gz", log, curl)
        with gzip.open(gz, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3 or fields[0] not in CANONICAL_CHROMS:
                    continue
                start, end = int(fields[1]), int(fields[2])
                if start >= 0 and end > start:
                    intervals.append((fields[0], start, end))
    order = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24}
    intervals.sort(key=lambda item: (order[item[0]], item[1], item[2]))
    merged: list[list[object]] = []
    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= int(merged[-1][2]):
            merged[-1][2] = max(int(merged[-1][2]), end)
        else:
            merged.append([chrom, start, end])
    peaks_dir = project / "peaks"
    peaks_dir.mkdir(parents=True, exist_ok=True)
    merged_path = peaks_dir / "merged_peaks.bed"
    with merged_path.open("w", encoding="utf-8") as handle:
        for chrom, start, end in merged:
            handle.write(f"{chrom}\t{start}\t{end}\n")
    analysis_path = peaks_dir / "merged_peaks_filtered.bed"
    shutil.copyfile(merged_path, analysis_path)
    count = len(merged)
    span = sum(int(end) - int(start) for _chrom, start, end in merged)
    observed_md5 = digest(merged_path, "md5")
    if spec.peak_universe is not None:
        expected = spec.peak_universe
        if (count, span, observed_md5) != (expected["regions"], expected["covered_bp"], expected["md5"]):
            raise ValueError(
                "Global conservative-IDR universe differs from the project specification: "
                f"observed count={count}, bp={span}, md5={observed_md5}"
            )
    (peaks_dir / "peak_universe_qc.json").write_text(
        json.dumps({"regions": count, "covered_bp": span, "md5": observed_md5, "merge_touching": True}, indent=2) + "\n",
        encoding="utf-8",
    )
    return merged_path, analysis_path


def prepare_blacklist(source: Path, project: Path) -> Path:
    if not source.is_file() or digest(source, "md5") != "393688b4f06c9ce26165d47433dd8c37":
        raise ValueError(f"ENCODE GRCh38 blacklist is missing or checksum-mismatched: {source}")
    output = project / "input" / "ENCFF356LFX_hg38_blacklist.bed"
    output.parent.mkdir(parents=True, exist_ok=True)
    if not output.is_file():
        with gzip.open(source, "rt", encoding="utf-8") as src, output.open("w", encoding="utf-8") as dst:
            shutil.copyfileobj(src, dst)
    return output


def write_project_tables(frame: pd.DataFrame, project: Path) -> Path:
    project.mkdir(parents=True, exist_ok=True)
    table = project / "samples.tsv"
    frame[["sample", "condition"]].to_csv(table, sep="\t", index=False)
    frame.to_csv(project / "input_manifest.tsv", sep="\t", index=False)
    return table


def validate_bigwig(path: Path) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty bigWig: {path}")
    values: list[float] = []
    with pyBigWig.open(str(path)) as handle:
        chroms = handle.chroms()
        header = handle.header()
        if not chroms or int(header.get("nBasesCovered", 0)) <= 0:
            raise ValueError(f"bigWig has no covered bases: {path}")
        for chrom in list(chroms)[:8]:
            entries = handle.intervals(chrom, 0, min(chroms[chrom], 5_000_000)) or []
            values.extend(float(entry[2]) for entry in entries[:250])
    if not values or not np.isfinite(values).all():
        raise ValueError(f"bigWig has no finite sampled values: {path}")
    return {"size": path.stat().st_size, "sha256": digest(path), "covered_bases": int(header["nBasesCovered"])}


def marker_valid(path: Path, signature: dict[str, object]) -> bool:
    if not path.is_file():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    if payload.get("signature") != signature:
        return False
    for output in payload.get("outputs", []):
        item = Path(output["path"])
        if not item.is_file() or item.stat().st_size != int(output["size"]) or digest(item) != output["sha256"]:
            return False
    return True


def write_marker(path: Path, signature: dict[str, object], outputs: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "completed_at": utc_now(),
        "signature": signature,
        "outputs": [{"path": str(item.resolve()), "size": item.stat().st_size, "sha256": digest(item)} for item in outputs],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def local_bam(row: object, hct_source: Path) -> Path | None:
    candidate = hct_source / f"{row.bam_accession}.bam"
    return candidate if candidate.is_file() else None


def ensure_bam(row: object, scratch: Path, hct_source: Path, curl: str, log: Path) -> tuple[Path, bool]:
    if pysam is None:
        raise RuntimeError("The ENCODE project runner requires pysam on Linux")
    existing = local_bam(row, hct_source)
    if existing is not None:
        bam = existing
        temporary = False
    else:
        bam = download(row.bam_accession, int(row.bam_size), row.bam_md5, scratch / row.sample / f"{row.bam_accession}.bam", log, curl)
        temporary = True
    if bam.stat().st_size != int(row.bam_size) or digest(bam, "md5") != row.bam_md5:
        raise ValueError(f"BAM failed size/MD5 validation: {bam}")
    result = pysam.quickcheck(str(bam))
    if result:
        raise ValueError(f"BAM quickcheck failed for {bam}: {result}")
    bai = Path(str(bam) + ".bai")
    if not bai.is_file() or bai.stat().st_size == 0:
        pysam.index("-@", "8", str(bam))
    return bam, temporary


def correction_phase(args: argparse.Namespace, frame: pd.DataFrame, project: Path, peaks: Path, blacklist: Path) -> None:
    genome = resolve(args.genome)
    if not genome.is_file() or not Path(str(genome) + ".fai").is_file():
        raise ValueError(f"Genome or FASTA index is missing: {genome}")
    peak_hash = digest(peaks)
    blacklist_hash = digest(blacklist)
    genome_hash = digest(genome)
    for row in frame.itertuples(index=False):
        output_dir = project / "samples" / row.sample / "atac_correct"
        corrected = output_dir / f"{row.sample}_corrected.bw"
        marker = project / "state" / f"{row.sample}.atac_correct.json"
        signature = {
            "step": "atac-correct", "sample": row.sample, "bam_accession": row.bam_accession,
            "bam_md5": row.bam_md5, "genome_sha256": genome_hash, "peaks_sha256": peak_hash,
            "blacklist_sha256": blacklist_hash, "write_tracks": ["corrected"], "defaults": True,
        }
        if marker_valid(marker, signature):
            continue
        bam, temporary = ensure_bam(
            row, resolve(args.scratch), resolve(args.hct116_source), args.curl,
            project / "logs" / f"{row.sample}.download.log",
        )
        run_logged(
            [
                str(ROOT / ".venv/bin/atac-correct"), "--bams", str(bam), "--genome", str(genome),
                "--peaks", str(peaks), "--blacklist", str(blacklist), "--outdir", str(output_dir),
                "--prefix", row.sample, "--write-tracks", "corrected", "--scale-corrected", "none",
                "--cores", str(args.cores), "--sample-workers", "1", "--verbosity", "3",
            ],
            project / "logs" / f"{row.sample}.atac_correct.log",
        )
        validate_bigwig(corrected)
        write_marker(marker, signature, [corrected])
        if temporary:
            bam.unlink(missing_ok=True)
            Path(str(bam) + ".bai").unlink(missing_ok=True)
            Path(str(bam) + ".part").unlink(missing_ok=True)


def downstream_phase(args: argparse.Namespace, frame: pd.DataFrame, project: Path, sample_table: Path, peaks: Path, spec: ProjectSpec) -> None:
    corrected = [project / "samples" / row.sample / "atac_correct" / f"{row.sample}_corrected.bw" for row in frame.itertuples()]
    if any(not path.is_file() for path in corrected):
        raise ValueError(f"All {len(frame)} corrected tracks are required before global normalization")
    normalized = [project / "samples" / row.sample / "normalize" / f"{row.sample}_corrected_q95_scaled.bw" for row in frame.itertuples()]
    norm_marker = project / "state" / "normalize_q95.json"
    norm_signature = {"step": "normalize-bigwig", "inputs": [digest(path) for path in corrected], "background": digest(peaks), "stat": "q95", "target": "median"}
    if not marker_valid(norm_marker, norm_signature):
        try:
            validate_normalization_outputs(frame, project, corrected, normalized)
        except (FileNotFoundError, OSError, KeyError, ValueError):
            run_logged(
                [
                    str(ROOT / ".venv/bin/normalize-bigwig"), "--sample-table", str(sample_table), "--layout", "project",
                    "--outdir", str(project), "--background", str(peaks), "--method", "background-scale",
                    "--stat", "q95", "--target", "median", "--workers", str(min(4, args.cores)),
                ],
                project / "logs" / "normalize_q95.log",
            )
            validate_normalization_outputs(frame, project, corrected, normalized)
        write_marker(norm_marker, norm_signature, normalized)

    footprints = [project / "samples" / row.sample / "footprints" / f"{row.sample}_footprints.bw" for row in frame.itertuples()]
    fp_marker = project / "state" / "call_footprints.json"
    fp_signature = {"step": "call-footprints", "inputs": [digest(path) for path in normalized], "regions": digest(peaks), "kernel": "fast", "defaults": True}
    if not marker_valid(fp_marker, fp_signature):
        run_logged(
            [
                str(ROOT / ".venv/bin/call-footprints"), "--sample-table", str(sample_table), "--layout", "project",
                "--outdir", str(project), "--regions", str(peaks), "--score", "footprint",
                "--footprint-kernel", "fast", "--cores", str(args.cores), "--sample-workers", "2", "--verbosity", "3",
            ],
            project / "logs" / "call_footprints.log",
        )
        for path in footprints:
            validate_bigwig(path)
        write_marker(fp_marker, fp_signature, footprints)

    match_marker = project / "state" / "match_motifs.json"
    motif_db = ROOT / "src/fp_tools/resources/motifs/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    match_signature = {
        "step": "match-motifs", "inputs": [digest(path) for path in footprints], "peaks": digest(peaks),
        "motifs": digest(motif_db), "motif_db": "jaspar2026_vertebrates", "motif_outputs": "summary",
        "normalization": "none", "motif_pvalue": 0.0001, "bound_pvalue": 0.001,
    }
    summary_paths = [project / "samples" / row.sample / "match_motifs" / "motif_matches_results.txt" for row in frame.itertuples()]
    cache_paths = [project / "samples" / row.sample / "match_motifs" / "cache" / name for row in frame.itertuples() for name in ("motif_sites.tsv.gz", "background_scores.tsv.gz", "manifest.json", "thresholds.json")]
    if not marker_valid(match_marker, match_signature):
        run_logged(
            [
                str(ROOT / ".venv/bin/match-motifs"), "--sample-table", str(sample_table), "--layout", "project",
                "--outdir", str(project), "--peaks", str(peaks), "--genome", str(resolve(args.genome)),
                "--motif-db", "jaspar2026_vertebrates", "--normalization", "none", "--plot-aggregate", "off",
                "--motif-outputs", "summary", "--skip-excel", "--match-scan-mode", "shared",
                "--cores", str(args.cores), "--verbosity", "3",
            ],
            project / "logs" / "match_motifs.log",
        )
        validate_motif_outputs(frame, project, spec.motifs_per_sample)
        write_marker(match_marker, match_signature, summary_paths + cache_paths)
    write_replicate_qc(frame, project, spec)


def validate_normalization_outputs(
    frame: pd.DataFrame,
    project: Path,
    corrected: list[Path],
    normalized: list[Path],
) -> None:
    for path in normalized:
        validate_bigwig(path)
    manifest = pd.read_csv(project / "normalize_bigwig_manifest.tsv", sep="\t")
    qc = pd.read_csv(project / "normalize_bigwig_qc.tsv", sep="\t")
    expected_samples = frame["sample"].tolist()
    expected_inputs = [str(path) for path in corrected]
    expected_outputs = [str(path) for path in normalized]
    for table, label in ((manifest, "manifest"), (qc, "QC")):
        if table["sample"].tolist() != expected_samples:
            raise ValueError(f"Global normalization {label} sample order is incomplete or mismatched")
        if table["input_bigwig"].tolist() != expected_inputs or table["output_bigwig"].tolist() != expected_outputs:
            raise ValueError(f"Global normalization {label} paths do not match the current project")
    numeric = qc.select_dtypes(include="number").to_numpy()
    if not np.isfinite(numeric).all() or (qc["scaling_value"] <= 0).any() or (qc["scale_factor"] <= 0).any():
        raise ValueError("Global q95 normalization QC is non-finite or non-positive")
    if set(qc["scaling_stat"]) != {"q95"} or qc["target_scaling_value"].nunique() != 1:
        raise ValueError("Global normalization QC does not describe one shared q95 target")


def score_column(table: pd.DataFrame) -> str:
    columns = [column for column in table if column.endswith("_mean_score")]
    if len(columns) != 1:
        raise ValueError(f"Expected one mean-score column, found {columns}")
    return columns[0]


def validate_motif_outputs(frame: pd.DataFrame, project: Path, motifs_per_sample: int = EXPECTED_MOTIFS) -> None:
    for row in frame.itertuples(index=False):
        match = project / "samples" / row.sample / "match_motifs"
        table = pd.read_csv(match / "motif_matches_results.txt", sep="\t")
        if len(table) != motifs_per_sample or table["output_prefix"].nunique() != motifs_per_sample:
            raise ValueError(f"{row.sample} motif summary does not contain exactly {motifs_per_sample} motifs")
        values = pd.to_numeric(table[score_column(table)], errors="coerce")
        if values.notna().sum() == 0 or not np.isfinite(values.dropna()).all():
            raise ValueError(f"{row.sample} motif scores are missing or non-finite")
        for name in ("motif_sites.tsv.gz", "background_scores.tsv.gz", "manifest.json", "thresholds.json"):
            path = match / "cache" / name
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing compact match-motifs cache: {path}")
        unexpected = [path for path in match.iterdir() if path.is_dir() and path.name != "cache"]
        if unexpected:
            raise ValueError(f"Summary-only output contains per-motif directories: {unexpected[:3]}")


def write_replicate_qc(frame: pd.DataFrame, project: Path, spec: ProjectSpec) -> None:
    series: dict[str, pd.Series] = {}
    conditions: dict[str, str] = {}
    for row in frame.itertuples(index=False):
        table = pd.read_csv(project / "samples" / row.sample / "match_motifs" / "motif_matches_results.txt", sep="\t")
        values = pd.to_numeric(table[score_column(table)], errors="coerce")
        series[row.sample] = pd.Series(values.to_numpy(), index=table["output_prefix"].astype(str), name=row.sample)
        conditions[row.sample] = row.condition
    matrix = pd.concat(series.values(), axis=1, join="outer")
    correlation = matrix.corr(method="spearman", min_periods=100)
    qc_dir = project / "reports" / "replicate_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    correlation.to_csv(qc_dir / "motif_score_spearman_matrix.tsv", sep="\t", index_label="sample")
    rows = []
    for condition in sorted(spec.conditions):
        samples = [sample for sample, value in conditions.items() if value == condition]
        for sample_1, sample_2 in combinations(samples, 2):
            paired = matrix[[sample_1, sample_2]].dropna()
            rho = spearmanr(paired.iloc[:, 0], paired.iloc[:, 1]).statistic if len(paired) >= 2 else math.nan
            difference = paired.iloc[:, 0] - paired.iloc[:, 1]
            rows.append({
                "condition": condition, "replicate_1": sample_1, "replicate_2": sample_2,
                "motifs_compared": len(paired), "spearman_rho": rho,
                "median_absolute_difference": float(difference.abs().median()),
                "rmse": float(np.sqrt(np.mean(np.square(difference)))) if len(difference) else math.nan,
            })
    pd.DataFrame(rows).to_csv(qc_dir / "within_condition_replicate_qc.tsv", sep="\t", index=False)


def verify_project(args: argparse.Namespace, frame: pd.DataFrame, project: Path, spec: ProjectSpec) -> dict[str, object]:
    observed_samples = {path.name for path in (project / "samples").iterdir() if path.is_dir()}
    expected_samples = set(frame["sample"])
    if observed_samples != expected_samples:
        raise ValueError(f"Project sample set mismatch: missing={expected_samples-observed_samples}, extra={observed_samples-expected_samples}")
    outputs: dict[str, dict[str, object]] = {}
    for row in frame.itertuples(index=False):
        sample_dir = project / "samples" / row.sample
        paths = {
            "corrected": sample_dir / "atac_correct" / f"{row.sample}_corrected.bw",
            "normalized": sample_dir / "normalize" / f"{row.sample}_corrected_q95_scaled.bw",
            "footprint": sample_dir / "footprints" / f"{row.sample}_footprints.bw",
        }
        outputs[row.sample] = {label: validate_bigwig(path) for label, path in paths.items()}
    validate_motif_outputs(frame, project, spec.motifs_per_sample)
    qc = pd.read_csv(project / "reports" / "replicate_qc" / "within_condition_replicate_qc.tsv", sep="\t")
    expected_pairs = sum(count * (count - 1) // 2 for count in spec.replicate_counts.values())
    if len(qc) != expected_pairs or set(qc["condition"]) != set(spec.conditions):
        raise ValueError("Replicate QC does not contain every within-condition replicate pair")
    scratch = resolve(args.scratch)
    retained_downloads = (
        list(scratch.glob("**/*.bam"))
        + list(scratch.glob("**/*.bam.bai"))
        + list(scratch.glob("**/*.bam.part"))
    )
    if retained_downloads:
        raise ValueError(f"Downloaded BAM artifacts remain after completion: {retained_downloads[:3]}")
    local_source = resolve(args.hct116_source)
    for row in frame.itertuples(index=False):
        candidate = local_source / f"{row.bam_accession}.bam"
        if candidate.is_file() and digest(candidate, "md5") != row.bam_md5:
            raise ValueError(f"Protected local source BAM changed: {row.bam_accession}")
    peak_universe = json.loads((project / "peaks" / "peak_universe_qc.json").read_text(encoding="utf-8"))
    summary = {
        "verified_at": utc_now(), "project": str(project), "project_id": spec.project_id,
        "samples": len(frame), "conditions": len(spec.conditions),
        "replicates_per_condition": spec.replicate_counts, "motifs_per_sample": spec.motifs_per_sample,
        "peak_universe": peak_universe,
        "outputs": outputs,
    }
    (project / "reports" / "verification.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def preflight(args: argparse.Namespace, frame: pd.DataFrame, project: Path, spec: ProjectSpec) -> tuple[Path, Path, Path]:
    free_gb = shutil.disk_usage(project.parent if project.parent.exists() else ROOT).free / (1024 ** 3)
    if free_gb < args.minimum_free_gb:
        raise ValueError(f"Only {free_gb:.1f} GiB free; at least {args.minimum_free_gb:.1f} GiB is required")
    for executable in ("curl", str(ROOT / ".venv/bin/atac-correct"), str(ROOT / ".venv/bin/normalize-bigwig"), str(ROOT / ".venv/bin/call-footprints"), str(ROOT / ".venv/bin/match-motifs")):
        if not (Path(executable).is_file() or shutil.which(executable)):
            raise ValueError(f"Required executable is unavailable: {executable}")
    project.mkdir(parents=True, exist_ok=True)
    audit_encode_files(frame, project / "input" / "encode_file_audit.tsv")
    audit_encode_replicates(frame, project / "input" / "encode_replicate_audit.tsv")
    merged, analysis = prepare_peaks(frame, project, args.curl, spec)
    blacklist = prepare_blacklist(resolve(args.blacklist), project)
    sample_table = write_project_tables(frame, project)
    return sample_table, analysis, blacklist


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("preflight", "run", "verify"))
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--spec", default=str(DEFAULT_SPEC))
    parser.add_argument("--project", default=str(DEFAULT_PROJECT))
    parser.add_argument("--scratch", default=str(DEFAULT_SCRATCH))
    parser.add_argument("--genome", default=str(DEFAULT_GENOME))
    parser.add_argument("--blacklist", default=str(DEFAULT_BLACKLIST))
    parser.add_argument("--hct116-source", default=str(DEFAULT_HCT116_SOURCE))
    parser.add_argument("--curl", default=shutil.which("curl") or "curl")
    parser.add_argument("--cores", type=int, default=min(28, os.cpu_count() or 1))
    parser.add_argument("--minimum-free-gb", type=float, default=150.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    spec = load_project_spec(resolve(args.spec))
    frame = read_manifest(resolve(args.manifest), spec)
    project = resolve(args.project)
    if args.mode == "verify":
        summary = verify_project(args, frame, project, spec)
        print(json.dumps({key: value for key, value in summary.items() if key != "outputs"}, indent=2))
        return 0
    sample_table, peaks, blacklist = preflight(args, frame, project, spec)
    if args.mode == "preflight":
        print(f"Preflight passed for {len(frame)} biological replicates; project={project}; peaks={peaks}")
        return 0
    correction_phase(args, frame, project, peaks, blacklist)
    downstream_phase(args, frame, project, sample_table, peaks, spec)
    summary = verify_project(args, frame, project, spec)
    print(json.dumps({key: value for key, value in summary.items() if key != "outputs"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
