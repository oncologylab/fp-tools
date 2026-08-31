from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import pandas as pd
import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from export_label_free_profile_sites import export_sites, parse_filter  # noqa: E402


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def test_export_sites_filters_without_leaking_labels(tmp_path: Path) -> None:
    sites_path = tmp_path / "source.tsv"
    pd.DataFrame(
        {
            "TFBS_chr": ["chr16", "chr17", "chr1"],
            "TFBS_start": [30, 10, 20],
            "TFBS_end": [31, 11, 21],
            "chromosome_split": ["validation", "validation", "train"],
            "chip_label": [1, 0, 1],
            "label_reason": ["positive", "negative", "positive"],
            "chip_accession": ["A", "A", "A"],
            "motif_score": [3.0, 2.0, 1.0],
        }
    ).to_csv(sites_path, sep="\t", index=False)
    artifact = tmp_path / "artifact.json"
    artifact.write_text(
        json.dumps({"sites": str(sites_path), "sites_sha256": _sha256(sites_path)}),
        encoding="utf-8",
    )
    output = tmp_path / "validation.tsv.gz"
    bed_output = tmp_path / "validation.bed"
    output_path, manifest_path = export_sites(
        artifact,
        output,
        [("chromosome_split", "validation")],
        bed_output=bed_output,
        flank=5,
    )
    exported = pd.read_csv(output_path, sep="\t")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert exported["TFBS_chr"].tolist() == ["chr16", "chr17"]
    assert not any("label" in column.lower() or "chip" in column.lower() for column in exported)
    assert manifest["labels_used"] is False
    assert manifest["rows"] == 2
    assert set(manifest["dropped_columns"]) == {
        "chip_label",
        "label_reason",
        "chip_accession",
    }
    assert manifest["output_sha256"] == _sha256(output)
    bed = pd.read_csv(
        bed_output,
        sep="\t",
        header=None,
        names=["chromosome", "start", "end"],
    )
    assert bed.to_dict("records") == [
        {"chromosome": "chr16", "start": 25, "end": 36},
        {"chromosome": "chr17", "start": 5, "end": 16},
    ]
    assert manifest["bed_output"]["sha256"] == _sha256(bed_output)


def test_export_sites_rejects_label_filter_and_bad_checksum(tmp_path: Path) -> None:
    with pytest.raises(argparse.ArgumentTypeError):
        parse_filter("chip_label=1")
    sites_path = tmp_path / "source.tsv"
    pd.DataFrame({"TFBS_chr": ["chr1"]}).to_csv(sites_path, sep="\t", index=False)
    artifact = tmp_path / "artifact.json"
    artifact.write_text(
        json.dumps({"sites": str(sites_path), "sites_sha256": "bad"}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="checksum"):
        export_sites(artifact, tmp_path / "out.tsv.gz", [])
