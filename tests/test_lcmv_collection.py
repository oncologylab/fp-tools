import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / "benchmarks/manifests/compact/lcmv_cd8_libraries.tsv"
RESOLVER = ROOT / "benchmarks/scripts/resolve_lcmv_cd8_collection.py"


def _load_resolver():
    spec = importlib.util.spec_from_file_location(
        "resolve_lcmv_cd8_collection", RESOLVER
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_lcmv_selection_has_expected_assay_and_pair_counts():
    with SELECTION.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    _load_resolver().validate_selection(rows)
    assert len(rows) == 58
    assert (
        sum(row["assay"] == "ATAC" and row["collection"] == "primary" for row in rows)
        == 27
    )
    assert (
        sum(
            row["assay"] == "ATAC" and row["collection"] == "supplemental"
            for row in rows
        )
        == 6
    )
    assert sum(row["assay"] == "RNA" for row in rows) == 25


def test_every_primary_pair_contains_atac_and_rna_without_cross_study_pairing():
    with SELECTION.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    groups = {}
    for row in rows:
        if row["collection"] != "primary":
            continue
        group = groups.setdefault(
            row["condition_pair_id"], {"assays": set(), "authors": set()}
        )
        group["assays"].add(row["assay"])
        group["authors"].add(row["author"])
    assert len(groups) == 9
    assert all(group["assays"] == {"ATAC", "RNA"} for group in groups.values())
    assert all(len(group["authors"]) == 1 for group in groups.values())


def test_supplemental_atac_is_never_marked_as_rna_paired():
    with SELECTION.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    supplemental = [row for row in rows if row["collection"] == "supplemental"]
    assert supplemental
    assert all(
        row["rna_match_status"] == "no_exact_within_study_match" for row in supplemental
    )
    assert all(
        row["include_in_primary_paired_analysis"] == "false" for row in supplemental
    )
