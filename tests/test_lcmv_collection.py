import csv
import importlib.util
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / "benchmarks/manifests/compact/lcmv_cd8_libraries.tsv"
SELECTION_V2 = ROOT / "benchmarks/manifests/compact/lcmv_cd8_libraries_v2.tsv"
COMPARISONS_V2 = ROOT / "benchmarks/manifests/compact/lcmv_cd8_comparisons_v2.tsv"
RESOLVER = ROOT / "benchmarks/scripts/resolve_lcmv_cd8_collection.py"
DOWNSTREAM = ROOT / "benchmarks/scripts/build_lcmv_cd8_downstream.py"
SUMMARIZER = ROOT / "benchmarks/scripts/summarize_lcmv_rna.py"
VALIDATOR = ROOT / "benchmarks/scripts/validate_lcmv_outputs.py"


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


def test_lcmv_v2_selection_and_comparison_contract():
    with SELECTION_V2.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    _load_resolver().validate_selection(rows)
    assert len(rows) == 87
    assert sum(row["assay"] == "ATAC" for row in rows) == 50
    assert sum(row["assay"] == "RNA" for row in rows) == 37
    assert sum(row["collection"] == "primary" and row["assay"] == "ATAC" for row in rows) == 39
    assert sum(row["collection"] == "primary" and row["assay"] == "RNA" for row in rows) == 35
    assert len({row["condition_pair_id"] for row in rows if row["collection"] == "primary"}) == 14

    with COMPARISONS_V2.open(encoding="utf-8", newline="") as handle:
        comparisons = list(csv.DictReader(handle, delimiter="\t"))
    assert len(comparisons) == 28
    assert sum(row["analysis_tier"].startswith("primary_") for row in comparisons) == 18
    assert sum(row["analysis_tier"] == "primary_matched_context" for row in comparisons) == 8
    assert sum(row["analysis_tier"] == "primary_contextual_trajectory" for row in comparisons) == 10
    assert sum(row["analysis_tier"] == "supporting_atac_only" for row in comparisons) == 7
    assert sum(row["analysis_tier"] == "supporting_rna_only" for row in comparisons) == 3

    spec = importlib.util.spec_from_file_location("build_lcmv_cd8_downstream", DOWNSTREAM)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    module.validate_comparisons(comparisons, rows)


def test_lcmv_v2_exact_added_accessions():
    with SELECTION.open(encoding="utf-8", newline="") as handle:
        v1 = {row["gsm_accession"] for row in csv.DictReader(handle, delimiter="\t")}
    with SELECTION_V2.open(encoding="utf-8", newline="") as handle:
        v2 = {row["gsm_accession"] for row in csv.DictReader(handle, delimiter="\t")}
    assert v1 < v2
    assert v2 - v1 == {
        "GSM2889446", "GSM2889447", "GSM2889448", "GSM2889449",
        "GSM2863680", "GSM2863681", "GSM2863682", "GSM2863683",
        "GSM2356780", "GSM2356781", "GSM2356782", "GSM2356784",
        "GSM2356785", "GSM2356786", "GSM2356787", "GSM2356788",
        "GSM2356818", "GSM2356819", "GSM2356820", "GSM2356821",
        "GSM2356822", "GSM2356823", "GSM2865601", "GSM2865602",
        "GSM2865605", "GSM2865606", "GSM2863678", "GSM2863679",
        "GSM2356783",
    }


def test_lcmv_tx2gene_keeps_versioned_transcript_ids():
    spec = importlib.util.spec_from_file_location("summarize_lcmv_rna", SUMMARIZER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory() as tmp:
        gtf = Path(tmp) / "genes.gtf"
        output = Path(tmp) / "tx2gene.tsv"
        gtf.write_text(
            'chr1\ttest\ttranscript\t1\t10\t.\t+\t.\tgene_id "ENSMUSG1.2"; '
            'transcript_id "ENSMUST1.3"; gene_name "Gene1";\n',
            encoding="utf-8",
        )
        module.write_tx2gene(gtf, output)
        with output.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert rows == [
            {
                "transcript_id": "ENSMUST1.3",
                "gene_id": "ENSMUSG1",
                "gene_symbol": "Gene1",
            }
        ]


def test_lcmv_analysis_units_merge_only_technical_partitions():
    spec = importlib.util.spec_from_file_location("build_lcmv_cd8_downstream", DOWNSTREAM)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    with SELECTION_V2.open(encoding="utf-8", newline="") as handle:
        atac = [row for row in csv.DictReader(handle, delimiter="\t") if row["assay"] == "ATAC"]
    units = module.analysis_units(atac)
    assert len(units) == 48
    merged = [{row["gsm_accession"] for row in group} for group in units if len(group) > 1]
    assert merged == [{"GSM2356780", "GSM2356781"}, {"GSM2356795", "GSM2356796"}]
