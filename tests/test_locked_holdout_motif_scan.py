import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "scan_locked_holdout_motif_sites.py"
spec = importlib.util.spec_from_file_location("scan_locked_holdout_motif_sites", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_peak_reader_merges_overlaps_and_filters_chromosomes(tmp_path):
    path = tmp_path / "peaks.bed"
    path.write_text(
        "chr1\t20\t30\nchr1\t10\t22\nchr1\t40\t50\nchrM\t1\t9\n",
        encoding="utf-8",
    )
    assert module.read_and_merge_peaks(path, {"chr1"}) == [
        ("chr1", 10, 30),
        ("chr1", 40, 50),
    ]


def test_scan_output_schema_contains_no_occupancy_labels():
    assert "chip_label" not in module.OUTPUT_COLUMNS
    assert not any("label" in column.lower() for column in module.OUTPUT_COLUMNS)
    assert {
        "tf",
        "motif_id",
        "motif_family",
        "TFBS_chr",
        "TFBS_start",
        "TFBS_end",
        "TFBS_strand",
        "chromosome_split",
    }.issubset(module.OUTPUT_COLUMNS)


def test_locked_holdout_motifs_exist_in_pinned_database():
    import json
    from fp_tools.utils.motifs import MotifList

    study = json.loads(
        (ROOT / "benchmarks" / "manifests" / "footprint_functional_v1.spec.json").read_text()
    )
    wanted = {
        task["motif_id"] for task in study["tasks"] if task["split"] == "locked_holdout"
    }
    motifs = MotifList().from_file(
        str(
            ROOT
            / "src/fp_tools/resources/motifs/"
            "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
        )
    )
    assert wanted.issubset({motif.id for motif in motifs})
