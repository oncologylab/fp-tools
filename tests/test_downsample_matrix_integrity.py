from pathlib import Path
import json
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import downsample_bam_depth_matrix as matrix
from downsample_bam_fragments import fragment_uniform


def source_bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "source.bam"
    with pysam.AlignmentFile(path, "wb", header={"HD": {"SO": "coordinate"},
                                                "SQ": [{"SN": "chr1", "LN": 10000}]}) as handle:
        for pair in range(40):
            for mate, flag in enumerate((99, 147)):
                read = pysam.AlignedSegment()
                read.query_name = f"pair_{pair}"
                read.query_sequence = "A" * 20
                read.flag = flag
                read.reference_id = read.next_reference_id = 0
                read.reference_start = pair * 50 + mate * 25
                read.next_reference_start = pair * 50 + (1-mate) * 25
                read.mapping_quality = 40
                read.cigarstring = "20M"
                handle.write(read)
    pysam.index(str(path))
    return path


@pytest.mark.parametrize("corrupt", ["bam", "index", "metadata", "source"])
def test_resume_rejects_changed_or_unproven_artifacts(tmp_path, corrupt):
    source = source_bam(tmp_path)
    targets = matrix.build_targets(tmp_path / "outputs", "sample", [20], [2026])
    matrix.downsample_matrix(source, targets, available_fragments=40)
    target = targets[0]
    path = {"bam": target.output_bam, "index": Path(str(target.output_bam) + ".bai"),
            "metadata": target.metadata_path, "source": source}[corrupt]
    if corrupt == "metadata":
        document = json.loads(path.read_text())
        document["available_fragments"] = 999
        path.write_text(json.dumps(document))
    else:
        with path.open("ab") as handle:
            handle.write(b"changed")
    with pytest.raises(ValueError, match="checksum|metadata|source"):
        matrix.downsample_matrix(source, targets, available_fragments=40)


def test_partial_output_is_preserved_and_never_silently_skipped(tmp_path):
    source = source_bam(tmp_path)
    targets = matrix.build_targets(tmp_path / "outputs", "sample", [20], [2026])
    target = targets[0]
    target.output_bam.parent.mkdir(parents=True)
    target.output_bam.write_bytes(b"interrupted")
    with pytest.raises(ValueError, match="incomplete"):
        matrix.downsample_matrix(source, targets, available_fragments=40)
    assert target.output_bam.read_bytes() == b"interrupted"


def test_valid_resume_keeps_exact_pairs_nested_depths_and_file_mtimes(tmp_path):
    pysam = pytest.importorskip("pysam")
    source = source_bam(tmp_path)
    targets = matrix.build_targets(tmp_path / "outputs", "sample", [10, 20, 40], [2026, 2027])
    first = matrix.downsample_matrix(source, targets, available_fragments=40)
    mtimes = {target.output_bam: target.output_bam.stat().st_mtime_ns for target in targets}
    second = matrix.downsample_matrix(source, targets, available_fragments=40)
    assert all(row["state"] == "skipped_existing" for row in second)
    assert mtimes == {target.output_bam: target.output_bam.stat().st_mtime_ns for target in targets}
    for target in targets:
        with pysam.AlignmentFile(target.output_bam, "rb") as bam:
            records = [(read.query_name, read.is_read1) for read in bam.fetch(until_eof=True)]
        expected = {f"pair_{i}" for i in range(40)
                    if fragment_uniform(f"pair_{i}", target.seed) < target.target_fragments / 40}
        assert {name for name, _ in records} == expected
        assert len(records) == 2 * len(expected)
        assert len({name for name, read1 in records if read1}) == len(expected)
        before = next(row for row in first if row["output_bam"] == str(target.output_bam))
        after = next(row for row in second if row["output_bam"] == str(target.output_bam))
        assert before["output_bam_sha256"] == after["output_bam_sha256"]


def test_frozen_source_checksum_is_required_before_any_output(tmp_path):
    source = source_bam(tmp_path)
    targets = matrix.build_targets(tmp_path / "outputs", "sample", [20], [2026])
    with pytest.raises(ValueError, match="source BAM checksum"):
        matrix.downsample_matrix(source, targets, available_fragments=40, expected_source_sha256="0" * 64)
    assert not (tmp_path / "outputs").exists()
