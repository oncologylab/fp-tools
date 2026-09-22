from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from build_strand_functional_profiles import extract_strand_cut_profiles


def bam_fixture(tmp_path, records):
    pysam = pytest.importorskip("pysam")
    bam = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                                    "SQ": [{"SN": "chr1", "LN": 1000}]}) as output:
        for index, (start, cigar, reverse) in enumerate(records):
            read = pysam.AlignedSegment()
            read.query_name = f"read{index}"
            read.flag = 16 if reverse else 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigartuples = cigar
            read.query_sequence = "A" * sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
            output.write(read)
    pysam.index(str(bam))
    return bam


def sites(centers):
    return pd.DataFrame({"TFBS_chr": ["chr1"] * len(centers),
                         "TFBS_start": centers, "TFBS_end": centers,
                         "TFBS_strand": ["+"] * len(centers)})


def extract(bam, centers):
    return extract_strand_cut_profiles(sites(centers), bam, flank=10,
                                       read_shift=(4, -5), minimum_mapq=30)


def test_long_alignment_crossing_disjoint_queries_is_counted_once(tmp_path):
    bam = bam_fixture(tmp_path, [(100, [(0, 150)], False)])
    plus, minus, valid = extract(bam, [104, 224])
    assert valid.all()
    assert plus[0, 10] == 1  # One record, even though its alignment overlaps both windows.
    assert plus.sum() == 1 and minus.sum() == 0


@pytest.mark.parametrize("start,cigar,reverse,center", [
    (140, [(4, 40), (0, 30)], False, 104),
    (30, [(0, 30), (4, 40)], True, 95),
])
def test_soft_clipped_cut_outside_alignment_is_not_missed(tmp_path, start, cigar, reverse, center):
    bam = bam_fixture(tmp_path, [(start, cigar, reverse)])
    plus, minus, valid = extract(bam, [center])
    assert valid.all()
    assert (minus if reverse else plus)[0, 10] == 1
    assert plus.sum() + minus.sum() == 1


def test_overlapping_sites_each_receive_one_cut_and_order_is_preserved(tmp_path):
    bam = bam_fixture(tmp_path, [(100, [(0, 150)], False)])
    plus, minus, valid = extract(bam, [109, 104, 99])
    assert valid.all() and plus.sum() == 3 and minus.sum() == 0
    assert [plus[i].argmax() for i in range(3)] == [5, 10, 15]


def test_corrected_extractor_matches_independent_full_width_oracle(tmp_path):
    import numpy as np
    from audit_strand_extraction import oracle

    bam = bam_fixture(tmp_path, [(100, [(0, 150)], False),
                                (300, [(4, 200), (0, 30)], False)])
    frame = sites([104, 324])
    actual = extract_strand_cut_profiles(frame, bam, flank=100, read_shift=(4, -5))
    expected = oracle(frame, bam)
    assert all(np.array_equal(a, b) for a, b in zip(actual, expected))
    assert actual[0][0, 100] == 2 and actual[0].sum() == 2


def test_nonconsecutive_dataframe_index_keeps_positional_order(tmp_path):
    bam = bam_fixture(tmp_path, [(100, [(0, 150)], False)])
    frame = sites([104, 224])
    frame.index = [7, 2]
    plus, minus, valid = extract_strand_cut_profiles(frame, bam, flank=10, read_shift=(4, -5))
    assert valid.all() and plus[0, 10] == 1 and plus[1].sum() == 0
