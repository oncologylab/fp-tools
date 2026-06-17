import tempfile
import unittest
from pathlib import Path

import numpy as np
import pyBigWig
import pysam

from fp_tools.tools.fixed_site_differential import (
    build_cutcount_matrix,
    build_score_matrix,
    load_motif_site_reference,
    moderated_footprint_score,
)


class FixedSiteDifferentialTest(unittest.TestCase):
    def test_bound_union_reference_deduplicates_sites(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bed_dir = root / "TF1_M1" / "beds"
            bed_dir.mkdir(parents=True)
            content = "chr1\t10\t20\tTF1\t5\t+\nchr1\t10\t20\tTF1\t5\t+\n"
            (bed_dir / "TF1_M1_A_bound.bed").write_text(content, encoding="utf-8")
            (bed_dir / "TF1_M1_B_bound.bed").write_text("chr1\t30\t40\tTF1\t4\t-\n", encoding="utf-8")

            ref, _meta = load_motif_site_reference([str(root)], site_set="bound-union", window=20)

        self.assertEqual(sorted(ref), ["TF1_M1"])
        self.assertEqual(len(ref["TF1_M1"]), 2)
        self.assertEqual((ref["TF1_M1"][0].start, ref["TF1_M1"][0].end), (5, 25))

    def test_build_cutcount_matrix_uses_shifted_tn5_insertions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "sample.bam"
            header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 200, "SN": "chr1"}]}
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
                read = pysam.AlignedSegment()
                read.query_name = "fwd"
                read.query_sequence = "A" * 20
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 46
                read.mapping_quality = 60
                read.cigar = ((0, 20),)
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                bam.write(read)

                rev = pysam.AlignedSegment()
                rev.query_name = "rev"
                rev.query_sequence = "A" * 20
                rev.flag = 16
                rev.reference_id = 0
                rev.reference_start = 70
                rev.mapping_quality = 60
                rev.cigar = ((0, 20),)
                rev.query_qualities = pysam.qualitystring_to_array("I" * 20)
                bam.write(rev)
            pysam.index(str(bam_path))

            ref, _ = load_motif_site_reference_from_lines("TF1_M1", ["chr1\t45\t55\tTF1\t0\t+"])
            matrix = build_cutcount_matrix(ref, [str(bam_path)], ["A_rep1"], read_shift=(4, -5), min_mapq=30)

        # forward cut: 46 + 4 = 50; reverse cut: 90 - 1 - 5 = 84, outside the 20 bp site window.
        self.assertEqual(int(matrix.loc["TF1_M1", "A_rep1"]), 1)

    def test_build_score_matrix_averages_footprint_bigwig_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bw_path = tmp / "score.bw"
            with pyBigWig.open(str(bw_path), "w") as bw:
                bw.addHeader([("chr1", 100)])
                bw.addEntries(["chr1"], [0], ends=[100], values=[2.0])

            ref, _ = load_motif_site_reference_from_lines("TF1_M1", ["chr1\t40\t50\tTF1\t0\t+"])
            matrix = build_score_matrix(ref, [str(bw_path)], ["A_rep1"])

        self.assertAlmostEqual(float(matrix.loc["TF1_M1", "A_rep1"]), 2.0)

    def test_moderated_footprint_score_returns_qvalues(self):
        import pandas as pd

        matrix = pd.DataFrame(
            {
                "A_rep1": [4.0, 1.0],
                "A_rep2": [4.2, 1.1],
                "A_rep3": [3.9, 1.0],
                "B_rep1": [1.0, 1.1],
                "B_rep2": [1.2, 1.0],
                "B_rep3": [0.9, 0.9],
            },
            index=["TF1_M1", "TF2_M2"],
        )
        result = moderated_footprint_score(
            matrix,
            list(matrix.columns),
            {"A": ["A_rep1", "A_rep2", "A_rep3"], "B": ["B_rep1", "B_rep2", "B_rep3"]},
            ("A", "B"),
        )
        row = result.set_index("output_prefix").loc["TF1_M1"]
        self.assertGreater(row["footprint_score_delta"], 0)
        self.assertLess(row["padj"], 0.05)


def load_motif_site_reference_from_lines(prefix, lines):
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        bed_dir = root / prefix / "beds"
        bed_dir.mkdir(parents=True)
        (bed_dir / f"{prefix}_A_bound.bed").write_text("\n".join(lines) + "\n", encoding="utf-8")
        return load_motif_site_reference([str(root)], site_set="bound-union", window=20)


if __name__ == "__main__":
    unittest.main()
