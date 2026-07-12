import tempfile
import unittest
from pathlib import Path

from benchmarks.compare_atac_profiles import hash_overlap, peak_metrics


class CompareAtacProfilesTest(unittest.TestCase):
    def test_peak_overlap_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left = root / "left.bed"
            right = root / "right.bed"
            left.write_text("chr1\t0\t10\nchr1\t20\t30\n", encoding="utf-8")
            right.write_text("chr1\t5\t15\nchr1\t20\t25\n", encoding="utf-8")
            metrics = peak_metrics(left, right)
            self.assertEqual(metrics["peak_overlap_bp"], 10)
            self.assertEqual(metrics["left_peak_count"], 2)
            self.assertAlmostEqual(metrics["peak_bp_jaccard"], 0.4)

    def test_streaming_hash_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            left = root / "left"
            right = root / "right"
            left.write_text("a\nb\nc\n", encoding="ascii")
            right.write_text("b\nc\nd\n", encoding="ascii")
            metrics = hash_overlap(left, right)
            self.assertEqual(metrics["read_name_intersection"], 2)
            self.assertEqual(metrics["read_name_union"], 4)
            self.assertEqual(metrics["read_name_jaccard"], 0.5)


if __name__ == "__main__":
    unittest.main()
