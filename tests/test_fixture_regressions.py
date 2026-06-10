import pathlib
import unittest

import pyBigWig


ROOT = pathlib.Path(__file__).resolve().parents[1]


class FixtureRegressionTest(unittest.TestCase):
    EXPECTED_BIGWIG_WINDOWS = {
        "test_data/Bcell_corrected.bw": (153, 4.732392, 0.030931),
        "test_data/Tcell_corrected.bw": (127, 43.777098, 0.344702),
        "test_data/Bcell_footprints.bw": (305, 401.060294, 1.314952),
        "test_data/Tcell_footprints.bw": (399, 1171.803251, 2.936850),
    }

    def test_bigwig_fixture_headers_and_signal_summaries_are_stable(self):
        for rel_path, expected in self.EXPECTED_BIGWIG_WINDOWS.items():
            with self.subTest(bigwig=rel_path):
                path = ROOT / rel_path
                self.assertTrue(path.exists(), rel_path)
                bw = pyBigWig.open(str(path))
                try:
                    self.assertEqual(bw.chroms(), {"chr4": 190214555})
                    intervals = bw.intervals("chr4", 74000, 75000) or []
                    count = len(intervals)
                    total = round(sum(float(item[2]) for item in intervals), 6)
                    mean = round(total / count, 6)
                    self.assertEqual((count, total, mean), expected)
                finally:
                    bw.close()

    def test_core_bed_fixtures_have_expected_minimum_shape(self):
        expected_min_lines = {
            "test_data/merged_peaks.bed": 1,
            "test_data/merged_peaks_annotated.bed": 1,
            "test_data/IRF1_all.bed": 1,
            "test_data/BATF_all.bed": 1,
        }
        for rel_path, min_lines in expected_min_lines.items():
            with self.subTest(bed=rel_path):
                rows = [line for line in (ROOT / rel_path).read_text(encoding="utf-8").splitlines() if line.strip()]
                self.assertGreaterEqual(len(rows), min_lines)
                first_fields = rows[0].split("	")
                self.assertGreaterEqual(len(first_fields), 3)
                self.assertLess(int(first_fields[1]), int(first_fields[2]))


if __name__ == "__main__":
    unittest.main()
