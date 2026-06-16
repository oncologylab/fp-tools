import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pyBigWig

from fp_tools.tools.normalize_bigwig import normalize_bigwigs
from fp_tools.tools.normalize_bigwig import _stat


class NormalizeBigwigTest(unittest.TestCase):
    def _write_bigwig(self, path: Path, values):
        bw = pyBigWig.open(str(path), "w")
        try:
            bw.addHeader([("chr1", len(values))])
            starts = list(range(len(values)))
            ends = [start + 1 for start in starts]
            bw.addEntries(["chr1"] * len(values), starts, ends=ends, values=[float(v) for v in values])
        finally:
            bw.close()

    def _read_values(self, path: Path, start=0, end=10):
        with pyBigWig.open(str(path)) as bw:
            return np.asarray(bw.values("chr1", start, end, numpy=True), dtype=float)

    def test_background_scale_writes_scaled_bigwigs_and_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bw1 = tmp / "sample1_corrected.bw"
            bw2 = tmp / "sample2_corrected.bw"
            self._write_bigwig(bw1, [1, 2, 3, 4, 5, 10, 10, 10, 10, 10])
            self._write_bigwig(bw2, [2, 4, 6, 8, 10, 20, 20, 20, 20, 20])
            background = tmp / "background.bed"
            background.write_text("".join(f"chr1\t{i}\t{i + 1}\n" for i in range(5)), encoding="utf-8")

            rows = normalize_bigwigs(
                [str(bw1), str(bw2)],
                background,
                tmp / "norm",
                method="background-scale",
                stat="q90",
                target="median",
            )

            self.assertEqual(len(rows), 2)
            self.assertAlmostEqual(rows[0].scale_factor, 1.5)
            self.assertAlmostEqual(rows[1].scale_factor, 0.75)
            np.testing.assert_allclose(self._read_values(Path(rows[0].output_bigwig), 0, 5), [1.5, 3.0, 4.5, 6.0, 7.5])
            np.testing.assert_allclose(self._read_values(Path(rows[1].output_bigwig), 0, 5), [1.5, 3.0, 4.5, 6.0, 7.5])

            qc = tmp / "norm" / "normalize_bigwig_qc.tsv"
            manifest = tmp / "norm" / "normalize_bigwig_manifest.tsv"
            self.assertTrue(qc.exists())
            self.assertTrue(manifest.exists())
            with qc.open(encoding="utf-8") as handle:
                data = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(data[0]["sample"], "sample1_corrected")
            self.assertEqual(data[0]["background_q90"], "4.6")
            self.assertEqual(data[0]["scale_factor"], "1.5")


    def test_background_scale_accepts_high_tail_quantiles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bw1 = tmp / "sample1_corrected.bw"
            bw2 = tmp / "sample2_corrected.bw"
            self._write_bigwig(bw1, [1, 2, 3, 4, 5, 10, 10, 10, 10, 10])
            self._write_bigwig(bw2, [2, 4, 6, 8, 10, 20, 20, 20, 20, 20])
            background = tmp / "background.bed"
            background.write_text("".join(f"chr1\t{i}\t{i + 1}\n" for i in range(5)), encoding="utf-8")

            rows = normalize_bigwigs(
                [str(bw1), str(bw2)],
                background,
                tmp / "norm",
                method="background-scale",
                stat="q95",
                target="median",
            )

            self.assertAlmostEqual(rows[0].background_q95, 4.8)
            self.assertEqual(rows[0].scaling_stat, "q95")
            self.assertAlmostEqual(rows[0].scaling_value, 4.8)
            self.assertTrue(rows[0].output_bigwig.endswith(".background_scale_q95.bw"))
            with (tmp / "norm" / "normalize_bigwig_qc.tsv").open(encoding="utf-8") as handle:
                data = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(data[0]["background_q95"], "4.8")
            self.assertEqual(data[0]["scaling_stat"], "q95")

    def test_quantile_stat_format_supports_decimal_filename(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bw = tmp / "sample_corrected.bw"
            self._write_bigwig(bw, [1, 2, 3, 4, 5])
            background = tmp / "background.bed"
            background.write_text("".join(f"chr1\t{i}\t{i + 1}\n" for i in range(5)), encoding="utf-8")

            rows = normalize_bigwigs([str(bw)], background, tmp / "norm", stat="q97.5")

            self.assertAlmostEqual(rows[0].scaling_value, np.quantile([1, 2, 3, 4, 5], 0.975))
            self.assertTrue(rows[0].output_bigwig.endswith(".background_scale_q97_5.bw"))

    def test_invalid_quantile_stats_fail_clearly(self):
        values = np.asarray([1, 2, 3], dtype=float)
        for stat in ("q0", "q100", "qabc"):
            with self.subTest(stat=stat):
                with self.assertRaisesRegex(ValueError, "Unsupported statistic|Quantile statistic"):
                    _stat(values, stat)

    def test_background_zscore_centers_by_background_median_and_mad(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bw = tmp / "footprints.bw"
            self._write_bigwig(bw, [1, 2, 3, 4, 100])
            background = tmp / "background.bed"
            background.write_text("".join(f"chr1\t{i}\t{i + 1}\n" for i in range(4)), encoding="utf-8")

            rows = normalize_bigwigs([str(bw)], background, tmp / "norm", method="background-zscore")

            self.assertAlmostEqual(rows[0].background_median, 2.5)
            self.assertAlmostEqual(rows[0].background_mad, 1.0)
            np.testing.assert_allclose(self._read_values(Path(rows[0].output_bigwig), 0, 4), [-1.5, -0.5, 0.5, 1.5])

    def test_missing_background_values_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bw = tmp / "sample.bw"
            self._write_bigwig(bw, [1, 2, 3])
            background = tmp / "background.bed"
            background.write_text("chr2\t0\t3\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "No background signal values"):
                normalize_bigwigs([str(bw)], background, tmp / "norm")


if __name__ == "__main__":
    unittest.main()
