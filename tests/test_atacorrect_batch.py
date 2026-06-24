import argparse
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from fp_tools.tools.atacorrect import _merge_peak_files, run_atacorrect


class AtacCorrectBatchTest(unittest.TestCase):
    def _args(self, tmp: Path, **overrides):
        values = {
            "bams": ["A.bam", "B.bam"],
            "genome": "genome.fa",
            "peaks": [str(tmp / "peaks.bed")],
            "sample_names": ["A", "B"],
            "prefix": None,
            "outdir": str(tmp / "out"),
            "track_off": [],
            "split_strands": False,
            "scale_corrected": "none",
            "scale_corrected_bigwigs": None,
            "scale_background": None,
            "scale_target": "median",
            "scale_chrom_sizes": None,
            "verbosity": 0,
        }
        values.update(overrides)
        return argparse.Namespace(**values)

    def test_merge_peak_files_writes_sorted_union(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            peaks_a = tmp / "a.bed"
            peaks_b = tmp / "b.bed"
            peaks_a.write_text("chr1\t10\t20\nchr1\t1\t5\n", encoding="utf-8")
            peaks_b.write_text("chr1\t18\t25\nchr2\t2\t4\n", encoding="utf-8")

            merged = Path(_merge_peak_files([peaks_a, peaks_b], tmp / "out"))

            self.assertEqual(merged.name, "merged_all.bed")
            self.assertEqual(
                merged.read_text(encoding="utf-8").splitlines(),
                ["chr1\t1\t5", "chr1\t10\t25", "chr2\t2\t4"],
            )

    def test_multi_bam_dispatch_uses_sample_subdirectories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            (tmp / "peaks.bed").write_text("chr1\t1\t5\n", encoding="utf-8")
            args = self._args(tmp)

            with mock.patch("fp_tools.tools.atacorrect._run_atacorrect_single") as run_one:
                run_atacorrect(args)

            self.assertEqual(run_one.call_count, 2)
            first = run_one.call_args_list[0].args[0]
            second = run_one.call_args_list[1].args[0]
            self.assertEqual(first.bam, "A.bam")
            self.assertEqual(first.prefix, "A")
            self.assertEqual(Path(first.outdir), tmp / "out" / "A")
            self.assertEqual(second.bam, "B.bam")
            self.assertEqual(second.prefix, "B")
            self.assertEqual(Path(second.outdir), tmp / "out" / "B")

    def test_multi_bam_shared_peak_does_not_warn(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            (tmp / "merged_peaks.bed").write_text("chr1\t1\t5\n", encoding="utf-8")
            args = self._args(tmp, peaks=[str(tmp / "merged_peaks.bed")], verbosity=1)

            stdout = io.StringIO()
            with mock.patch("fp_tools.tools.atacorrect._run_atacorrect_single"), contextlib.redirect_stdout(stdout):
                run_atacorrect(args)

            self.assertNotIn("WARNING", stdout.getvalue())

    def test_multi_bam_matching_peak_names_do_not_warn(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            peaks_a = tmp / "A_peaks.bed"
            peaks_b = tmp / "B_peaks.bed"
            peaks_a.write_text("chr1\t1\t5\n", encoding="utf-8")
            peaks_b.write_text("chr1\t10\t15\n", encoding="utf-8")
            args = self._args(tmp, peaks=[str(peaks_a), str(peaks_b)], verbosity=1)

            stdout = io.StringIO()
            with mock.patch("fp_tools.tools.atacorrect._run_atacorrect_single"), contextlib.redirect_stdout(stdout):
                run_atacorrect(args)

            self.assertNotIn("WARNING", stdout.getvalue())

    def test_multi_bam_mismatched_peak_count_warns_and_continues(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            peaks_a = tmp / "A_peaks.bed"
            peaks_b = tmp / "B_peaks.bed"
            peaks_c = tmp / "C_peaks.bed"
            for i, path in enumerate([peaks_a, peaks_b, peaks_c], start=1):
                path.write_text(f"chr1\t{i}\t{i + 5}\n", encoding="utf-8")
            args = self._args(tmp, peaks=[str(peaks_a), str(peaks_b), str(peaks_c)], verbosity=1)

            stdout = io.StringIO()
            with mock.patch("fp_tools.tools.atacorrect._run_atacorrect_single") as run_one, contextlib.redirect_stdout(stdout):
                run_atacorrect(args)

            self.assertEqual(run_one.call_count, 2)
            self.assertIn("counts differ", stdout.getvalue())

    def test_multi_bam_mismatched_peak_names_warn_and_continue(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            peaks_a = tmp / "B_peaks.bed"
            peaks_b = tmp / "A_peaks.bed"
            peaks_a.write_text("chr1\t1\t5\n", encoding="utf-8")
            peaks_b.write_text("chr1\t10\t15\n", encoding="utf-8")
            args = self._args(tmp, peaks=[str(peaks_a), str(peaks_b)], verbosity=1)

            stdout = io.StringIO()
            with mock.patch("fp_tools.tools.atacorrect._run_atacorrect_single") as run_one, contextlib.redirect_stdout(stdout):
                run_atacorrect(args)

            self.assertEqual(run_one.call_count, 2)
            self.assertIn("do not appear to match", stdout.getvalue())

    def test_multi_bam_rejects_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            args = self._args(tmp, prefix="bad")
            with self.assertRaises(SystemExit):
                run_atacorrect(args)


if __name__ == "__main__":
    unittest.main()
