import argparse
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from fp_tools.tools.diff_footprints import _run_project_comparison_table
from fp_tools.utils.project_layout import (
    analysis_peaks_path,
    comparison_dir,
    read_comparison_table,
    read_sample_table,
    samples_for_condition,
)
from fp_tools.tools.normalize_bigwig import project_scaled_output_path
from fp_tools.tools.score_bigwig import _scorebigwig_batch_items


class ProjectLayoutPathTest(unittest.TestCase):
    def test_normalize_project_scaled_output_path(self):
        path = project_scaled_output_path("project/samples", "cy249")
        self.assertEqual(path, Path("project/samples/cy249/normalize/cy249_corrected_q95_scaled.bw"))

    def test_call_footprints_project_layout_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "samples"
            args = argparse.Namespace(
                signals=["A_corrected_q95_scaled.bw", "B_corrected_q95_scaled.bw"],
                signal=None,
                output=None,
                outputs=None,
                output_bed=None,
                output_beds=None,
                output_bed_dir=None,
                output_multiscale_npz=None,
                output_multiscale_npzs=None,
                outdir=None,
                sample_names=["A", "B"],
                sample_output_root=str(root),
            )

            items = _scorebigwig_batch_items(args)

            self.assertEqual(items[0][1], str(root / "A" / "footprints" / "A_footprints.bw"))
            self.assertEqual(items[1][1], str(root / "B" / "footprints" / "B_footprints.bw"))
            self.assertEqual(items[0][2], str(root / "A" / "footprints" / "A_candidate_footprints.bed"))
            self.assertEqual(items[1][2], str(root / "B" / "footprints" / "B_candidate_footprints.bed"))

    def test_generic_sample_and_comparison_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            sample_table = tmp / "samples.tsv"
            sample_table.write_text(
                "sample\tcondition\tbam\tpeaks\n"
                "A\tCtrl\tA.bam\tA.bed\n"
                "B\tTreat\tB.bam\tB.bed\n"
                "C\tCtrl\tC.bam\tC.bed\n",
                encoding="utf-8",
            )
            comparison_table = tmp / "comparisons.tsv"
            comparison_table.write_text(
                "comparison\tcond1\tcond2\n"
                "Treat_vs_Ctrl\tTreat\tCtrl\n",
                encoding="utf-8",
            )

            samples = read_sample_table(sample_table)
            comparisons = read_comparison_table(comparison_table)

            self.assertEqual([row.sample for row in samples], ["A", "B", "C"])
            self.assertEqual(comparisons[0].comparison, "Treat_vs_Ctrl")
            self.assertEqual([row.sample for row in samples_for_condition(samples, "Ctrl")], ["A", "C"])
            self.assertEqual(comparison_dir(tmp, comparisons[0].comparison), tmp / "comparisons" / "Treat_vs_Ctrl")
            self.assertEqual(analysis_peaks_path(tmp), tmp / "peaks" / "merged_peaks.analysis.bed")

    def test_diff_footprints_project_comparison_table_expands_replicates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            sample_table = tmp / "samples.tsv"
            sample_table.write_text(
                "sample\tcondition\tbam\tpeaks\n"
                "A1\tCtrl\tA1.bam\tA1.bed\n"
                "A2\tCtrl\tA2.bam\tA2.bed\n"
                "B1\tTreat\tB1.bam\tB1.bed\n",
                encoding="utf-8",
            )
            comparison_table = tmp / "comparisons.tsv"
            comparison_table.write_text(
                "comparison\tcond1\tcond2\n"
                "Treat_vs_Ctrl\tTreat\tCtrl\n",
                encoding="utf-8",
            )
            for sample in ["A1", "A2", "B1"]:
                normalize_dir = tmp / "samples" / sample / "normalize"
                normalize_dir.mkdir(parents=True)
                (normalize_dir / f"{sample}_corrected_q95_scaled.bw").write_text("placeholder", encoding="utf-8")
            args = argparse.Namespace(
                outdir=str(tmp),
                sample_table=str(sample_table),
                comparison_table=str(comparison_table),
                peaks=None,
                project_dir=None,
                sample_dirs=None,
                sample_names=None,
                cond_names=None,
                plot_aggregate="sig",
                motif_outputs="auto",
                skip_excel=False,
            )

            with mock.patch("fp_tools.tools.diff_footprints.run_diff_footprints") as run_one:
                _run_project_comparison_table(args, argparse.ArgumentParser())

            self.assertEqual(run_one.call_count, 1)
            call_args = run_one.call_args.args[0]
            self.assertEqual(call_args.sample_names, ["B1", "A1", "A2"])
            self.assertEqual(call_args.cond_names, ["Treat", "Ctrl", "Ctrl"])
            self.assertEqual(Path(call_args.outdir), tmp / "comparisons" / "Treat_vs_Ctrl")
            self.assertEqual(call_args.plot_aggregate, "sig")
            self.assertEqual(
                [Path(path) for path in call_args.aggregate_signals],
                [
                    tmp / "samples" / "B1" / "normalize" / "B1_corrected_q95_scaled.bw",
                    tmp / "samples" / "A1" / "normalize" / "A1_corrected_q95_scaled.bw",
                    tmp / "samples" / "A2" / "normalize" / "A2_corrected_q95_scaled.bw",
                ],
            )
            self.assertEqual(call_args.motif_outputs, "summary")
            self.assertTrue(call_args.skip_excel)


if __name__ == "__main__":
    unittest.main()
