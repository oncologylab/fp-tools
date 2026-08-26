import argparse
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from fp_tools.tools.diff_footprints import _apply_match_motifs_project_layout, _run_project_comparison_table
from fp_tools.utils.project_layout import (
    analysis_peaks_path,
    comparison_dir,
    merged_peaks_path,
    project_analysis_peaks,
    read_comparison_table,
    read_sample_table,
    samples_for_condition,
    write_analysis_peaks,
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
                call_candidates=False,
                outdir=None,
                sample_names=["A", "B"],
                sample_output_root=str(root),
            )

            items = _scorebigwig_batch_items(args)

            self.assertEqual(items[0][1], str(root / "A" / "footprints" / "A_footprints.bw"))
            self.assertEqual(items[1][1], str(root / "B" / "footprints" / "B_footprints.bw"))
            self.assertIsNone(items[0][2])
            self.assertIsNone(items[1][2])

            args.call_candidates = True
            items = _scorebigwig_batch_items(args)
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
            self.assertEqual(analysis_peaks_path(tmp), tmp / "peaks" / "merged_peaks_filtered.bed")

    def test_project_analysis_peaks_excludes_mitochondrial_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            merged = merged_peaks_path(tmp)
            merged.parent.mkdir(parents=True)
            merged.write_text(
                "chr1\t10\t20\n"
                "chrM\t1\t10\n"
                "MT\t2\t12\n"
                "chr2\t30\t40\n",
                encoding="utf-8",
            )
            analysis = write_analysis_peaks(merged, analysis_peaks_path(tmp))

            self.assertTrue(os.path.samefile(project_analysis_peaks(tmp), analysis))
            self.assertTrue(os.path.samefile(project_analysis_peaks(tmp, merged), analysis))
            self.assertTrue(os.path.samefile(project_analysis_peaks(tmp, analysis), analysis))
            self.assertEqual(analysis.read_text(encoding="utf-8"), "chr1\t10\t20\nchr2\t30\t40\n")

    def test_match_motifs_project_layout_uses_analysis_peaks_when_raw_project_peaks_are_passed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            sample_table = tmp / "samples.tsv"
            sample_table.write_text("sample\tcondition\nA\tA\n", encoding="utf-8")
            raw_peaks = merged_peaks_path(tmp)
            raw_peaks.parent.mkdir(parents=True)
            raw_peaks.write_text("chr1\t1\t10\nchrM\t1\t10\n", encoding="utf-8")
            write_analysis_peaks(raw_peaks, analysis_peaks_path(tmp))
            footprint_dir = tmp / "samples" / "A" / "footprints"
            footprint_dir.mkdir(parents=True)
            (footprint_dir / "A_footprints.bw").write_text("placeholder", encoding="utf-8")
            args = argparse.Namespace(
                layout="project",
                sample_table=str(sample_table),
                outdir=str(tmp),
                peaks=str(raw_peaks),
                signals=None,
                sample_names=None,
                cond_names=None,
                sample_output_root=None,
            )

            _apply_match_motifs_project_layout(args, argparse.ArgumentParser())

            self.assertTrue(os.path.samefile(args.peaks, analysis_peaks_path(tmp)))
            self.assertEqual(args.signals, [str(footprint_dir / "A_footprints.bw")])

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
            self.assertEqual(
                Path(call_args.outdir).parts[-2:],
                ("comparisons", "Treat_vs_Ctrl"),
            )
            self.assertEqual(call_args.plot_aggregate, "sig")
            expected_signals = [
                tmp / "samples" / "B1" / "normalize" / "B1_corrected_q95_scaled.bw",
                tmp / "samples" / "A1" / "normalize" / "A1_corrected_q95_scaled.bw",
                tmp / "samples" / "A2" / "normalize" / "A2_corrected_q95_scaled.bw",
            ]
            self.assertTrue(
                all(os.path.samefile(observed, expected) for observed, expected in zip(call_args.aggregate_signals, expected_signals))
            )
            self.assertEqual(call_args.motif_outputs, "summary")
            self.assertTrue(call_args.skip_excel)


if __name__ == "__main__":
    unittest.main()
