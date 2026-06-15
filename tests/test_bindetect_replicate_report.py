import tempfile
import unittest
from pathlib import Path

import pandas as pd

from fp_tools.tools.bindetect_replicate_report import (
    build_replicate_report,
    infer_comparisons,
    infer_conditions,
    read_replicate_map,
    replicate_uncertainty,
)


RESULTS = (
    "name\tmotif_id\tcluster\t"
    "treated_mean_score\tcontrol_mean_score\t"
    "treated_bound\tcontrol_bound\t"
    "treated_control_change\ttreated_control_pvalue\ttreated_control_qvalue_bh\ttreated_control_significant_fdr05\n"
    "MA0001\tm1\tc1\t0.8\t0.2\t120\t40\t0.5\t0.0001\t0.0003\tTrue\n"
    "MA0002\tm2\tc1\t0.3\t0.4\t60\t70\t-0.2\t0.2\t0.3\tFalse\n"
    "MA0003\tm3\tc2\t0.6\t0.6\t90\t90\t0.0\t0.6\t0.6\tFalse\n"
)


def _write(tmp: Path, text: str) -> Path:
    path = tmp / "Atest_results.txt"
    path.write_text(text, encoding="utf-8")
    return path


class ReplicateMapParsingTest(unittest.TestCase):
    def test_n_replicates_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            path = tmp / "map.tsv"
            path.write_text(
                "condition\tn_replicates\ntreated\t3\ncontrol\t1\n", encoding="utf-8"
            )
            self.assertEqual(read_replicate_map(path), {"treated": 3, "control": 1})

    def test_condition_replicate_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            path = tmp / "map.tsv"
            path.write_text(
                "condition\treplicate\ntreated\trep1\ntreated\trep2\ncontrol\trepA\n",
                encoding="utf-8",
            )
            self.assertEqual(read_replicate_map(path), {"treated": 2, "control": 1})

    def test_headerless_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            path = tmp / "map.tsv"
            path.write_text("treated\trep1\ntreated\trep2\ncontrol\trepA\n", encoding="utf-8")
            self.assertEqual(read_replicate_map(path), {"treated": 2, "control": 1})

    def test_none_returns_empty(self):
        self.assertEqual(read_replicate_map(None), {})


class ComparisonInferenceTest(unittest.TestCase):
    def test_infer_conditions_and_comparisons(self):
        frame = pd.read_csv(pd.io.common.StringIO(RESULTS), sep="\t")
        conditions = infer_conditions(frame)
        self.assertEqual(set(conditions), {"treated", "control"})
        comparisons = infer_comparisons(frame, conditions)
        self.assertEqual(comparisons, [("treated_control", "treated", "control")])

    def test_no_comparison_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            path = tmp / "bad.txt"
            path.write_text("name\ttreated_mean_score\nMA0001\t0.5\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                build_replicate_report(path, tmp / "out.tsv")


class ReplicateUncertaintyTest(unittest.TestCase):
    def test_se_and_shrinkage(self):
        result = replicate_uncertainty(change=0.5, pvalue=0.0001, n_min=3)
        self.assertGreater(result["z_score"], 0)
        self.assertGreater(result["effect_se"], 0)
        self.assertAlmostEqual(result["shrinkage_weight"], 0.75)
        self.assertAlmostEqual(result["shrunk_change"], 0.5 * 0.75)
        self.assertLess(result["ci_lower"], 0.5)
        self.assertGreater(result["ci_upper"], 0.5)

    def test_more_replicates_weight_closer_to_one(self):
        low = replicate_uncertainty(0.5, 0.01, 1)["shrinkage_weight"]
        high = replicate_uncertainty(0.5, 0.01, 10)["shrinkage_weight"]
        self.assertLess(low, high)

    def test_missing_values_return_nan(self):
        result = replicate_uncertainty(change=float("nan"), pvalue=0.01, n_min=2)
        self.assertNotEqual(result["shrinkage_weight"], 0.0)
        import math
        self.assertTrue(math.isnan(result["effect_se"]))


class ReportOutputTest(unittest.TestCase):
    def test_report_summary_and_figure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            results = _write(tmp, RESULTS)
            rep_map = tmp / "map.tsv"
            rep_map.write_text(
                "condition\tn_replicates\ntreated\t3\ncontrol\t2\n", encoding="utf-8"
            )
            out = tmp / "report.tsv"
            summary_out = tmp / "summary.tsv"
            figure_out = tmp / "report.png"

            report, summary = build_replicate_report(
                results,
                out,
                summary_output=summary_out,
                figure_output=figure_out,
                replicate_map=rep_map,
                alpha=0.05,
            )

            self.assertTrue(out.exists())
            self.assertTrue(summary_out.exists())
            self.assertTrue(figure_out.exists())

            self.assertEqual(len(report), 3)
            self.assertTrue((report["replicate_support"] == "replicate-supported").all())
            sig = report.loc[report["name"] == "MA0001"].iloc[0]
            self.assertTrue(bool(sig["significant"]))
            self.assertTrue(bool(sig["significant_fdr05"]))
            self.assertAlmostEqual(float(sig["qvalue_bh"]), 0.0003)
            self.assertEqual(sig["direction"], "treated")
            nonsig = report.loc[report["name"] == "MA0002"].iloc[0]
            self.assertFalse(bool(nonsig["significant"]))
            self.assertEqual(nonsig["direction"], "control")

            for column in ("effect_se", "z_score", "shrinkage_weight", "shrunk_change", "ci_lower", "ci_upper"):
                self.assertIn(column, report.columns)
            self.assertTrue((report["shrinkage_weight"] > 0).all())

            on_disk = pd.read_csv(out, sep="\t")
            self.assertEqual(len(on_disk), 3)
            row = summary.loc[summary["replicate_support"] == "replicate-supported"].iloc[0]
            self.assertEqual(int(row["n_motifs"]), 3)
            self.assertEqual(int(row["significant"]), 1)
            self.assertEqual(int(row["significant_fdr05"]), 1)

    def test_single_replicate_fallback(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            results = _write(tmp, RESULTS)
            out = tmp / "report.tsv"
            report, _ = build_replicate_report(results, out)
            self.assertTrue((report["replicate_support"] == "single-replicate").all())


if __name__ == "__main__":
    unittest.main()
