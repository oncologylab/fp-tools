import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from fp_tools.tools.bindetect_functions import build_bindetect_aggregate_payload


def load_review_module():
    script = Path(__file__).resolve().parents[1] / "scripts" / "build_four_method_gsea_review.py"
    spec = importlib.util.spec_from_file_location("build_four_method_gsea_review", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class NormalizationReviewTest(unittest.TestCase):
    def test_method1_uses_bonferroni_delta_and_bound_gate(self):
        review = load_review_module()
        df = pd.DataFrame(
            {
                "output_prefix": ["strong", "small_delta", "few_sites", "weak_p"],
                "name": ["strong", "small_delta", "few_sites", "weak_p"],
                "motif_id": ["M1", "M2", "M3", "M4"],
                "cluster": ["C1", "C2", "C3", "C4"],
                "total_tfbs": [1000, 1000, 100, 1000],
                "K562_bound": [700, 700, 400, 700],
                "HepG2_bound": [650, 650, 300, 650],
                "K562_HepG2_change": [0.2, 0.05, 0.2, 0.2],
                "K562_HepG2_pvalue": [1e-7, 1e-7, 1e-7, 1e-2],
                "K562_HepG2_qvalue_bh": [1e-7, 1e-7, 1e-7, 1e-2],
            }
        )
        out = review.apply_method1_significance(df).set_index("output_prefix")
        self.assertTrue(bool(out.loc["strong", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["small_delta", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["few_sites", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["weak_p", "K562_HepG2_significant_fdr05"]))
        self.assertAlmostEqual(float(out.loc["strong", "K562_HepG2_bonferroni_pvalue"]), 4e-7)

    def test_fdr_methods_use_delta_and_bound_gate(self):
        review = load_review_module()
        df = pd.DataFrame(
            {
                "output_prefix": ["strong", "small_delta", "few_sites", "weak_fdr"],
                "name": ["strong", "small_delta", "few_sites", "weak_fdr"],
                "motif_id": ["M1", "M2", "M3", "M4"],
                "cluster": ["C1", "C2", "C3", "C4"],
                "total_tfbs": [1000, 1000, 100, 1000],
                "K562_bound": [700, 700, 400, 700],
                "HepG2_bound": [650, 650, 300, 650],
                "K562_HepG2_change": [-0.2, -0.05, -0.2, -0.2],
                "K562_HepG2_pvalue": [1e-8, 1e-8, 1e-8, 1e-8],
                "K562_HepG2_qvalue_bh": [1e-4, 1e-4, 1e-4, 0.002],
            }
        )
        out = review.apply_fdr_significance(df).set_index("output_prefix")
        self.assertTrue(bool(out.loc["strong", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["small_delta", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["few_sites", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["weak_fdr", "K562_HepG2_significant_fdr05"]))

    def test_empirical_bayes_log_matrix_uses_log_score_effect_scale(self):
        review = load_review_module()
        matrix = pd.DataFrame(
            [
                [0.20, 0.20, 0.20, 0.10, 0.10, 0.10],
                [0.05, 0.06, 0.07, 0.04, 0.05, 0.06],
            ],
            index=["M1", "M2"],
            columns=review.SAMPLES,
        )
        metadata = pd.DataFrame(
            {
                "output_prefix": ["M1", "M2"],
                "name": ["M1", "M2"],
                "motif_id": ["motif1", "motif2"],
                "cluster": ["C1", "C2"],
                "total_tfbs": [1000, 1000],
                "K562_bound": [700, 700],
                "HepG2_bound": [700, 700],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            native = review.empirical_bayes_log_matrix(matrix, metadata, Path(tmpdir) / "native.tsv").set_index("output_prefix")
        expected = np.log2(0.20 * review.PSEUDOCOUNT_SCALE + 1.0) - np.log2(0.10 * review.PSEUDOCOUNT_SCALE + 1.0)
        raw_delta = 0.20 - 0.10
        self.assertAlmostEqual(float(native.loc["M1", "footprint_score_delta"]), expected)
        self.assertAlmostEqual(float(native.loc["M1", "raw_score_delta"]), raw_delta)
        self.assertGreater(abs(float(native.loc["M1", "footprint_score_delta"])), abs(raw_delta))

    def test_native_bindetect_to_diff_uses_native_change_and_keeps_matrix_audit(self):
        review = load_review_module()
        native = pd.DataFrame(
            {
                "output_prefix": ["M1", "M2"],
                "name": ["M1", "M2"],
                "motif_id": ["motif1", "motif2"],
                "cluster": ["C1", "C2"],
                "total_tfbs": [1000, 1000],
                "K562_bound": [700, 700],
                "HepG2_bound": [700, 700],
                "K562_HepG2_change": [0.2, 0.05],
                "K562_HepG2_pvalue": [1e-8, 1e-8],
                "K562_HepG2_qvalue_bh": [1e-5, 1e-5],
                "K562_HepG2_mean_delta_fp": [0.01, 0.02],
                "K562_HepG2_mean_log2fc": [0.11, 0.12],
            }
        )
        audit = pd.DataFrame(
            {
                "output_prefix": ["M1", "M2"],
                "footprint_score_delta": [0.9, 0.8],
                "raw_score_delta": [0.09, 0.08],
                "pvalue": [0.1, 0.2],
                "padj": [0.3, 0.4],
                "moderated_se": [0.05, 0.06],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            out = review.native_bindetect_to_diff(native, audit, Path(tmpdir) / "results.tsv", "Q95 native").set_index("output_prefix")
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_change"]), 0.2)
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_matrix_log_score_delta"]), 0.9)
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_raw_score_delta"]), 0.09)
        self.assertTrue(bool(out.loc["M1", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["M2", "K562_HepG2_significant_fdr05"]))

    def test_method5_view_uses_limma_effect_and_preserves_native_audit(self):
        review = load_review_module()
        method4 = pd.DataFrame(
            {
                "output_prefix": ["M1", "M2"],
                "name": ["M1", "M2"],
                "motif_id": ["motif1", "motif2"],
                "cluster": ["C1", "C2"],
                "total_tfbs": [1000, 1000],
                "K562_bound": [700, 700],
                "HepG2_bound": [700, 700],
                "K562_HepG2_change": [0.2, 0.2],
                "K562_HepG2_pvalue": [0.5, 1e-8],
                "K562_HepG2_qvalue_bh": [0.5, 1e-5],
                "K562_HepG2_matrix_log_score_delta": [0.35, -0.05],
                "K562_HepG2_matrix_log_score_pvalue": [1e-8, 0.5],
                "K562_HepG2_matrix_log_score_qvalue_bh": [1e-5, 0.5],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            out = review.method5_limma_ebayes_view(method4, Path(tmpdir) / "results.tsv").set_index("output_prefix")
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_change"]), 0.35)
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_limma_effect"]), 0.35)
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_native_bindetect_change"]), 0.2)
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_pvalue"]), 1e-8)
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_native_bindetect_pvalue"]), 0.5)
        self.assertAlmostEqual(float(out.loc["M1", "K562_HepG2_limma_pvalue"]), 1e-8)
        self.assertAlmostEqual(float(out.loc["M2", "K562_HepG2_native_bindetect_pvalue"]), 1e-8)
        self.assertAlmostEqual(float(out.loc["M2", "K562_HepG2_limma_pvalue"]), 0.5)
        self.assertTrue(bool(out.loc["M1", "K562_HepG2_significant_fdr05"]))
        self.assertFalse(bool(out.loc["M2", "K562_HepG2_significant_fdr05"]))

    def test_append_method5_to_old_summary_preserves_existing_methods(self):
        review = load_review_module()
        old = pd.DataFrame(
            {
                "output_prefix": ["M1", "M2"],
                "name": ["M1", "M2"],
                "motif_id": ["motif1", "motif2"],
                "cluster": ["C1", "C2"],
                "total_tfbs": [1000, 1000],
                "method1_tobias_qnorm_delta_score": [0.2, -0.2],
                "method1_tobias_qnorm_significant": [True, False],
                "method1_tobias_qnorm_direction": ["K562_up", "HepG2_up"],
                "method2_qnorm_limma_delta_score": [0.1, -0.1],
                "method2_qnorm_limma_significant": [False, True],
                "method2_qnorm_limma_direction": ["K562_up", "HepG2_up"],
                "n_methods_significant": [1, 1],
                "methods_significant": ["method1_tobias_qnorm", "method2_qnorm_limma"],
                "direction_agreement": [True, True],
            }
        )
        method5 = pd.DataFrame(
            {
                "output_prefix": ["M1", "M2"],
                "name": ["M1", "M2"],
                "motif_id": ["motif1", "motif2"],
                "cluster": ["C1", "C2"],
                "total_tfbs": [1000, 1000],
                "K562_bound": [800, 800],
                "HepG2_bound": [700, 700],
                "K562_HepG2_change": [0.3, -0.3],
                "K562_HepG2_pvalue": [1e-6, 1e-6],
                "K562_HepG2_qvalue_bh": [1e-5, 1e-5],
                "K562_HepG2_significant_fdr05": [True, True],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            old_path = Path(tmpdir) / "old.csv"
            out_path = Path(tmpdir) / "combined.csv"
            old.to_csv(old_path, index=False)
            combined = review.append_method5_to_old_summary(old_path, method5, out_path).set_index("output_prefix")
        self.assertIn("method2_qnorm_limma_significant", combined.columns)
        self.assertIn("method5_q95_limma_ebayes_significant", combined.columns)
        self.assertEqual(int(combined.loc["M1", "n_methods_significant"]), 2)
        self.assertEqual(int(combined.loc["M2", "n_methods_significant"]), 2)
        self.assertEqual(
            combined.loc["M1", "methods_significant"],
            "method1_tobias_qnorm;method5_q95_limma_ebayes",
        )
        self.assertEqual(
            combined.loc["M2", "methods_significant"],
            "method2_qnorm_limma;method5_q95_limma_ebayes",
        )

    def test_aggregate_sig_mode_can_disable_fallback(self):
        review = load_review_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            rows = pd.DataFrame(
                {
                    "output_prefix": ["sig", "nonsig"],
                    "name": ["sig", "nonsig"],
                    "motif_id": ["M1", "M2"],
                    "K562_HepG2_change": [0.2, 0.3],
                    "K562_HepG2_pvalue": [1e-5, 1.0],
                    "K562_HepG2_highlighted": [True, False],
                }
            )
            args = SimpleNamespace(
                outdir=str(root),
                signals=["a.bw", "b.bw"],
                aggregate_signals=["a.bw", "b.bw"],
                cond_groups={"K562": [0], "HepG2": [1]},
                sample_names=["a", "b"],
                normalization="none",
                aggregate_normalization="none",
                aggregate_site_set="all",
                plot_aggregate="sig",
                plot_aggregate_top_n=500,
                aggregate_pvalue_threshold=review.SIG_P_CUTOFF,
                aggregate_sig_only=True,
                aggregate_sig_no_fallback=True,
                aggregate_flank=1,
                aggregate_max_sites=500,
                cores=1,
            )
            motifs = [SimpleNamespace(prefix="sig"), SimpleNamespace(prefix="nonsig")]
            with patch("fp_tools.tools.bindetect_functions._fit_aggregate_normalizers", return_value={}):
                with patch("fp_tools.tools.bindetect_functions._aggregate_payload_for_row", side_effect=lambda task: {"prefix": task[0]["output_prefix"]}):
                    payload = build_bindetect_aggregate_payload(motifs, rows, ("K562", "HepG2"), args)
            self.assertEqual([m["prefix"] for m in payload["motifs"]], ["sig"])

            rows["K562_HepG2_pvalue"] = 1.0
            with patch("fp_tools.tools.bindetect_functions._fit_aggregate_normalizers", return_value={}):
                with patch("fp_tools.tools.bindetect_functions._aggregate_payload_for_row", side_effect=lambda task: {"prefix": task[0]["output_prefix"]}):
                    empty = build_bindetect_aggregate_payload(motifs, rows, ("K562", "HepG2"), args)
            self.assertEqual(empty["motifs"], [])

    def test_make_aggregate_args_uses_all_significant_motifs_from_all_beds(self):
        review = load_review_module()
        args = review.make_aggregate_args(Path("."), ["a.bw", "b.bw"], cores=1)
        self.assertEqual(args.plot_aggregate, "sig")
        self.assertEqual(args.plot_aggregate_top_n, 500)
        self.assertTrue(args.aggregate_sig_no_fallback)
        self.assertTrue(args.aggregate_sig_only)
        self.assertEqual(args.aggregate_pvalue_threshold, 0.001)
        self.assertEqual(args.aggregate_site_set, "all")
        self.assertEqual(args.aggregate_max_sites, None)


if __name__ == "__main__":
    unittest.main()
