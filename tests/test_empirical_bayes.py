import unittest
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import norm

from fp_tools.utils.empirical_bayes import (
    benjamini_hochberg,
    estimate_variance_prior,
    fit_moderated_contrast,
    fit_moderated_paired_contrast,
)
from fp_tools.tools.diff_footprints import _apply_replicate_empirical_bayes


class EmpiricalBayesTest(unittest.TestCase):

    def test_paired_contrast_uses_within_replicate_effects(self):
        matrix = pd.DataFrame(
            {
                "rep1": [2.0, -1.0, 0.1, 0.0],
                "rep2": [2.2, -1.1, -0.1, 0.1],
                "rep3": [1.8, -0.9, 0.0, -0.1],
            },
            index=["positive", "negative", "null1", "null2"],
        )
        result = fit_moderated_paired_contrast(matrix)
        self.assertGreater(result.loc["positive", "effect"], 1.5)
        self.assertLess(result.loc["negative", "effect"], -0.8)
        self.assertLess(result.loc["positive", "pvalue"], result.loc["null1", "pvalue"])
        self.assertTrue(np.isfinite(result.loc["positive", "ci_lower"]))

    def setUp(self):
        self.conditions = {
            "treated_1": "treated",
            "treated_2": "treated",
            "treated_3": "treated",
            "control_1": "control",
            "control_2": "control",
            "control_3": "control",
        }

    def test_bh_preserves_missing_values_and_monotonicity(self):
        adjusted = benjamini_hochberg([0.01, np.nan, 0.04, 0.03])
        self.assertTrue(np.isnan(adjusted[1]))
        self.assertTrue(np.allclose(adjusted[[0, 2, 3]], [0.03, 0.04, 0.04]))

    def test_prior_and_moderated_contrast_are_finite(self):
        matrix = pd.DataFrame(
            {
                "treated_1": [3.0, 1.0, 2.0, 4.0],
                "treated_2": [3.2, 1.2, 1.9, 4.4],
                "treated_3": [2.8, 0.8, 2.1, 3.7],
                "control_1": [1.0, 1.1, 2.0, 4.1],
                "control_2": [1.2, 0.9, 2.2, 3.9],
                "control_3": [0.8, 1.0, 1.8, 4.0],
            },
            index=["strong", "null", "flat", "variable"],
        )
        result = fit_moderated_contrast(
            matrix,
            self.conditions,
            "treated",
            "control",
        )
        self.assertAlmostEqual(result.loc["strong", "effect"], 2.0)
        self.assertGreater(result.loc["strong", "moderated_t"], 0)
        self.assertLess(result.loc["strong", "pvalue"], result.loc["null", "pvalue"])
        self.assertTrue((result["moderated_se"] > 0).all())
        self.assertTrue((result["posterior_variance"] > 0).all())
        self.assertTrue((result["moderated_df"] > 4).all())

    def test_effect_direction_follows_requested_contrast(self):
        rng = np.random.default_rng(11)
        base = rng.normal(size=(30, 3))
        matrix = pd.DataFrame(
            np.column_stack([base + 0.7, base]),
            columns=list(self.conditions),
        )
        forward = fit_moderated_contrast(matrix, self.conditions, "treated", "control")
        reverse = fit_moderated_contrast(matrix, self.conditions, "control", "treated")
        self.assertTrue(np.allclose(forward["effect"], -reverse["effect"]))
        self.assertTrue(np.allclose(forward["pvalue"], reverse["pvalue"]))

    def test_unbalanced_three_vs_two_design_uses_all_replicates(self):
        matrix = pd.DataFrame(
            {
                "a_1": [1.0, 2.0, 4.0],
                "a_2": [1.2, 2.2, 4.1],
                "a_3": [0.8, 1.8, 3.9],
                "b_1": [0.0, 2.1, 3.0],
                "b_2": [0.2, 1.9, 3.2],
            },
            index=["m1", "m2", "m3"],
        )
        conditions = {"a_1": "A", "a_2": "A", "a_3": "A", "b_1": "B", "b_2": "B"}
        result = fit_moderated_contrast(matrix, conditions, "A", "B")
        self.assertAlmostEqual(result.loc["m1", "effect"], 0.9)
        self.assertAlmostEqual(result.loc["m3", "effect"], 0.9)
        self.assertTrue(np.isfinite(result["moderated_t"]).all())
        self.assertTrue((result["moderated_df"] > 3).all())

    def test_equal_residual_dispersion_uses_normal_limit(self):
        offsets = np.array([0.0, 0.1, -0.1])
        matrix = pd.DataFrame(
            {
                "treated_1": 1.5 + offsets[0] + np.arange(4),
                "treated_2": 1.5 + offsets[1] + np.arange(4),
                "treated_3": 1.5 + offsets[2] + np.arange(4),
                "control_1": 1.0 + offsets[0] + np.arange(4),
                "control_2": 1.0 + offsets[1] + np.arange(4),
                "control_3": 1.0 + offsets[2] + np.arange(4),
            },
            index=["m1", "m2", "m3", "m4"],
        )
        result = fit_moderated_contrast(
            matrix,
            self.conditions,
            "treated",
            "control",
        )
        self.assertTrue(np.isposinf(result["moderated_df"]).all())
        expected = 2.0 * norm.sf(np.abs(result["moderated_t"]))
        self.assertTrue(np.allclose(result["pvalue"], expected))
        self.assertTrue(np.isfinite(result["ci_lower"]).all())
        self.assertTrue(np.isfinite(result["ci_upper"]).all())

    def test_requires_biological_replicates(self):
        matrix = pd.DataFrame({"treated_1": [1.0, 2.0], "control_1": [0.0, 1.0]})
        with self.assertRaises(ValueError):
            fit_moderated_contrast(
                matrix,
                {"treated_1": "treated", "control_1": "control"},
                "treated",
                "control",
            )

    def test_variance_prior_rejects_no_positive_variance(self):
        with self.assertRaises(ValueError):
            estimate_variance_prior([0.0, 0.0], 4)

    def test_diff_footprints_integration_writes_matrix_and_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            samples = list(self.conditions)
            frame = pd.DataFrame(
                {
                    "total_tfbs": [200, 220, 240, 260],
                    **{
                        f"{sample}_mean_score": values
                        for sample, values in {
                            "treated_1": [0.8, 0.2, 0.5, 0.4],
                            "treated_2": [0.9, 0.3, 0.4, 0.6],
                            "treated_3": [0.7, 0.1, 0.6, 0.5],
                            "control_1": [0.2, 0.2, 0.5, 0.5],
                            "control_2": [0.3, 0.1, 0.6, 0.4],
                            "control_3": [0.1, 0.3, 0.4, 0.6],
                        }.items()
                    },
                },
                index=["m1", "m2", "m3", "m4"],
            )
            args = SimpleNamespace(
                sample_names=samples,
                sample_to_condition=self.conditions,
                condition_replicates={"treated": 3, "control": 3},
                comparisons=[("treated", "control")],
                outdir=tmpdir,
                prefix="probe",
            )
            result = _apply_replicate_empirical_bayes(frame, args)
            self.assertIn("treated_control_ebayes_moderated_t", result.columns)
            self.assertIn("treated_control_ebayes_qvalue_bh", result.columns)
            matrix = Path(tmpdir) / "probe_replicate_motif_score_matrix.tsv"
            self.assertTrue(matrix.exists())
            self.assertEqual(len(pd.read_csv(matrix, sep="\t")), 4)


if __name__ == "__main__":
    unittest.main()
