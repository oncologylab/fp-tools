import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from fp_tools.tools import diff_footprint_helpers


class DifferentialFootprintStatisticsTest(unittest.TestCase):
    @staticmethod
    def _args(outdir):
        return SimpleNamespace(
            verbosity=0,
            log_q=None,
            write_motif_outputs=False,
            write_cache_motif_all=False,
            aggregate_site_set="all",
            aggregate_signals=None,
            plot_aggregate="off",
            tmp_tfbs_root=outdir,
            outdir=outdir,
            output_peaks=None,
            cond_names=["A", "B"],
            comparisons=[("A", "B")],
            peak_header_list=[],
            sample_names=["A_rep1", "B_rep1"],
            normalization="none",
            thresholds={"A": 0.0, "B": 0.0},
            condition_samples={"A": ["A_rep1"], "B": ["B_rep1"]},
            condition_replicates={"A": 1, "B": 1},
            pseudo=1.0,
            per_motif_plots=False,
            skip_excel=True,
            keep_tmp_tfbs_for_cache=False,
        )

    @staticmethod
    def _row(start, a_score, b_score):
        return ["chr1", start, start + 10, "TF1", 1, "+", a_score, b_score]

    def test_equal_observed_and_background_means_keep_pvalue_one(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [self._row(0, 2.0, 2.0), self._row(20, 4.0, 4.0)]
            result = diff_footprint_helpers.process_tfbs(
                "TF1",
                self._args(tmpdir),
                {("A", "B"): (0.0, 0.5)},
                bed_rows=rows,
            )

        self.assertEqual(float(result.at["TF1", "A_B_change"]), 0.0)
        self.assertEqual(float(result.at["TF1", "A_B_pvalue"]), 1.0)

    def test_zero_variance_unequal_means_produce_finite_statistics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [self._row(0, 2.0, 1.0)]
            result = diff_footprint_helpers.process_tfbs(
                "TF1",
                self._args(tmpdir),
                {("A", "B"): (0.0, 0.0)},
                bed_rows=rows,
            )

        self.assertTrue(np.isfinite(float(result.at["TF1", "A_B_change"])))
        self.assertTrue(np.isfinite(float(result.at["TF1", "A_B_pvalue"])))


if __name__ == "__main__":
    unittest.main()
