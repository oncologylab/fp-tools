import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.compare_diff_footprints_normalization import compare_results, write_outputs


class CompareDiffFootprintsNormalizationTest(unittest.TestCase):
    def test_compare_results_writes_method_specific_sets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            a = tmp / "corrected_q95.tsv"
            b = tmp / "none.tsv"
            cols = [
                "output_prefix",
                "name",
                "motif_id",
                "K562_HepG2_change",
                "K562_HepG2_pvalue",
                "K562_HepG2_qvalue_bh",
                "K562_HepG2_mean_delta_fp",
                "K562_HepG2_mean_log2fc",
                "K562_HepG2_significant_fdr05",
                "K562_bound",
                "HepG2_bound",
            ]
            pd.DataFrame(
                [
                    ["TF1_M1", "TF1", "M1", 1.0, 0.001, 0.01, 0.2, 0.3, True, 10, 5],
                    ["TF2_M2", "TF2", "M2", -1.0, 0.2, 0.5, -0.1, -0.2, False, 3, 9],
                    ["TF3_M3", "TF3", "M3", -0.5, 0.01, 0.02, -0.1, -0.4, True, 4, 8],
                ],
                columns=cols,
            ).to_csv(a, sep="\t", index=False)
            pd.DataFrame(
                [
                    ["TF1_M1", "TF1", "M1", 0.8, 0.001, 0.02, 0.1, 0.2, True, 11, 6],
                    ["TF2_M2", "TF2", "M2", -0.9, 0.001, 0.03, -0.2, -0.3, True, 2, 10],
                    ["TF3_M3", "TF3", "M3", 0.7, 0.01, 0.04, 0.1, 0.5, True, 7, 4],
                ],
                columns=cols,
            ).to_csv(b, sep="\t", index=False)

            merged = compare_results(a, b, "corrected_q95", "none", "K562_HepG2")
            self.assertEqual(set(merged["significance_class"]), {"shared", "none_only"})
            self.assertEqual(int(merged["direction_flip"].sum()), 1)

            outdir = tmp / "out"
            write_outputs(merged, outdir, "corrected_q95", "none", "K562_HepG2")
            self.assertTrue((outdir / "corrected_q95_vs_none_all_motifs.csv").exists())
            self.assertTrue((outdir / "corrected_q95_only_significant.csv").exists())
            self.assertTrue((outdir / "none_only_significant.csv").exists())
            self.assertTrue((outdir / "shared_significant.csv").exists())
            summary = pd.read_csv(outdir / "significance_summary.csv")
            self.assertIn("direction_flip_shared", set(summary["metric"]))


if __name__ == "__main__":
    unittest.main()
