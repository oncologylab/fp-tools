import importlib.util
from pathlib import Path
import unittest

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts" / "plot_pbmc5k_pseudobulk_markers.py"
SPEC = importlib.util.spec_from_file_location("plot_pbmc5k_pseudobulk_markers", SCRIPT)
plot_pbmc5k = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(plot_pbmc5k)


class DirectionalMarkerRowsTest(unittest.TestCase):
    def test_labels_only_markers_from_compared_groups_in_expected_direction(self):
        plot_df = pd.DataFrame(
            [
                {"name": "PAX5", "change": 0.4, "qvalue": 0.001, "pvalue": 0.001, "neg_log10_p": 3.0, "status": "higher in first"},
                {"name": "EBF1", "change": -0.4, "qvalue": 0.001, "pvalue": 0.001, "neg_log10_p": 3.0, "status": "higher in second"},
                {"name": "CEBPB", "change": -0.5, "qvalue": 0.001, "pvalue": 0.001, "neg_log10_p": 3.0, "status": "higher in second"},
                {"name": "TCF7", "change": -0.6, "qvalue": 0.001, "pvalue": 0.001, "neg_log10_p": 3.0, "status": "higher in second"},
                {"name": "SPIB", "change": 0.6, "qvalue": 0.2, "pvalue": 0.2, "neg_log10_p": 0.7, "status": "not significant"},
            ]
        )

        labels = plot_pbmc5k.directional_marker_rows(plot_df, "B_cell", "Monocyte")

        self.assertEqual(["PAX5", "CEBPB"], labels["name"].tolist())
        self.assertEqual(["B_cell", "Monocyte"], labels["marker_group"].tolist())

    def test_labels_tnk_markers_when_tnk_is_second_condition(self):
        plot_df = pd.DataFrame(
            [
                {"name": "RUNX3", "change": -0.2, "qvalue": 0.001, "pvalue": 0.001, "neg_log10_p": 3.0, "status": "higher in second"},
                {"name": "GATA3", "change": 0.2, "qvalue": 0.001, "pvalue": 0.001, "neg_log10_p": 3.0, "status": "higher in first"},
                {"name": "CEBPA", "change": -0.2, "qvalue": 0.001, "pvalue": 0.001, "neg_log10_p": 3.0, "status": "higher in second"},
            ]
        )

        labels = plot_pbmc5k.directional_marker_rows(plot_df, "B_cell", "T_NK_cell")

        self.assertEqual(["RUNX3"], labels["name"].tolist())
        self.assertEqual(["T_NK_cell"], labels["marker_group"].tolist())


if __name__ == "__main__":
    unittest.main()
