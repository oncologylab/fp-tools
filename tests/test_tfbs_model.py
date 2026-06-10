import tempfile
import unittest
from pathlib import Path

import pandas as pd

from fp_tools.tools.tfbs_model import infer_feature_columns, predict_tabular_model, train_tabular_model


class TfbsModelTest(unittest.TestCase):
    def test_train_and_predict_tabular_model_roundtrip(self):
        frame = pd.DataFrame(
            {
                "site_id": [f"site_{idx}" for idx in range(8)],
                "label": [0, 0, 0, 0, 1, 1, 1, 1],
                "motif_score": [0.1, 0.2, 0.2, 0.3, 1.1, 1.2, 1.3, 1.4],
                "footprint_score": [0.0, 0.1, 0.2, 0.1, 2.0, 2.1, 2.2, 2.3],
                "ms_16": [0.1, 0.1, 0.2, 0.2, 1.0, 1.1, 1.2, 1.3],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            train = tmp / "train.tsv"
            model = tmp / "model.pkl"
            metrics = tmp / "metrics.tsv"
            predictions = tmp / "predictions.tsv"
            frame.to_csv(train, sep="\t", index=False)

            metric_values = train_tabular_model(train, model, metrics_out=metrics, seed=7)
            out = predict_tabular_model(model, train, predictions)

            self.assertTrue(model.exists())
            self.assertTrue(metrics.exists())
            self.assertIn("auroc", metric_values)
            self.assertIn("binding_probability", out.columns)
            self.assertGreater(out.loc[out["label"] == 1, "binding_probability"].mean(), out.loc[out["label"] == 0, "binding_probability"].mean())
            self.assertEqual(len(pd.read_csv(predictions, sep="\t")), len(frame))

    def test_infer_feature_columns_uses_numeric_non_label_columns(self):
        frame = pd.DataFrame({"label": [0, 1], "name": ["a", "b"], "motif_score": [0.1, 1.0]})
        self.assertEqual(infer_feature_columns(frame, label_column="label"), ["motif_score"])


if __name__ == "__main__":
    unittest.main()
