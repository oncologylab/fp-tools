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

    def _toy_frame(self):
        return pd.DataFrame(
            {
                "site_id": [f"site_{idx}" for idx in range(8)],
                "label": [0, 0, 0, 0, 1, 1, 1, 1],
                "motif_score": [0.1, 0.2, 0.2, 0.3, 1.1, 1.2, 1.3, 1.4],
                "footprint_score": [0.0, 0.1, 0.2, 0.1, 2.0, 2.1, 2.2, 2.3],
                "ms_16": [0.1, 0.1, 0.2, 0.2, 1.0, 1.1, 1.2, 1.3],
            }
        )

    def test_training_is_reproducible_with_fixed_seed(self):
        frame = self._toy_frame()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            train = tmp / "train.tsv"
            frame.to_csv(train, sep="\t", index=False)

            preds = []
            metrics = []
            for run in range(2):
                model = tmp / f"model_{run}.pkl"
                pred_path = tmp / f"pred_{run}.tsv"
                metrics.append(train_tabular_model(train, model, seed=11))
                out = predict_tabular_model(model, train, pred_path)
                preds.append(out["binding_probability"].to_numpy())

            # Identical seed must yield identical out-of-sample probabilities and metrics.
            self.assertTrue((preds[0] == preds[1]).all())
            self.assertEqual(metrics[0]["auroc"], metrics[1]["auroc"])
            self.assertEqual(metrics[0]["brier"], metrics[1]["brier"])

    def test_metrics_include_calibration_and_ranking_scores(self):
        frame = self._toy_frame()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            train = tmp / "train.tsv"
            model = tmp / "model.pkl"
            frame.to_csv(train, sep="\t", index=False)

            metrics = train_tabular_model(train, model, seed=11)
            for key in ("auroc", "auprc", "brier"):
                self.assertIn(key, metrics)
            # Brier score is a proper-scoring calibration metric bounded in [0, 1].
            self.assertGreaterEqual(metrics["brier"], 0.0)
            self.assertLessEqual(metrics["brier"], 1.0)
            self.assertGreaterEqual(metrics["auroc"], 0.0)
            self.assertLessEqual(metrics["auroc"], 1.0)


if __name__ == "__main__":
    unittest.main()
