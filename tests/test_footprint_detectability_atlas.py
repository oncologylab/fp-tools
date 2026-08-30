import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

import pandas as pd


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "build_footprint_detectability_atlas.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_footprint_detectability_atlas", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


atlas = load_module()


def write_result(path, conditions, cluster="C_TF1", first_name="TF1"):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(
        {
            "motif_id": ["M1", "M2", "M3"],
            "name": [first_name, "TF2", "TF3"],
            "output_prefix": ["TF1_M1", "TF2_M2", "TF3_M3"],
            "cluster": [cluster, "C_TF2", "C_TF3"],
            **{
                f"{condition}_mean_score": [scores[0], scores[1], scores[2]]
                for condition, scores in conditions.items()
            },
        }
    )
    frame.to_csv(path, sep="\t", index=False)


class DetectabilityAtlasTest(unittest.TestCase):
    def test_build_collapses_repeated_biological_contexts_and_selects_candidate(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            encode = root / "encode"
            nutrient = root / "nutrient"
            write_result(
                encode / "pairs/A_vs_B/results/diff_footprints_results.txt",
                {"A": [0, 1, 2], "B": [0, 1, 2]},
                cluster="C_TF1A",
            )
            write_result(
                encode / "pairs/A_vs_C/results/diff_footprints_results.txt",
                {"A": [0, 2, 1], "C": [0, 1, 2]},
                cluster="C_TF1B",
            )
            write_result(
                nutrient / "comparisons/stress_vs_control/diff_footprints_results.txt",
                {"stress": [0, 1, 2], "control": [0, 1, 2]},
                cluster="C_TF1A",
            )
            write_result(
                nutrient / "comparisons/other_vs_control/diff_footprints_results.txt",
                {"other": [0, 2, 1], "control": [0, 1, 2]},
                cluster="C_TF1A",
            )
            expression = root / "expression.tsv"
            pd.DataFrame(
                {"gene_key": ["TF1", "TF2", "TF3"], "sample1": [5.0, 5.0, 5.0]}
            ).to_csv(expression, sep="\t", index=False)

            artifacts = atlas.build_atlas(
                encode,
                [("CELL", nutrient)],
                expression,
                root / "out",
                thresholds=atlas.AtlasThresholds(
                    low_percentile=0.34,
                    minimum_low_context_fraction=0.70,
                    minimum_expression=4.0,
                    minimum_encode_contexts=2,
                    minimum_nutrient_contexts=2,
                ),
                path_root=root,
            )

            contexts = pd.read_csv(artifacts.context_scores, sep="\t")
            self.assertEqual(
                contexts.groupby("cohort")["biological_context"].nunique().to_dict(),
                {"ENCODE": 3, "NUTRIENT": 3},
            )
            repeated_a = contexts[
                (contexts["cohort"] == "ENCODE")
                & (contexts["biological_context"] == "A")
                & (contexts["motif_id"] == "M1")
            ].iloc[0]
            self.assertEqual(repeated_a["n_source_analyses"], 2)
            self.assertEqual(repeated_a["cluster"], "C_TF1A;C_TF1B")

            candidates = pd.read_csv(artifacts.candidates, sep="\t")
            self.assertEqual(candidates["motif_id"].tolist(), ["M1"])
            self.assertEqual(candidates.loc[0, "candidate_status"], "weak_shape_expressed")
            self.assertIn("orthogonal occupancy", candidates.loc[0, "interpretation"])

            metadata = json.loads(artifacts.metadata.read_text(encoding="utf-8"))
            self.assertEqual(metadata["contexts"], {"ENCODE": 3, "NUTRIENT": 3})
            self.assertEqual(metadata["weak_shape_candidates"], 1)
            manifest = pd.read_csv(artifacts.input_manifest, sep="\t")
            self.assertTrue(manifest["path"].str.startswith(("encode/", "nutrient/", "expression.tsv")).all())
            report = artifacts.report.read_text(encoding="utf-8")
            self.assertIn("weak aggregate-shape hypotheses", report)
            self.assertIn("TF1", report)

    def test_inconsistent_stable_motif_metadata_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            first = root / "one.tsv"
            second = root / "two.tsv"
            write_result(first, {"A": [0, 1, 2]}, first_name="TF1")
            write_result(second, {"B": [0, 1, 2]}, first_name="DIFFERENT")
            records = pd.concat(
                [
                    atlas.read_result_contexts(first, "ENCODE", "", "one"),
                    atlas.read_result_contexts(second, "ENCODE", "", "two"),
                ],
                ignore_index=True,
            )
            with self.assertRaisesRegex(ValueError, "Inconsistent name metadata for M1"):
                atlas.collapse_biological_contexts(records, 0.10)

    def test_project_spec_requires_cell_and_directory(self):
        self.assertEqual(atlas.parse_project_spec("PANC1=project"), ("PANC1", pathlib.Path("project")))
        with self.assertRaisesRegex(ValueError, "CELL=DIR"):
            atlas.parse_project_spec("project")


if __name__ == "__main__":
    unittest.main()
