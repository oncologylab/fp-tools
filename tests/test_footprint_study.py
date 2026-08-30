import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd


ROOT = pathlib.Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


study = load_module(
    "validate_footprint_study",
    "benchmarks/scripts/validate_footprint_study.py",
)
labels = load_module(
    "build_footprint_site_labels",
    "benchmarks/scripts/build_footprint_site_labels.py",
)


class FootprintStudySpecTest(unittest.TestCase):
    def test_committed_study_spec_is_valid_and_locked(self):
        path = ROOT / "benchmarks/manifests/footprint_detectability_v1.spec.json"
        spec = study.load_spec(path)
        self.assertEqual(study.validate_spec(spec), [])
        self.assertEqual(len(spec["tasks"]), 35)
        self.assertEqual(
            {task["cell"] for task in spec["tasks"] if task["split"] == "development"},
            {"K562", "HepG2"},
        )
        self.assertEqual(spec["nutrient_application"]["external_pdac_accession"], "GSE144833")

    def test_spec_rejects_chromosome_leakage_and_duplicate_tasks(self):
        path = ROOT / "benchmarks/manifests/footprint_detectability_v1.spec.json"
        spec = json.loads(path.read_text(encoding="utf-8"))
        spec["chromosome_split"]["validation"].append("chr1")
        spec["tasks"].append(dict(spec["tasks"][0]))
        errors = study.validate_spec(spec)
        self.assertTrue(any("multiple splits" in error for error in errors))
        self.assertTrue(any("duplicate task" in error for error in errors))


class FootprintSiteLabelTest(unittest.TestCase):
    def test_summit_supported_positive_far_negative_and_indeterminate(self):
        peaks = {"chr1": [labels.Peak(100, 200, 150)]}
        sites = pd.DataFrame(
            [
                {"chrom": "chr1", "start": 145, "end": 155, "strand": "+", "site_id": "p", "motif_score": 8.0},
                {"chrom": "chr1", "start": 105, "end": 115, "strand": "+", "site_id": "inside", "motif_score": 7.0},
                {"chrom": "chr1", "start": 250, "end": 260, "strand": "+", "site_id": "near", "motif_score": 6.0},
                {"chrom": "chr1", "start": 800, "end": 810, "strand": "+", "site_id": "n", "motif_score": 5.0},
                {"chrom": "chr2", "start": 10, "end": 20, "strand": "+", "site_id": "empty", "motif_score": 4.0},
            ]
        )
        output = labels.label_sites(
            sites, peaks, positive_summit_distance=20, negative_peak_distance=100
        ).set_index("site_id")
        self.assertEqual(int(output.loc["p", "label"]), 1)
        self.assertEqual(int(output.loc["inside", "label"]), -1)
        self.assertEqual(int(output.loc["near", "label"]), -1)
        self.assertEqual(int(output.loc["n", "label"]), 0)
        self.assertEqual(int(output.loc["empty", "label"]), 0)

    def test_reads_narrowpeak_summit_offset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = pathlib.Path(tmpdir) / "labels.narrowPeak"
            path.write_text(
                "chr1\t100\t200\tpeak\t1000\t.\t5\t10\t8\t45\n",
                encoding="utf-8",
            )
            peaks = labels.read_peaks(path)
            self.assertEqual(peaks["chr1"], [labels.Peak(100, 200, 145)])

    def test_propensity_matching_is_deterministic_and_balanced(self):
        rng = np.random.default_rng(4)
        frame = pd.DataFrame(
            {
                "chrom": ["chr1"] * 30,
                "start": np.arange(30) * 10,
                "end": np.arange(30) * 10 + 5,
                "strand": ["+"] * 30,
                "site_id": [f"s{i}" for i in range(30)],
                "motif_score": np.r_[rng.normal(2, 0.2, 10), rng.normal(2, 0.5, 20)],
                "accessibility": np.r_[rng.normal(5, 0.2, 10), rng.normal(5, 0.5, 20)],
                "label": [1] * 10 + [0] * 20,
                "label_reason": ["positive"] * 10 + ["negative"] * 20,
                "nearest_peak_distance": [0] * 10 + [1000] * 20,
                "nearest_summit_distance": [0] * 10 + [1000] * 20,
            }
        )
        first = labels.propensity_match(
            frame, ["motif_score", "accessibility"], negative_ratio=1, seed=7
        )
        second = labels.propensity_match(
            frame, ["motif_score", "accessibility"], negative_ratio=1, seed=7
        )
        self.assertEqual(first["label"].value_counts().to_dict(), {1: 10, 0: 10})
        pd.testing.assert_frame_equal(first, second)


if __name__ == "__main__":
    unittest.main()
