import importlib.util
import pathlib
import sys
import unittest

import pandas as pd


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "classify_footprint_failure_modes.py"


def load_module():
    spec = importlib.util.spec_from_file_location("classify_footprint_failure_modes", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


diagnostics = load_module()


def task_rows(tf, current, raw, alternative, **metadata):
    base = {
        "cell": "CELL",
        "tf": tf,
        "motif_id": f"M_{tf}",
        "positive_sites": 1000,
        "coverage_pass": True,
        "depth_plateau": False,
        "protein_supported": True,
        "motif_ambiguous": False,
        "bias_residual": False,
        **metadata,
    }
    return [
        {**base, "method": "current", "auroc": current, "auprc": current - 0.1},
        {**base, "method": "raw", "auroc": raw, "auprc": raw - 0.1},
        {**base, "method": "alternative", "auroc": alternative, "auprc": alternative - 0.1},
    ]


class FootprintFailureModesTest(unittest.TestCase):
    def test_classifies_prespecified_diagnostic_modes(self):
        rows = []
        rows += task_rows("OVERCORRECTED", current=0.65, raw=0.72, alternative=0.68)
        rows += task_rows("SCORER", current=0.64, raw=0.63, alternative=0.75)
        rows += task_rows(
            "NOINFO", current=0.60, raw=0.61, alternative=0.62, depth_plateau=True
        )
        rows += task_rows("DETECTABLE", current=0.82, raw=0.75, alternative=0.80)
        rows += task_rows("AMBIGUOUS", current=0.80, raw=0.75, alternative=0.81, motif_ambiguous=True)
        rows += task_rows("LOWLABELS", current=0.80, raw=0.75, alternative=0.81, positive_sites=20)
        output = diagnostics.classify_failure_modes(
            pd.DataFrame(rows), current_method="current", raw_method="raw"
        ).set_index("tf")
        self.assertEqual(output.loc["OVERCORRECTED", "diagnostic_status"], "correction_sensitive")
        self.assertEqual(output.loc["SCORER", "diagnostic_status"], "scorer_limited")
        self.assertEqual(output.loc["NOINFO", "diagnostic_status"], "atac_information_limited")
        self.assertEqual(output.loc["DETECTABLE", "diagnostic_status"], "detectable")
        self.assertEqual(output.loc["AMBIGUOUS", "diagnostic_status"], "not_callable_motif_ambiguous")
        self.assertEqual(output.loc["LOWLABELS", "diagnostic_status"], "insufficient_orthogonal_labels")

    def test_information_limit_requires_protein_support_and_depth_plateau(self):
        frame = pd.DataFrame(
            task_rows(
                "TF",
                current=0.60,
                raw=0.61,
                alternative=0.62,
                depth_plateau=False,
                protein_supported=True,
            )
        )
        output = diagnostics.classify_failure_modes(frame, "current", "raw")
        self.assertEqual(output.loc[0, "diagnostic_status"], "weak_site_discrimination_unresolved")

    def test_rejects_inconsistent_task_metadata(self):
        frame = pd.DataFrame(task_rows("TF", current=0.7, raw=0.7, alternative=0.7))
        frame.loc[1, "positive_sites"] = 999
        with self.assertRaisesRegex(ValueError, "Inconsistent positive_sites"):
            diagnostics.classify_failure_modes(frame, "current", "raw")

    def test_current_method_requires_an_auroc(self):
        frame = pd.DataFrame(task_rows("TF", current=float("nan"), raw=0.7, alternative=0.8))
        with self.assertRaisesRegex(ValueError, "has no AUROC"):
            diagnostics.classify_failure_modes(frame, "current", "raw")


if __name__ == "__main__":
    unittest.main()
