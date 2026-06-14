import argparse
import base64
import gzip
import json
import re
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from fp_tools.parsers import add_bindetect_arguments
import pandas as pd

from fp_tools.tools import bindetect
from fp_tools.tools import bindetect_functions
from fp_tools.tools.bindetect_functions import plot_interactive_bindetect


class InteractiveBindetectHtmlTest(unittest.TestCase):
    def test_aggregate_payload_is_serialized_as_compressed_json(self):
        motif = SimpleNamespace(name="TF1", group="Bcell_up", change=1.2, pvalue=0.001, base="")
        aggregate_data = {
            "x": [-1, 0, 1],
            "motifs": [
                {
                    "prefix": "TF1_MA0001.1",
                    "name": "TF1",
                    "motif_id": "MA0001.1",
                    "conditions": [
                        {"name": "Bcell", "profile": [0.2, 0.1, 0.2], "samples": [{"name": "Bcell_rep1", "profile": [0.18, 0.09, 0.19]}]},
                        {"name": "Tcell", "profile": [0.1, 0.3, 0.1], "samples": [{"name": "Tcell_rep1", "profile": [0.11, 0.29, 0.12]}]},
                    ],
                }
            ],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report.html"
            plot_interactive_bindetect([motif], ["Bcell", "Tcell"], str(out), aggregate_data=aggregate_data)
            html = out.read_text()
        self.assertIn("const reportPayloadB64=", html)
        self.assertIn("DecompressionStream", html)
        self.assertIn("Download volcano SVG", html)
        self.assertIn("Download aggregate SVG", html)
        self.assertIn("aggregate-search", html)
        self.assertIn("combo-option", html)
        self.assertIn("function setSelectedMotif", html)
        self.assertIn("function setupAggregateSearch", html)
        self.assertIn("function motifLabel", html)
        self.assertIn("function setAggregateLayout", html)
        self.assertIn('id="aggregate-width"', html)
        self.assertIn("aggregate-wide", html)
        self.assertIn("aggregate-full", html)
        self.assertIn('id="logo-title"', html)
        self.assertIn("logoTitle.textContent", html)
        self.assertIn('stroke-width="0.9" stroke-opacity="0.30"', html)
        self.assertIn('stroke-width="2.2"', html)
        self.assertIn("mean</text>", html)
        self.assertIn("class=\"pt", html)
        self.assertIn("Distance from motif center (bp)", html)
        match = re.search(r'const reportPayloadB64="([^"]+)"', html)
        self.assertIsNotNone(match)
        payload = json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))
        self.assertEqual(payload["aggregate"]["motifs"][0]["prefix"], "TF1_MA0001.1")
        self.assertEqual(payload["aggregate"]["motifs"][0]["motif_id"], "MA0001.1")
        self.assertEqual(payload["aggregate"]["motifs"][0]["conditions"][0]["samples"][0]["name"], "Bcell_rep1")
        self.assertNotIn("json.dumps", html)


    def test_aggregate_payload_uses_parallel_executor_when_multiple_cores(self):
        class FakeExecutor:
            used = False

            def __init__(self, max_workers):
                self.max_workers = max_workers
                FakeExecutor.used = True

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, func, tasks):
                return [func(task) for task in tasks]

        info_table = pd.DataFrame(
            [
                {"output_prefix": "TF1", "name": "TF1", "motif_id": "M1", "Bcell_Tcell_change": 1.0, "Bcell_Tcell_pvalue": 0.001},
                {"output_prefix": "TF2", "name": "TF2", "motif_id": "M2", "Bcell_Tcell_change": -1.0, "Bcell_Tcell_pvalue": 0.002},
            ]
        )
        args = SimpleNamespace(
            aggregate_signals=["B.bw", "T.bw"],
            signals=["B.fp.bw", "T.fp.bw"],
            plot_aggregate="sig",
            plot_aggregate_top_n=20,
            aggregate_pvalue_threshold=0.05,
            aggregate_flank=1,
            outdir=".",
            cond_groups={"Bcell": [0], "Tcell": [1]},
            cores=2,
        )

        def fake_row_worker(task):
            row = task[0]
            return {"prefix": row["output_prefix"], "name": row["name"], "conditions": []}

        with patch.object(bindetect_functions, "ProcessPoolExecutor", FakeExecutor):
            with patch.object(bindetect_functions, "_aggregate_payload_for_row", side_effect=fake_row_worker):
                payload = bindetect_functions.build_bindetect_aggregate_payload([], info_table, ["Bcell", "Tcell"], args)
        self.assertTrue(FakeExecutor.used)
        self.assertEqual([motif["prefix"] for motif in payload["motifs"]], ["TF1", "TF2"])



    def test_aggregate_profile_normalization_uses_report_level_normalizers(self):
        raw_profiles = [[-2.0, 0.0, 2.0], [10.0, 20.0, 30.0]]
        names = ["sample_1", "sample_2"]
        cond_groups = {"Bcell": [0], "Tcell": [1]}
        norm_spec = {
            "sample": {
                "sample_1": bindetect_functions.AggregateAffineNorm(0.0, 2.0, 1.0),
                "sample_2": bindetect_functions.AggregateAffineNorm(20.0, 0.5, 1.0),
            }
        }
        normalized = bindetect_functions._normalize_aggregate_profiles(
            raw_profiles, names, ("Bcell", "Tcell"), cond_groups, "sample-quantile", norm_spec
        )
        self.assertEqual(set(normalized), set(names))
        self.assertEqual(normalized["sample_1"].tolist(), [-3.0, 1.0, 5.0])
        self.assertEqual(normalized["sample_2"].tolist(), [-4.0, 1.0, 6.0])

        unchanged = bindetect_functions._normalize_aggregate_profiles(
            raw_profiles, names, ("Bcell", "Tcell"), cond_groups, "sample-quantile"
        )
        self.assertEqual(unchanged["sample_1"].tolist(), raw_profiles[0])

    def test_aggregate_affine_normalizers_preserve_shape_and_align_scale(self):
        arrays = [pd.Series([-2.0, 0.0, 2.0]).to_numpy(), pd.Series([10.0, 20.0, 30.0]).to_numpy()]
        norms = bindetect_functions._robust_affine_normalizers(arrays, ["a", "b"])
        self.assertEqual(set(norms), {"a", "b"})
        a = norms["a"].normalize(arrays[0])
        b = norms["b"].normalize(arrays[1])
        self.assertLess(a[0], a[1])
        self.assertLess(a[1], a[2])
        self.assertLess(b[0], b[1])
        self.assertLess(b[1], b[2])
        self.assertAlmostEqual(float(pd.Series(a).median()), float(pd.Series(b).median()))

    def test_aggregate_payload_for_row_keeps_replicate_profiles(self):
        row = {"output_prefix": "TF1_MA0001.1", "name": "TF1", "motif_id": "MA0001.1", "Bcell_Tcell_change": 1.0, "Bcell_Tcell_pvalue_numeric": 0.001}
        cond_groups = {"Bcell": [0, 1], "Tcell": [2, 3]}
        task = (row, ("Bcell", "Tcell"), ".", ["B1.bw", "B2.bw", "T1.bw", "T2.bw"], cond_groups, 1, 2, "Bcell_Tcell", "none", {}, ["Bcell_rep1", "Bcell_rep2", "Tcell_rep1", "Tcell_rep2"])
        profiles = {
            "B1.bw": [1.0, 3.0],
            "B2.bw": [3.0, 5.0],
            "T1.bw": [10.0, 20.0],
            "T2.bw": [30.0, 40.0],
        }
        with patch.object(bindetect_functions, "_read_bed_centers", return_value=[("chr1", 10)]):
            with patch.object(bindetect_functions, "_mean_profile", side_effect=lambda path, centers, flank, norm=None: profiles[path]):
                payload = bindetect_functions._aggregate_payload_for_row(task)
        self.assertEqual(payload["motif_id"], "MA0001.1")
        self.assertEqual(payload["conditions"][0]["profile"], [2.0, 4.0])
        self.assertEqual(payload["conditions"][1]["profile"], [20.0, 30.0])
        self.assertEqual([s["name"] for s in payload["conditions"][0]["samples"]], ["Bcell_rep1", "Bcell_rep2"])
        self.assertEqual(payload["conditions"][1]["samples"][1]["profile"], [30.0, 40.0])

    def test_reuse_existing_results_regenerates_html_without_scanning(self):
        aggregate_data = {
            "x": [-1, 1],
            "motifs": [
                {
                    "prefix": "TF1_MA0001.1",
                    "name": "TF1",
                    "n_sites": 1,
                    "conditions": [
                        {"name": "Bcell", "profile": [0.2, 0.2]},
                        {"name": "Tcell", "profile": [0.1, 0.1]},
                    ],
                }
            ],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            motif_dir = outdir / "TF1_MA0001.1"
            beds_dir = motif_dir / "beds"
            beds_dir.mkdir(parents=True)
            (beds_dir / "TF1_MA0001.1_all.bed").write_text("chr1\t10\t20\n")
            (motif_dir / "TF1_MA0001.1.png").write_bytes(b"not-a-real-png")
            (outdir / "diff_footprints_results.txt").write_text(
                "output_prefix\tname\tmotif_id\tcluster\ttotal_tfbs\tBcell_mean_score\tTcell_mean_score\t"
                "Bcell_n_replicates\tTcell_n_replicates\tBcell_Tcell_change\tBcell_Tcell_pvalue\tBcell_Tcell_highlighted\n"
                "TF1_MA0001.1\tTF1\tMA0001.1\tTF1\t1\t1.0\t0.5\t1\t1\t1.25\t1.0E-04\tTrue\n"
            )
            parser = add_bindetect_arguments(argparse.ArgumentParser())
            args = parser.parse_args([
                "--signals", "B1.bw", "T1.bw",
                "--cond-names", "Bcell", "Tcell",
                "--outdir", str(outdir),
                "--prefix", "diff_footprints",
                "--reuse-existing-results",
                "--aggregate-signals", "B1_corrected.bw", "T1_corrected.bw",
                "--plot-aggregate", "top",
                "--replicate-report", "off",
            ])
            with patch.object(bindetect, "scan_and_score", side_effect=AssertionError("scan should not run")):
                with patch.object(bindetect, "build_bindetect_aggregate_payload", return_value=aggregate_data):
                    bindetect.run_bindetect(args)
            html = (outdir / "diff_footprints_Bcell_Tcell.html").read_text()
        self.assertIn("const reportPayloadB64=", html)
        match = re.search(r'const reportPayloadB64="([^"]+)"', html)
        payload = json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))
        self.assertEqual(payload["aggregate"]["motifs"][0]["prefix"], "TF1_MA0001.1")

    def test_reuse_existing_results_requires_results_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            parser = add_bindetect_arguments(argparse.ArgumentParser())
            args = parser.parse_args([
                "--signals", "B1.bw", "T1.bw",
                "--cond-names", "Bcell", "Tcell",
                "--outdir", tmpdir,
                "--prefix", "diff_footprints",
                "--reuse-existing-results",
            ])
            with self.assertRaises(FileNotFoundError):
                bindetect.run_bindetect(args)


if __name__ == "__main__":
    unittest.main()
