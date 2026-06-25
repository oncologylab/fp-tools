import base64
import gzip
import json
import re
import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.plot_aggregate_batch import (
    _compressed_json_b64,
    _discover_motifs,
    _motif_bed_path,
    build_payload,
    merge_payloads,
    read_embedded_payload,
    write_html,
)


def _payload_from_html(path: Path) -> dict:
    html = path.read_text(encoding="utf-8")
    match = re.search(r'reportPayloadB64="([^"]+)"', html)
    assert match is not None
    return json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))


class PlotAggregateBatchTest(unittest.TestCase):
    def test_discovers_motifs_from_match_dir_and_writes_compressed_grid_html(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bed_dir = root / "CTCF_MA0139.1" / "beds"
            bed_dir.mkdir(parents=True)
            (bed_dir / "CTCF_MA0139.1_all.bed").write_text("chr1\t10\t20\tsite1\n", encoding="utf-8")
            motifs = _discover_motifs(root)
            out = root / "report.html"
            write_html({"x": [-1, 1], "samples": [{"sample": "S1", "label": "S1", "condition": "Bcell", "motifs": [{"prefix": "CTCF_MA0139.1", "name": "CTCF", "score": 1, "sites": 1, "profile": [0.1, 0.2]}]}]}, out, "Demo")
            html = out.read_text(encoding="utf-8")
            payload = _payload_from_html(out)

        self.assertEqual(motifs[0]["prefix"], "CTCF_MA0139.1")
        self.assertIn("Demo", html)
        self.assertIn("reportPayloadB64", html)
        self.assertIn("DecompressionStream", html)
        self.assertNotIn("motif-search", html)
        self.assertNotIn("sample-search", html)
        self.assertNotIn("slot-controls", html)
        self.assertNotIn("subplots-panel", html)
        self.assertIn("motif-columns", html)
        self.assertIn("motif-col", html)
        self.assertIn("motif-select", html)
        self.assertIn("sample-panel", html)
        self.assertIn("legend-panel", html)
        self.assertNotIn("Advanced line styles", html)
        self.assertNotIn("sample-picker", html)
        self.assertNotIn("data-slot-sample", html)
        self.assertIn("data-download-panel", html)
        self.assertIn("data-visible", html)
        self.assertIn("data-alpha", html)
        self.assertIn("data-width", html)
        self.assertIn("data-type", html)
        self.assertIn("Group autoscale", html)
        self.assertIn("Download motif logo panel", html)
        self.assertIn("Download motif aggregate panel", html)
        self.assertIn("Download SVG", html)
        self.assertNotIn("Download selected panel SVG", html)
        self.assertNotIn("TF site summary", html)
        self.assertNotIn("summary-sort", html)
        self.assertNotIn("summary-rows", html)
        self.assertNotIn("Union sites", html)
        self.assertNotIn("site counts use plotted BED", html)
        self.assertNotIn("bound sites", html)
        self.assertNotIn("FP score", html)
        self.assertNotIn("footprint score", html)
        self.assertNotIn("Visible TF condition summary", html)
        self.assertNotIn("function drawSummary", html)
        self.assertIn("plot-count", html)
        self.assertNotIn("plot-cols", html)
        self.assertIn("width=340,height=340", html)
        self.assertIn("payload.groups_defined", html)
        self.assertIn("'Samples'", html)
        self.assertNotIn("const data =", html)
        self.assertEqual(payload["schema"], "fp-tools.aggregate.batch.v2")
        self.assertEqual(payload["motifs"][0]["prefix"], "CTCF_MA0139.1")
        self.assertEqual(payload["motifs"][0]["series"][0]["condition"], "Bcell")
        self.assertFalse(payload["groups_defined"])

    def test_build_payload_adds_sample_and_condition_series(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            for sample in ["B1", "B2"]:
                bed_dir = root / sample / "TF1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_all.bed").write_text("chr1\t10\t20\n", encoding="utf-8")
            rows = [
                {"sample": "B1", "label": "B rep1", "condition": "Bcell", "signal": "B1.bw", "match_dir": str(root / "B1")},
                {"sample": "B2", "label": "B rep2", "condition": "Bcell", "signal": "B2.bw", "match_dir": str(root / "B2")},
            ]
            import fp_tools.tools.plot_aggregate_batch as mod
            old = mod._mean_profile
            try:
                mod._mean_profile = lambda signal, centers, flank: [1.0, 3.0] if str(signal) == "B1.bw" else [3.0, 5.0]
                payload = build_payload(rows, flank=1, top_n=1)
            finally:
                mod._mean_profile = old
        series = payload["motifs"][0]["series"]
        self.assertEqual([s["kind"] for s in series].count("sample"), 2)
        self.assertEqual([s for s in series if s["kind"] == "sample"][0]["bed_source"], "all.bed")
        self.assertEqual(payload["motifs"][0]["site_set"], "all.bed")
        self.assertEqual(payload["motifs"][0]["sites"], 1)
        self.assertEqual([s for s in series if s["kind"] == "condition"][0]["sites"], 1)
        self.assertEqual([s for s in series if s["kind"] == "condition"][0]["profile"], [2.0, 4.0])
        self.assertFalse(payload["groups_defined"])

    def test_build_payload_marks_explicit_groups_and_supports_more_than_two_groups(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            rows = []
            for sample, condition in [("S1", "A"), ("S2", "B"), ("S3", "C")]:
                bed_dir = root / sample / "TF1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_all.bed").write_text("chr1\t10\t20\n", encoding="utf-8")
                rows.append({
                    "sample": sample,
                    "label": sample,
                    "condition": condition,
                    "_groups_defined": True,
                    "signal": f"{sample}.bw",
                    "match_dir": str(root / sample),
                })
            import fp_tools.tools.plot_aggregate_batch as mod
            old = mod._mean_profile
            try:
                mod._mean_profile = lambda signal, centers, flank: [1.0, 2.0]
                payload = build_payload(rows, flank=1, top_n=1)
            finally:
                mod._mean_profile = old

        self.assertTrue(payload["groups_defined"])
        self.assertEqual(payload["conditions"], ["A", "B", "C"])

    def test_build_payload_prefers_condition_bound_bed_for_profiles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bed_dir = root / "match" / "TF1" / "beds"
            bed_dir.mkdir(parents=True)
            (bed_dir / "TF1_all.bed").write_text("chr1\t10\t20\nchr1\t30\t40\n", encoding="utf-8")
            (bed_dir / "TF1_Bcell_bound.bed").write_text("chr1\t10\t20\n", encoding="utf-8")
            rows = [{"sample": "B1", "label": "B rep1", "condition": "Bcell", "signal": "B1.bw", "match_dir": str(root / "match")}]
            import fp_tools.tools.plot_aggregate_batch as mod
            seen = []
            old = mod._mean_profile
            try:
                mod._mean_profile = lambda signal, centers, flank: seen.append(len(centers)) or [1.0, 2.0]
                payload = build_payload(rows, flank=1, top_n=1)
                self.assertEqual(_motif_bed_path(root / "match", "TF1", condition="Bcell").name, "TF1_Bcell_bound.bed")
            finally:
                mod._mean_profile = old

        self.assertEqual(seen, [1])
        self.assertEqual(payload["motifs"][0]["series"][0]["sites"], 1)
        self.assertEqual(payload["motifs"][0]["series"][0]["bed_source"], "bound.bed")
        self.assertEqual(payload["motifs"][0]["sites"], 1)
        self.assertEqual(payload["motifs"][0]["site_set"], "bound.bed")

    def test_build_payload_deduplicates_condition_and_union_sites_across_replicates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            for sample in ["B1", "B2"]:
                bed_dir = root / sample / "TF1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_all.bed").write_text("chr1\t10\t20\nchr1\t30\t40\n", encoding="utf-8")
                (bed_dir / "TF1_Bcell_bound.bed").write_text("chr1\t10\t20\n", encoding="utf-8")
            rows = [
                {"sample": "B1", "label": "B rep1", "condition": "Bcell", "signal": "B1.bw", "match_dir": str(root / "B1")},
                {"sample": "B2", "label": "B rep2", "condition": "Bcell", "signal": "B2.bw", "match_dir": str(root / "B2")},
            ]
            import fp_tools.tools.plot_aggregate_batch as mod
            old = mod._mean_profile
            try:
                mod._mean_profile = lambda signal, centers, flank: [1.0, 2.0]
                payload = build_payload(rows, flank=1, top_n=1)
            finally:
                mod._mean_profile = old

        motif = payload["motifs"][0]
        condition = [s for s in motif["series"] if s["kind"] == "condition"][0]
        self.assertEqual([s["sites"] for s in motif["series"] if s["kind"] == "sample"], [1, 1])
        self.assertEqual(condition["sites"], 1)
        self.assertEqual(motif["sites"], 1)

    def test_reads_diff_footprints_html_and_merges_payload(self):
        diff_payload = {
            "colors": {"Bcell_up": "#2563eb", "Tcell_up": "#dc2626"},
            "aggregate": {
                "x": [-1, 1],
                "normalization": "sample-quantile",
                "site_set": "bound",
                "x_label": "Distance from motif center (bp)",
                "y_label": "Quantile-scaled corrected cut-site signal (a.u.)",
                "motifs": [
                    {
                        "prefix": "TF1_MA0001.1",
                        "name": "TF1",
                        "motif_id": "MA0001.1",
                        "n_sites": 2,
                        "conditions": [
                            {"name": "Bcell", "profile": [0.1, 0.2], "fp_score": -0.4, "samples": [{"name": "B1", "profile": [0.1, 0.2], "fp_score": -0.3}]},
                            {"name": "Tcell", "profile": [0.3, 0.4], "fp_score": 0.6, "samples": [{"name": "T1", "profile": [0.3, 0.4], "fp_score": 0.5}]},
                        ],
                    }
                ],
            },
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            src = root / "diff.html"
            src.write_text(f'<script>const reportPayloadB64="{_compressed_json_b64(diff_payload)}";</script>', encoding="utf-8")
            embedded = read_embedded_payload(src)
            merged = merge_payloads([embedded])
            out = root / "batch.html"
            write_html(merged, out, "Batch", default_layout="2x3")
            html = out.read_text(encoding="utf-8")
            payload = _payload_from_html(out)

        self.assertEqual(payload["motifs"][0]["prefix"], "TF1_MA0001.1")
        self.assertEqual(payload["motifs"][0]["series"][0]["label"], "B1")
        self.assertEqual(payload["default_layout"], "2x3")
        self.assertEqual(payload["colors"]["Bcell"], "#2563eb")
        self.assertEqual(payload["colors"]["Tcell"], "#dc2626")
        self.assertEqual(payload["site_set"], "bound.bed")
        self.assertEqual(payload["motifs"][0]["site_set"], "bound.bed")
        self.assertEqual([s for s in payload["motifs"][0]["series"] if s["kind"] == "condition"][0]["fp_score"], -0.4)
        self.assertNotIn("sample-search", html)
        self.assertIn("motif-columns", html)
        self.assertIn("sample-panel", html)
        self.assertNotIn("data-panel-tf", html)
        self.assertNotIn("sample-picker", html)
        self.assertNotIn("waterfall-chart", html)
        self.assertNotIn("TF site summary", html)
        self.assertNotIn("summary-sort", html)
        self.assertNotIn("summary-rows", html)
        self.assertNotIn("Union sites", html)
        self.assertNotIn("bound sites", html)
        self.assertNotIn("FP score", html)
        self.assertIn("plot_aggregate_grid.svg", html)
        self.assertIn("plot_aggregate_motif_logo_panel.svg", html)

    def test_write_html_uses_column_layout_without_summary_sidebar(self):
        payload = {
            "schema": "fp-tools.aggregate.batch.v2",
            "x": [-1, 1],
            "motifs": [
                {
                    "prefix": "TF1",
                    "name": "TF1",
                    "motif_id": "",
                    "score": 1,
                    "sites": 1,
                    "site_set": "bound.bed",
                    "series": [
                        {
                            "id": "sample::S1",
                            "label": "S1",
                            "kind": "sample",
                            "condition": "S1",
                            "profile": [0.1, 0.2],
                            "avg_score": 0.15,
                            "sites": 1,
                            "bed_source": "bound.bed",
                        }
                    ],
                }
            ],
            "conditions": ["S1"],
            "colors": {"S1": "#2563eb"},
            "normalization": "none",
            "site_set": "bound.bed",
            "x_label": "Distance from motif center (bp)",
            "y_label": "Corrected cut-site signal (a.u.)",
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "aggregate.html"
            write_html(payload, out, "Selected TFs", show_summary=False)
            html = out.read_text(encoding="utf-8")
        self.assertNotIn("waterfall-card", html)
        self.assertNotIn("TF site summary", html)
        self.assertIn("motif-columns", html)
        self.assertIn("sample-panel", html)


if __name__ == "__main__":
    unittest.main()
