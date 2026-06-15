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
        self.assertIn("motif-search", html)
        self.assertIn("slot-controls", html)
        self.assertIn("Download grid SVG", html)
        self.assertIn("Download selected panel SVG", html)
        self.assertIn("1x1", html)
        self.assertIn("2x3", html)
        self.assertNotIn("const data =", html)
        self.assertEqual(payload["schema"], "fp-tools.aggregate.batch.v2")
        self.assertEqual(payload["motifs"][0]["prefix"], "CTCF_MA0139.1")
        self.assertEqual(payload["motifs"][0]["series"][0]["condition"], "Bcell")

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
        self.assertEqual([s for s in series if s["kind"] == "condition"][0]["profile"], [2.0, 4.0])

    def test_reads_diff_footprints_html_and_merges_payload(self):
        diff_payload = {
            "colors": {"Bcell_up": "#2563eb", "Tcell_up": "#dc2626"},
            "aggregate": {
                "x": [-1, 1],
                "normalization": "sample-quantile",
                "x_label": "Distance from motif center (bp)",
                "y_label": "Quantile-scaled corrected cut-site signal (a.u.)",
                "motifs": [
                    {
                        "prefix": "TF1_MA0001.1",
                        "name": "TF1",
                        "motif_id": "MA0001.1",
                        "n_sites": 2,
                        "conditions": [
                            {"name": "Bcell", "profile": [0.1, 0.2], "samples": [{"name": "B1", "profile": [0.1, 0.2]}]},
                            {"name": "Tcell", "profile": [0.3, 0.4], "samples": [{"name": "T1", "profile": [0.3, 0.4]}]},
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
        self.assertIn("condition_samples::", html)
        self.assertIn("plot_aggregate_batch_grid.svg", html)


if __name__ == "__main__":
    unittest.main()
