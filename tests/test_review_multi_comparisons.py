import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.review_multi_comparisons import (
    _compressed_json_b64,
    build_review_payload,
    discover_input_htmls,
    read_diff_html_payload,
    write_review_html,
)


def _diff_payload(label="A vs B"):
    return {
        "title": f"Differential footprint report ({label})",
        "report_label": f"Method: {label}",
        "conditions": ["A", "B"],
        "groups": ["A_up", "B_up", "n.s."],
        "colors": {"A_up": "#dc2626", "B_up": "#2563eb", "n.s.": "#8a94a6"},
        "change_label": "Differential footprint score",
        "points": [
            {"prefix": "TF1", "name": "TF1", "motif_id": "M1", "group": "A_up", "change": 0.4, "pvalue": 1e-6, "fdr": 1e-4, "neglog10p": 6.0},
            {"prefix": "TF2", "name": "TF2", "motif_id": "M2", "group": "B_up", "change": -0.3, "pvalue": 1e-5, "fdr": 1e-3, "neglog10p": 5.0},
        ],
        "motif_matrices": {
            "TF1": [
                [10, 0, 0, 1],
                [0, 10, 1, 0],
                [0, 0, 10, 0],
                [0, 0, 0, 10],
            ]
        },
        "logos": {},
        "aggregate": {
            "x": [-1, 1],
            "x_label": "Distance from motif center (bp)",
            "y_label": "Corrected cut-site signal",
            "motifs": [
                {
                    "prefix": "TF1",
                    "name": "TF1",
                    "motif_id": "M1",
                    "n_sites": 10,
                    "conditions": [
                        {"name": "A", "samples": [{"name": "A_rep1", "profile": [0.1, 0.2]}]},
                        {"name": "B", "samples": [{"name": "B_rep1", "profile": [0.2, 0.1]}]},
                    ],
                }
            ],
        },
    }


def _write_diff_html(path: Path, payload: dict):
    path.write_text(f'<script>const reportPayloadB64="{_compressed_json_b64(payload)}";</script>', encoding="utf-8")


class ReviewMultiComparisonsTest(unittest.TestCase):
    def test_discovers_direct_child_diff_htmls_from_directories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            first = root / "diff_footprints_A_B.html"
            second = root / "diff_footprints_C_D.html"
            ignored = root / "other.html"
            _write_diff_html(first, _diff_payload("A vs B"))
            _write_diff_html(second, _diff_payload("C vs D"))
            _write_diff_html(ignored, _diff_payload("ignored"))

            self.assertEqual(discover_input_htmls([root]), [first, second])

    def test_discovers_nested_diff_htmls_from_directories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            first = root / "A_vs_B" / "diff_footprints_A_B.html"
            second = root / "C_vs_D" / "nested" / "diff_footprints_C_D.html"
            first.parent.mkdir(parents=True)
            second.parent.mkdir(parents=True)
            _write_diff_html(first, _diff_payload("A vs B"))
            _write_diff_html(second, _diff_payload("C vs D"))

            self.assertEqual(discover_input_htmls([root]), [first, second])

    def test_builds_payload_and_writes_standalone_html(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            first = root / "diff_footprints_A_B.html"
            second = root / "diff_footprints_C_D.html"
            _write_diff_html(first, _diff_payload("A vs B"))
            _write_diff_html(second, _diff_payload("C vs D"))

            payload = build_review_payload([first, second], labels=["one", "two"], title="Review")
            out = root / "review.html"
            write_review_html(payload, out)
            html = out.read_text(encoding="utf-8")
            first_payload = read_diff_html_payload(first)

        self.assertEqual(first_payload["points"][0]["prefix"], "TF1")
        self.assertEqual(payload["schema"], "fp-tools.review-multi-comparisons.v1")
        self.assertEqual([c["label"] for c in payload["comparisons"]], ["one", "two"])
        self.assertIn("reportPayloadB64", html)
        self.assertIn("Comparison ${slot+1}", html)
        self.assertIn("drawRank", html)
        self.assertIn("drawVolcano", html)
        self.assertIn("activePrefix", html)
        self.assertIn('id="rank-rows"', html)
        self.assertIn('id="rank-rows-slider"', html)
        self.assertIn("Top rows", html)
        self.assertIn("syncRankRows", html)
        self.assertIn("niceStep", html)
        self.assertIn("aggregateLegendHtml", html)
        self.assertIn("agg-legend-row", html)
        self.assertIn("function sampleDisplayName", html)
        self.assertIn("label=sampleDisplayName(s,s.condition)", html)
        self.assertIn("position", html)
        self.assertIn("bits", html)
        self.assertIn("1 selected motif", html)
        self.assertIn("Group autoscale", html)
        self.assertIn("downloadLogoPanel", html)
        self.assertIn("review_multi_comparisons_motif_logo_panel.svg", html)
        self.assertIn("review_multi_comparisons_panel.svg", html)
        self.assertIn("motifLogoSvgFromCounts", html)
        self.assertNotIn("Motif matrix embedded", html)
        self.assertNotIn("motif-detail", html)
        self.assertNotIn("motifDetail", html)
        self.assertNotIn("&#916;FP", html)
        self.assertNotIn("FDR =", html)


if __name__ == "__main__":
    unittest.main()
