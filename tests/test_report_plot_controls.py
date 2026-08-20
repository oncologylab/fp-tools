import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.review_multi_comparisons import write_review_html
from fp_tools.tools.static_comparison_browser import build_static_browser


ROOT = Path(__file__).resolve().parents[1]
RESOURCE_ROOT = ROOT / "src" / "fp_tools" / "resources" / "static_browser"


def _payload(include_aggregates: bool) -> dict:
    aggregate_motifs = []
    if include_aggregates:
        aggregate_motifs = [
            {
                "prefix": "JUN_M1",
                "name": "JUN",
                "motif_id": "M1",
                "conditions": [
                    {"name": "KD", "samples": [{"name": "KD_1", "profile": [0.1, 0.2]}]},
                    {"name": "P", "samples": [{"name": "P_1", "profile": [0.2, 0.1]}]},
                ],
            }
        ]
    return {
        "conditions": ["KD", "P"],
        "groups": ["KD_up", "P_up", "n.s."],
        "colors": {"KD_up": "#dc2626", "P_up": "#2563eb", "n.s.": "#8a94a6"},
        "change_label": "Differential footprint score",
        "points": [
            {"prefix": "JUN_M1", "name": "JUN", "motif_id": "M1", "group": "KD_up", "change": 0.2, "pvalue": 1e-10, "fdr": 1e-8, "neglog10p": 10.0},
            {"prefix": "FOS_M2", "name": "FOS", "motif_id": "M2", "group": "P_up", "change": -0.8, "pvalue": 1e-2, "fdr": 0.05, "neglog10p": 2.0},
        ],
        "motif_matrices": {},
        "logos": {},
        "aggregate": {"x": [-1, 1], "motifs": aggregate_motifs},
    }


class ReportPlotControlsTest(unittest.TestCase):
    def test_shared_browser_exposes_requested_controls(self):
        document = (RESOURCE_ROOT / "index.html").read_text(encoding="utf-8")

        self.assertIn('id="rank-sort-toggle"', document)
        self.assertIn('role="switch"', document)
        self.assertIn('id="volcano-highlight"', document)
        self.assertIn('<option value="none">(none)</option>', document)
        self.assertIn('id="volcano-labels"', document)
        self.assertIn('placeholder="e.g. JUN, HIF2A"', document)
        self.assertIn('<script src="plot_controls.js" defer></script>', document)
        app = (RESOURCE_ROOT / "app.js").read_text(encoding="utf-8")
        self.assertIn("function comparisonTitle()", app)
        self.assertIn("Comparison: ${esc(comparisonTitle())}", app)
        self.assertIn(
            'font-family="Helvetica,Arial,sans-serif" font-size="12" font-weight="900" fill="#111827" stroke="none"',
            app,
        )
        self.assertNotIn("paint-order:stroke", app)
        self.assertNotIn(
            '<rect width="${width}" height="${height}" fill="#fff"/><text x="${width / 2}" y="24"',
            app,
        )

    def test_bundle_and_embedded_reports_share_plot_control_code(self):
        review = {
            "schema": "fp-tools.review-multi-comparisons.v1",
            "title": "Review",
            "comparisons": [{"label": "KD vs P", "payload": _payload(True)}],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bundle = root / "bundle"
            build_static_browser([_payload(True)], bundle, "Review")
            embedded = root / "with-aggregates.html"
            write_review_html(review, embedded)
            aggregate_free = root / "aggregate-free.html"
            review["comparisons"][0]["payload"] = _payload(False)
            write_review_html(review, aggregate_free)

            self.assertTrue((bundle / "plot_controls.js").is_file())
            shared = (bundle / "plot_controls.js").read_text(encoding="utf-8")
            with_html = embedded.read_text(encoding="utf-8")
            without_html = aggregate_free.read_text(encoding="utf-8")

        self.assertIn("rankMotifs", shared)
        self.assertIn("matchingMotifs", shared)
        for document in (with_html, without_html):
            self.assertIn('id="rank-sort-toggle"', document)
            self.assertIn('id="volcano-highlight"', document)
            self.assertIn('id="volcano-labels"', document)
            self.assertIn("rankMotifs", document)
            self.assertNotIn('<script src="plot_controls.js"', document)

    @unittest.skipUnless(shutil.which("node"), "Node.js is required for JavaScript behavior checks")
    def test_rank_modes_and_tf_name_matching(self):
        script = RESOURCE_ROOT / "plot_controls.js"
        program = r"""
const controls = require(process.argv[1]);
const motifs = [
  {prefix:'A', name:'JUNB', motif_id:'M1', effect:0.2, neglog10p:10, pvalue:1e-10},
  {prefix:'B', name:'FOS', motif_id:'M2', effect:0.8, neglog10p:2, pvalue:1e-2},
  {prefix:'C', name:'HIF2A', motif_id:'M3', effect:-0.3, neglog10p:12, pvalue:1e-12},
  {prefix:'D', name:'TEAD1', motif_id:'M4', effect:-0.9, neglog10p:1, pvalue:1e-1},
  {prefix:'E', name:'JUN', motif_id:'M5', effect:0.1, neglog10p:3, pvalue:1e-3},
];
const effect = controls.rankMotifs(motifs, 'effect', 2);
const significance = controls.rankMotifs(motifs, 'significance', 2);
const effectPositiveLow = controls.rankColor(motifs[4], 'effect', {max: 10}, {first:'#dc2626', second:'#2563eb', neutralCenter:'#f8fafc'});
const effectPositiveHigh = controls.rankColor(motifs[0], 'effect', {max: 10}, {first:'#dc2626', second:'#2563eb', neutralCenter:'#f8fafc'});
const effectNegativeHigh = controls.rankColor(motifs[2], 'effect', {max: 12}, {first:'#dc2626', second:'#2563eb', neutralCenter:'#f8fafc'});
process.stdout.write(JSON.stringify({
  effect: [...effect.negative, ...effect.positive].map(item => item.prefix),
  significance: [...significance.negative, ...significance.positive].map(item => item.prefix),
  signedPositive: controls.rankMetric(motifs[0], 'significance'),
  signedNegative: controls.rankMetric(motifs[2], 'significance'),
  oppositeEffect: controls.oppositeMetric(motifs[0], 'effect'),
  oppositeSignificance: controls.oppositeMetric(motifs[0], 'significance'),
  effectPositiveLow,
  effectPositiveHigh,
  effectNegativeHigh,
  matches: controls.matchingMotifs(motifs, 'jun, hif2').map(item => item.prefix),
}));
"""
        completed = subprocess.run(
            ["node", "-e", program, str(script)],
            check=True,
            capture_output=True,
            text=True,
        )
        observed = json.loads(completed.stdout)

        self.assertEqual(observed["effect"], ["D", "B"])
        self.assertEqual(observed["significance"], ["C", "A"])
        self.assertEqual(observed["signedPositive"], 10)
        self.assertEqual(observed["signedNegative"], -12)
        self.assertEqual(observed["oppositeEffect"], 10)
        self.assertEqual(observed["oppositeSignificance"], 0.2)
        self.assertNotEqual(observed["effectPositiveLow"], observed["effectPositiveHigh"])
        self.assertEqual(observed["effectPositiveHigh"], "#dc2626")
        self.assertEqual(observed["effectNegativeHigh"], "#2563eb")
        self.assertEqual(observed["matches"], ["C", "E"])


if __name__ == "__main__":
    unittest.main()
