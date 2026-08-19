import tempfile
import unittest
import json
from pathlib import Path
from unittest.mock import patch

from fp_tools.tools.review_multi_comparisons import (
    _compressed_json_b64,
    build_review_payload,
    build_parser,
    count_missing_aggregate_profiles,
    discover_input_htmls,
    fill_missing_aggregate_profiles,
    read_diff_html_payload,
    write_review_html,
)
from fp_tools.tools.static_comparison_browser import build_static_browser


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
        self.assertIn("aggregate-tile-label", html)
        self.assertIn("Comparison ${slot+1}", html)
        self.assertIn("downloadLogoPanel", html)
        self.assertIn("review_multi_comparisons_motif_logo_panel.svg", html)
        self.assertIn("review_multi_comparisons_panel.svg", html)
        self.assertIn("motifLogoSvgFromCounts", html)
        self.assertNotIn("Motif matrix embedded", html)
        self.assertNotIn("motif-detail", html)
        self.assertNotIn("motifDetail", html)
        self.assertNotIn("&#916;FP", html)
        self.assertNotIn("FDR =", html)

    def test_writes_scalable_static_browser_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "browser"
            index = build_static_browser([_diff_payload()], output, "Review")
            self.assertEqual(index, output / "index.html")
            self.assertTrue((output / "app.js").is_file())
            self.assertTrue((output / "styles.css").is_file())
            self.assertTrue((output / "data" / "metadata.json").is_file())
            self.assertTrue((output / "data" / "reports" / "A_vs_B.json.gz").is_file())
            self.assertTrue((output / "data" / "profiles" / "A_vs_B" / "00.json.gz").is_file())
            self.assertFalse((output / "data" / "review_payload.json.gz").exists())
            app = (output / "app.js").read_text(encoding="utf-8")
            self.assertIn("condition-1", app)
            self.assertIn("condition-2", app)
            self.assertIn("profile_shards", app)

    def test_static_browser_records_explicit_defaults(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            payload = _diff_payload()
            payload["aggregate"]["motifs"].append({
                "prefix": "TF2", "name": "TF2", "motif_id": "M2", "n_sites": 8,
                "conditions": [
                    {"name": "A", "n_sites": 5, "samples": [{"name": "A::rep1", "display_name": "rep1", "profile": [0.2, 0.3]}]},
                    {"name": "B", "n_sites": 3, "samples": [{"name": "B::rep1", "display_name": "rep1", "profile": [0.1, 0.2]}]},
                ],
            })
            build_static_browser(
                [payload], root / "browser", "Review",
                default_comparison=("B", "A"),
                default_motifs=["M2", "M1"],
                default_aggregate_plots=2,
                documentation_url="../../../",
            )
            metadata = json.loads((root / "browser/data/metadata.json").read_text())
            app = (root / "browser/app.js").read_text()

        self.assertEqual(metadata["default_comparison"], {"condition1": "B", "condition2": "A"})
        self.assertEqual(metadata["default_aggregate_motifs"], ["TF2", "TF1"])
        self.assertEqual(metadata["default_aggregate_plots"], 2)
        self.assertEqual(metadata["documentation_url"], "../../../")
        self.assertIn("metadata.default_aggregate_motifs", app)
        self.assertIn("regions:", app)
        self.assertIn("mean</title>", app)

    def test_can_fill_missing_aggregate_profiles_before_writing_review(self):
        payload = {
            "schema": "fp-tools.review-multi-comparisons.v1",
            "title": "Review",
            "comparisons": [{"label": "A vs B", "payload": _diff_payload("A vs B")}],
        }
        self.assertEqual(count_missing_aggregate_profiles(payload), (1, 2))
        filled_tf2 = {
            "prefix": "TF2",
            "name": "TF2",
            "motif_id": "M2",
            "n_sites": 6,
            "x": [-1, 1],
            "conditions": [
                {"name": "A", "samples": [{"name": "A_rep1", "profile": [0.3, 0.4]}]},
                {"name": "B", "samples": [{"name": "B_rep1", "profile": [0.1, 0.2]}]},
            ],
        }
        with patch("fp_tools.tools.motif_aggregate_grid.prepare_aggregate_maps", return_value={(0, "TF2"): (filled_tf2, "assembled")}):
            stats = fill_missing_aggregate_profiles(payload, fill_missing=True, recompute_missing=False, aggregate_flank="auto")

        self.assertEqual(stats["filled"], 1)
        self.assertEqual(stats["after_missing"], 0)
        aggregate_prefixes = {m["prefix"] for m in payload["comparisons"][0]["payload"]["aggregate"]["motifs"]}
        self.assertEqual(aggregate_prefixes, {"TF1", "TF2"})
        tf2 = next(m for m in payload["comparisons"][0]["payload"]["aggregate"]["motifs"] if m["prefix"] == "TF2")
        self.assertEqual(tf2["profile_source"], "assembled")
        self.assertEqual(tf2["conditions"][0]["samples"][0]["profile"], [0.3, 0.4])

    def test_writes_eight_panel_review_html(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            inputs = []
            for idx in range(8):
                path = root / f"diff_footprints_C{idx}_D{idx}.html"
                _write_diff_html(path, _diff_payload(f"C{idx} vs D{idx}"))
                inputs.append(path)

            payload = build_review_payload(inputs, title="Review")
            out = root / "review_8.html"
            write_review_html(payload, out, display_panels=8, aggregate_legends="hide")
            html = out.read_text(encoding="utf-8")

        self.assertIn("initialDisplayPanels=8", html)
        self.assertIn('initialAggregateLegends="hide"', html)
        self.assertIn('id="panel-count"', html)
        self.assertIn('id="aggregate-legends"', html)
        self.assertIn("function setPanelCount", html)
        self.assertIn("--comparison-cols", html)
        self.assertIn("--aggregate-cols", html)
        self.assertIn("hide-legends", html)
        self.assertIn("showAggregateLegends", html)
        self.assertIn("function panelColumnCount", html)
        self.assertIn("Math.ceil(n/2)", html)
        self.assertIn("aggregateGrid.style.setProperty('--aggregate-cols',n)", html)
        self.assertIn("cols=tileData.length", html)
        self.assertIn("compact-panels", html)
        self.assertIn("${slotComparisons.length} displayed panels", html)
        self.assertNotIn("slotComparisons=[0,1,2,3]", html)

    def test_parser_accepts_display_panels_option(self):
        parser = build_parser()
        args = parser.parse_args(["--inputs", "a.html", "--output-dir", "review", "--display-panels", "8", "--aggregate-legends", "hide"])
        self.assertEqual(args.output_dir, "review")
        self.assertEqual(args.display_panels, 8)
        self.assertEqual(args.aggregate_legends, "hide")

    def test_parser_accepts_static_browser_defaults(self):
        parser = build_parser()
        args = parser.parse_args([
            "--inputs", "a.html", "--output-dir", "review",
            "--default-comparison", "HNF4A + FOXA2", "No HNF4A/FOXA2",
            "--default-aggregate-motifs", "MA1494.2", "MA0047.4",
            "--default-aggregate-plots", "8",
            "--documentation-url", "../../../",
        ])
        self.assertEqual(args.default_comparison, ["HNF4A + FOXA2", "No HNF4A/FOXA2"])
        self.assertEqual(args.default_aggregate_motifs, ["MA1494.2", "MA0047.4"])
        self.assertEqual(args.default_aggregate_plots, 8)
        self.assertEqual(args.documentation_url, "../../../")

    def test_static_payload_keeps_only_matrices_with_aggregate_profiles(self):
        from fp_tools.tools.static_comparison_browser import split_browser_payload

        payload = _diff_payload("A vs B")
        payload["motif_matrices"]["TF2"] = [[1], [1], [1], [1]]
        compact, _shards = split_browser_payload(payload)
        self.assertEqual(set(compact["motif_matrices"]), {"TF1"})

    def test_parser_accepts_missing_aggregate_profile_options(self):
        parser = build_parser()
        args = parser.parse_args(
            [
                "--outdir",
                "project",
                "--recompute-missing-aggregate-profiles",
                "--fill-missing-aggregate-profiles",
                "--aggregate-flank",
                "100",
                "--cores",
                "4",
            ]
        )
        self.assertTrue(args.recompute_missing_aggregate_profiles)
        self.assertTrue(args.fill_missing_aggregate_profiles)
        self.assertEqual(args.aggregate_flank, "100")
        self.assertEqual(args.cores, 4)


if __name__ == "__main__":
    unittest.main()
