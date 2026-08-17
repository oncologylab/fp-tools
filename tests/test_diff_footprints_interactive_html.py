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

from fp_tools.parsers import add_diff_footprints_arguments
import pandas as pd

from fp_tools.tools import diff_footprints
from fp_tools.tools import diff_footprint_helpers
from fp_tools.tools.diff_footprint_helpers import plot_interactive_diff_footprints


class InteractiveDiffFootprintsHtmlTest(unittest.TestCase):
    def test_aggregate_payload_is_serialized_as_compressed_json(self):
        motif = SimpleNamespace(
            name="TF1",
            prefix="TF1_MA0001.1",
            id="MA0001.1",
            group="Bcell_up",
            change=1.2,
            pvalue=0.001,
            base="old-png-should-not-be-embedded-when-counts-exist",
            counts=[
                [10.0, 0.0, 0.0],
                [0.0, 10.0, 0.0],
                [0.0, 0.0, 10.0],
                [0.0, 0.0, 0.0],
            ],
        )
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
            plot_interactive_diff_footprints(
                [motif],
                ["Bcell", "Tcell"],
                str(out),
                aggregate_data=aggregate_data,
                report_label="Method: test label",
            )
            html = out.read_text()
        self.assertIn("const reportPayloadB64=", html)
        self.assertIn("DecompressionStream", html)
        self.assertIn("Download motif logo panel", html)
        self.assertNotIn("currently this is not working", html)
        self.assertIn("Download bar plot panel", html)
        self.assertIn("Download volcano plot", html)
        self.assertIn("Download motif aggregate panel", html)
        self.assertIn("Download combined panel", html)
        self.assertIn("<summary>Show/Hide</summary>", html)
        self.assertNotIn("Sort A-Z", html)
        self.assertIn("function downloadMotifLogoPanel", html)
        self.assertIn("function motifLogoPanelSvg", html)
        self.assertIn("function motifLogoSvg", html)
        self.assertIn("function motifLogoHtml", html)
        self.assertIn("payload.motif_matrices", html)
        self.assertIn("motifLogoHtml(prefix)", html)
        self.assertIn("motifLogoSvg(prefix,", html)
        self.assertIn("diff_footprints_motif_logo_panel.svg", html)
        self.assertNotIn("function downloadSelectedMotifLogo", html)
        self.assertIn('id="download-rank"', html)
        self.assertIn('id="report-method"', html)
        self.assertIn("Method: test label", html)
        self.assertIn("Top differential motifs", html)
        self.assertNotIn("Scale is differential score", html)
        self.assertNotIn("Differential footprint evidence", html)
        self.assertIn("Selected motifs", html)
        self.assertLess(html.index('class="option-col options-samples"'), html.index('class="selected-card option-col"'))
        self.assertIn('id="rank-chart"', html)
        self.assertIn('id="rank-rows"', html)
        self.assertIn('id="rank-rows-slider"', html)
        self.assertIn('type="range"', html)
        self.assertIn('step="1"', html)
        self.assertIn('id="selected-grid"', html)
        self.assertIn('id="aggregate-grid"', html)
        self.assertIn('id="aggregate-legend"', html)
        self.assertNotIn("Panel ${idx+1}", html)
        self.assertIn("Motif for aggregate plot ${idx+1}", html)
        self.assertIn("data-panel-tf", html)
        self.assertIn("selected-head", html)
        self.assertIn("data-sample-visible", html)
        self.assertIn('class="sample-visible"', html)
        self.assertNotIn("data-panel-sample", html)
        self.assertNotIn("samplePickerHtml", html)
        self.assertNotIn("panel-tools", html)
        self.assertNotIn("data-download-panel", html)
        self.assertNotIn(">Download SVG</button>", html)
        self.assertNotIn("diff_footprints_panel_", html)
        self.assertIn("panelPrefixes", html)
        self.assertIn("visibleSamples", html)
        self.assertNotIn("panelSamples", html)
        self.assertIn("function drawTopMotifs", html)
        self.assertIn("perDir", html)
        self.assertIn("positive=points.filter(p=>p.change>0)", html)
        self.assertIn("negative=points.filter(p=>p.change<0)", html)
        self.assertIn("function drawAggregatePanel", html)
        self.assertIn("function setPanelMotif", html)
        self.assertIn("function setSelectedMotif", html)
        self.assertIn("function motifLabel", html)
        self.assertIn("function reportSummary", html)
        self.assertIn("payload.report_label", html)
        self.assertIn("Aggregate normalization: ${norm}; input beds: ${bed}", html)
        self.assertIn("function renderHeader", html)
        self.assertIn("titleCond1.style.color", html)
        self.assertIn("titleCond2.style.color", html)
        self.assertNotIn('id="aggregate-width"', html)
        self.assertNotIn('id="aggregate-show-mean"', html)
        self.assertNotIn("Aggregate profile", html)
        self.assertNotIn("Show mean", html)
        self.assertIn('id="aggregate-sample-styles"', html)
        self.assertIn('data-sample-color=', html)
        self.assertIn('data-sample-alpha=', html)
        self.assertIn('aria-label="Alpha for', html)
        self.assertIn('data-sample-width=', html)
        self.assertIn('aria-label="Line width for', html)
        self.assertIn('data-sample-type=', html)
        self.assertNotIn('id="aggregate-mean-width"', html)
        self.assertNotIn('id="aggregate-mean-type"', html)
        self.assertIn('sample-style-group', html)
        self.assertIn('sample-style-group-title', html)
        self.assertIn('data-sample-group=', html)
        self.assertIn('sample-style-head', html)
        self.assertIn('${escText(cond)} samples', html)
        self.assertIn('<span>Show</span><span>Sample</span>', html)
        self.assertIn("function sampleDisplayName", html)
        self.assertIn("label=sampleDisplayName(name,cond)", html)
        self.assertIn('Adjust ${escText(label)}', html)
        self.assertIn('Alpha for ${escText(label)}', html)
        self.assertIn('<span>Alpha</span>', html)
        self.assertIn('<span>Width</span>', html)
        self.assertIn('function alphaValue', html)
        self.assertIn("aggregateDisplayBp=60", html)
        self.assertIn("p.v>=-aggregateDisplayBp&&p.v<=aggregateDisplayBp", html)
        self.assertIn("xTicks=[-aggregateDisplayBp,0,aggregateDisplayBp]", html)
        self.assertIn("640px", html)
        self.assertIn('grid-template-areas:"actions samples selected"', html)
        self.assertIn('grid-template-areas:"actions samples" "selected selected"', html)
        self.assertIn('grid-template-areas:"actions" "samples" "selected"', html)
        self.assertIn("@media(max-width:700px)", html)
        self.assertIn("grid-template-columns:300px 520px minmax(0,1fr)", html)
        self.assertIn("grid-template-columns:500px 640px minmax(0,1fr)", html)
        self.assertIn(".volcano-card #chart{max-width:628px}", html)
        self.assertIn('viewBox="0 0 760 760"', html)
        self.assertIn("tickStyle='font-size:15px;font-weight:900;font-family:Arial,Helvetica,sans-serif'", html)
        self.assertIn("axisStyle='font-size:17px;font-weight:900;font-family:Arial,Helvetica,sans-serif'", html)
        self.assertIn('style="${tickStyle}"', html)
        self.assertIn('style="${axisStyle}"', html)
        self.assertIn('x="16"', html)
        self.assertIn("rotate(-90 16", html)
        self.assertIn('font-size="14" font-weight="900"', html)
        self.assertIn('font-size="24" font-weight="900"', html)
        self.assertIn("margin={top:64,bottom:68,left:128,right:14}", html)
        self.assertIn("axisY=height-60", html)
        self.assertIn('y1="${margin.top-36}"', html)
        self.assertIn('y="${margin.top-22}"', html)
        self.assertIn('id="plot-count"', html)
        self.assertIn('<option value="12">12</option>', html)
        self.assertIn("function gridShape", html)
        self.assertIn("--aggregate-cols", html)
        self.assertIn("--aggregate-rows", html)
        self.assertIn("plotSvgStyle", html)
        self.assertIn("function styledSvgClone", html)
        self.assertIn("svg,text{font-family:Arial,Helvetica,sans-serif}", html)
        self.assertIn("clone.setAttribute('font-family','Arial,Helvetica,sans-serif')", html)
        self.assertIn("function downloadStandalonePlot", html)
        self.assertIn("downloadStandalonePlot(renderVolcano,chart,'diff_footprints_volcano.svg')", html)
        self.assertIn("downloadStandalonePlot(drawTopMotifs,rankChart,'diff_footprints_barplot.svg')", html)
        self.assertIn("drawFn(false);downloadBlob(svgBlob(node),filename);renderAll(false)", html)
        self.assertNotIn("function downloadCurrentSvg", html)
        self.assertNotIn("function downloadPlainPlot", html)
        self.assertNotIn("function safeFilename", html)
        self.assertIn("function downloadDashboardPanel", html)
        self.assertIn("function renderAggregateLegend", html)
        self.assertIn("function exportLegendSvg", html)
        self.assertIn("aggregate-export-legend", html)
        self.assertIn("gridW+gap", html)
        self.assertIn("diff_footprints_aggregate_grid.svg", html)
        self.assertIn("diff_footprints_panel.svg", html)
        self.assertIn('function visibleSiteCount', html)
        self.assertIn('function visibleConditionNames', html)
        self.assertIn('cond.n_sites!==undefined', html)
        self.assertIn('siteCount=visibleSiteCount(motif,samples)', html)
        self.assertIn('${siteCount}</text>', html)
        self.assertNotIn('${samples.length} sample', html)
        self.assertIn('f<=2.5?2.5', html)
        self.assertIn('Math.abs(value)/5', html)
        self.assertIn('stroke-width="${lineWidthValue(style.width,2)}"${dash}', html)
        self.assertIn("width:stored.width||2", html)
        self.assertIn("Math.max(2,lineWidthValue(row.style.width,2))", html)
        self.assertIn('stroke-opacity="${alphaValue(style.alpha,.9)}"', html)
        self.assertIn('stroke="${style.color}"', html)
        self.assertNotIn('aggregateMeanWidth', html)
        self.assertIn('width=300,height=300', html)
        self.assertIn('viewBox="0 0 ${width} ${height}"', html)
        self.assertLess(html.index('class="motif-logo"'), html.index('class="detail-grid"'))
        self.assertIn("class=\"pt", html)
        self.assertIn("payload.conditions[1]+'_up'", html)
        self.assertIn("payload.conditions[0]+'_up'", html)
        self.assertIn("function allSelectableMotifs", html)
        self.assertIn("localeCompare(motifLabel(b)", html)
        self.assertIn("function defaultPanelPrefixes", html)
        self.assertIn("function targetPlotCount", html)
        self.assertIn("function visiblePanelPrefixes", html)
        self.assertIn("hasAgg=new Set(allAggregateMotifs().map(m=>m.prefix))", html)
        self.assertIn("negN=Math.floor(target/2)", html)
        self.assertIn("posN=target-negN", html)
        self.assertIn("panelPrefixes=panelPrefixes.filter(Boolean).slice(0,12)", html)
        self.assertNotIn("panelPrefixes=panelPrefixes.filter(Boolean).slice(0,target)", html)
        self.assertNotIn("panelPrefixes=allAggregateMotifs().slice(0,target)", html)
        self.assertIn("No aggregate profile", html)
        self.assertIn("Use --plot-aggregate all", html)
        self.assertIn("or increase --plot-aggregate-top-n", html)
        self.assertIn("no aggregate", html)
        self.assertIn("groupColors", html)
        self.assertIn("drawRows(negative)", html)
        self.assertIn("drawRows(positive)", html)
        self.assertNotIn("function drawSection", html)
        self.assertNotIn("sites - ${bedLabel", html)
        self.assertNotIn("motif-site set", html)
        self.assertIn("Distance from motif center (bp)", html)
        self.assertIn("&#916;FP = ${fmtDelta(change)}", html)
        self.assertIn("FDR = ${fmtSci(fdr)}", html)
        self.assertNotIn("<strong>Group:</strong>", html)
        match = re.search(r'const reportPayloadB64="([^"]+)"', html)
        self.assertIsNotNone(match)
        payload = json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))
        self.assertEqual(payload["title"], "Differential footprint report (Bcell vs Tcell)")
        self.assertEqual(payload["report_label"], "Method: test label")
        self.assertEqual(payload["change_label"], "Differential footprint score")
        self.assertEqual(
            payload["motif_matrices"]["TF1_MA0001.1"],
            [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0], [0.0, 0.0, 0.0]],
        )
        self.assertNotIn("TF1_MA0001.1", payload["logos"])
        self.assertNotIn("old-png-should-not-be-embedded-when-counts-exist", html)
        self.assertIn("fdr", payload["points"][0])
        self.assertEqual(payload["colors"]["Bcell_up"], "#dc2626")
        self.assertEqual(payload["colors"]["Tcell_up"], "#2563eb")
        self.assertEqual(payload["aggregate"]["motifs"][0]["prefix"], "TF1_MA0001.1")
        self.assertEqual(payload["aggregate"]["motifs"][0]["motif_id"], "MA0001.1")
        self.assertEqual(payload["aggregate"]["motifs"][0]["conditions"][0]["samples"][0]["name"], "Bcell_rep1")
        self.assertNotIn("json.dumps", html)

    def test_interactive_report_keeps_png_fallback_when_matrix_is_missing(self):
        motif = SimpleNamespace(name="TF1", group="Bcell_up", change=1.2, pvalue=0.001, base="ZmFrZS1wbmc=")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report.html"
            plot_interactive_diff_footprints([motif], ["Bcell", "Tcell"], str(out))
            html = out.read_text()
        match = re.search(r'const reportPayloadB64="([^"]+)"', html)
        self.assertIsNotNone(match)
        payload = json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))
        self.assertEqual(payload["motif_matrices"], {})
        self.assertEqual(payload["logos"]["TF1"]["png"], "data:image/png;base64,ZmFrZS1wbmc=")

    def test_interactive_report_accepts_custom_change_label(self):
        motif = SimpleNamespace(name="TF1", group="Bcell_up", change=0.25, pvalue=0.001, base="")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report.html"
            plot_interactive_diff_footprints(
                [motif],
                ["Bcell", "Tcell"],
                str(out),
                change_label="Mean unique-footprint log2FC",
            )
            html = out.read_text()
        match = re.search(r'const reportPayloadB64="([^"]+)"', html)
        self.assertIsNotNone(match)
        payload = json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))
        self.assertEqual(payload["change_label"], "Mean unique-footprint log2FC")
        self.assertIn("${escText(changeLabel)}", html)
        self.assertIn("&#916;FP = ${fmtDelta(change)}", html)

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

        with patch.object(diff_footprint_helpers, "ProcessPoolExecutor", FakeExecutor):
            with patch.object(diff_footprint_helpers, "_aggregate_payload_for_row", side_effect=fake_row_worker):
                payload = diff_footprint_helpers.build_diff_footprint_aggregate_payload([], info_table, ["Bcell", "Tcell"], args)
        self.assertTrue(FakeExecutor.used)
        self.assertEqual([motif["prefix"] for motif in payload["motifs"]], ["TF1", "TF2"])

    def test_sig_aggregate_mode_caps_significant_rows_by_top_n(self):
        info_table = pd.DataFrame(
            [
                {"output_prefix": "TF1", "name": "TF1", "motif_id": "M1", "Bcell_Tcell_change": 1.0, "Bcell_Tcell_pvalue": 0.001},
                {"output_prefix": "TF2", "name": "TF2", "motif_id": "M2", "Bcell_Tcell_change": 0.8, "Bcell_Tcell_pvalue": 0.002},
                {"output_prefix": "TF3", "name": "TF3", "motif_id": "M3", "Bcell_Tcell_change": -0.6, "Bcell_Tcell_pvalue": 0.003},
                {"output_prefix": "TF4", "name": "TF4", "motif_id": "M4", "Bcell_Tcell_change": 0.1, "Bcell_Tcell_pvalue": 0.5},
            ]
        )
        args = SimpleNamespace(
            aggregate_signals=["B.bw", "T.bw"],
            signals=["B.fp.bw", "T.fp.bw"],
            plot_aggregate="sig",
            plot_aggregate_top_n=2,
            aggregate_pvalue_threshold=0.05,
            aggregate_flank=1,
            outdir=".",
            cond_groups={"Bcell": [0], "Tcell": [1]},
            cores=1,
        )

        def fake_row_worker(task):
            row = task[0]
            return {"prefix": row["output_prefix"], "name": row["name"], "conditions": []}

        with patch.object(diff_footprint_helpers, "_aggregate_payload_for_row", side_effect=fake_row_worker):
            payload = diff_footprint_helpers.build_diff_footprint_aggregate_payload([], info_table, ["Bcell", "Tcell"], args)
        self.assertEqual([motif["prefix"] for motif in payload["motifs"]], ["TF1", "TF2"])

    def test_summary_mode_aggregate_reads_temporary_motif_beds(self):
        info_table = pd.DataFrame(
            [{
                "output_prefix": "TF1",
                "name": "TF1",
                "motif_id": "M1",
                "Bcell_Tcell_change": 1.0,
                "Bcell_Tcell_pvalue": 0.001,
            }]
        )
        args = SimpleNamespace(
            aggregate_signals=["B.bw", "T.bw"],
            signals=["B.fp.bw", "T.fp.bw"],
            plot_aggregate="all",
            aggregate_flank=1,
            aggregate_normalization="none",
            aggregate_site_set="bound",
            outdir="/final/results",
            tmp_tfbs_root="/tmp/summary-sites",
            cond_groups={"Bcell": [0], "Tcell": [1]},
            sample_names=["B_rep1", "T_rep1"],
            cores=1,
        )
        observed_roots = []

        def fake_row_worker(task):
            observed_roots.append(task[2])
            return {"prefix": task[0]["output_prefix"], "conditions": []}

        with patch.object(diff_footprint_helpers, "_aggregate_payload_for_row", side_effect=fake_row_worker):
            payload = diff_footprint_helpers.build_diff_footprint_aggregate_payload(
                [], info_table, ["Bcell", "Tcell"], args
            )
        self.assertEqual(observed_roots, ["/tmp/summary-sites"])
        self.assertEqual(payload["motifs"][0]["prefix"], "TF1")

    def test_summary_mode_writes_only_temporary_aggregate_site_beds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_dir = Path(tmpdir) / "TF1" / "beds"
            bed_dir.mkdir(parents=True)
            args = SimpleNamespace(
                verbosity=0,
                log_q=None,
                write_motif_outputs=False,
                write_cache_motif_all=False,
                aggregate_site_set="bound",
                aggregate_signals=["B.bw"],
                plot_aggregate="all",
                tmp_tfbs_root=tmpdir,
                outdir=tmpdir,
                output_peaks=None,
                cond_names=["Bcell"],
                comparisons=[],
                peak_header_list=[],
                sample_names=["B1"],
                normalization="none",
                thresholds={"Bcell": 0.3},
                condition_samples={"Bcell": ["B1"]},
                condition_replicates={"Bcell": 1},
                per_motif_plots=False,
                skip_excel=True,
                keep_tmp_tfbs_for_cache=False,
            )
            diff_footprint_helpers.process_tfbs("TF1", args, {}, bed_rows=[])
            self.assertTrue((bed_dir / "TF1_all.bed").is_file())
            self.assertTrue((bed_dir / "TF1_Bcell_bound.bed").is_file())
            self.assertFalse((bed_dir / "TF1_Bcell_unbound.bed").exists())
            self.assertFalse((Path(tmpdir) / "TF1" / "TF1_overview.txt").exists())

    def test_cached_match_dirs_provide_aggregate_centers_without_comparison_beds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            all_bed = tmp / "TF1_all.bed"
            b_bound = tmp / "TF1_Bcell_bound.bed"
            t_bound = tmp / "TF1_Tcell_bound.bed"
            all_bed.write_text("chr1\t10\t20\nchr1\t30\t40\n", encoding="utf-8")
            b_bound.write_text("chr1\t10\t20\n", encoding="utf-8")
            t_bound.write_text("chr1\t30\t40\n", encoding="utf-8")

            maps = [
                {"TF1": {"all": str(all_bed), "bound": str(b_bound)}},
                {"TF1": {"all": str(all_bed), "bound": str(t_bound)}},
            ]

            centers_by_condition, all_centers = diff_footprint_helpers._aggregate_centers_for_row(
                str(tmp / "comparison"),
                "TF1",
                ("Bcell", "Tcell"),
                "all",
                aggregate_site_maps=maps,
                cond_groups={"Bcell": [0], "Tcell": [1]},
            )
            self.assertEqual(centers_by_condition["Bcell"], [("chr1", 15), ("chr1", 35)])
            self.assertEqual(centers_by_condition["Tcell"], [("chr1", 15), ("chr1", 35)])
            self.assertEqual(all_centers, [("chr1", 15), ("chr1", 35)])

            centers_by_condition, all_centers = diff_footprint_helpers._aggregate_centers_for_row(
                str(tmp / "comparison"),
                "TF1",
                ("Bcell", "Tcell"),
                "bound",
                aggregate_site_maps=maps,
                cond_groups={"Bcell": [0], "Tcell": [1]},
            )
            self.assertEqual(centers_by_condition["Bcell"], [("chr1", 15)])
            self.assertEqual(centers_by_condition["Tcell"], [("chr1", 35)])
            self.assertEqual(all_centers, [("chr1", 15), ("chr1", 35)])



    def test_aggregate_profile_normalization_uses_report_level_normalizers(self):
        raw_profiles = [[-2.0, 0.0, 2.0], [10.0, 20.0, 30.0]]
        names = ["sample_1", "sample_2"]
        cond_groups = {"Bcell": [0], "Tcell": [1]}
        norm_spec = {
            "sample": {
                "sample_1": diff_footprint_helpers.AggregateAffineNorm(0.0, 2.0, 1.0),
                "sample_2": diff_footprint_helpers.AggregateAffineNorm(20.0, 0.5, 1.0),
            }
        }
        normalized = diff_footprint_helpers._normalize_aggregate_profiles(
            raw_profiles, names, ("Bcell", "Tcell"), cond_groups, "sample-quantile", norm_spec
        )
        self.assertEqual(set(normalized), set(names))
        self.assertEqual(normalized["sample_1"].tolist(), [-3.0, 1.0, 5.0])
        self.assertEqual(normalized["sample_2"].tolist(), [-4.0, 1.0, 6.0])

        unchanged = diff_footprint_helpers._normalize_aggregate_profiles(
            raw_profiles, names, ("Bcell", "Tcell"), cond_groups, "sample-quantile"
        )
        self.assertEqual(unchanged["sample_1"].tolist(), raw_profiles[0])

    def test_aggregate_affine_normalizers_preserve_shape_and_align_scale(self):
        arrays = [pd.Series([-2.0, 0.0, 2.0]).to_numpy(), pd.Series([10.0, 20.0, 30.0]).to_numpy()]
        norms = diff_footprint_helpers._robust_affine_normalizers(arrays, ["a", "b"])
        self.assertEqual(set(norms), {"a", "b"})
        a = norms["a"].normalize(arrays[0])
        b = norms["b"].normalize(arrays[1])
        self.assertLess(a[0], a[1])
        self.assertLess(a[1], a[2])
        self.assertLess(b[0], b[1])
        self.assertLess(b[1], b[2])
        self.assertAlmostEqual(float(pd.Series(a).median()), float(pd.Series(b).median()))

    def test_aggregate_size_factor_normalizers_divide_by_sample_factor(self):
        arrays = [pd.Series([1.0, 2.0, 3.0]).to_numpy(), pd.Series([10.0, 20.0, 30.0]).to_numpy()]
        norms = diff_footprint_helpers._size_factor_normalizers(arrays, ["low", "high"])
        low = norms["low"].normalize(arrays[0])
        high = norms["high"].normalize(arrays[1])
        self.assertEqual(set(norms), {"low", "high"})
        self.assertAlmostEqual(float(low.mean()), float(high.mean()))
        self.assertAlmostEqual(float(low[1] / low[0]), 2.0)
        self.assertAlmostEqual(float(high[1] / high[0]), 2.0)
        self.assertGreater(norms["high"].size_factor, norms["low"].size_factor)

    def test_bound_site_set_uses_condition_specific_bound_beds(self):
        paths = diff_footprint_helpers._aggregate_bed_paths("out", "TF1", ("Bcell", "Tcell"), "bound")
        self.assertEqual(paths["Bcell"], "out/TF1/beds/TF1_Bcell_bound.bed")
        self.assertEqual(paths["Tcell"], "out/TF1/beds/TF1_Tcell_bound.bed")
        all_paths = diff_footprint_helpers._aggregate_bed_paths("out", "TF1", ("Bcell", "Tcell"), "all")
        self.assertEqual(set(all_paths.values()), {"out/TF1/beds/TF1_all.bed"})

    def test_aggregate_center_limit_is_deterministic(self):
        centers = [("chr1", idx) for idx in range(10)]
        limited = diff_footprint_helpers._limit_aggregate_centers(centers, 4)
        self.assertEqual(limited, [("chr1", 0), ("chr1", 3), ("chr1", 6), ("chr1", 9)])
        self.assertEqual(diff_footprint_helpers._limit_aggregate_centers(centers, None), centers)

    def test_aggregate_payload_for_row_keeps_replicate_profiles(self):
        row = {"output_prefix": "TF1_MA0001.1", "name": "TF1", "motif_id": "MA0001.1", "Bcell_Tcell_change": 1.0, "Bcell_Tcell_pvalue_numeric": 0.001}
        cond_groups = {"Bcell": [0, 1], "Tcell": [2, 3]}
        task = (row, ("Bcell", "Tcell"), ".", ["B1.bw", "B2.bw", "T1.bw", "T2.bw"], cond_groups, 1, 2, "Bcell_Tcell", "none", {}, ["Bcell_rep1", "Bcell_rep2", "Tcell_rep1", "Tcell_rep2"], "all")
        profiles = {
            "B1.bw": [1.0, 3.0],
            "B2.bw": [3.0, 5.0],
            "T1.bw": [10.0, 20.0],
            "T2.bw": [30.0, 40.0],
        }
        with patch.object(diff_footprint_helpers, "_read_bed_centers", return_value=[("chr1", 10)]):
            with patch.object(diff_footprint_helpers, "_mean_profile", side_effect=lambda path, centers, flank, norm=None: profiles[path]):
                payload = diff_footprint_helpers._aggregate_payload_for_row(task)
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
            parser = add_diff_footprints_arguments(argparse.ArgumentParser())
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
            with patch.object(diff_footprints, "scan_and_score", side_effect=AssertionError("scan should not run")):
                with patch.object(diff_footprints, "build_diff_footprint_aggregate_payload", return_value=aggregate_data):
                    diff_footprints.run_diff_footprints(args)
            html = (outdir / "diff_footprints_Bcell_Tcell.html").read_text()
        self.assertIn("const reportPayloadB64=", html)
        match = re.search(r'const reportPayloadB64="([^"]+)"', html)
        payload = json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))
        self.assertEqual(payload["aggregate"]["motifs"][0]["prefix"], "TF1_MA0001.1")

    def test_reuse_existing_results_requires_results_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            parser = add_diff_footprints_arguments(argparse.ArgumentParser())
            args = parser.parse_args([
                "--signals", "B1.bw", "T1.bw",
                "--cond-names", "Bcell", "Tcell",
                "--outdir", tmpdir,
                "--prefix", "diff_footprints",
                "--reuse-existing-results",
            ])
            with self.assertRaises(FileNotFoundError):
                diff_footprints.run_diff_footprints(args)


if __name__ == "__main__":
    unittest.main()
