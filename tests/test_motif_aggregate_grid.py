import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fp_tools.tools.motif_aggregate_grid import (
    SampleView,
    collect_motifs,
    collect_motifs_from_payloads,
    compute_rna_fc_map,
    nutrient_sort_key,
    ordered_comparisons,
    plot_grid_pdf,
    prepare_aggregate_maps,
    write_source_table,
)
from fp_tools.tools import motif_aggregate_grid as motif_grid_module
from fp_tools.tools.review_multi_comparisons import write_review_html
from fp_tools.tools.plot_aggregate_batch import read_embedded_payload


def _payload(label, condition, change, include_aggregate=True):
    aggregate_motifs = []
    if include_aggregate:
        aggregate_motifs.append(
            {
                "prefix": "TF1_M1",
                "name": "TF1",
                "motif_id": "M1",
                "n_sites": 7,
                "conditions": [
                    {"name": condition, "samples": [{"name": f"{condition}_rep1", "profile": [0.2, 0.0, 0.3]}]},
                    {"name": "10_FBS_Ctrl", "samples": [{"name": "10_FBS_Ctrl_rep1", "profile": [0.1, 0.2, 0.1]}]},
                ],
            }
        )
    return {
        "title": f"Differential footprint report ({label})",
        "conditions": [condition, "10_FBS_Ctrl"],
        "colors": {f"{condition}_up": "#dc2626", "10_FBS_Ctrl_up": "#2563eb"},
        "change_label": "Differential footprint score",
        "points": [
            {"prefix": "TF1_M1", "name": "TF1", "motif_id": "M1", "change": change, "pvalue": 0.001, "fdr": 0.01, "group": f"{condition}_up"},
            {"prefix": "TF2_M2", "name": "TF2", "motif_id": "M2", "change": -0.1, "pvalue": 0.2, "fdr": 0.3, "group": "n.s."},
        ],
        "aggregate": {"x": [-60, 0, 60], "motifs": aggregate_motifs},
    }


def _review_payload():
    return {
        "schema": "fp-tools.review-multi-comparisons.v1",
        "title": "Review",
        "comparisons": [
            {"label": "0_Lys vs 10_FBS_Ctrl", "payload": _payload("0_Lys vs 10_FBS_Ctrl", "0_Lys", 0.2, True)},
            {"label": "0.4_FBS vs 10_FBS_Ctrl", "payload": _payload("0.4_FBS vs 10_FBS_Ctrl", "0.4_FBS", -0.4, True)},
            {"label": "0_FBS vs 10_FBS_Ctrl", "payload": _payload("0_FBS vs 10_FBS_Ctrl", "0_FBS", -0.3, False)},
            {"label": "5_Gln.Arg vs 10_FBS_Ctrl", "payload": _payload("5_Gln.Arg vs 10_FBS_Ctrl", "5_Gln.Arg", 0.1, True)},
            {"label": "10_Gln.Arg vs 10_FBS_Ctrl", "payload": _payload("10_Gln.Arg vs 10_FBS_Ctrl", "10_Gln.Arg", 0.5, True)},
        ],
    }


class MotifAggregateGridTest(unittest.TestCase):
    def test_nutrient_order_uses_group_order_and_high_to_low_concentration(self):
        values = ["0_Lys", "0_FBS", "0.4_FBS", "5_Gln.Arg", "10_Gln.Arg", "0.05_Glc"]
        ordered = sorted(values, key=nutrient_sort_key)
        self.assertEqual(ordered, ["0.4_FBS", "0_FBS", "0.05_Glc", "10_Gln.Arg", "5_Gln.Arg", "0_Lys"])

    def test_collects_available_aggregate_motifs_and_orders_comparisons(self):
        comparisons = ordered_comparisons(_review_payload())
        self.assertEqual([c.condition for c in comparisons], ["0.4_FBS", "0_FBS", "10_Gln.Arg", "5_Gln.Arg", "0_Lys"])
        motifs = collect_motifs(comparisons)
        self.assertEqual([m.prefix for m in motifs], ["TF1_M1"])
        self.assertEqual(motifs[0].sort_score, 0.5)

    def test_writes_pdf_and_source_table_from_review_html_payload(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            html = root / "review.html"
            pdf = root / "grid.pdf"
            tsv = root / "grid_source.tsv"
            write_review_html(_review_payload(), html)
            payload = read_embedded_payload(html)
            motif_count, page_count = plot_grid_pdf(payload, pdf, rows_per_page=1, flank=60, title="Test")
            row_count = write_source_table(payload, tsv)

            self.assertEqual(motif_count, 2)
            self.assertEqual(page_count, 2)
            self.assertGreater(pdf.stat().st_size, 1000)
            self.assertEqual(row_count, 10)
            text = tsv.read_text(encoding="utf-8")
            self.assertIn("motif_prefix\tmotif_name\tmotif_id", text.splitlines()[0])
            self.assertIn("profile_source", text.splitlines()[0])
            self.assertIn("delta_fp", text)
            self.assertIn("0_FBS vs 10_FBS_Ctrl", text)
            self.assertIn("False", text)

    def test_row_column_labels_repeat_comparison_names_inside_panels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            pdf = root / "grid_row_labels.pdf"
            labels = []
            original_add_text = motif_grid_module._add_fig_text

            def capture_text(fig, x, y, text, **kwargs):
                labels.append(str(text))
                return original_add_text(fig, x, y, text, **kwargs)

            with patch("fp_tools.tools.motif_aggregate_grid._add_fig_text", side_effect=capture_text):
                motif_count, page_count = plot_grid_pdf(_review_payload(), pdf, rows_per_page=1, flank=60, title="Test", repeat_column_labels="row")

            self.assertEqual(motif_count, 2)
            self.assertEqual(page_count, 2)
            self.assertGreater(pdf.stat().st_size, 1000)
            self.assertGreaterEqual(labels.count("0_Lys"), 4)

    def test_collects_global_motif_order_from_multiple_review_payloads(self):
        payload1 = _review_payload()
        payload2 = _review_payload()
        payload2["comparisons"][0]["payload"]["points"][1]["change"] = -2.5
        motifs = collect_motifs_from_payloads([payload1, payload2])
        self.assertEqual([motif.prefix for motif in motifs[:2]], ["TF2_M2", "TF1_M1"])
        self.assertEqual(motifs[0].sort_score, 2.5)

    def test_prepare_aggregate_maps_can_fill_missing_profiles_from_project_samples(self):
        payload = _review_payload()
        sample_views = [
            SampleView(sample="case1", condition="0_FBS", corrected_bigwig=Path("/tmp/case.bw"), match_dir=Path("/tmp/case_match")),
            SampleView(sample="ctrl1", condition="10_FBS_Ctrl", corrected_bigwig=Path("/tmp/ctrl.bw"), match_dir=Path("/tmp/ctrl_match")),
        ]
        with (
            patch("fp_tools.tools.motif_aggregate_grid._read_project_samples", return_value=sample_views),
            patch("fp_tools.tools.motif_aggregate_grid._motif_all_bed", return_value=Path("/tmp/motif.bed")),
            patch("pathlib.Path.exists", return_value=True),
            patch("fp_tools.tools.motif_aggregate_grid._read_bed_centers", return_value=[("chr1", 100)]),
            patch("fp_tools.tools.motif_aggregate_grid._mean_profile", return_value=[0.2, 0.1, 0.3, 0.2]),
        ):
            maps = prepare_aggregate_maps(payload, project="/tmp/project", fill_missing=True, recompute_missing=True, flank=2, cores=1)
        comparison = next(comp for comp in ordered_comparisons(payload) if comp.condition == "0_FBS")
        aggregate, source = maps[(comparison.index, "TF1_M1")]
        self.assertEqual(source, "recomputed")
        self.assertEqual(aggregate["x"], [-2, -1, 0, 1])
        self.assertEqual([condition["name"] for condition in aggregate["conditions"]], ["0_FBS", "10_FBS_Ctrl"])

    def test_compute_rna_fc_map_uses_motif_gene_map_and_expression_filter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            project = root / "project"
            (project / "metadata").mkdir(parents=True)
            (project / "metadata" / "samples.tsv").write_text(
                "sample\tcondition\ncase1\t0_FBS\ncase2\t0_FBS\nctrl1\t10_FBS_Ctrl\nctrl2\t10_FBS_Ctrl\n",
                encoding="utf-8",
            )
            log2norm = root / "rna.tsv"
            log2norm.write_text(
                "gene_key\tcase1\tcase2\tctrl1\tctrl2\n"
                "TF1\t5\t7\t3\t3\n"
                "TF2\t4\t4\t6\t8\n"
                "LOW\t9\t9\t2\t2\n",
                encoding="utf-8",
            )
            raw = root / "raw.tsv"
            raw.write_text(
                "gene_key\tensembl_gene_id\tHGNC\tcase1\tcase2\tctrl1\tctrl2\n"
                "TF1\tENSG1\tTF1\t10\t12\t8\t9\n"
                "TF2\tENSG2\tTF2\t20\t18\t30\t32\n"
                "LOW\tENSG3\tLOW\t0\t0\t0\t0\n",
                encoding="utf-8",
            )
            motif_map = root / "motifs.tsv"
            motif_map.write_text("motif\tgene_symbol\nTF1_M1\tTF1::LOW\nTF2_M2\tTF2\n", encoding="utf-8")
            payload = _review_payload()
            rna = compute_rna_fc_map(payload, project, None, log2norm, raw, motif_map, min_raw_mean=1.0)
            comparison = next(comp for comp in ordered_comparisons(payload) if comp.condition == "0_FBS")
            tf1 = rna[(comparison.index, "TF1_M1")]
            tf2 = rna[(comparison.index, "TF2_M2")]
            self.assertEqual(tf1.label, "RNA TF1=+3.00")
            self.assertEqual(tf2.label, "RNA TF2=-3.00")

    def test_project_sample_reader_accepts_sample_table_tsv_fallback(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            project = root / "project"
            (project / "metadata").mkdir(parents=True)
            (project / "metadata" / "sample_table.tsv").write_text(
                "sample\tcondition\tbam\tpeaks\n"
                "case1\tA\tcase.bam\tcase.bed\n"
                "ctrl1\tB\tctrl.bam\tctrl.bed\n",
                encoding="utf-8",
            )
            samples = motif_grid_module._read_project_samples(project)
            self.assertEqual([(s.sample, s.condition) for s in samples], [("case1", "A"), ("ctrl1", "B")])
            self.assertEqual(motif_grid_module._read_sample_conditions(project), {"case1": "A", "ctrl1": "B"})


if __name__ == "__main__":
    unittest.main()
