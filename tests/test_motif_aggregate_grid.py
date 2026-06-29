import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fp_tools.tools.motif_aggregate_grid import (
    SampleView,
    collect_motifs,
    collect_motifs_from_payloads,
    nutrient_sort_key,
    ordered_comparisons,
    plot_grid_pdf,
    prepare_aggregate_maps,
    write_source_table,
)
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
            maps = prepare_aggregate_maps(payload, project="/tmp/project", fill_missing=True, recompute_missing=True, flank=2)
        comparison = next(comp for comp in ordered_comparisons(payload) if comp.condition == "0_FBS")
        aggregate, source = maps[(comparison.index, "TF1_M1")]
        self.assertEqual(source, "recomputed")
        self.assertEqual(aggregate["x"], [-2, -1, 0, 1])
        self.assertEqual([condition["name"] for condition in aggregate["conditions"]], ["0_FBS", "10_FBS_Ctrl"])


if __name__ == "__main__":
    unittest.main()
