import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fp_tools.tools import motif_aggregate_grid as motif_grid_module
from fp_tools.tools.motif_aggregate_grid import (
    SampleView,
    collect_motifs,
    collect_motifs_from_payloads,
    ordered_comparisons,
    plot_grid_pdf,
    prepare_aggregate_maps,
    write_source_table,
)


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
                    {"name": "Control", "samples": [{"name": "Control_rep1", "profile": [0.1, 0.2, 0.1]}]},
                ],
            }
        )
    return {
        "title": f"Differential footprint report ({label})",
        "conditions": [condition, "Control"],
        "colors": {f"{condition}_up": "#dc2626", "Control_up": "#2563eb"},
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
            {"label": "CellA vs Control", "payload": _payload("CellA vs Control", "CellA", 0.2, True)},
            {"label": "CellB vs Control", "payload": _payload("CellB vs Control", "CellB", -0.4, True)},
            {"label": "CellC vs Control", "payload": _payload("CellC vs Control", "CellC", -0.3, False)},
        ],
    }


class MotifAggregateGridTest(unittest.TestCase):
    def test_preserves_comparison_input_order(self):
        comparisons = ordered_comparisons(_review_payload())
        self.assertEqual([item.condition for item in comparisons], ["CellA", "CellB", "CellC"])
        motifs = collect_motifs(comparisons)
        self.assertEqual([motif.prefix for motif in motifs], ["TF1_M1"])

    def test_writes_pdf_and_source_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            pdf = root / "grid.pdf"
            tsv = root / "grid_source.tsv"
            payload = _review_payload()
            motif_count, page_count = plot_grid_pdf(payload, pdf, rows_per_page=1, flank=60, title="Test")
            row_count = write_source_table(payload, tsv)
            self.assertEqual((motif_count, page_count, row_count), (2, 2, 6))
            self.assertGreater(pdf.stat().st_size, 1000)
            header = tsv.read_text(encoding="utf-8").splitlines()[0]
            self.assertIn("profile_source", header)
            self.assertNotIn("rna", header.lower())

    def test_collects_global_motif_order(self):
        payload2 = _review_payload()
        payload2["comparisons"][0]["payload"]["points"][1]["change"] = -2.5
        motifs = collect_motifs_from_payloads([_review_payload(), payload2])
        self.assertEqual([motif.prefix for motif in motifs[:2]], ["TF2_M2", "TF1_M1"])

    def test_prepare_maps_can_recompute_missing_profiles(self):
        payload = _review_payload()
        sample_views = [
            SampleView("case1", "CellC", Path("/tmp/case.bw"), Path("/tmp/case_match")),
            SampleView("ctrl1", "Control", Path("/tmp/ctrl.bw"), Path("/tmp/ctrl_match")),
        ]
        with (
            patch("fp_tools.tools.motif_aggregate_grid._read_project_samples", return_value=sample_views),
            patch("fp_tools.tools.motif_aggregate_grid._motif_all_bed", return_value=Path("/tmp/motif.bed")),
            patch("pathlib.Path.exists", return_value=True),
            patch("fp_tools.tools.motif_aggregate_grid._read_bed_centers", return_value=[("chr1", 100)]),
            patch("fp_tools.tools.motif_aggregate_grid._mean_profile", return_value=[0.2, 0.1, 0.3, 0.2]),
        ):
            maps = prepare_aggregate_maps(payload, project="/tmp/project", fill_missing=True, recompute_missing=True, flank=2, cores=1)
        comparison = ordered_comparisons(payload)[2]
        aggregate, source = maps[(comparison.index, "TF1_M1")]
        self.assertEqual(source, "recomputed")
        self.assertEqual(aggregate["x"], [-2, -1, 0, 1])

    def test_project_sample_reader_accepts_fallback_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project = Path(tmpdir) / "project"
            (project / "metadata").mkdir(parents=True)
            (project / "metadata" / "sample_table.tsv").write_text(
                "sample\tcondition\tbam\tpeaks\ncase1\tA\tcase.bam\tcase.bed\nctrl1\tB\tctrl.bam\tctrl.bed\n",
                encoding="utf-8",
            )
            samples = motif_grid_module._read_project_samples(project)
            self.assertEqual([(sample.sample, sample.condition) for sample in samples], [("case1", "A"), ("ctrl1", "B")])


if __name__ == "__main__":
    unittest.main()
