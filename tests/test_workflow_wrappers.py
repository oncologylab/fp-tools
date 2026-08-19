import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.bulk_footprinting import build_commands, build_parser, run_bulk_footprinting
from fp_tools.tools.find_signature_fp import read_marker_sites_from_diff


class WorkflowWrapperTest(unittest.TestCase):
    def test_bulk_wrapper_plans_complete_command_chain(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            samples = root / "samples.tsv"
            comparisons = root / "comparisons.tsv"
            samples.write_text(
                "sample\tcondition\tbam\tpeaks\nA1\tA\tA1.bam\tA1.bed\nB1\tB\tB1.bam\tB1.bed\n",
                encoding="utf-8",
            )
            comparisons.write_text(
                "comparison\tcond1\tcond2\nA_vs_B\tA\tB\n",
                encoding="utf-8",
            )
            parser = build_parser()
            args = parser.parse_args(
                [
                    "--sample-table", str(samples),
                    "--comparison-table", str(comparisons),
                    "--genome", "genome.fa",
                    "--outdir", str(root / "project"),
                    "--dry-run",
                ]
            )
            self.assertEqual(run_bulk_footprinting(args), 0)
            commands = build_commands(args)
            self.assertEqual(
                [label for label, _ in commands],
                ["atac-correct", "call-footprints", "match-motifs", "diff-footprints", "review-multi-comparisons"],
            )
            review = commands[-1][1]
            self.assertIn("--output-dir", review)
            differential = commands[-2][1]
            self.assertIn("--comparison-table", differential)
            self.assertIn("--aggregate-site-set", differential)

    def test_single_cell_marker_sites_come_from_own_diff_results(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results = root / "results.tsv"
            results.write_text(
                "output_prefix\tname\tmotif_id\ttotal_tfbs\n"
                "STAT6_M1\tSTAT6\tM1\t2\n"
                "CEBPA_M2\tCEBPA\tM2\t1\n",
                encoding="utf-8",
            )
            for prefix, lines in {
                "STAT6_M1": "chr1\t10\t14\t.\t8\nchr2\t20\t24\t.\t7\n",
                "CEBPA_M2": "chr3\t30\t34\t.\t6\n",
            }.items():
                bed_dir = root / prefix / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / f"{prefix}_all.bed").write_text(lines, encoding="utf-8")
            sites = read_marker_sites_from_diff(results, root, ["STAT6", "CEBPA"], 10)
            self.assertEqual(sites["STAT6"], [("chr1", 12), ("chr2", 22)])
            self.assertEqual(sites["CEBPA"], [("chr3", 32)])
            all_sites = read_marker_sites_from_diff(results, root, ["STAT6"], 0)
            self.assertEqual(all_sites["STAT6"], [("chr1", 12), ("chr2", 22)])


if __name__ == "__main__":
    unittest.main()
