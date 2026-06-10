import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.rerank import read_table, rerank_sites


class CandidateRerankTest(unittest.TestCase):
    def test_rerank_sites_combines_scores_and_family_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            sites = tmp / "sites.tsv"
            families = tmp / "families.tsv"
            out = tmp / "ranked.tsv"
            sites.write_text(
                "chrom\tstart\tend\tname\tmotif_id\tcandidate_score\tbinding_probability\n"
                "chr1\t10\t20\ta\tM1\t1.0\t0.90\n"
                "chr1\t30\t40\tb\tM2\t5.0\t0.20\n"
                "chr1\t50\t60\tc\tM3\t4.0\t0.80\n",
                encoding="utf-8",
            )
            families.write_text("M1\tETS\nM2\tAP1\nM3\tETS\n", encoding="utf-8")

            ranked = rerank_sites(
                sites,
                out,
                score_columns=["binding_probability", "candidate_score"],
                weights=[2.0, 1.0],
                family_map=families,
                family_bonus=0.1,
            )
            self.assertTrue(out.exists())

        self.assertEqual(ranked.iloc[0]["name"], "c")
        self.assertEqual(ranked.iloc[0]["motif_family"], "ETS")
        self.assertEqual(ranked["rank"].tolist(), [1, 2, 3])

    def test_rerank_reads_fp_tools_candidate_header_and_limits_family(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            sites = tmp / "candidates.bed"
            out = tmp / "ranked.tsv"
            sites.write_text(
                "#chrom\tstart\tend\tname\tscore\tstrand\tgenerator\trank\tsource_region\tmotif_family\n"
                "chr1\t1\t5\ta\t2.0\t.\tmotif-free\t1\tchr1:1-5\tF1\n"
                "chr1\t6\t9\tb\t4.0\t.\tmotif-free\t2\tchr1:6-9\tF1\n"
                "chr1\t10\t12\tc\t1.0\t.\tmotif-free\t3\tchr1:10-12\tF2\n",
                encoding="utf-8",
            )
            frame = read_table(sites)
            ranked = rerank_sites(sites, out, score_columns=["score"], top_per_family=1)

        self.assertEqual(frame.columns[0], "chrom")
        self.assertEqual(ranked["name"].tolist(), ["b", "c"])
        self.assertIn("input_rank", ranked.columns)


if __name__ == "__main__":
    unittest.main()
