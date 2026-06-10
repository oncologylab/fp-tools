import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.motif_discovery import (
    export_candidate_fasta,
    meme_command,
    parse_meme_txt,
    parse_tomtom_tsv,
    read_candidate_sites,
    summarize_motif_outputs,
    write_motif_discovery_plan,
    write_motif_summary_html,
    write_motif_summary_tsv,
)


class MotifDiscoveryPrepTest(unittest.TestCase):
    def test_export_candidate_fasta_centered_on_sites(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            genome = tmp / "genome.fa"
            candidates = tmp / "candidates.bed"
            out = tmp / "candidates.fa"
            genome.write_text(">chr1\nAAAACCCCGGGGTTTT\n", encoding="utf-8")
            candidates.write_text("chr1\t4\t8\tcand1\t2.5\nchr2\t0\t2\tmissing\t1\n", encoding="utf-8")

            written = export_candidate_fasta(candidates, genome, out, flank=3)
            text = out.read_text(encoding="utf-8")

        self.assertEqual(written, 1)
        self.assertIn(">cand1|chr1:3-9|score=2.5", text)
        self.assertIn("ACCCCG", text)

    def test_read_candidate_sites_and_meme_command(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = Path(tmpdir) / "sites.bed"
            bed.write_text("#header\nchr1\t1\t5\tsiteA\t4.2\n", encoding="utf-8")
            sites = read_candidate_sites(bed)
        self.assertEqual(sites[0].name, "siteA")
        self.assertEqual(meme_command("sites.fa", "motifs", method="dreme", extra_args=["-dna"]), ["dreme", "sites.fa", "-oc", "motifs", "-dna"])

    def test_write_motif_discovery_plan_includes_tomtom_and_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "sites.fa"
            known = tmp / "known.meme"
            script = tmp / "run.sh"
            fasta.write_text(">site1\nACGT\n", encoding="utf-8")
            known.write_text("MEME version 4\n", encoding="utf-8")

            path = write_motif_discovery_plan(
                fasta,
                tmp / "motifs",
                script,
                method="meme",
                known_motifs=known,
                extra_args=["-dna", "-nmotifs", "5"],
            )
            text = path.read_text(encoding="utf-8")

        self.assertIn("meme", text)
        self.assertIn("-nmotifs 5", text)
        self.assertIn("tomtom", text)
        self.assertIn("fp-tools-summarize-motifs", text)
        self.assertTrue(text.startswith("#!/usr/bin/env bash"))

    def test_parse_meme_tomtom_and_write_reports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            meme = tmp / "meme.txt"
            tomtom = tmp / "tomtom.tsv"
            out_tsv = tmp / "summary.tsv"
            out_html = tmp / "summary.html"
            meme.write_text(
                "MEME version 5\n"
                "MOTIF motif_1 CTCF_like\n"
                "letter-probability matrix: alength= 4 w= 4 nsites= 12 E= 1.2e-05\n"
                "0.90 0.05 0.03 0.02\n"
                "0.01 0.92 0.04 0.03\n"
                "0.02 0.03 0.91 0.04\n"
                "0.03 0.02 0.05 0.90\n",
                encoding="utf-8",
            )
            tomtom.write_text(
                "Query_ID\tTarget_ID\tOptimal_offset\tp-value\tE-value\tq-value\tOverlap\tQuery_consensus\tTarget_consensus\tOrientation\n"
                "motif_1\tMA0139.1_CTCF\t0\t1e-6\t1e-4\t1e-3\t4\tACGT\tACGT\t+\n",
                encoding="utf-8",
            )

            meme_rows = parse_meme_txt(meme)
            tomtom_rows = parse_tomtom_tsv(tomtom)
            rows = summarize_motif_outputs(meme, tomtom)
            write_motif_summary_tsv(rows, out_tsv)
            write_motif_summary_html(rows, out_html, title="Motif Report")

            self.assertEqual(meme_rows[0]["motif_id"], "motif_1")
            self.assertEqual(meme_rows[0]["consensus"], "ACGT")
            self.assertEqual(tomtom_rows[0]["target_id"], "MA0139.1_CTCF")
            self.assertEqual(len(rows), 2)
            self.assertIn("motif_1", out_tsv.read_text(encoding="utf-8"))
            html = out_html.read_text(encoding="utf-8")
            self.assertIn("Motif Report", html)
            self.assertIn("<svg", html)
            self.assertIn("ACGT", html)


if __name__ == "__main__":
    unittest.main()
