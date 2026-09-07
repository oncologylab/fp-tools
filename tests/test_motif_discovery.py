import contextlib
import io
import shlex
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from fp_tools.tools import motif_discovery

from fp_tools.tools.motif_discovery import (
    discovery_motif_text,
    export_candidate_fasta,
    meme_command,
    motif_discovery_plan_main,
    parse_meme_txt,
    parse_tomtom_tsv,
    prepare_known_motifs_for_tomtom,
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

    def test_streme_command_and_output_path(self):
        self.assertEqual(
            meme_command("sites.fa", "motifs", method="streme", extra_args=["--dna", "--nmotifs", "8"]),
            ["streme", "--p", "sites.fa", "--oc", "motifs", "--dna", "--nmotifs", "8"],
        )
        self.assertEqual(discovery_motif_text("streme", "motifs"), Path("motifs/streme.txt"))

    def test_write_motif_discovery_plan_includes_tomtom_and_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "sites.fa"
            known = tmp / "known.meme"
            script = tmp / "run.sh"
            fasta.write_text(">site1\nACGT\n", encoding="utf-8")
            known.write_text(
                "MEME version 4\n\n"
                "ALPHABET= ACGT\n\n"
                "strands: + -\n\n"
                "Background letter frequencies\n"
                "A 0.25 C 0.25 G 0.25 T 0.25\n\n"
                "MOTIF M1 known\n"
                "letter-probability matrix: alength= 4 w= 2 nsites= 10 E= 0\n"
                "0.7 0.1 0.1 0.1\n"
                "0.1 0.1 0.1 0.7\n",
                encoding="utf-8",
            )

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
        self.assertIn("summarize-motifs", text)
        self.assertTrue(text.startswith("#!/usr/bin/env bash"))

    def test_jaspar_inputs_are_converted_to_valid_ordered_meme(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            jaspar = root / "known.jaspar"
            jaspar.write_text(
                ">M1 First motif name\n"
                "A [ 8 1 ]\nC [ 1 1 ]\nG [ 1 1 ]\nT [ 0 7 ]\n"
                ">M2 Second motif\n"
                "A [ 0 5 ]\nC [ 5 0 ]\nG [ 0 0 ]\nT [ 0 0 ]\n",
                encoding="utf-8",
            )
            outputs = prepare_known_motifs_for_tomtom(
                [jaspar], root / "results"
            )
            converted = outputs[0]
            motif_format, motifs = motif_discovery._validated_motif_list(converted)

            self.assertEqual(
                converted, root / "results" / "known_motifs" / "known.meme"
            )
            self.assertEqual(motif_format, "meme")
            self.assertEqual(
                [(motif.id, motif.name) for motif in motifs],
                [("M1", "First motif name"), ("M2", "Second motif")],
            )
            self.assertEqual(motifs[0].counts, [[8, 1], [1, 1], [1, 1], [0, 7]])

    def test_valid_meme_input_is_passed_through_without_changes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "known.meme"
            source.write_text(
                "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
                "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
                "MOTIF M1 Known\n"
                "letter-probability matrix: alength= 4 w= 1 nsites= 10 E= 0\n"
                "0.25 0.25 0.25 0.25\n",
                encoding="utf-8",
            )
            before = source.read_bytes()
            outputs = prepare_known_motifs_for_tomtom(source, root / "results")
            self.assertEqual(outputs, [source])
            self.assertEqual(source.read_bytes(), before)
            self.assertFalse((root / "results" / "known_motifs").exists())

    def test_motif_discovery_plan_accepts_builtin_known_motif_db(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "sites.fa"
            outdir = tmp / "motifs"
            fasta.write_text(">site1\nACGT\n", encoding="utf-8")

            code = motif_discovery_plan_main(
                [
                    "--fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                    "--known-motif-db",
                    "jaspar2026",
                ]
            )
            script = outdir / "run_motif_discovery.sh"
            text = script.read_text(encoding="utf-8")
            converted = outdir / "known_motifs" / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.meme"
            self.assertTrue(converted.is_file())
            converted_count = len(
                motif_discovery._validated_motif_list(converted)[1]
            )

        self.assertEqual(code, 0)
        self.assertIn("tomtom", text)
        self.assertIn("known_motifs", text)
        self.assertIn("JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.meme", text)
        self.assertEqual(converted_count, 1019)

    def test_summary_command_uses_frozen_dispatch_and_quotes_every_token(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "directory with spaces"
            root.mkdir()
            fasta = root / "sites file.fa"
            fasta.write_text(">site1\nACGT\n", encoding="utf-8")
            with mock.patch(
                "fp_tools.utils.subprocess_commands.is_frozen", return_value=True
            ), mock.patch.object(sys, "executable", "/Applications/fp tools"):
                script = write_motif_discovery_plan(
                    fasta, root / "motif results", root / "run plan.sh"
                )
            summary_tokens = shlex.split(script.read_text(encoding="utf-8").splitlines()[-1])

        self.assertEqual(
            summary_tokens[:3],
            [
                "/Applications/fp tools",
                "--fp-tools-internal-command",
                "summarize-motifs",
            ],
        )
        self.assertIn(str(root / "motif results" / "meme" / "meme.txt"), summary_tokens)

    def test_execute_failure_returns_actionable_nonzero_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "sites.fa"
            fasta.write_text(">site1\nACGT\n", encoding="utf-8")
            stderr = io.StringIO()
            with mock.patch.object(
                motif_discovery, "prepare_command_runtime", return_value=None
            ), mock.patch.object(
                motif_discovery.subprocess,
                "run",
                return_value=SimpleNamespace(returncode=7),
            ), contextlib.redirect_stderr(stderr):
                with self.assertRaises(SystemExit) as raised:
                    motif_discovery_plan_main(
                        [
                            "--fasta",
                            str(fasta),
                            "--outdir",
                            str(root / "results"),
                            "--execute",
                            "--runtime",
                            "system",
                        ]
                    )
            self.assertEqual(raised.exception.code, 7)
            self.assertIn("generated workflow failed with exit code 7", stderr.getvalue())
            self.assertIn("run_motif_discovery.sh", stderr.getvalue())

    def test_frozen_execute_resets_pyinstaller_environment_for_shell_child(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "sites.fa"
            fasta.write_text(">site1\nACGT\n", encoding="utf-8")
            with mock.patch.object(
                motif_discovery, "prepare_command_runtime", return_value=None
            ), mock.patch.object(
                motif_discovery.subprocess,
                "run",
                return_value=SimpleNamespace(returncode=0),
            ) as runner, mock.patch.object(sys, "frozen", True, create=True):
                code = motif_discovery_plan_main(
                    [
                        "--fasta",
                        str(fasta),
                        "--outdir",
                        str(root / "results"),
                        "--execute",
                        "--runtime",
                        "system",
                    ]
                )

            self.assertEqual(code, 0)
            environment = runner.call_args.kwargs["env"]
            self.assertEqual(environment["PYINSTALLER_RESET_ENVIRONMENT"], "1")

    def test_write_streme_plan_uses_streme_txt(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "sites.fa"
            script = tmp / "run.sh"
            fasta.write_text(">site1\nACGT\n", encoding="utf-8")

            path = write_motif_discovery_plan(
                fasta,
                tmp / "motifs",
                script,
                method="streme",
                extra_args=["--dna", "--nmotifs", "5"],
            )
            text = path.read_text(encoding="utf-8")

        self.assertIn("streme --p", text)
        self.assertIn("streme/streme.txt", text.replace("\\", "/"))
        self.assertIn("--nmotifs 5", text)

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
