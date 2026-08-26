import argparse
import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from fp_tools.utils import bigwig as pyBigWig

try:
    import pysam
except ImportError:
    pysam = None

from fp_tools.tools.prepare_atac import (
    DEFAULTS,
    PROFILE_DEFAULTS,
    ReferenceBundle,
    _bedgraph_to_bigwig,
    _fragment_metrics,
    _relative_link,
    _samtools_sort_memory,
    _tool_version,
    _tss_enrichment,
    _write_chromosome_subset,
    _write_project_metadata,
    build_parser,
    download_file,
    load_settings,
    prepare_reference,
    read_preprocess_metadata,
    resolve_ena_fastqs,
    write_default_config,
)
from fp_tools.tools.prepare_atac_legacy import (
    _bowtie_command,
    _ensure_bigwig,
    _remove_xs,
    _tag_directory_command,
    _trim_command,
)


class PrepareAtacMetadataTest(unittest.TestCase):
    def test_tool_version_replaces_invalid_utf8(self):
        with tempfile.TemporaryDirectory() as tmp:
            tool = Path(tmp) / "invalid-version"
            tool.write_bytes(b"#!/bin/sh\nprintf '\\377tool 1.0\\n'\n")
            tool.chmod(0o755)
            self.assertEqual(_tool_version(str(tool)), "�tool 1.0")

    def test_reads_existing_gse_style_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "samples.txt"
            table.write_text(
                "ID\tAntibody\tCell\tSample\n"
                "SRR1\tATAC\tNIH3T3\tDox_BATF_1\n"
                "SRR2\tATAC\tNIH3T3\tDox_BATF_2\n",
                encoding="utf-8",
            )
            samples = read_preprocess_metadata(table)
            self.assertEqual(
                [row.sample for row in samples], ["Dox_BATF_1", "Dox_BATF_2"]
            )
            self.assertEqual(samples[0].condition, "Dox_BATF_1")
            self.assertEqual(samples[0].runs[0].accession, "SRR1")

    def test_groups_technical_runs_and_accepts_csv_fastqs(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "samples.csv"
            table.write_text(
                "run_accession,sample,condition,replicate,fastq_1,fastq_2\n"
                "lane1,A,treated,2,/x/a1.fq.gz,/x/a2.fq.gz\n"
                "lane2,A,treated,2,/x/b1.fq.gz,/x/b2.fq.gz\n",
                encoding="utf-8",
            )
            samples = read_preprocess_metadata(table)
            self.assertEqual(len(samples), 1)
            self.assertEqual(len(samples[0].runs), 2)
            self.assertEqual(samples[0].replicate, "2")

    def test_rejects_duplicate_accessions_and_conflicting_groups(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "samples.tsv"
            table.write_text(
                "ID\tSample\tCondition\nSRR1\tA\tx\nSRR1\tB\ty\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "Duplicate run accession"):
                read_preprocess_metadata(table)
            table.write_text(
                "ID\tSample\tCondition\nSRR1\tA\tx\nSRR2\tA\ty\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "disagree"):
                read_preprocess_metadata(table)


class PrepareAtacConfigTest(unittest.TestCase):
    def test_defaults_config_and_override_precedence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "defaults.yml"
            write_default_config(path)
            loaded = load_settings(
                path, {"resources": {"cores": 3}, "peaks": {"qvalue": 0.05}}
            )
            self.assertEqual(loaded["resources"]["cores"], 3)
            self.assertEqual(loaded["peaks"]["qvalue"], 0.05)
            self.assertEqual(loaded["align"]["mapq"], DEFAULTS["align"]["mapq"])

    def test_unknown_config_section_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.yml"
            path.write_text("mystery: true\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unknown"):
                load_settings(path)

    def test_homer_profile_has_expected_atac_defaults_and_all_cores(self):
        settings = load_settings(overrides={"profile": "homer-atac"})
        self.assertEqual(settings["profile"], "homer-atac")
        self.assertEqual(settings["align"]["max_insert"], 1000)
        self.assertEqual(settings["align"]["mapq"], 0)
        self.assertIn("--very-sensitive-local", settings["align"]["extra_args"])
        self.assertFalse(settings["filter"]["remove_mito"])
        self.assertEqual(settings["peaks"]["homer_local_fold"], 15)
        self.assertEqual(
            settings["resources"]["cores"],
            PROFILE_DEFAULTS["homer-atac"]["resources"]["cores"],
        )

    def test_sort_memory_is_bounded_by_total_run_budget(self):
        settings = load_settings(
            overrides={"resources": {"cores": 32, "memory_gb": 24}}
        )
        self.assertEqual(_samtools_sort_memory(settings), "384M")
        settings = load_settings(overrides={"resources": {"cores": 8, "memory_gb": 8}})
        self.assertEqual(_samtools_sort_memory(settings), "512M")

    def test_yaml_profile_and_cli_override_precedence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "legacy.yml"
            path.write_text(
                "profile: legacy-atac\nalign:\n  max_insert: 750\n", encoding="utf-8"
            )
            settings = load_settings(path, {"align": {"max_insert": 900}})
            self.assertEqual(settings["profile"], "homer-atac")
            self.assertEqual(settings["align"]["max_insert"], 900)

    def test_writes_profile_specific_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "homer.yml"
            write_default_config(path, "legacy-atac")
            text = path.read_text(encoding="utf-8")
            self.assertIn("profile: homer-atac", text)
            self.assertIn("sample_memory_gb: 16", text)
            self.assertNotIn("legacy", text.lower())
            self.assertIn("homer_local_size: 150000", text)

    def test_old_memory_setting_is_normalized(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "old.yml"
            path.write_text(
                "profile: legacy-atac\nresources:\n  legacy_sample_memory_gb: 20\n",
                encoding="utf-8",
            )
            settings = load_settings(path)
            self.assertEqual(settings["profile"], "homer-atac")
            self.assertEqual(settings["resources"]["sample_memory_gb"], 20)
            self.assertNotIn("legacy_sample_memory_gb", settings["resources"])

    def test_parser_has_simple_run_and_management_options(self):
        parser = build_parser()
        args = parser.parse_args(
            ["--samples", "s.tsv", "--genome", "mm10", "--outdir", "out", "--dry-run"]
        )
        self.assertTrue(args.dry_run)
        self.assertEqual(args.genome, "mm10")
        homer = parser.parse_args(["--profile", "homer-atac", "--doctor"])
        self.assertEqual(homer.profile, "homer-atac")
        compatibility = parser.parse_args(["--profile", "legacy-atac", "--doctor"])
        self.assertEqual(compatibility.profile, "homer-atac")
        help_text = parser.format_help().lower()
        self.assertIn("{modern,homer-atac}", help_text)
        self.assertNotIn("legacy", help_text)
        self.assertIsInstance(parser, argparse.ArgumentParser)


class PrepareAtacDownloadTest(unittest.TestCase):
    @mock.patch("urllib.request.urlopen")
    def test_resolves_ena_paired_fastqs_and_checksums(self, urlopen):
        text = (
            "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\tlibrary_layout\n"
            "SRR1\tftp.sra.ebi.ac.uk/a_1.fastq.gz;ftp.sra.ebi.ac.uk/a_2.fastq.gz\tmd51;md52\t10;11\tPAIRED\n"
        )
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = text.encode()
        urlopen.return_value = response
        result = resolve_ena_fastqs("SRR1")
        self.assertEqual(len(result), 2)
        self.assertEqual(
            result[0], ("https://ftp.sra.ebi.ac.uk/a_1.fastq.gz", "md51", 10)
        )

    def test_existing_download_is_reused_only_when_md5_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "reads.fastq.gz"
            output.write_bytes(b"reads")
            digest = hashlib.md5(b"reads").hexdigest()
            with (
                mock.patch("shutil.which", return_value=None),
                mock.patch("urllib.request.urlopen") as urlopen,
            ):
                result = download_file("https://example.invalid/reads", output, digest)
            self.assertEqual(result, output)
            urlopen.assert_not_called()

    def test_existing_download_is_reused_when_size_matches_without_md5(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "reads.fastq.gz"
            output.write_bytes(b"reads")
            with (
                mock.patch("shutil.which", return_value=None),
                mock.patch("urllib.request.urlopen") as urlopen,
            ):
                result = download_file(
                    "https://example.invalid/reads",
                    output,
                    expected_size=5,
                )
            self.assertEqual(result, output)
            urlopen.assert_not_called()

    def test_truncated_download_is_replaced_when_expected_size_is_known(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "reads.fastq.gz"
            output.write_bytes(b"bad")
            response = mock.MagicMock()
            response.__enter__.return_value.read.side_effect = [b"reads", b""]
            with (
                mock.patch("shutil.which", return_value=None),
                mock.patch("urllib.request.urlopen", return_value=response),
            ):
                result = download_file(
                    "https://example.invalid/reads",
                    output,
                    expected_size=5,
                )
            self.assertEqual(result.read_bytes(), b"reads")
            self.assertFalse(output.with_suffix(".gz.partial").exists())

    def test_unvalidated_existing_download_is_not_trusted(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "reads.fastq.gz"
            output.write_bytes(b"possibly-partial")
            response = mock.MagicMock()
            response.__enter__.return_value.read.side_effect = [b"complete", b""]
            with (
                mock.patch("shutil.which", return_value=None),
                mock.patch("urllib.request.urlopen", return_value=response) as urlopen,
            ):
                result = download_file("https://example.invalid/reads", output)
            self.assertEqual(result.read_bytes(), b"complete")
            urlopen.assert_called_once()

    def test_size_mismatch_removes_invalid_partial_and_preserves_no_final(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "reads.fastq.gz"
            response = mock.MagicMock()
            response.__enter__.return_value.read.side_effect = [b"short", b""]
            with (
                mock.patch("shutil.which", return_value=None),
                mock.patch("urllib.request.urlopen", return_value=response),
            ):
                with self.assertRaisesRegex(RuntimeError, "Size mismatch"):
                    download_file(
                        "https://example.invalid/reads",
                        output,
                        expected_size=10,
                    )
            self.assertFalse(output.exists())
            self.assertFalse(output.with_suffix(".gz.partial").exists())


class PrepareAtacOutputTest(unittest.TestCase):
    def test_legacy_commands_adapt_to_single_end_without_paired_flags(self):
        settings = load_settings(overrides={"profile": "legacy-atac"})
        trim = _trim_command("SRR1", Path("reads.fq.gz"), None, Path("work"), "32", [])
        align = _bowtie_command(
            Path("index/mm10"), Path("trimmed.fq.gz"), None, "32", settings
        )
        tags = _tag_directory_command(Path("tags"), Path("reads.bam"), False)
        self.assertNotIn("--paired", trim)
        self.assertIn("-U", align)
        self.assertIn("--very-sensitive-local", align)
        for flag in ("--no-mixed", "--no-discordant", "--dovetail", "-X"):
            self.assertNotIn(flag, align)
        self.assertNotIn("-sspe", tags)

    def test_legacy_paired_commands_retain_historical_flags(self):
        settings = load_settings(overrides={"profile": "legacy-atac"})
        trim = _trim_command(
            "SRR1",
            Path("r1.fq.gz"),
            Path("r2.fq.gz"),
            Path("work"),
            "32",
            [],
        )
        align = _bowtie_command(
            Path("index/mm10"),
            Path("r1.fq.gz"),
            Path("r2.fq.gz"),
            "32",
            settings,
        )
        tags = _tag_directory_command(Path("tags"), Path("reads.bam"), True)
        self.assertIn("--paired", trim)
        for flag in ("-1", "-2", "--no-mixed", "--dovetail", "-X", "--interleaved"):
            self.assertIn(flag, align)
        self.assertIn("-sspe", tags)

    def test_chromosome_filter_indexes_source_bam_before_region_view(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.bam"
            output = root / "subset.bam"
            log = root / "commands.log"
            with mock.patch("fp_tools.tools.prepare_atac._run") as run:
                _write_chromosome_subset(source, output, ["chr1", "chr2"], "4", log)
            self.assertEqual(run.call_count, 2)
            self.assertEqual(
                run.call_args_list[0].args[0],
                ["samtools", "index", "-@", "4", str(source)],
            )
            self.assertEqual(run.call_args_list[1].args[0][-2:], ["chr1", "chr2"])

    def test_legacy_homer_bedgraph_fallback_becomes_bigwig(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sizes = root / "sizes"
            sizes.write_text("chr1\t20\n", encoding="utf-8")
            output = root / "signal.bw"
            output.write_text("chr1\t2\t5\t3.5\n", encoding="utf-8")

            def fake_sort(command, check, stdout):
                stdout.write(Path(command[-1]).read_text(encoding="utf-8"))
                return mock.MagicMock(returncode=0)

            with mock.patch(
                "fp_tools.tools.prepare_atac_legacy.subprocess.run",
                side_effect=fake_sort,
            ):
                _ensure_bigwig(output, sizes)
            with pyBigWig.open(str(output)) as bw:
                self.assertEqual(bw.values("chr1", 2, 5), [3.5, 3.5, 3.5])

    @unittest.skipUnless(pysam is not None, "pysam is required to write BAM fixtures")
    def test_legacy_xs_filter_retains_only_unique_alignments(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.bam"
            output = root / "unique.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
            with pysam.AlignmentFile(source, "wb", header=header) as bam:
                for idx in range(2):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{idx}"
                    read.query_sequence = "A" * 20
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = idx
                    read.mapping_quality = 30
                    read.cigar = [(0, 20)]
                    read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                    if idx:
                        read.set_tag("XS", 10)
                    bam.write(read)
            self.assertEqual(_remove_xs(source, output), 1)
            with pysam.AlignmentFile(output, "rb") as bam:
                self.assertEqual(sum(1 for _ in bam.fetch(until_eof=True)), 1)

    @unittest.skipUnless(pysam is not None, "pysam is required to write BAM fixtures")
    def test_tss_enrichment_uses_shifted_cut_sites(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam_path = root / "reads.bam"
            header = {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "chr1", "LN": 5000}],
            }
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                starts = [546] + [2496] * 10 + [4446]
                for idx, start in enumerate(starts):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{idx}"
                    read.query_sequence = "A" * 30
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 60
                    read.cigar = [(0, 30)]
                    read.query_qualities = pysam.qualitystring_to_array("I" * 30)
                    bam.write(read)
            pysam.index(str(bam_path))
            tss = root / "tss.bed"
            tss.write_text("chr1\t2500\t2501\tgene\t0\t+\n", encoding="utf-8")
            score = _tss_enrichment(bam_path, tss)
            self.assertIsNotNone(score)
            self.assertGreater(score, 100)
            fragment_metrics = _fragment_metrics(
                bam_path, root / "fragment_lengths.tsv"
            )
            self.assertEqual(fragment_metrics["median_fragment_length"], 30)
            self.assertEqual(fragment_metrics["nucleosome_free_fraction"], 1.0)

    def test_custom_reference_dry_run_does_not_require_download_or_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "custom.fa"
            fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
            bundle = prepare_reference(
                "custom1",
                Path(tmp) / "refs",
                fasta=fasta,
                bowtie2_index=Path(tmp) / "idx" / "custom1",
                macs_genome_size="1000",
                dry_run=True,
            )
            self.assertEqual(bundle.assembly, "custom1")
            self.assertEqual(bundle.macs_genome_size, "1000")

    def test_bedgraph_conversion_and_relative_compatibility_link(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sizes = root / "sizes"
            sizes.write_text("chr1\t10\n", encoding="utf-8")
            bedgraph = root / "signal.bedgraph"
            bedgraph.write_text("chr1\t1\t3\t5\n", encoding="utf-8")
            bigwig = root / "signal.bw"
            _bedgraph_to_bigwig(bedgraph, sizes, bigwig)
            with pyBigWig.open(str(bigwig)) as bw:
                self.assertEqual(bw.values("chr1", 1, 3), [5.0, 5.0])
            link = root / "legacy" / "sample.rp10m.bw"
            _relative_link(bigwig, link)
            self.assertEqual(link.resolve(), bigwig.resolve())

    def test_writes_downstream_sample_table_and_runnable_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.tsv"
            source.write_text("ID\nSRR1\n", encoding="utf-8")
            ref = ReferenceBundle(
                "mm10", root / "mm10.fa", root / "idx", None, root / "sizes", None, "mm"
            )
            result = {
                "sample": "S1",
                "condition": "C1",
                "bam": "/b.bam",
                "peaks": "/p.bed",
                "bigwig": "/s.bw",
            }
            _write_project_metadata([result], source, root, ref, load_settings())
            table = (root / "metadata" / "samples.tsv").read_text(encoding="utf-8")
            self.assertIn("S1\tC1\t/b.bam\t/p.bed\t/s.bw", table)
            config = (root / "metadata" / "atac_correct.yml").read_text(
                encoding="utf-8"
            )
            self.assertIn("tool: atac-correct", config)


if __name__ == "__main__":
    unittest.main()
