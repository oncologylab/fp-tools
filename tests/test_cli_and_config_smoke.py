import contextlib
import io
import argparse
import pathlib
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from fp_tools.cli_batch import run_config_file
from fp_tools.gui_config import (
    build_cli_command,
    expand_jobs,
    load_yaml_config,
    normalize_config,
    validate_config,
    validate_gui_config,
)
from fp_tools.parsers import add_aggregate_arguments, add_atacorrect_arguments, add_diff_footprints_arguments, add_scorebigwig_arguments
from fp_tools.tools import diff_footprints
from fp_tools.tools.diff_footprints import _prepare_condition_metadata


ROOT = pathlib.Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "examples" / "gui_configs"


class CliAndConfigSmokeTest(unittest.TestCase):
    def test_bulk_config_accepts_exactly_one_input_level(self):
        base = {
            "tool": "bulk-footprinting",
            "comparison_table": "comparisons.tsv",
            "genome": "hg38",
            "outdir": "results",
        }
        raw = {"samples": [{**base, "reads_table": "reads.tsv"}], "comparisons": []}
        aligned = {"samples": [{**base, "sample_table": "samples.tsv"}], "comparisons": []}
        both = {
            "samples": [{**base, "reads_table": "reads.tsv", "sample_table": "samples.tsv"}],
            "comparisons": [],
        }
        self.assertEqual(validate_config(raw), [])
        self.assertEqual(validate_config(aligned), [])
        self.assertTrue(any("exactly one" in value for value in validate_config(both)))

    def test_gui_accepts_bam_bulk_config_and_rejects_raw_reads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bam = root / "sample.bam"
            bai = root / "sample.bam.bai"
            peaks = root / "peaks.bed"
            genome = root / "genome.fa"
            comparisons = root / "comparisons.tsv"
            samples = root / "samples.tsv"
            for path in (bam, bai, peaks, genome):
                path.write_bytes(b"fixture")
            comparisons.write_text("comparison\tcond1\tcond2\n", encoding="utf-8")
            samples.write_text(
                f"sample\tcondition\tbam\tpeaks\nsample\tA\t{bam}\t{peaks}\n",
                encoding="utf-8",
            )
            base = {
                "tool": "bulk-footprinting",
                "comparison_table": str(comparisons),
                "genome": str(genome),
                "outdir": str(root / "results"),
            }
            aligned = {"samples": [{**base, "sample_table": str(samples)}], "comparisons": []}
            raw = {"samples": [{**base, "reads_table": "reads.tsv"}], "comparisons": []}
            raw_extra = {
                "samples": [
                    {**base, "sample_table": str(samples), "extra_args": ["--reads-table=reads.tsv"]}
                ],
                "comparisons": [],
            }
            self.assertEqual(validate_gui_config(aligned), [])
            self.assertTrue(any("GUI bulk workflows" in error for error in validate_gui_config(raw)))
            self.assertTrue(any("GUI bulk workflows" in error for error in validate_gui_config(raw_extra)))

    def test_gui_rejects_prepare_atac_even_on_linux(self):
        config = {
            "samples": [
                {
                    "tool": "prepare-atac",
                    "samples": "reads.tsv",
                    "genome": "hg38",
                    "outdir": "project",
                }
            ],
            "comparisons": [],
        }
        self.assertTrue(any("GUI starts from BAM/BAI" in error for error in validate_gui_config(config)))

    def test_gui_rejects_invalid_enum_values_loaded_from_yaml(self):
        config = {
            "samples": [
                {
                    "tool": "bulk-footprinting",
                    "sample_table": "missing.tsv",
                    "comparison_table": "missing-comparisons.tsv",
                    "genome": "missing.fa",
                    "outdir": "results",
                    "review_format": "not-a-format",
                    "plot_aggregate": "not-a-mode",
                    "normalization": "not-a-normalization",
                }
            ],
            "comparisons": [],
        }
        errors = validate_gui_config(config)
        self.assertTrue(any("unsupported 'review_format'" in error for error in errors))
        self.assertTrue(any("unsupported 'plot_aggregate'" in error for error in errors))
        self.assertTrue(any("unsupported 'normalization'" in error for error in errors))

    def test_gui_requires_motif_summary_input(self):
        config = {
            "samples": [
                {
                    "tool": "summarize-motifs",
                    "meme_txt": "",
                    "tomtom_tsv": "",
                    "out_tsv": "summary.tsv",
                }
            ],
            "comparisons": [],
        }
        self.assertTrue(any("provide 'meme_txt' or 'tomtom_tsv'" in error for error in validate_gui_config(config)))

    def test_region_comparison_yaml_does_not_require_peaks(self):
        config = {
            "version": 1,
            "run_mode": "single",
            "samples": [{
                "sample_id": "regions",
                "tool": "diff-footprints",
                "comparison_axis": "regions",
                "signals": ["rep1.bw", "rep2.bw"],
                "sample_names": ["rep1", "rep2"],
                "regions": ["bound.bed", "control.bed"],
                "region_labels": ["bound", "control"],
                "genome": "hg38.fa",
            }],
            "comparisons": [],
        }
        self.assertEqual(validate_config(config), [])
        job = expand_jobs(config)[0]
        self.assertIn("--comparison-axis", job.command)
        self.assertIn("--regions", job.command)
        self.assertNotIn("--peaks", job.command)

    def test_review_multi_comparisons_yaml_supports_standalone_output(self):
        item = {
            "tool": "review-multi-comparisons",
            "inputs": ["comparison_a.html", "comparison_b.html"],
            "labels": ["A", "B"],
            "output_html": "review.html",
        }
        config = {"samples": [], "comparisons": [item]}
        self.assertEqual(validate_config(config), [])
        command = build_cli_command(item["tool"], item)
        self.assertIn("--output-html", command)
        self.assertNotIn("--output-dir", command)

    def test_review_multi_comparisons_yaml_rejects_two_output_modes(self):
        item = {
            "tool": "review-multi-comparisons",
            "inputs": ["comparison.html"],
            "output_dir": "review",
            "output_html": "review.html",
        }
        errors = validate_config({"samples": [], "comparisons": [item]})
        self.assertTrue(any("mutually exclusive" in error for error in errors))

    def test_all_example_yaml_configs_expand_to_jobs(self):
        config_paths = sorted(CONFIG_DIR.glob("*.yml"))
        self.assertGreaterEqual(len(config_paths), 1)

        for path in config_paths:
            with self.subTest(config=path.name):
                config = normalize_config(load_yaml_config(path))
                jobs = expand_jobs(config)
                self.assertGreaterEqual(len(jobs), 1)
                for job in jobs:
                    self.assertTrue(job.job_id)
                    self.assertIn(
                        job.tool,
                        {
                            "atac-correct",
                            "call-footprints",
                            "diff-footprints",
                            "normalize-bigwig",
                            "plot-aggregate",
                            "review-multi-comparisons",
                            "bulk-footprinting",
                            "discover-motifs",
                            "summarize-motifs",
                            "pseudobulk-fragments",
                            "find-signature-fp",
                            "sc-footprinting",
                        },
                    )
                    self.assertEqual(job.command[0], job.tool)

    def test_all_example_yaml_configs_support_dry_run(self):
        for path in sorted(CONFIG_DIR.glob("*.yml")):
            with self.subTest(config=path.name):
                stdout = io.StringIO()
                with contextlib.redirect_stdout(stdout):
                    code = run_config_file(path, dry_run=True)
                self.assertEqual(code, 0)
                self.assertIn("[", stdout.getvalue())

    def test_packaged_entry_points_print_help(self):
        commands = [
            "atac-correct",
            "call-footprints",
            "match-motifs",
            "diff-footprints",
            "normalize-bigwig",
            "plot-aggregate",
            "review-multi-comparisons",
            "bulk-footprinting",
            "fp-tools-runtime",
            "run-yaml-workflow",
            "discover-motifs",
            "summarize-motifs",
            "pseudobulk-fragments",
            "find-signature-fp",
            "sc-footprinting",
        ]
        for command in commands:
            exe = ROOT / ".venv" / "bin" / command
            if not exe.exists():
                self.skipTest(f"{exe} is not available in this checkout")
            with self.subTest(command=command):
                result = subprocess.run(
                    [str(exe), "--help"],
                    cwd=ROOT,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=30,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout.lower())

    def test_call_footprints_accepts_batch_signal_arguments(self):
        parser = add_scorebigwig_arguments(argparse.ArgumentParser())
        args = parser.parse_args(
            [
                "--signals",
                "A_corrected.bw",
                "B_corrected.bw",
                "--outputs",
                "A_footprints.bw",
                "B_footprints.bw",
                "--regions",
                "peaks.bed",
                "--output-beds",
                "A_candidates.bed",
                "B_candidates.bed",
            ]
        )
        self.assertEqual(args.signals, ["A_corrected.bw", "B_corrected.bw"])
        self.assertEqual(args.outputs, ["A_footprints.bw", "B_footprints.bw"])
        self.assertEqual(args.output_beds, ["A_candidates.bed", "B_candidates.bed"])

    def test_call_footprints_uses_public_reference_kernel_name(self):
        parser = add_scorebigwig_arguments(argparse.ArgumentParser())
        reference = parser.parse_args(["--footprint-kernel", "reference"])
        compatibility = parser.parse_args(["--footprint-kernel", "legacy"])
        self.assertEqual(reference.footprint_kernel, "reference")
        self.assertEqual(compatibility.footprint_kernel, "reference")
        help_text = parser.format_help().lower()
        self.assertIn("{fast,reference}", help_text)
        self.assertNotIn("legacy", help_text)

    def test_atac_correct_accepts_batch_bam_and_peak_arguments(self):
        parser = add_atacorrect_arguments(argparse.ArgumentParser())
        args = parser.parse_args(
            [
                "--bams",
                "A.bam",
                "B.bam",
                "--sample-names",
                "A",
                "B",
                "--genome",
                "genome.fa",
                "--peaks",
                "A_peaks.bed",
                "B_peaks.bed",
            ]
        )
        self.assertEqual(args.bams, ["A.bam", "B.bam"])
        self.assertEqual(args.sample_names, ["A", "B"])
        self.assertEqual(args.peaks, ["A_peaks.bed", "B_peaks.bed"])

    def test_plot_aggregate_accepts_match_dir_html_arguments(self):
        parser = add_aggregate_arguments(argparse.ArgumentParser())
        args = parser.parse_args(
            [
                "--match-dir",
                "results/motif_matches/sample",
                "--signals",
                "sample_corrected.bw",
                "--motifs",
                "SPIB",
                "CEBPB",
                "--format",
                "html",
                "--output",
                "aggregate.html",
            ]
        )
        self.assertEqual(args.match_dir, ["results/motif_matches/sample"])
        self.assertEqual(args.motifs, ["SPIB", "CEBPB"])
        self.assertEqual(args.format, "html")
        self.assertNotIn("legacy", parser.format_help().lower())

    def test_diff_footprints_accepts_sample_names_separate_from_conditions(self):
        parser = add_diff_footprints_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args(
            [
                "--signals",
                "K562_rep1_footprints.bw",
                "K562_rep2_footprints.bw",
                "HepG2_rep1_footprints.bw",
                "HepG2_rep2_footprints.bw",
                "--sample-names",
                "K562_R1",
                "K562_R2",
                "HepG2_R1",
                "HepG2_R2",
                "--cond-names",
                "K562",
                "K562",
                "HepG2",
                "HepG2",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        prepared = _prepare_condition_metadata(args)
        self.assertEqual(prepared.cond_names, ["K562", "HepG2"])
        self.assertEqual(prepared.sample_names, ["K562_R1", "K562_R2", "HepG2_R1", "HepG2_R2"])
        self.assertEqual(prepared.condition_samples["K562"], ["K562_R1", "K562_R2"])
        self.assertEqual(prepared.condition_samples["HepG2"], ["HepG2_R1", "HepG2_R2"])
        self.assertEqual(prepared.sample_to_condition["K562_R2"], "K562")

    def test_diff_footprints_accepts_folder_inputs(self):
        parser = add_diff_footprints_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args(
            [
                "--sample-dirs",
                "results/samples/A",
                "results/samples/B",
                "--cond-names",
                "A",
                "B",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        self.assertEqual(args.sample_dirs, ["results/samples/A", "results/samples/B"])
        self.assertIsNone(args.project_dir)

    def test_diff_footprints_disambiguates_sample_condition_name_collisions(self):
        parser = add_diff_footprints_arguments(argparse.ArgumentParser(prog="diff-footprints"))
        args = parser.parse_args(
            [
                "--signals",
                "Bcell_footprints.bw",
                "Tcell_footprints.bw",
                "--sample-names",
                "Bcell",
                "Tcell",
                "--cond-names",
                "Bcell",
                "Tcell",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        prepared = _prepare_condition_metadata(args)
        self.assertEqual(prepared.sample_names, ["Bcell_sample", "Tcell_sample"])
        self.assertEqual(prepared.condition_samples["Bcell"], ["Bcell_sample"])
        self.assertEqual(prepared.condition_samples["Tcell"], ["Tcell_sample"])

    def test_cached_match_dirs_merge_to_tfbs_tmp_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            for sample, score in [("A", "1.25000"), ("B", "2.50000")]:
                bed_dir = root / sample / "match_motifs" / "TF1_MA0001.1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_MA0001.1_all.bed").write_text(
                    "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t" + score + "\n",
                    encoding="utf-8",
                )
            outdir = root / "out"
            (outdir / "TF1_MA0001.1" / "beds").mkdir(parents=True)
            args = SimpleNamespace(
                cached_match_dirs=[
                    str(root / "A" / "match_motifs"),
                    str(root / "B" / "match_motifs"),
                ],
                sample_names=["A_sample", "B_sample"],
                outdir=str(outdir),
                peak_header_list=["peak_chr", "peak_start", "peak_end"],
            )
            logger = diff_footprints.FpToolsLogger("test", 0)
            diff_footprints._write_cached_tfbs_tmp_files(args, ["TF1_MA0001.1"], logger)
            tmp_file = outdir / "TF1_MA0001.1" / "beds" / "TF1_MA0001.1.tmp"
            self.assertEqual(
                tmp_file.read_text(encoding="utf-8").strip(),
                "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t1.25000\t2.50000",
            )

    def test_cached_match_dirs_summary_mode_defers_tmp_files_to_workers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            for sample, score in [("A", "1.25000"), ("B", "2.50000")]:
                bed_dir = root / sample / "match_motifs" / "TF1_MA0001.1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_MA0001.1_all.bed").write_text(
                    "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t" + score + "\n",
                    encoding="utf-8",
                )
            outdir = root / "out"
            (outdir / "TF1_MA0001.1" / "beds").mkdir(parents=True)
            args = SimpleNamespace(
                cached_match_dirs=[
                    str(root / "A" / "match_motifs"),
                    str(root / "B" / "match_motifs"),
                ],
                sample_names=["A_sample", "B_sample"],
                outdir=str(outdir),
                peak_header_list=["peak_chr", "peak_start", "peak_end"],
                static_plots=False,
                write_motif_outputs=False,
            )
            logger = diff_footprints.FpToolsLogger("test", 0)
            diff_footprints._write_cached_tfbs_tmp_files(args, ["TF1_MA0001.1"], logger)
            tmp_file = outdir / "TF1_MA0001.1" / "beds" / "TF1_MA0001.1.tmp"
            self.assertFalse(tmp_file.exists())
            self.assertEqual(len(args.cached_motif_bed_maps), 2)
            self.assertIn("TF1_MA0001.1", args.cached_motif_bed_maps[0])

    def test_cached_match_dirs_use_compact_motif_site_cache_when_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            logger = diff_footprints.FpToolsLogger("test", 0)
            for sample, score in [("A", "1.25000"), ("B", "2.50000")]:
                match_dir = root / sample / "match_motifs"
                bed_dir = match_dir / "TF1_MA0001.1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_MA0001.1_all.bed").write_text(
                    "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t" + score + "\t" + score + "\tNA\n",
                    encoding="utf-8",
                )
                args = SimpleNamespace(
                    outdir=str(match_dir),
                    peak_header_list=["peak_chr", "peak_start", "peak_end"],
                )
                diff_footprints._write_match_motif_site_cache(args, ["TF1_MA0001.1"], logger)
                (bed_dir / "TF1_MA0001.1_all.bed").unlink()

            outdir = root / "out"
            (outdir / "TF1_MA0001.1" / "beds").mkdir(parents=True)
            args = SimpleNamespace(
                cached_match_dirs=[
                    str(root / "A" / "match_motifs"),
                    str(root / "B" / "match_motifs"),
                ],
                sample_names=["A_sample", "B_sample"],
                outdir=str(outdir),
                peak_header_list=["peak_chr", "peak_start", "peak_end"],
            )
            diff_footprints._write_cached_tfbs_tmp_files(args, ["TF1_MA0001.1"], logger)
            tmp_file = outdir / "TF1_MA0001.1" / "beds" / "TF1_MA0001.1.tmp"
            self.assertEqual(
                tmp_file.read_text(encoding="utf-8").strip(),
                "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t1.25000\t2.50000",
            )

    def test_cached_match_dirs_build_shards_from_compact_cache_in_summary_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            logger = diff_footprints.FpToolsLogger("test", 0)
            for sample, score in [("A", "1.25000"), ("B", "2.50000")]:
                match_dir = root / sample / "match_motifs"
                bed_dir = match_dir / "TF1_MA0001.1" / "beds"
                bed_dir.mkdir(parents=True)
                (bed_dir / "TF1_MA0001.1_all.bed").write_text(
                    "chr1\t10\t20\tTF1\t8.0\t+\tchr1\t1\t100\t" + score + "\t" + score + "\tNA\n",
                    encoding="utf-8",
                )
                args = SimpleNamespace(
                    outdir=str(match_dir),
                    peak_header_list=["peak_chr", "peak_start", "peak_end"],
                )
                diff_footprints._write_match_motif_site_cache(args, ["TF1_MA0001.1"], logger)
                shutil.rmtree(match_dir / "TF1_MA0001.1")

            outdir = root / "out"
            args = SimpleNamespace(
                cached_match_dirs=[
                    str(root / "A" / "match_motifs"),
                    str(root / "B" / "match_motifs"),
                ],
                sample_names=["A_sample", "B_sample"],
                outdir=str(outdir),
                peak_header_list=["peak_chr", "peak_start", "peak_end"],
                static_plots=False,
                write_motif_outputs=False,
                cores=2,
            )
            diff_footprints._write_cached_tfbs_tmp_files(args, ["TF1_MA0001.1"], logger)
            self.assertEqual(len(args.cached_motif_bed_maps), 2)
            self.assertTrue((root / "A" / "match_motifs" / "cache" / "motif_sites_by_motif" / "TF1_MA0001.1.bed").exists())
            self.assertTrue((root / "B" / "match_motifs" / "cache" / "motif_sites_by_motif" / "TF1_MA0001.1.bed").exists())
            self.assertFalse((outdir / "TF1_MA0001.1" / "beds" / "TF1_MA0001.1.tmp").exists())

    def test_cached_match_dirs_load_background_without_bigwigs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            logger = diff_footprints.FpToolsLogger("test", 0)
            for sample, scores in [("A", [1.0, 3.0]), ("B", [2.0, 4.0])]:
                match_dir = root / sample / "match_motifs"
                match_dir.mkdir(parents=True)
                args = SimpleNamespace(
                    outdir=str(match_dir),
                    sample_names=[sample],
                    cond_names=[sample],
                    peak_header_list=["peak_chr", "peak_start", "peak_end"],
                    normalization="none",
                )
                background = {
                    "keys": [["chr1", "1", "100", "5"], ["chr1", "1", "100", "25"]],
                    "sample_signal": {sample: scores},
                }
                diff_footprints._write_match_motifs_cache(args, background, logger)

            args = SimpleNamespace(
                cached_without_bigwigs=True,
                cached_match_dirs=[str(root / "A" / "match_motifs"), str(root / "B" / "match_motifs")],
                sample_names=["A", "B"],
                cond_names=["condition"],
                condition_samples={"condition": ["A", "B"]},
            )
            background = diff_footprints._collect_cached_background([], args, None, 1)
            self.assertEqual(background["keys"], [["chr1", "1", "100", "5"], ["chr1", "1", "100", "25"]])
            self.assertEqual(background["sample_signal"]["A"].tolist(), [1.0, 3.0])
            self.assertEqual(background["sample_signal"]["B"].tolist(), [2.0, 4.0])
            self.assertEqual(background["signal"]["condition"].tolist(), [1.5, 3.5])

    def test_match_motifs_sample_names_label_match_only_outputs(self):
        parser = add_diff_footprints_arguments(argparse.ArgumentParser(prog="match-motifs"), command_name="match-motifs")
        args = parser.parse_args(
            [
                "--signals",
                "file_stem_A.bw",
                "file_stem_B.bw",
                "--sample-names",
                "SampleA",
                "SampleB",
                "--genome",
                "hg38.fa.gz",
                "--peaks",
                "peaks.bed",
            ]
        )
        args.match_only = True
        prepared = _prepare_condition_metadata(args)
        self.assertEqual(prepared.cond_names, ["SampleA", "SampleB"])
        self.assertEqual(prepared.sample_names, ["SampleA", "SampleB"])
        self.assertEqual(prepared.condition_samples["SampleA"], ["SampleA"])


if __name__ == "__main__":
    unittest.main()
