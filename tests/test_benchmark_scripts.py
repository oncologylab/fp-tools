import csv
import importlib.util
import pathlib
import shutil
import sys
import sysconfig
import tempfile
import unittest

import numpy as np
import pandas as pd

from fp_tools.utils.multiscale import write_multiscale_npz
from fp_tools.utils import bigwig as pyBigWig

try:
    import pysam
except ImportError:
    pysam = None


ROOT = pathlib.Path(__file__).resolve().parents[1]


def console_script(name):
    scripts_dir = pathlib.Path(sysconfig.get_path("scripts"))
    return shutil.which(name, path=str(scripts_dir)) or str(scripts_dir / name)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_optional_module(name, path):
    if not path.exists():
        return None
    try:
        return load_module(name, path)
    except ModuleNotFoundError as exc:
        missing = getattr(exc, "name", "")
        if missing in {"plot_benchmark_panels", "plot_calibration_panels", "plot_multiscale_npz", "plot_method_comparison"}:
            return None
        raise


download_manifest = load_module("download_manifest", ROOT / "benchmarks" / "scripts" / "download_manifest.py")
compute_binary_metrics = load_module("compute_binary_metrics", ROOT / "benchmarks" / "scripts" / "compute_binary_metrics.py")
compute_calibration = load_module("compute_calibration", ROOT / "benchmarks" / "scripts" / "compute_calibration.py")
plot_benchmark_panels = load_optional_module("plot_benchmark_panels", ROOT / "manuscript" / "scripts" / "plot_benchmark_panels.py")
plot_calibration_panels = load_optional_module("plot_calibration_panels", ROOT / "manuscript" / "scripts" / "plot_calibration_panels.py")
plot_multiscale_npz = load_optional_module("plot_multiscale_npz", ROOT / "manuscript" / "scripts" / "plot_multiscale_npz.py")
build_encode_manifest = load_module("build_encode_manifest", ROOT / "benchmarks" / "scripts" / "build_encode_manifest.py")
build_motif_removal_benchmark = load_module("build_motif_removal_benchmark", ROOT / "benchmarks" / "scripts" / "build_motif_removal_benchmark.py")
build_label_overlap_benchmark = load_module("build_label_overlap_benchmark", ROOT / "benchmarks" / "scripts" / "build_label_overlap_benchmark.py")
run_benchmark_pipeline = load_optional_module("run_benchmark_pipeline", ROOT / "benchmarks" / "scripts" / "run_benchmark_pipeline.py")
score_peaks_with_pwm = load_module("score_peaks_with_pwm", ROOT / "benchmarks" / "scripts" / "score_peaks_with_pwm.py")
footprint_from_bam = load_module("footprint_from_bam", ROOT / "benchmarks" / "scripts" / "footprint_from_bam.py")
footprint_occupancy_score = load_module("footprint_occupancy_score", ROOT / "benchmarks" / "scripts" / "footprint_occupancy_score.py")
benchmark_footprint_kernel = load_module("benchmark_footprint_kernel", ROOT / "benchmarks" / "scripts" / "benchmark_footprint_kernel.py")
build_tf_feature_table = load_module("build_tf_feature_table", ROOT / "benchmarks" / "scripts" / "build_tf_feature_table.py")
evaluate_methods = load_module("evaluate_methods", ROOT / "benchmarks" / "scripts" / "evaluate_methods.py")
plot_method_comparison = load_optional_module("plot_method_comparison", ROOT / "manuscript" / "scripts" / "plot_method_comparison.py")
validate_manifests = load_module("validate_manifests", ROOT / "benchmarks" / "scripts" / "validate_manifests.py")
run_engineering_benchmark = load_module("run_engineering_benchmark", ROOT / "benchmarks" / "scripts" / "run_engineering_benchmark.py")
evaluate_bigwig_site_scores = load_module("evaluate_bigwig_site_scores", ROOT / "benchmarks" / "scripts" / "evaluate_bigwig_site_scores.py")
search_tf_footprint_models = load_module("search_tf_footprint_models", ROOT / "benchmarks" / "scripts" / "search_tf_footprint_models.py")
compare_frozen_tf_candidates = load_module("compare_frozen_tf_candidates", ROOT / "benchmarks" / "scripts" / "compare_frozen_tf_candidates.py")
discover_encode_chip_peaks = load_module("discover_encode_chip_peaks", ROOT / "benchmarks" / "scripts" / "discover_encode_chip_peaks.py")
build_footprint_site_labels = load_module("build_footprint_site_labels", ROOT / "benchmarks" / "scripts" / "build_footprint_site_labels.py")
build_encode_tf_site_matrix = load_module("build_encode_tf_site_matrix", ROOT / "benchmarks" / "scripts" / "build_encode_tf_site_matrix.py")
match_tf_sites_on_accessibility = load_module("match_tf_sites_on_accessibility", ROOT / "benchmarks" / "scripts" / "match_tf_sites_on_accessibility.py")
plot_frozen_tf_profiles = load_module("plot_frozen_tf_profiles", ROOT / "benchmarks" / "scripts" / "plot_frozen_tf_profiles.py")
summarize_tf_footprint_search = load_module("summarize_tf_footprint_search", ROOT / "benchmarks" / "scripts" / "summarize_tf_footprint_search.py")
evaluate_tf_correction_transfer = load_module("evaluate_tf_correction_transfer", ROOT / "benchmarks" / "scripts" / "evaluate_tf_correction_transfer.py")


class TfFootprintModelSearchTest(unittest.TestCase):
    def test_geometry_search_recovers_center_depletion(self):
        rng = np.random.default_rng(12)
        profiles = rng.normal(2.0, 0.15, size=(120, 81))
        labels = np.repeat([0, 1], 60)
        profiles[labels == 1, 35:46] -= 1.0
        candidate = search_tf_footprint_models.Candidate(
            "raw", center_width=11, flank_width=12, gap=2
        )
        scores = search_tf_footprint_models.score_candidate(profiles, candidate)
        metrics = search_tf_footprint_models.binary_metrics(labels, scores)
        self.assertGreater(metrics["auroc"], 0.99)
        self.assertGreater(metrics["auprc"], 0.99)

    def test_asymmetry_penalty_rejects_one_sided_artifact(self):
        profiles = np.ones((2, 81), dtype=float)
        profiles[:, 35:46] = 0.5
        profiles[1, 21:33] = 4.0
        base = search_tf_footprint_models.Candidate(
            "DWM", center_width=11, flank_width=12, gap=2
        )
        penalized = search_tf_footprint_models.replace(base, asymmetry_penalty=1.0)
        base_scores = search_tf_footprint_models.score_candidate(profiles, base)
        penalized_scores = search_tf_footprint_models.score_candidate(profiles, penalized)
        self.assertAlmostEqual(base_scores[0], penalized_scores[0])
        self.assertLess(penalized_scores[1], base_scores[1])

    def test_sampler_caps_each_label_and_split(self):
        rows = []
        for split in ("train", "validation"):
            for label in (0, 1):
                for index in range(9):
                    rows.append(
                        {
                            "cell": "K562", "tf": "CTCF", "chromosome_split": split,
                            "chip_label": label, "TFBS_chr": "chr1",
                            "TFBS_start": index + label * 100, "TFBS_end": index + label * 100 + 1,
                        }
                    )
        sampled = search_tf_footprint_models.deterministic_class_sample(
            pd.DataFrame(rows), maximum_per_class=3, seed=4
        )
        counts = sampled.groupby(["chromosome_split", "chip_label"]).size()
        self.assertTrue((counts == 3).all())

        pooled = search_tf_footprint_models.deterministic_class_sample(
            pd.DataFrame(rows), maximum_per_class=2, seed=4, negative_pool_multiplier=3
        )
        pooled_counts = pooled.groupby(["chromosome_split", "chip_label"]).size()
        self.assertEqual(int(pooled_counts.loc[("train", 0)]), 6)
        self.assertEqual(int(pooled_counts.loc[("train", 1)]), 2)

    def test_evaluation_regions_are_padded_and_merged(self):
        sites = pd.DataFrame(
            {
                "TFBS_chr": ["chr1", "chr1", "chr2"],
                "TFBS_start": [100, 105, 50],
                "TFBS_end": [101, 106, 51],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output = pathlib.Path(tmpdir) / "regions.bed"
            count = search_tf_footprint_models.write_merged_regions(sites, output, padding=10)
            self.assertEqual(count, 2)
            self.assertEqual(
                output.read_text(encoding="utf-8").splitlines()[0].split("\t")[:3],
                ["chr1", "90", "116"],
            )

    def test_frozen_comparison_uses_common_finite_sites(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            sites = pd.DataFrame(
                {
                    "cell": ["K562"] * 4,
                    "tf": ["CTCF"] * 4,
                    "TFBS_chr": ["chr1"] * 4,
                    "TFBS_start": [40, 50, 60, 70],
                    "TFBS_end": [41, 51, 61, 71],
                    "chromosome_split": ["validation"] * 4,
                    "chip_label": [0, 0, 1, 1],
                }
            )
            signal = root / "baseline.bw"
            handle = pyBigWig.open(str(signal), "w")
            handle.addHeader([("chr1", 100)])
            handle.addEntries("chr1", 40, values=[0.0, 0.0, 1.0], span=10, step=10)
            handle.close()
            profiles = np.ones((4, 81), dtype=np.float32)
            profiles[2:, 35:46] = 0.0
            np.savez(root / "K562.raw.flank40.npz", profiles=profiles, valid=np.ones(4, bool))
            winners = pd.DataFrame(
                [
                    {
                        "cell": "K562", "tf": "CTCF", "correction": "raw",
                        "center_width": 11, "flank_width": 12, "gap": 2,
                        "shoulder": "mean", "center": "mean", "normalization": "none",
                        "asymmetry_penalty": 0.0,
                    }
                ]
            )
            baselines = pd.DataFrame(
                [{"cell": "K562", "method": "legacy", "signal": str(signal)}]
            )
            result = compare_frozen_tf_candidates.compare(
                sites, winners, baselines, root, flank=40, split="validation"
            )
            self.assertEqual(int(result.loc[0, "n_sites"]), 3)
            self.assertEqual(float(result.loc[0, "candidate_auroc"]), 1.0)

    def test_encode_selector_prefers_unperturbed_optimal_idr(self):
        frame = pd.DataFrame(
            [
                {
                    "cell": "K562", "tf": "CTCF", "file_accession": "NEW",
                    "output_type": "IDR thresholded peaks", "perturbed": True,
                    "preferred_default": True, "biological_replicate_count": 2,
                    "date_created": "2026", "experiment_accession": "E1",
                },
                {
                    "cell": "K562", "tf": "CTCF", "file_accession": "OPT",
                    "output_type": "optimal IDR thresholded peaks", "perturbed": False,
                    "preferred_default": True, "biological_replicate_count": 2,
                    "date_created": "2020", "experiment_accession": "E2",
                },
            ]
        )
        selected = discover_encode_chip_peaks.select_candidates(frame)
        self.assertEqual(selected.loc[0, "file_accession"], "OPT")

    def test_extracts_only_study_motifs(self):
        study = {
            "tasks": [
                {"motif_id": "MA0001.1", "split": "development"},
                {"motif_id": "MA0002.1", "split": "locked_holdout"},
            ]
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            database = root / "all.jaspar"
            database.write_text(
                ">MA0001.1 ONE\nA [ 1 2 ]\n>MA0002.1 TWO\nA [ 3 4 ]\n",
                encoding="utf-8",
            )
            output = root / "subset.jaspar"
            count = discover_encode_chip_peaks.extract_task_motifs(
                study, database, output, "development"
            )
            self.assertEqual(count, 1)
            self.assertIn("MA0001.1", output.read_text(encoding="utf-8"))
            self.assertNotIn("MA0002.1", output.read_text(encoding="utf-8"))

    def test_finds_one_motif_site_bed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            bed = root / "TF_MA0001.1" / "beds" / "TF_MA0001.1_all.bed"
            bed.parent.mkdir(parents=True)
            bed.write_text("chr1\t1\t2\n", encoding="utf-8")
            self.assertEqual(
                build_encode_tf_site_matrix.motif_site_file(root, "MA0001.1"), bed
            )
            self.assertEqual(build_encode_tf_site_matrix.cell_motif_root(root, "K562"), root)
            (root / "K562").mkdir()
            self.assertEqual(
                build_encode_tf_site_matrix.cell_motif_root(root, "K562"), root / "K562"
            )

    def test_accessibility_matching_reduces_feature_imbalance(self):
        rows = []
        for label, values in ((1, [2.0, 3.0, 8.0]), (0, [1.9, 3.1, 7.9, 20.0])):
            for index, value in enumerate(values):
                rows.append(
                    {
                        "cell": "K562", "tf": "CTCF", "chromosome_split": "train",
                        "chip_label": label, "motif_score": value,
                        "accessibility": value, "TFBS_chr": "chr1",
                        "TFBS_start": index + label * 100, "TFBS_end": index + label * 100 + 1,
                    }
                )
        sites = pd.DataFrame(rows)
        matched, diagnostics = match_tf_sites_on_accessibility.match_sites(
            sites, ["motif_score", "accessibility"], negative_ratio=1, seed=4
        )
        self.assertEqual(len(matched), 6)
        self.assertLess(
            abs(diagnostics.loc[0, "after_smd_accessibility"]),
            abs(diagnostics.loc[0, "before_smd_accessibility"]),
        )

    def test_display_normalization_removes_outer_flank_level(self):
        profiles = np.tile(np.linspace(2.0, 4.0, 81), (3, 1))
        normalized = plot_frozen_tf_profiles.normalize_profiles_for_display(
            profiles, outer_width=10
        )
        outer = np.concatenate([normalized[:, :10], normalized[:, -10:]], axis=1)
        np.testing.assert_allclose(np.mean(outer, axis=1), 0.0, atol=1e-12)

    def test_detectability_summary_prioritizes_abstention_reasons(self):
        row = pd.Series(
            {
                "positive_sites": 800, "after_smd_accessibility": 0.5,
                "after_smd_motif_score": 0.0, "candidate_auroc": 0.9,
            }
        )
        status = summarize_tf_footprint_search.classify_row(
            row, 500, 0.25, 0.5, 0.65, 0.70, 0.75
        )
        self.assertEqual(status, "accessibility_confounded")
        row["positive_sites"] = 20
        self.assertEqual(
            summarize_tf_footprint_search.classify_row(
                row, 500, 0.25, 0.5, 0.65, 0.70, 0.75
            ),
            "underpowered",
        )

    def test_correction_transfer_keeps_geometry_fixed(self):
        labels = np.array([0, 0, 1, 1])
        raw = np.ones((4, 81), dtype=float)
        raw[labels == 1, 35:46] = 0.0
        flat = np.ones_like(raw)
        candidate = search_tf_footprint_models.Candidate(
            "raw", center_width=11, flank_width=12, gap=2
        )
        metrics = evaluate_tf_correction_transfer.evaluate_corrections(
            {"raw": raw, "DWM": flat}, labels, candidate
        )
        self.assertEqual(metrics.loc[metrics["auroc"].idxmax(), "correction"], "raw")


class BigwigSiteScoreTest(unittest.TestCase):
    def test_scores_fixed_site_centers_and_computes_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            signal = root / "scores.bw"
            handle = pyBigWig.open(str(signal), "w")
            handle.addHeader([("chr1", 20)])
            handle.addEntries("chr1", 0, values=[float(value) for value in range(20)], span=1, step=1)
            handle.close()
            sites = pd.DataFrame(
                {
                    "cell": ["K562"] * 4,
                    "tf": ["CTCF"] * 4,
                    "TFBS_chr": ["chr1"] * 4,
                    "TFBS_start": [0, 2, 10, 12],
                    "TFBS_end": [2, 4, 12, 14],
                    "chip_label": [0, 0, 1, 1],
                }
            )
            signals = pd.DataFrame(
                {"cell": ["K562"], "method": ["candidate"], "signal": [str(signal)]}
            )
            predictions, metrics = evaluate_bigwig_site_scores.evaluate(sites, signals)
            self.assertEqual(predictions["score"].tolist(), [1.0, 3.0, 11.0, 13.0])
            self.assertEqual(float(metrics.loc[0, "auroc"]), 1.0)
            self.assertEqual(float(metrics.loc[0, "auprc"]), 1.0)

    def test_uses_common_finite_sites_across_methods(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = pathlib.Path(tmpdir)
            full_signal = root / "full.bw"
            short_signal = root / "short.bw"
            for path, length in ((full_signal, 20), (short_signal, 12)):
                handle = pyBigWig.open(str(path), "w")
                handle.addHeader([("chr1", length)])
                handle.addEntries(
                    "chr1", 0, values=[float(value) for value in range(length)], span=1, step=1
                )
                handle.close()
            sites = pd.DataFrame(
                {
                    "cell": ["K562"] * 4,
                    "tf": ["CTCF"] * 4,
                    "TFBS_chr": ["chr1"] * 4,
                    "TFBS_start": [0, 2, 10, 12],
                    "TFBS_end": [2, 4, 12, 14],
                    "chip_label": [0, 0, 1, 1],
                }
            )
            signals = pd.DataFrame(
                {
                    "cell": ["K562", "K562"],
                    "method": ["full", "short"],
                    "signal": [str(full_signal), str(short_signal)],
                }
            )
            predictions, metrics = evaluate_bigwig_site_scores.evaluate(sites, signals)
            self.assertEqual(len(predictions), 6)
            self.assertEqual(set(metrics["n_sites"]), {3})

    def test_paired_chromosome_bootstrap_uses_matched_blocks(self):
        rows = []
        for chrom in ("chr17", "chr18"):
            for index, label in enumerate([0, 0, 1, 1]):
                for method, score in (
                    ("raw", float(index % 2)),
                    ("candidate", float(index)),
                ):
                    rows.append(
                        {
                            "cell": "K562",
                            "tf": "CTCF",
                            "TFBS_chr": chrom,
                            "TFBS_start": index,
                            "TFBS_end": index + 1,
                            "chip_label": label,
                            "method": method,
                            "score": score,
                        }
                    )
        result = evaluate_bigwig_site_scores.paired_chromosome_bootstrap(
            pd.DataFrame(rows), baseline_method="raw", n_bootstrap=20, seed=4
        )
        self.assertEqual(set(result["metric"]), {"auroc", "auprc"})
        self.assertTrue((result["successful_bootstraps"] == 20).all())
        self.assertTrue((result["probability_delta_gt_zero"] == 1.0).all())


class ManifestValidationTest(unittest.TestCase):
    def test_committed_manifests_validate(self):
        errors = validate_manifests.validate_manifests(ROOT / "benchmarks" / "manifests")
        self.assertEqual(errors, [])

    def test_full_manifest_reports_missing_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = pathlib.Path(tmpdir) / "bad.tsv"
            path.write_text("source\turl\nENCODE\thttps://example.org/file.bam\n", encoding="utf-8")
            errors = validate_manifests.validate_manifest(path)
            self.assertTrue(errors)
            self.assertIn("missing full-manifest columns", errors[0])


class EngineeringBenchmarkHelperTest(unittest.TestCase):
    def test_records_command_runtime_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = pathlib.Path(tmpdir) / "runtime.tsv"
            row = run_engineering_benchmark.run_benchmark(
                [sys.executable, "-c", "print('ok')"],
                out,
                "python-smoke",
                cores=1,
            )
            self.assertEqual(row["exit_code"], 0)
            table = pd.read_csv(out, sep="\t")
            self.assertEqual(table.loc[0, "label"], "python-smoke")
            self.assertIn("wall_seconds", table.columns)
            self.assertIn("peak_rss_kb", table.columns)

    def test_footprint_kernel_benchmark_compares_fast_and_legacy_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            regions = tmp / "small_regions.bed"
            rows = [
                line
                for line in (ROOT / "test_data" / "merged_peaks.bed").read_text(encoding="utf-8").splitlines()
                if line.strip() and not line.startswith("#")
            ][:8]
            regions.write_text("\n".join(rows) + "\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "signal": ROOT / "test_data" / "Bcell_corrected.bw",
                    "regions": regions,
                    "outdir": tmp / "kernel_benchmark",
                    "call_footprints": console_script("call-footprints"),
                    "cores": 1,
                    "chunk_size": 1_000_000,
                    "verbosity": 1,
                    "workflow_second_signal": None,
                    "workflow_first_name": "Bcell",
                    "workflow_second_name": "Tcell",
                    "workflow_first_condition": "Bcell",
                    "workflow_second_condition": "Tcell",
                    "genome": None,
                    "motifs": None,
                    "motif_db": None,
                    "match_motifs": console_script("match-motifs"),
                    "diff_footprints": console_script("diff-footprints"),
                },
            )()
            summary = benchmark_footprint_kernel.run_kernel_benchmark(args)
            benchmark_footprint_kernel.write_summary(summary, args.outdir)

            self.assertGreater(summary["legacy_seconds"], 0)
            self.assertGreater(summary["fast_seconds"], 0)
            self.assertGreater(summary["speedup"], 0)
            self.assertLessEqual(summary["bigwig_max_abs_diff"], 2e-5)
            self.assertLessEqual(summary["bigwig_mean_abs_diff"], 1e-6)
            self.assertGreaterEqual(summary["bed_coordinate_jaccard"], 0.999)
            self.assertTrue((args.outdir / "kernel_benchmark_summary.tsv").exists())
            self.assertTrue((args.outdir / "kernel_benchmark_summary.json").exists())


class FeatureTableAndMethodsTest(unittest.TestCase):

    def test_gc_content(self):
        self.assertAlmostEqual(build_tf_feature_table.gc_content("GGCC"), 1.0)
        self.assertAlmostEqual(build_tf_feature_table.gc_content("ATAT"), 0.0)
        self.assertAlmostEqual(build_tf_feature_table.gc_content("ACGT"), 0.5)

    def test_read_accessibility(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed = pathlib.Path(tmpdir) / "a.bed"
            bed.write_text("#chrom\tstart\tend\tname\tscore\nchr1\t0\t10\tp1\t3.5\nchr1\t20\t30\tp2\t7.0\n", encoding="utf-8")
            acc = build_tf_feature_table.read_accessibility(bed)
            self.assertEqual(acc, {"p1": 3.5, "p2": 7.0})

    def test_integrated_oof_scores_separable(self):
        import numpy as np

        rng = np.random.default_rng(0)
        n = 200
        motif = np.concatenate([rng.normal(2, 1, n // 2), rng.normal(-2, 1, n // 2)])
        frame = pd.DataFrame({
            "accessibility": rng.normal(0, 1, n),
            "motif": motif,
            "gc": rng.normal(0.5, 0.1, n),
            "label": np.array([1] * (n // 2) + [0] * (n // 2)),
        })
        scores = evaluate_methods.integrated_oof_scores(frame, ["accessibility", "motif", "gc"], seed=1)
        self.assertEqual(len(scores), n)
        # Positive class should get higher mean score than negative class.
        self.assertGreater(scores[frame["label"] == 1].mean(), scores[frame["label"] == 0].mean())


class FootprintFromBamTest(unittest.TestCase):

    def test_footprint_score_positive_when_center_depleted(self):
        counts = np.array([5, 5, 5, 0, 0, 0, 5, 5, 5], dtype=float)
        # center = indices 3:6 (all zero), flanks of width 3 on each side (all 5s)
        score = footprint_from_bam.footprint_score(counts, 3, 6, 3)
        self.assertGreater(score, 0)

    def test_footprint_score_negative_when_center_enriched(self):
        counts = np.array([0, 0, 0, 9, 9, 9, 0, 0, 0], dtype=float)
        self.assertLess(footprint_from_bam.footprint_score(counts, 3, 6, 3), 0)

    @unittest.skipUnless(pysam is not None, "pysam is required to write BAM fixtures")
    def test_build_cutsites_counts_tn5_insertions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            bam_path = tmp / "mini.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": 200}]}
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
                a = pysam.AlignedSegment()
                a.query_name = "r1"
                a.query_sequence = "A" * 20
                a.flag = 0  # forward
                a.reference_id = 0
                a.reference_start = 50
                a.mapping_quality = 60
                a.cigarstring = "20M"
                a.query_qualities = pysam.qualitystring_to_array("I" * 20)
                bam.write(a)
            pysam.index(str(bam_path))
            counts = footprint_from_bam.build_cutsites(str(bam_path), "chrT", 200)
            # forward read start 50 + shift 4 = insertion at 54
            self.assertEqual(int(counts[54]), 1)
            self.assertEqual(int(counts.sum()), 1)


class FootprintOccupancyScoreTest(unittest.TestCase):

    def test_best_match_finds_motif_offset(self):
        from fp_tools.tools.variants import read_pwm_motifs

        with tempfile.TemporaryDirectory() as tmpdir:
            meme = pathlib.Path(tmpdir) / "m.meme"
            meme.write_text(
                "MEME version 4\n\nALPHABET= ACGT\n\nMOTIF m M\n"
                "letter-probability matrix: alength= 4 w= 4 nsites= 1 E= 0\n"
                "0.91 0.03 0.03 0.03\n0.91 0.03 0.03 0.03\n0.91 0.03 0.03 0.03\n0.91 0.03 0.03 0.03\n",
                encoding="utf-8",
            )
            motif = read_pwm_motifs(meme)[0]
            score, off = footprint_occupancy_score.best_match("CCCCAAAACCCC", motif)
            self.assertEqual(off, 4)
            self.assertGreater(score, 0)

    def test_footprint_score_sign(self):
        import numpy as np

        depleted = np.array([4.0, 4.0, 0.0, 0.0, 4.0, 4.0])
        self.assertGreater(footprint_occupancy_score.footprint_score(depleted, 2, 4, 2), 0)


class PlotMethodComparisonTest(unittest.TestCase):

    @unittest.skipIf(plot_method_comparison is None, "manuscript plotting scripts are not published")
    def test_plot_writes_three_formats(self):
        metrics = pd.DataFrame(
            {
                "group": ["CTCF/accessibility", "CTCF/fp-tools-motif", "global"],
                "auroc": [0.53, 0.87, 0.70],
                "auprc": [0.36, 0.81, 0.58],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "cmp"
            outputs = plot_method_comparison.plot_method_comparison(metrics, prefix)
            self.assertEqual(len(outputs), 3)
            for path in outputs:
                self.assertTrue(path.exists())


class ScorePeaksWithPwmTest(unittest.TestCase):

    def test_best_pwm_match_ranks_motif_bearing_peak_higher(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            # Synthetic contig: a strong A-run motif planted in peak 1, absent in peak 2.
            seq = ("A" * 40) + ("AAAAAAAA") + ("C" * 40) + ("CGCGCGCG") + ("C" * 40)
            fa = tmp / "mini.fa"
            fa.write_text(f">chrT\n{seq}\n", encoding="utf-8")
            if pysam is not None:
                pysam.faidx(str(fa))

            # MEME motif favouring poly-A (length 8).
            meme = tmp / "polyA.meme"
            rows = "\n".join("0.97 0.01 0.01 0.01" for _ in range(8))
            meme.write_text(
                "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
                "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
                "MOTIF polyA TEST\nletter-probability matrix: alength= 4 w= 8 nsites= 1 E= 0\n"
                + rows + "\n",
                encoding="utf-8",
            )

            peaks = tmp / "peaks.bed"
            peaks.write_text("chrT\t0\t48\tpeak_A\nchrT\t88\t136\tpeak_C\n", encoding="utf-8")
            out = tmp / "scored.bed"
            n = score_peaks_with_pwm.score_peaks_with_pwm(peaks, fa, meme, out)
            self.assertEqual(n, 2)

            scored = pd.read_csv(out, sep="\t")
            by_name = dict(zip(scored["name"], scored["score"]))
            self.assertGreater(by_name["peak_A"], by_name["peak_C"])

    def test_chroms_not_in_genome_are_skipped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            fa = tmp / "mini.fa"
            fa.write_text(">chrT\n" + "ACGT" * 20 + "\n", encoding="utf-8")
            if pysam is not None:
                pysam.faidx(str(fa))
            meme = tmp / "m.meme"
            meme.write_text(
                "MEME version 4\n\nALPHABET= ACGT\n\nMOTIF m M\n"
                "letter-probability matrix: alength= 4 w= 4 nsites= 1 E= 0\n"
                "0.7 0.1 0.1 0.1\n0.1 0.7 0.1 0.1\n0.1 0.1 0.7 0.1\n0.1 0.1 0.1 0.7\n",
                encoding="utf-8",
            )
            peaks = tmp / "peaks.bed"
            peaks.write_text("chrT\t0\t40\tkeep\nchrZ\t0\t40\tdrop\n", encoding="utf-8")
            out = tmp / "scored.bed"
            n = score_peaks_with_pwm.score_peaks_with_pwm(peaks, fa, meme, out)
            self.assertEqual(n, 1)
            scored = pd.read_csv(out, sep="\t")
            self.assertEqual(list(scored["name"]), ["keep"])


class BenchmarkScriptsTest(unittest.TestCase):

    def test_encode_manifest_choose_file_prefers_requested_assembly(self):
        files = [
            {"status": "released", "output_type": "alignments", "file_format": "bam", "assembly": "hg19", "accession": "OLD"},
            {"status": "released", "output_type": "alignments", "file_format": "bam", "assembly": "GRCh38", "accession": "NEW"},
        ]
        selected = build_encode_manifest.choose_file(files, "alignments", "bam", assembly="GRCh38")
        self.assertEqual(selected["accession"], "NEW")


    def test_encode_smoke_manifest_is_versioned_and_downloadable(self):
        manifest = ROOT / "benchmarks" / "manifests" / "encode_k562_ctcf_smoke.tsv"
        rows = download_manifest.read_manifest(manifest)
        self.assertEqual(len(rows), 3)
        self.assertEqual({row["source"] for row in rows}, {"ENCODE"})
        self.assertEqual({row["assembly"] for row in rows}, {"GRCh38"})
        self.assertTrue(any(row["file_format"] == "bam" for row in rows))
        self.assertTrue(all(row["url"].startswith("https://www.encodeproject.org/") for row in rows))

    def test_download_manifest_dry_run_writes_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            manifest = tmp / "manifest.tsv"
            report = tmp / "download_report.tsv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["file_accession", "url", "local_path", "checksum"], delimiter="	")
                writer.writeheader()
                writer.writerow({"file_accession": "FILE1", "url": "https://example.org/file.bam", "local_path": str(tmp / "file.bam"), "checksum": ""})
            rows = download_manifest.read_manifest(manifest)
            results = [download_manifest.download_one(row, dry_run=True) for row in rows]
            download_manifest.write_report(results, report)
            text = report.read_text(encoding="utf-8")
            self.assertIn("dry_run", text)
            self.assertIn("FILE1", text)

    def test_download_manifest_uses_validated_atomic_partial(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            source = tmp / "source.fastq.gz"
            source.write_bytes(b"validated-fastq")
            output = tmp / "downloads" / "sample.fastq.gz"
            row = {
                "file_accession": "SRR_TEST_1",
                "url": source.as_uri(),
                "local_path": str(output),
                "checksum": download_manifest.md5sum(source),
                "expected_bytes": str(source.stat().st_size),
            }
            result = download_manifest.download_one(row, downloader="urllib")
            self.assertEqual(result.status, "downloaded")
            self.assertEqual(output.read_bytes(), source.read_bytes())
            self.assertFalse(output.with_name(output.name + ".partial").exists())

    def test_compute_binary_metrics_global_and_grouped(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0, 1, 0],
                "score": [0.95, 0.05, 0.80, 0.10, 0.70, 0.20],
                "tf": ["A", "A", "A", "A", "B", "B"],
                "cell": ["C", "C", "C", "C", "C", "C"],
                "method": ["fp", "fp", "fp", "fp", "fp", "fp"],
            }
        )
        metrics = compute_binary_metrics.compute_metrics(df, "label", "score", ["tf", "cell", "method"])
        self.assertIn("global", set(metrics["group"]))
        self.assertGreaterEqual(metrics.loc[metrics["group"] == "global", "auprc"].iloc[0], 0.99)

    def test_build_motif_removal_benchmark_zeroes_baseline_and_summarizes(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0],
                "tf": ["CTCF", "CTCF", "GATA1", "GATA1"],
                "cell": ["K562", "K562", "K562", "K562"],
                "motif_family": ["CTCF", "CTCF", "GATA", "GATA"],
                "motif_score": [12.0, 8.0, 10.0, 7.0],
                "rank_score": [0.91, 0.15, 0.77, 0.20],
                "candidate_score": [5.0, 1.0, 4.0, 1.5],
            }
        )
        table = build_motif_removal_benchmark.build_motif_removal_table(
            df,
            remove_col="motif_family",
            remove_values=["CTCF"],
            baseline_score_col="motif_score",
            recovery_score_cols=["rank_score", "candidate_score"],
        )
        self.assertEqual(set(table["method"]), {"motif_removed_baseline", "rank_score", "candidate_score"})
        self.assertEqual(len(table), 6)
        baseline = table[table["method"] == "motif_removed_baseline"]
        self.assertTrue((baseline["score"] == 0.0).all())
        self.assertTrue(table["removed"].all())
        self.assertEqual(set(table["removal_target"]), {"CTCF"})

        summary = build_motif_removal_benchmark.summarize_removal_table(table)
        self.assertEqual(set(summary["method"]), {"motif_removed_baseline", "rank_score", "candidate_score"})
        self.assertTrue((summary["n"] == 2).all())

    def test_motif_removal_table_can_include_controls(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1],
                "motif_family": ["ETS", "AP1", "ETS"],
                "rank_score": [0.8, 0.2, 0.7],
            }
        )
        table = build_motif_removal_benchmark.build_motif_removal_table(
            df,
            remove_col="motif_family",
            remove_values=["ETS"],
            recovery_score_cols=["rank_score"],
            include_controls=True,
        )
        self.assertEqual(len(table), 3)
        self.assertEqual(table["removed"].tolist(), [True, False, True])

    def test_build_label_overlap_benchmark_outputs_metrics_ready_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            predictions = tmp / "predictions.bed"
            labels = tmp / "labels.bed"
            out = tmp / "benchmark.tsv"
            predictions.write_text(
                "#chrom\tstart\tend\tname\tscore\tmotif_family\n"
                "chr1\t10\t20\tsite1\t0.9\tCTCF\n"
                "chr1\t30\t40\tsite2\t0.2\tCTCF\n"
                "chr2\t5\t15\tsite3\t0.7\tGATA\n",
                encoding="utf-8",
            )
            labels.write_text("chr1\t15\t25\nchr2\t0\t6\n", encoding="utf-8")

            table = build_label_overlap_benchmark.build_label_overlap_table(
                predictions,
                labels,
                out,
                min_overlap_bp=2,
                method="reranked",
                tf="CTCF",
                cell="K562",
                metadata_cols=["name", "motif_family"],
            )
            self.assertTrue(out.exists())

        self.assertEqual(table["label"].tolist(), [1, 0, 0])
        self.assertEqual(table["method"].unique().tolist(), ["reranked"])
        self.assertEqual(table["tf"].unique().tolist(), ["CTCF"])
        self.assertIn("motif_family", table.columns)

    def test_compute_binary_metrics_bootstrap_confidence_intervals(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0, 1, 0],
                "score": [0.95, 0.05, 0.80, 0.10, 0.70, 0.20],
                "method": ["fp", "fp", "fp", "fp", "fp", "fp"],
            }
        )
        boot = compute_binary_metrics.bootstrap_confidence_intervals(
            df,
            "label",
            "score",
            ["method"],
            n_bootstrap=25,
            seed=7,
        )
        auprc = boot[(boot["group"] == "global") & (boot["metric"] == "auprc")].iloc[0]
        self.assertEqual(int(auprc["n_bootstrap"]), 25)
        self.assertGreaterEqual(float(auprc["ci_high"]), float(auprc["ci_low"]))
        self.assertGreater(int(auprc["successful_bootstraps"]), 0)

    def test_compute_binary_metrics_resamples_complete_genomic_blocks(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2", "chr2", "chr3", "chr3"],
                "label": [1, 0, 1, 0, 1, 0],
                "score": [0.95, 0.05, 0.80, 0.10, 0.70, 0.20],
            }
        )
        first = compute_binary_metrics.bootstrap_confidence_intervals(
            df,
            "label",
            "score",
            [],
            n_bootstrap=20,
            seed=11,
            block_cols=["chrom"],
        )
        second = compute_binary_metrics.bootstrap_confidence_intervals(
            df,
            "label",
            "score",
            [],
            n_bootstrap=20,
            seed=11,
            block_cols=["chrom"],
        )
        pd.testing.assert_frame_equal(first, second)
        self.assertEqual(set(first["resampling_unit"]), {"chrom"})
        self.assertEqual(set(first["n_blocks"]), {3})

        with self.assertRaisesRegex(ValueError, "block columns are missing"):
            compute_binary_metrics.bootstrap_confidence_intervals(
                df,
                "label",
                "score",
                [],
                n_bootstrap=2,
                block_cols=["peak_id"],
            )

    def test_compute_binary_metrics_allows_raw_non_probability_scores(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0],
                "score": [5.0, 1.0, 4.0, 2.0],
                "method": ["raw", "raw", "raw", "raw"],
            }
        )
        metrics = compute_binary_metrics.compute_metrics(df, "label", "score", ["method"])
        global_row = metrics[metrics["group"] == "global"].iloc[0]
        self.assertGreaterEqual(global_row["auprc"], 0.99)
        self.assertTrue(np.isnan(global_row["brier"]))

    @unittest.skipIf(run_benchmark_pipeline is None, "benchmark figure pipeline requires unpublished plotting helpers")
    def test_run_benchmark_pipeline_writes_tables_figures_and_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = pathlib.Path(tmpdir)
            predictions = tmp / "predictions.tsv"
            outdir = tmp / "benchmark_run"
            pd.DataFrame(
                {
                    "label": [1, 0, 1, 0],
                    "score": [0.9, 0.1, 0.8, 0.2],
                    "tf": ["CTCF", "CTCF", "IRF1", "IRF1"],
                    "cell": ["K562", "K562", "K562", "K562"],
                    "method": ["fp-tools", "fp-tools", "fp-tools", "fp-tools"],
                }
            ).to_csv(predictions, sep="\t", index=False)

            outputs = run_benchmark_pipeline.run_benchmark_pipeline(
                [predictions],
                outdir,
                bins=5,
                bootstrap=5,
                title="synthetic benchmark",
            )

            self.assertTrue(outputs["combined_predictions"].exists())
            self.assertTrue(outputs["binary_metrics"].exists())
            self.assertTrue(outputs["binary_metrics_bootstrap"].exists())
            self.assertTrue(outputs["run_summary"].exists())
            for key in ("benchmark_figures", "calibration_figures"):
                self.assertEqual(len(outputs[key]), 3)
                self.assertTrue(all(path.exists() for path in outputs[key]))

    @unittest.skipIf(plot_calibration_panels is None, "manuscript plotting scripts are not published")
    def test_compute_and_plot_calibration_reports(self):
        df = pd.DataFrame(
            {
                "label": [1, 0, 1, 0, 1, 0, 1, 0],
                "score": [0.95, 0.10, 0.80, 0.15, 0.70, 0.30, 0.60, 0.40],
                "tf": ["A", "A", "A", "A", "B", "B", "B", "B"],
                "cell": ["K562"] * 8,
                "method": ["model"] * 8,
            }
        )
        bins, summary = compute_calibration.compute_calibration(df, "label", "score", ["tf", "cell", "method"], bins=5)
        self.assertIn("global", set(summary["group"]))
        self.assertTrue((summary["ece"] >= 0).all())
        self.assertEqual(len(bins[bins["group"] == "global"]), 5)

        with tempfile.TemporaryDirectory() as tmpdir:
            out_prefix = pathlib.Path(tmpdir) / "figure_calibration"
            outputs = plot_calibration_panels.plot_calibration(bins, summary, out_prefix)
            self.assertEqual({path.suffix for path in outputs}, {".pdf", ".svg", ".png"})
            for output in outputs:
                self.assertTrue(output.exists())

    @unittest.skipIf(plot_multiscale_npz is None, "manuscript plotting scripts are not published")
    def test_plot_multiscale_npz_writes_all_formats(self):
        records = [
            (("chr1", 0, 5), {8: np.array([1, 2, 3, 2, 1]), 16: np.array([0, 1, 2, 1, 0])}),
            (("chr1", 10, 15), {8: np.array([2, 3, 4, 3, 2]), 16: np.array([1, 2, 3, 2, 1])}),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            npz_path = pathlib.Path(tmpdir) / "multiscale.npz"
            out_prefix = pathlib.Path(tmpdir) / "figure_multiscale"
            write_multiscale_npz(str(npz_path), records, [8, 16], "max")
            outputs = plot_multiscale_npz.plot_multiscale_npz(npz_path, out_prefix)
            self.assertEqual([path.suffix for path in outputs], [".pdf", ".svg", ".png"])
            for output in outputs:
                self.assertTrue(output.exists())
                self.assertGreater(output.stat().st_size, 0)

    @unittest.skipIf(plot_benchmark_panels is None, "manuscript plotting scripts are not published")
    def test_plot_benchmark_panels_writes_all_formats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_prefix = pathlib.Path(tmpdir) / "figure_benchmark_summary"
            metrics = pd.DataFrame(
                {
                    "group": ["global", "A/C/fp", "B/C/fp"],
                    "n": [6, 4, 2],
                    "positives": [3, 2, 1],
                    "auroc": [1.0, 1.0, 1.0],
                    "auprc": [1.0, 1.0, 1.0],
                    "recall_at_1pct_fdr": [1.0, 1.0, 1.0],
                    "recall_at_5pct_fdr": [1.0, 1.0, 1.0],
                    "recall_at_10pct_fdr": [1.0, 1.0, 1.0],
                    "brier": [0.05, 0.04, 0.06],
                }
            )
            outputs = plot_benchmark_panels.plot_metrics(metrics, out_prefix)
            self.assertEqual({path.suffix for path in outputs}, {".pdf", ".svg", ".png"})
            for path in outputs:
                self.assertTrue(path.exists())


if __name__ == "__main__":
    unittest.main()
