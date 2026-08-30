import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd


ROOT = pathlib.Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


study = load_module(
    "validate_footprint_study",
    "benchmarks/scripts/validate_footprint_study.py",
)
labels = load_module(
    "build_footprint_site_labels",
    "benchmarks/scripts/build_footprint_site_labels.py",
)
downsample = load_module(
    "downsample_bam_fragments",
    "benchmarks/scripts/downsample_bam_fragments.py",
)
ablation = load_module(
    "build_footprint_ablation_plan",
    "benchmarks/scripts/build_footprint_ablation_plan.py",
)
ablation_summary = load_module(
    "summarize_footprint_ablation",
    "benchmarks/scripts/summarize_footprint_ablation.py",
)
promotion = load_module(
    "evaluate_footprint_promotion",
    "benchmarks/scripts/evaluate_footprint_promotion.py",
)
nutrient = load_module(
    "evaluate_nutrient_footprint_replication",
    "benchmarks/scripts/evaluate_nutrient_footprint_replication.py",
)
ablation_runner = load_module(
    "run_footprint_ablation_plan",
    "benchmarks/scripts/run_footprint_ablation_plan.py",
)
evidence_fusion = load_module(
    "evaluate_site_evidence_fusion",
    "benchmarks/scripts/evaluate_site_evidence_fusion.py",
)


class FootprintStudySpecTest(unittest.TestCase):
    def test_committed_study_spec_is_valid_and_locked(self):
        path = ROOT / "benchmarks/manifests/footprint_detectability_v1.spec.json"
        spec = study.load_spec(path)
        self.assertEqual(study.validate_spec(spec), [])
        self.assertEqual(len(spec["tasks"]), 35)
        self.assertEqual(
            {task["cell"] for task in spec["tasks"] if task["split"] == "development"},
            {"K562", "HepG2"},
        )
        self.assertEqual(spec["nutrient_application"]["external_pdac_accession"], "GSE144833")

    def test_spec_rejects_chromosome_leakage_and_duplicate_tasks(self):
        path = ROOT / "benchmarks/manifests/footprint_detectability_v1.spec.json"
        spec = json.loads(path.read_text(encoding="utf-8"))
        spec["chromosome_split"]["validation"].append("chr1")
        spec["tasks"].append(dict(spec["tasks"][0]))
        errors = study.validate_spec(spec)
        self.assertTrue(any("multiple splits" in error for error in errors))
        self.assertTrue(any("duplicate task" in error for error in errors))


class SiteEvidenceFusionStudyTest(unittest.TestCase):
    def test_locked_holdout_requires_explicit_unlock(self):
        spec = {
            "chromosome_split": {
                "train": ["chr1"],
                "validation": ["chr17"],
                "test": ["chr19"],
            },
            "tasks": [
                {"cell": "K562", "tf": "CTCF", "split": "development", "role": "positive_control", "motif_family": "CTCF"},
                {"cell": "A549", "tf": "CTCF", "split": "locked_holdout", "role": "positive_control", "motif_family": "CTCF"},
            ],
        }
        rows = []
        for cell in ("K562", "A549"):
            for chrom in ("chr17", "chr19"):
                for index, label in enumerate((0, 1)):
                    rows.append(
                        {
                            "cell": cell,
                            "tf": "CTCF",
                            "TFBS_chr": chrom,
                            "TFBS_start": index,
                            "TFBS_end": index + 1,
                            "chip_label": label,
                            "footprint_score": float(index),
                            "pwm_score": float(index),
                        }
                    )
        frame = evidence_fusion.add_candidate_scores(
            evidence_fusion.attach_design(pd.DataFrame(rows), spec)
        )
        validation = evidence_fusion.select_evaluation_rows(
            frame,
            unlock_development_test=False,
            unlock_holdout=False,
        )
        self.assertEqual(set(validation["cell"]), {"K562"})
        self.assertEqual(set(validation["TFBS_chr"]), {"chr17"})
        unlocked = evidence_fusion.select_evaluation_rows(
            frame,
            unlock_development_test=True,
            unlock_holdout=True,
        )
        self.assertEqual(set(unlocked["cell"]), {"K562", "A549"})
        self.assertEqual(
            set(unlocked.loc[unlocked["cell"] == "A549", "TFBS_chr"]),
            {"chr19"},
        )

    def test_paired_chromosome_bootstrap_reports_candidate_probability(self):
        rows = []
        for chrom in ("chr17", "chr18"):
            for index in range(20):
                label = int(index >= 10)
                rows.append(
                    {
                        "cell_split": "development",
                        "chromosome_split": "validation",
                        "cell": "K562",
                        "tf": "CTCF",
                        "role": "positive_control",
                        "motif_family": "CTCF",
                        "TFBS_chr": chrom,
                        "chip_label": label,
                        "footprint_score": float(index % 3),
                        "evidence_fusion_score": float(index),
                    }
                )
        result = evidence_fusion.paired_chromosome_bootstrap(
            pd.DataFrame(rows), n_bootstrap=20, seed=4
        )
        self.assertEqual(set(result["metric"]), {"auroc", "auprc"})
        self.assertTrue((result["successful_bootstraps"] == 20).all())
        self.assertTrue((result["probability_delta_gt_zero"] == 1.0).all())


class FootprintSiteLabelTest(unittest.TestCase):
    def test_summit_supported_positive_far_negative_and_indeterminate(self):
        peaks = {"chr1": [labels.Peak(100, 200, 150)]}
        sites = pd.DataFrame(
            [
                {"chrom": "chr1", "start": 145, "end": 155, "strand": "+", "site_id": "p", "motif_score": 8.0},
                {"chrom": "chr1", "start": 105, "end": 115, "strand": "+", "site_id": "inside", "motif_score": 7.0},
                {"chrom": "chr1", "start": 250, "end": 260, "strand": "+", "site_id": "near", "motif_score": 6.0},
                {"chrom": "chr1", "start": 800, "end": 810, "strand": "+", "site_id": "n", "motif_score": 5.0},
                {"chrom": "chr2", "start": 10, "end": 20, "strand": "+", "site_id": "empty", "motif_score": 4.0},
            ]
        )
        output = labels.label_sites(
            sites, peaks, positive_summit_distance=20, negative_peak_distance=100
        ).set_index("site_id")
        self.assertEqual(int(output.loc["p", "label"]), 1)
        self.assertEqual(int(output.loc["inside", "label"]), -1)
        self.assertEqual(int(output.loc["near", "label"]), -1)
        self.assertEqual(int(output.loc["n", "label"]), 0)
        self.assertEqual(int(output.loc["empty", "label"]), 0)

    def test_reads_narrowpeak_summit_offset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = pathlib.Path(tmpdir) / "labels.narrowPeak"
            path.write_text(
                "chr1\t100\t200\tpeak\t1000\t.\t5\t10\t8\t45\n",
                encoding="utf-8",
            )
            peaks = labels.read_peaks(path)
            self.assertEqual(peaks["chr1"], [labels.Peak(100, 200, 145)])

    def test_propensity_matching_is_deterministic_and_balanced(self):
        rng = np.random.default_rng(4)
        frame = pd.DataFrame(
            {
                "chrom": ["chr1"] * 30,
                "start": np.arange(30) * 10,
                "end": np.arange(30) * 10 + 5,
                "strand": ["+"] * 30,
                "site_id": [f"s{i}" for i in range(30)],
                "motif_score": np.r_[rng.normal(2, 0.2, 10), rng.normal(2, 0.5, 20)],
                "accessibility": np.r_[rng.normal(5, 0.2, 10), rng.normal(5, 0.5, 20)],
                "label": [1] * 10 + [0] * 20,
                "label_reason": ["positive"] * 10 + ["negative"] * 20,
                "nearest_peak_distance": [0] * 10 + [1000] * 20,
                "nearest_summit_distance": [0] * 10 + [1000] * 20,
            }
        )
        first = labels.propensity_match(
            frame, ["motif_score", "accessibility"], negative_ratio=1, seed=7
        )
        second = labels.propensity_match(
            frame, ["motif_score", "accessibility"], negative_ratio=1, seed=7
        )
        self.assertEqual(first["label"].value_counts().to_dict(), {1: 10, 0: 10})
        pd.testing.assert_frame_equal(first, second)


class FootprintAblationPlanTest(unittest.TestCase):
    def test_fragment_hash_is_deterministic_and_depth_subsets_are_nested(self):
        first = [
            name
            for name in (f"fragment-{index}" for index in range(1000))
            if downsample.fragment_uniform(name, 9) < 0.2
        ]
        second = [
            name
            for name in (f"fragment-{index}" for index in range(1000))
            if downsample.fragment_uniform(name, 9) < 0.5
        ]
        repeated = [
            name
            for name in (f"fragment-{index}" for index in range(1000))
            if downsample.fragment_uniform(name, 9) < 0.2
        ]
        self.assertEqual(first, repeated)
        self.assertTrue(set(first).issubset(second))
        self.assertNotEqual(
            downsample.fragment_uniform("fragment-1", 9),
            downsample.fragment_uniform("fragment-1", 10),
        )

    def test_ablation_plan_reuses_each_depth_subset_across_corrections(self):
        spec = {
            "random_seed": 20,
            "depth_randomizations": 2,
            "depth_fragments": [10_000_000, "full"],
            "native_corrections": [
                "raw",
                "fp_tools_pwm",
                "fp_tools_dwm",
                "fp_tools_reused_bias",
            ],
            "whole_methods": ["fp_tools", "TOBIAS"],
            "tasks": [
                {
                    "cell": "K562",
                    "tf": "CTCF",
                    "motif_id": "MA0139.2",
                    "motif_family": "CTCF",
                    "role": "positive_control",
                    "split": "development",
                }
            ],
        }
        samples = pd.DataFrame(
            [
                {
                    "sample": "K562_rep1",
                    "cell": "K562",
                    "bam": "/inputs/a sample.bam",
                    "peaks": "/inputs/peaks.bed",
                    "fragments": 100_000_000,
                }
            ]
        )
        signal_plan = ablation.build_signal_plan(
            spec,
            samples,
            pathlib.Path("/ref/genome.fa"),
            pathlib.Path("/ref/blacklist.bed"),
            pathlib.Path("/results"),
            python="python",
            cores=2,
        )
        self.assertEqual(len(signal_plan), 13)
        downsample_jobs = signal_plan[signal_plan["stage"] == "downsample"]
        self.assertEqual(len(downsample_jobs), 2)
        subset = signal_plan[(signal_plan["depth"] == 10_000_000) & (signal_plan["seed"] == 20)]
        self.assertEqual(set(subset["correction"]), {"", "raw", "fp_tools_pwm", "fp_tools_dwm", "fp_tools_reused_bias"})
        correction_jobs = subset[subset["stage"] == "correction"]
        self.assertTrue(all("downsample:K562_rep1.10m.s20" in value for value in correction_jobs["depends_on"]))
        self.assertIn("'/inputs/a sample.bam'", downsample_jobs.iloc[0]["command"])

        evaluation = ablation.build_evaluation_plan(spec, signal_plan)
        self.assertEqual(len(evaluation), 22)
        self.assertEqual(set(evaluation["split"]), {"development"})

    def test_ablation_summary_detects_plateau_and_correction_gain(self):
        rows = []
        for correction, values in {
            "raw": [0.60, 0.61, 0.61],
            "fp_tools_dwm": [0.68, 0.70, 0.705],
        }.items():
            for depth, auroc in zip([10_000_000, 50_000_000, "full"], values):
                for seed in [1, 2]:
                    rows.append(
                        {
                            "cell": "K562",
                            "tf": "CTCF",
                            "motif_id": "MA0139.2",
                            "correction": correction,
                            "method": "fp_tools",
                            "depth": depth,
                            "seed": seed,
                            "auroc": auroc,
                            "auprc": auroc - 0.1,
                        }
                    )
        aggregated = ablation_summary.aggregate_randomizations(pd.DataFrame(rows))
        depth = ablation_summary.depth_diagnostics(aggregated, plateau_delta=0.01)
        dwm = depth[depth["correction"] == "fp_tools_dwm"].iloc[0]
        self.assertTrue(bool(dwm["depth_plateau"]))
        correction = ablation_summary.correction_diagnostics(aggregated)
        self.assertAlmostEqual(float(correction.iloc[0]["correction_gain_auroc"]), 0.095)

    def test_promotion_gate_uses_locked_roles_and_clustered_bootstrap(self):
        tasks = []
        rows = []
        for index, (cell, family, role) in enumerate(
            [
                ("A549", "CTCF", "positive_control"),
                ("A549", "FOX", "difficult"),
                ("MCF-7", "GATA", "difficult"),
                ("MCF-7", "TCF7", "transfer"),
            ]
        ):
            tf = f"TF{index}"
            motif_id = f"MA{index}"
            tasks.append(
                {
                    "cell": cell,
                    "tf": tf,
                    "motif_id": motif_id,
                    "motif_family": family,
                    "role": role,
                    "split": "locked_holdout",
                }
            )
            rows.extend(
                [
                    {"cell": cell, "tf": tf, "motif_id": motif_id, "method": "baseline", "auroc": 0.70, "auprc": 0.40},
                    {"cell": cell, "tf": tf, "motif_id": motif_id, "method": "candidate", "auroc": 0.75, "auprc": 0.46},
                ]
            )
        gates = {
            "minimum_relative_auprc_gain": 0.10,
            "minimum_auroc_gain": 0.03,
            "maximum_strong_control_auroc_loss": 0.02,
            "minimum_difficult_tf_families_improved": 2,
            "minimum_holdout_cells_improved": 2,
            "minimum_detectability_probability": 0.80,
            "maximum_naked_dna_false_positive_rate": 0.05,
            "maximum_naked_dna_false_positive_rate_increase": 0.01,
        }
        study_spec = {"tasks": tasks, "promotion_gates": gates}
        pairs = promotion.prepare_pairs(
            pd.DataFrame(rows), study_spec, "candidate", "baseline", "locked_holdout"
        )
        paired, summary = promotion.evaluate_gates(
            pairs,
            gates,
            n_bootstrap=100,
            seed=12,
            negative_controls=pd.DataFrame(
                {
                    "method": ["baseline", "candidate"],
                    "false_positive_rate": [0.04, 0.045],
                }
            ),
        )
        self.assertTrue(summary["passed"])
        self.assertTrue(paired["improved"].all())
        self.assertEqual(summary["improved_cells"], 2)

    def test_nutrient_candidate_requires_local_rna_recovery_and_occupancy(self):
        rows = [
            {
                "cohort": "local",
                "cell": cell,
                "contrast_class": "stress",
                "motif_id": "MA0833.3",
                "tf": "ATF4",
                "delta_footprint": delta,
                "fdr": 0.01,
                "rna_log2fc": 1.2,
                "rna_fdr": 0.01,
            }
            for cell, delta in [("AsPC-1", 1.0), ("HPAF-II", 0.8), ("Panc1", 0.7)]
        ]
        rows.extend(
            [
                {"cohort": "external_pdac", "cell": "SUIT-2", "contrast_class": "stress", "motif_id": "MA0833.3", "tf": "ATF4", "delta_footprint": 0.6, "fdr": 0.01},
                {"cohort": "external_pdac", "cell": "SUIT-2", "contrast_class": "recovery", "motif_id": "MA0833.3", "tf": "ATF4", "delta_footprint": -0.5, "fdr": 0.01},
                {"cohort": "external_mechanistic", "cell": "THP1", "contrast_class": "occupancy", "motif_id": "MA0833.3", "tf": "ATF4", "delta_footprint": np.nan, "fdr": np.nan, "occupancy_log2fc": 1.1, "occupancy_fdr": 0.01},
            ]
        )
        rules = {
            "required_directional_cell_lines": 3,
            "top_fraction_absolute_change": 1.0,
            "minimum_top_cell_lines": 2,
            "minimum_rna_concordant_cell_lines": 2,
            "minimum_external_reversal_fraction": 0.5,
            "external_fdr": 0.05,
        }
        result = nutrient.evaluate_replication(pd.DataFrame(rows), rules)
        self.assertEqual(result.iloc[0]["evidence_tier"], "mechanism_supported")
        self.assertTrue(bool(result.iloc[0]["mechanism_pass"]))

    def test_ablation_runner_resolves_dependencies_in_dry_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            plan = pd.DataFrame(
                [
                    {"job_id": "second", "stage": "correction", "depends_on": "first", "expected_output": f"{tmpdir}/second", "command": "tool second"},
                    {"job_id": "first", "stage": "downsample", "depends_on": "", "expected_output": f"{tmpdir}/first", "command": "tool first"},
                ]
            )
            status = pathlib.Path(tmpdir) / "status.tsv"
            result = ablation_runner.execute_plan(plan, status, dry_run=True)
            self.assertEqual(result["job_id"].tolist(), ["first", "second"])
            self.assertTrue(status.exists())


if __name__ == "__main__":
    unittest.main()
