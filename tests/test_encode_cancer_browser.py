import importlib.util
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks/scripts/build_encode_cancer_browser.py"
SPEC = importlib.util.spec_from_file_location("build_encode_cancer_browser", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class EncodeCancerBrowserTest(unittest.TestCase):
    def test_static_browser_is_dependency_free_and_uses_requested_font_stack(self):
        browser = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
        for name in ("index.html", "styles.css", "app.js"):
            self.assertTrue((browser / name).is_file(), name)
        source = "\n".join((browser / name).read_text(encoding="utf-8") for name in ("index.html", "styles.css", "app.js"))
        self.assertIsNone(re.search(r"(?:src|href)=[\"']https?://", source))
        self.assertIn("Arial,Helvetica,sans-serif", source)
        self.assertIn("Download all results", source)
        self.assertIn('value="svg"', source)
        self.assertIn('value="png"', source)
        self.assertIn('value="pdf"', source)
        self.assertIn('id="condition-1"', source)
        self.assertIn('id="condition-2"', source)
        self.assertIn('id="selected-grid"', source)
        self.assertIn('id="rank-chart"', source)
        self.assertIn('id="aggregate-grid"', source)
        self.assertNotIn('id="motif-table"', source)

    def test_browser_motif_matrices_match_every_report_motif(self):
        browser = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
        matrices = json.loads((browser / "data/motif_matrices.json").read_text(encoding="utf-8"))
        comparison = json.loads(
            (browser / "data/comparisons/K562_vs_MCF-7.json").read_text(encoding="utf-8")
        )
        matrix_map = matrices["motifs"]
        prefixes = {motif["prefix"] for motif in comparison["motifs"]}
        self.assertEqual(matrices["schema"], "fp-tools.motif-matrices.v1")
        self.assertEqual(len(matrix_map), 1019)
        self.assertEqual(set(matrix_map), prefixes)
        self.assertTrue(all(len(matrix) == 4 for matrix in matrix_map.values()))

    def test_static_browser_stays_within_documented_size_budget(self):
        browser = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
        size = sum(path.stat().st_size for path in browser.rglob("*") if path.is_file())
        self.assertLessEqual(size, 50 * 1024 * 1024)

    def test_design_is_exactly_seven_lines_15_replicates_and_21_pairs(self):
        manifest, comparisons, spec = MODULE.read_design(
            MODULE.DEFAULT_MANIFEST,
            MODULE.DEFAULT_SPEC,
            MODULE.DEFAULT_COMPARISONS,
        )
        self.assertEqual(len(manifest), 15)
        self.assertEqual(len(comparisons), 21)
        self.assertEqual(len(spec["conditions"]), 7)
        self.assertEqual(sum(len(samples) for samples in spec["conditions"].values()), 15)
        self.assertEqual(set(manifest.loc[manifest.condition.eq("A549"), "peak_accession"]), {"ENCFF876UEM"})

    def test_result_validation_requires_complete_empirical_bayes_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "result.tsv"
            table = self._result_table("A549", "HCT116")
            table.to_csv(path, sep="\t", index=False)
            observed = MODULE.validate_result(path, "A549", "HCT116")
            self.assertEqual(len(observed), 1019)
            table = table.drop(columns=["A549_HCT116_ebayes_qvalue_bh"])
            table.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "missing empirical-Bayes"):
                MODULE.validate_result(path, "A549", "HCT116")

    def test_result_validation_rejects_nonfinite_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "result.tsv"
            table = self._result_table("A549", "HCT116")
            table.loc[0, "A549_HCT116_ebayes_effect"] = np.nan
            table.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "non-finite effect"):
                MODULE.validate_result(path, "A549", "HCT116")

    def test_result_validation_accepts_positive_infinite_normal_limit_df(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "result.tsv"
            table = self._result_table("A549", "HCT116")
            table["A549_HCT116_ebayes_moderated_df"] = np.inf
            table.to_csv(path, sep="\t", index=False)
            observed = MODULE.validate_result(path, "A549", "HCT116")
            self.assertTrue(np.isposinf(observed["A549_HCT116_ebayes_moderated_df"]).all())
            table.loc[0, "A549_HCT116_ebayes_moderated_df"] = -np.inf
            table.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "invalid moderated_df"):
                MODULE.validate_result(path, "A549", "HCT116")

    def test_completed_replicate_matrix_recovers_result_without_site_rerun(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project = Path(tmpdir)
            output_dir = project / "comparisons/all_21_pairwise"
            summary_dir = project / "samples/a1/match_motifs"
            output_dir.mkdir(parents=True)
            summary_dir.mkdir(parents=True)
            matrix = pd.DataFrame({
                "motif": ["TF1_MA0001.1", "TF2_MA0002.1", "TF3_MA0003.1"],
                "n_sites": [100, 120, 140],
                "a1": [1.0, 2.0, 3.0],
                "a2": [1.2, 2.1, 2.8],
                "b1": [0.2, 2.0, 2.0],
                "b2": [0.4, 1.9, 2.2],
            })
            matrix.to_csv(
                output_dir / "diff_footprints_replicate_motif_score_matrix.tsv",
                sep="\t",
                index=False,
            )
            pd.DataFrame({
                "output_prefix": matrix["motif"],
                "name": ["TF1", "TF2", "TF3"],
                "motif_id": ["MA0001.1", "MA0002.1", "MA0003.1"],
                "cluster": ["C1", "C2", "C3"],
                "total_tfbs": matrix["n_sites"],
            }).to_csv(summary_dir / "motif_matches_results.txt", sep="\t", index=False)
            manifest = pd.DataFrame({
                "condition": ["A", "A", "B", "B"],
                "sample": ["a1", "a2", "b1", "b2"],
                "biological_replicate": ["1", "2", "1", "2"],
            })
            comparisons = pd.DataFrame({"comparison": ["A_vs_B"], "cond1": ["A"], "cond2": ["B"]})
            with mock.patch.object(MODULE, "EXPECTED_MOTIFS", 3):
                path = MODULE.recover_result_from_replicate_matrix(project, manifest, comparisons)
                recovered = MODULE.validate_result(path, "A", "B")
            self.assertEqual(len(recovered), 3)
            self.assertAlmostEqual(recovered.loc[0, "A_B_ebayes_effect"], 0.8)

    def test_comparison_column_contract_contains_directional_interval(self):
        columns = MODULE._comparison_columns("A549", "HCT116")
        self.assertEqual(columns["effect"], "A549_HCT116_ebayes_effect")
        self.assertEqual(columns["ci_lower"], "A549_HCT116_ebayes_ci_lower")
        self.assertEqual(columns["ci_upper"], "A549_HCT116_ebayes_ci_upper")

    def test_browser_record_labels_positive_infinite_df_as_normal_limit(self):
        row = SimpleNamespace(
            prefix="TF1_MA0001.1", name="TF1", motif_id="MA0001.1", cluster="C1",
            n_sites=100, mean1=0.2, sd1=0.01, mean2=0.1, sd2=0.02,
            effect=0.1, ci_lower=0.05, ci_upper=0.15, moderated_t=4.0,
            moderated_df=np.inf, pvalue=0.001, qvalue=0.01, significant=True,
        )
        record = MODULE._motif_record(row)
        self.assertTrue(record["normal_limit"])
        self.assertIsNone(record["moderated_df"])

    def test_tiled_profiles_reverse_negative_strand_and_preserve_total_counts(self):
        class FakeBigWig:
            def chroms(self):
                return {"chr1": 100}

            def values(self, chrom, start, end, numpy=True):
                self.assert_request = (chrom, start, end, numpy)
                return np.arange(start, end, dtype=float)

        sites = pd.DataFrame({
            "motif": ["TF_A", "TF_A", "TF_B", "TF_B", "TF_B"],
            "TFBS_chr": ["chr1"] * 5,
            "TFBS_start": [29, 69, 9, 19, 89],
            "TFBS_end": [31, 71, 11, 21, 91],
            "TFBS_strand": ["+", "-", "+", "+", "+"],
        })
        selected, totals = MODULE._select_profile_sites(sites, max_sites_per_motif=2)
        self.assertEqual(totals, {"TF_A": 2, "TF_B": 3})
        self.assertEqual(selected.groupby("motif").size().to_dict(), {"TF_A": 2, "TF_B": 2})
        sums, used = MODULE._profile_sums_by_motif(FakeBigWig(), selected, flank=2, tile_bp=100)
        np.testing.assert_array_equal(sums["TF_A"], np.array([99.0, 99.0, 99.0, 99.0]))
        self.assertEqual(used["TF_A"], 2)

    def test_shared_profile_site_selection_streams_and_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project = Path(tmpdir)
            cache = project / "samples/sample_1/match_motifs/cache/motif_sites.tsv.gz"
            summary = project / "samples/sample_1/match_motifs/motif_matches_results.txt"
            marker = project / "state/match_motifs.json"
            cache.parent.mkdir(parents=True)
            marker.parent.mkdir(parents=True)
            sites = pd.DataFrame({
                "motif": ["TF_A"] * 5 + ["TF_B"] * 3,
                "TFBS_chr": ["chr1"] * 8,
                "TFBS_start": list(range(10, 18)),
                "TFBS_end": list(range(11, 19)),
                "TFBS_name": ["site"] * 8,
                "TFBS_score": [1.0] * 8,
                "TFBS_strand": ["+"] * 8,
                "peak_chr": ["chr1"] * 8,
                "peak_start": [0] * 8,
                "peak_end": [100] * 8,
                "score": [0.5] * 8,
            })
            sites.to_csv(cache, sep="\t", index=False, compression="gzip")
            pd.DataFrame({
                "output_prefix": ["TF_A", "TF_B"],
                "total_tfbs": [5, 3],
            }).to_csv(summary, sep="\t", index=False)
            marker.write_text(json.dumps({"verified": True}), encoding="utf-8")
            manifest = pd.DataFrame({"sample": ["sample_1"]})
            with mock.patch.object(MODULE, "EXPECTED_MOTIFS", 2):
                selected_path, counts = MODULE.prepare_shared_profile_sites(project, manifest, 2)
                selected = pd.read_csv(selected_path, sep="\t")
                cached_path, cached_counts = MODULE.prepare_shared_profile_sites(project, manifest, 2)
            self.assertEqual(counts, {"TF_A": 5, "TF_B": 3})
            self.assertEqual(cached_counts, counts)
            self.assertEqual(cached_path, selected_path)
            self.assertEqual(selected.groupby("motif").size().to_dict(), {"TF_A": 2, "TF_B": 2})
            self.assertEqual(selected.loc[selected.motif.eq("TF_A"), "TFBS_start"].tolist(), [10, 14])
            self.assertEqual(selected.loc[selected.motif.eq("TF_B"), "TFBS_start"].tolist(), [15, 17])

    def test_all_pairs_run_in_one_15_replicate_analysis(self):
        manifest, comparisons, _spec = MODULE.read_design(
            MODULE.DEFAULT_MANIFEST,
            MODULE.DEFAULT_SPEC,
            MODULE.DEFAULT_COMPARISONS,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            project = Path(tmpdir) / "project"
            peaks = project / "peaks/merged_peaks_filtered.bed"
            genome = Path(tmpdir) / "hg38.fa"
            peaks.parent.mkdir(parents=True)
            peaks.touch()
            genome.touch()
            with (
                mock.patch.object(MODULE.subprocess, "run", return_value=mock.Mock(returncode=0)) as run,
                mock.patch.object(MODULE, "validate_result", return_value=pd.DataFrame()) as validate,
            ):
                MODULE.run_comparisons(project, manifest, comparisons, genome, cores=8)
            self.assertEqual(run.call_count, 1)
            command = run.call_args.args[0]
            sample_start = command.index("--sample-dirs") + 1
            sample_end = command.index("--sample-names")
            condition_start = command.index("--cond-names") + 1
            condition_end = command.index("--peaks")
            self.assertEqual(len(command[sample_start:sample_end]), 15)
            self.assertEqual(len(command[condition_start:condition_end]), 15)
            self.assertEqual(command[command.index("--replicate-report") + 1], "off")
            self.assertEqual(command[command.index("--normalization") + 1], "sample-quantile")
            self.assertEqual(validate.call_count, 21)

    @staticmethod
    def _result_table(cond1, cond2):
        size = MODULE.EXPECTED_MOTIFS
        base = f"{cond1}_{cond2}_ebayes"
        return pd.DataFrame({
            "output_prefix": [f"TF{i}_MA{i:04d}.1" for i in range(size)],
            "name": [f"TF{i}" for i in range(size)],
            "motif_id": [f"MA{i:04d}.1" for i in range(size)],
            "cluster": ["C_TEST"] * size,
            "total_tfbs": [100] * size,
            f"{cond1}_mean_score": np.linspace(0.1, 1.0, size),
            f"{cond1}_score_sd": [0.1] * size,
            f"{cond2}_mean_score": np.linspace(0.2, 1.1, size),
            f"{cond2}_score_sd": [0.1] * size,
            f"{base}_effect": np.linspace(-0.5, 0.5, size),
            f"{base}_ci_lower": np.linspace(-0.6, 0.4, size),
            f"{base}_ci_upper": np.linspace(-0.4, 0.6, size),
            f"{base}_moderated_t": np.linspace(-4, 4, size),
            f"{base}_moderated_df": [10.0] * size,
            f"{base}_pvalue": np.linspace(0.001, 0.9, size),
            f"{base}_qvalue_bh": np.linspace(0.01, 0.95, size),
            f"{base}_significant_fdr05": [True] * 44 + [False] * (size - 44),
        })


if __name__ == "__main__":
    unittest.main()
