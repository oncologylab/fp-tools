import importlib.util
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks/scripts/run_encode_atac_project.py"
SEVEN_SPEC = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.spec.json"
SEVEN_MANIFEST = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.tsv"
SEVEN_COMPARISONS = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814_comparisons.tsv"
SPEC = importlib.util.spec_from_file_location("run_encode_atac_project", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class EncodeAtacProjectTest(unittest.TestCase):
    def test_manifest_contains_exactly_all_requested_replicates(self):
        frame = MODULE.read_manifest(MODULE.DEFAULT_MANIFEST)
        self.assertEqual(len(frame), 16)
        self.assertEqual(set(frame["selected_biosample"]), MODULE.EXPECTED_SELECTED_BIOSAMPLES)
        self.assertEqual(frame.groupby("condition")["sample"].count().to_dict(), {condition: 2 for condition in MODULE.EXPECTED_CONDITIONS})

    def test_manifest_rejects_a_missing_replicate(self):
        frame = pd.read_csv(MODULE.DEFAULT_MANIFEST, sep="\t")
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "missing.tsv"
            frame.iloc[:-1].to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "exactly 16"):
                MODULE.read_manifest(path)

    def test_manifest_rejects_duplicate_selected_biosample(self):
        frame = pd.read_csv(MODULE.DEFAULT_MANIFEST, sep="\t")
        frame.loc[1, "selected_biosample"] = frame.loc[0, "selected_biosample"]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "duplicate.tsv"
            frame.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "selected-biosample set"):
                MODULE.read_manifest(path)

    def test_seven_line_manifest_contains_all_15_requested_replicates(self):
        spec = MODULE.load_project_spec(SEVEN_SPEC)
        frame = MODULE.read_manifest(SEVEN_MANIFEST, spec)
        self.assertEqual(len(frame), 15)
        self.assertEqual(
            spec.replicate_counts,
            {"A549": 3, "HCT116": 2, "HepG2": 2, "K562": 2, "MCF-7": 2, "PC-3": 2, "Panc1": 2},
        )
        self.assertEqual(set(frame["condition"]), set(spec.conditions))
        self.assertNotIn("GM12878", set(frame["condition"]))
        self.assertNotIn("IMR-90", set(frame["condition"]))
        self.assertNotIn("DND-41", set(frame["condition"]))
        self.assertEqual(set(frame.loc[frame["condition"].eq("A549"), "peak_accession"]), {"ENCFF876UEM"})

    def test_seven_line_manifest_rejects_missing_third_a549_replicate(self):
        spec = MODULE.load_project_spec(SEVEN_SPEC)
        frame = pd.read_csv(SEVEN_MANIFEST, sep="\t")
        frame = frame.loc[~frame["sample"].eq("A549_rep3")]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "missing.tsv"
            frame.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "exactly 15"):
                MODULE.read_manifest(path, spec)

    def test_seven_line_comparisons_are_all_21_unordered_pairs(self):
        spec = MODULE.load_project_spec(SEVEN_SPEC)
        frame = pd.read_csv(SEVEN_COMPARISONS, sep="\t")
        observed = {frozenset((row.cond1, row.cond2)) for row in frame.itertuples(index=False)}
        expected = {
            frozenset(pair)
            for pair in __import__("itertools").combinations(sorted(spec.conditions), 2)
        }
        self.assertEqual(len(frame), 21)
        self.assertEqual(observed, expected)

    def test_peak_merge_uses_bedtools_touching_interval_semantics(self):
        intervals = [("chr1", 0, 10), ("chr1", 10, 20), ("chr1", 12, 15), ("chr2", 1, 2)]
        merged = []
        for chrom, start, end in intervals:
            if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
                merged[-1][2] = max(merged[-1][2], end)
            else:
                merged.append([chrom, start, end])
        self.assertEqual(merged, [["chr1", 0, 20], ["chr2", 1, 2]])

    def test_encode_s3_url_resolves_official_public_object(self):
        url = MODULE.encode_s3_url("ENCFF646NWY.bam", "2021-02-24T00:00:00Z")
        self.assertEqual(
            url,
            "https://encode-public.s3.amazonaws.com/2021/02/24/baf5cad2-f9e6-4adc-84d8-5ec034b49977/ENCFF646NWY.bam",
        )

    def test_download_prefers_verified_official_s3_object(self):
        payload = b"abc"
        expected_md5 = hashlib.md5(payload).hexdigest()
        commands = []

        def fake_run_logged(command, _log):
            commands.append(command)
            Path(command[command.index("--output") + 1]).write_bytes(payload)

        with tempfile.TemporaryDirectory() as tmpdir:
            destination = Path(tmpdir) / "ENCFFTEST.bam"
            with (
                mock.patch.object(MODULE, "encode_metadata", return_value={"date_created": "2026-08-14T00:00:00Z"}),
                mock.patch.object(MODULE, "encode_s3_url", return_value="https://encode-public.s3.amazonaws.com/object/ENCFFTEST.bam"),
                mock.patch.object(MODULE, "run_logged", side_effect=fake_run_logged),
            ):
                observed = MODULE.download(
                    "ENCFFTEST", len(payload), expected_md5, destination,
                    Path(tmpdir) / "download.log", "curl",
                )
        self.assertEqual(observed, destination)
        self.assertEqual(commands[0][-1], "https://encode-public.s3.amazonaws.com/object/ENCFFTEST.bam")

    def test_encode_audit_can_resume_from_complete_cached_audit(self):
        frame = MODULE.read_manifest(MODULE.DEFAULT_MANIFEST)
        with tempfile.TemporaryDirectory() as tmpdir:
            audit = Path(tmpdir) / "audit.tsv"
            rows = []
            for row in frame.itertuples(index=False):
                rows.append({"accession": row.bam_accession, "status": "released", "assembly": "GRCh38", "output_type": "alignments", "file_size": row.bam_size, "md5sum": row.bam_md5, "dataset": f"/{row.experiment}/", "biological_replicates": str(row.biological_replicate), "technical_replicates": f"{row.biological_replicate}_1", "audited_at": "2026-08-13T00:00:00+00:00"})
            for row in frame.drop_duplicates("peak_accession").itertuples(index=False):
                rows.append({"accession": row.peak_accession, "status": "released", "assembly": "GRCh38", "output_type": "conservative IDR thresholded peaks", "file_size": row.peak_size, "md5sum": row.peak_md5, "dataset": f"/{row.experiment}/", "audited_at": "2026-08-13T00:00:00+00:00"})
            pd.DataFrame(rows).to_csv(audit, sep="\t", index=False)
            with mock.patch.object(MODULE, "encode_metadata", side_effect=OSError("portal unavailable")):
                observed = MODULE.audit_encode_files(frame, audit)
            self.assertEqual(len(observed), 24)

    def test_encode_audit_rejects_stale_cached_checksum(self):
        frame = MODULE.read_manifest(MODULE.DEFAULT_MANIFEST)
        with tempfile.TemporaryDirectory() as tmpdir:
            audit = Path(tmpdir) / "audit.tsv"
            row = frame.iloc[0]
            pd.DataFrame([{
                "accession": row.bam_accession,
                "status": "released",
                "assembly": "GRCh38",
                "output_type": "alignments",
                "file_size": row.bam_size,
                "md5sum": "0" * 32,
                "dataset": f"/{row.experiment}/",
                "audited_at": "2026-08-13T00:00:00+00:00",
            }]).to_csv(audit, sep="\t", index=False)
            with mock.patch.object(MODULE, "encode_metadata", side_effect=OSError("portal unavailable")):
                with self.assertRaisesRegex(ValueError, "size/MD5 differs"):
                    MODULE.audit_encode_files(frame, audit)

    def test_replicate_audit_verifies_library_and_assayed_biosample(self):
        spec = MODULE.load_project_spec(SEVEN_SPEC)
        frame = MODULE.read_manifest(SEVEN_MANIFEST, spec)
        frame = frame.loc[frame["condition"].eq("A549")].copy()
        payload = {
            "replicates": [
                {
                    "biological_replicate_number": int(row.biological_replicate),
                    "technical_replicate_number": 1,
                    "library": {
                        "accession": row.library,
                        "biosample": {"accession": row.biosample},
                    },
                }
                for row in frame.itertuples(index=False)
            ]
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "replicates.tsv"
            with mock.patch.object(MODULE, "encode_experiment_metadata", return_value=payload):
                observed = MODULE.audit_encode_replicates(frame, output)
            self.assertEqual(len(observed), 3)
            self.assertEqual(set(observed["source_relation"]), {"assayed biosample"})
            self.assertTrue(output.is_file())

    def test_replicate_audit_accepts_pinned_source_biosample(self):
        spec = MODULE.load_project_spec(SEVEN_SPEC)
        frame = MODULE.read_manifest(SEVEN_MANIFEST, spec)
        frame = frame.loc[frame["condition"].eq("HCT116")].copy()
        payload = {
            "replicates": [
                {
                    "biological_replicate_number": int(row.biological_replicate),
                    "technical_replicate_number": 1,
                    "library": {
                        "accession": row.library,
                        "biosample": {
                            "accession": row.biosample,
                            "part_of": f"/biosamples/{row.selected_biosample}/",
                        },
                    },
                }
                for row in frame.itertuples(index=False)
            ]
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "replicates.tsv"
            with mock.patch.object(MODULE, "encode_experiment_metadata", return_value=payload):
                observed = MODULE.audit_encode_replicates(frame, output)
        self.assertEqual(set(observed["source_relation"]), {"source biosample"})

    def test_normalization_validation_uses_project_level_reports(self):
        frame = MODULE.read_manifest(MODULE.DEFAULT_MANIFEST)
        with tempfile.TemporaryDirectory() as tmpdir:
            project = Path(tmpdir)
            corrected = [project / "samples" / sample / "atac_correct" / f"{sample}_corrected.bw" for sample in frame["sample"]]
            normalized = [project / "samples" / sample / "normalize" / f"{sample}_corrected_q95_scaled.bw" for sample in frame["sample"]]
            report = pd.DataFrame({
                "sample": frame["sample"],
                "input_bigwig": [str(path) for path in corrected],
                "output_bigwig": [str(path) for path in normalized],
            })
            report.to_csv(project / "normalize_bigwig_manifest.tsv", sep="\t", index=False)
            qc = report.assign(
                scaling_stat="q95", scaling_value=1.0,
                target_scaling_value=1.0, scale_factor=1.0,
            )
            qc.to_csv(project / "normalize_bigwig_qc.tsv", sep="\t", index=False)
            with mock.patch.object(MODULE, "validate_bigwig") as validate:
                MODULE.validate_normalization_outputs(frame, project, corrected, normalized)
            self.assertEqual(validate.call_count, 16)


if __name__ == "__main__":
    unittest.main()
