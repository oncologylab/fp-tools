import gzip
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = ROOT / "benchmarks/scripts/run_encode_cancer_pairwise_q95.py"
BUILDER_PATH = ROOT / "benchmarks/scripts/build_encode_cancer_q95_site.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUNNER = load_module("run_encode_cancer_pairwise_q95", RUNNER_PATH)
BUILDER = load_module("build_encode_cancer_q95_site", BUILDER_PATH)


class EncodeCancerBrowserTest(unittest.TestCase):
    def test_locked_design_has_17_replicates_30_peak_files_and_21_pairs(self):
        samples, peaks, comparisons = RUNNER.load_design()
        self.assertEqual(len(samples), 17)
        self.assertEqual(samples.groupby("condition")["sample"].count().to_dict(), RUNNER.EXPECTED_REPLICATES)
        self.assertEqual(
            samples.loc[samples.condition.eq("HepG2"), "biological_replicate"].tolist(),
            ["2", "3", "1"],
        )
        self.assertEqual(
            samples.loc[samples.condition.eq("K562"), "biological_replicate"].tolist(),
            ["3", "2", "1"],
        )
        self.assertEqual(len(peaks), 30)
        self.assertEqual(len(comparisons), 21)
        self.assertEqual(
            samples.loc[samples.condition.eq("K562"), "bam_accession"].tolist(),
            ["ENCFF077FBI", "ENCFF128WZG", "ENCFF534DCE"],
        )
        self.assertEqual(
            samples.loc[samples.condition.eq("HepG2"), "bam_accession"].tolist(),
            ["ENCFF624SON", "ENCFF926KFU", "ENCFF990VCP"],
        )
        self.assertEqual(set(peaks["output_type"]), {"IDR thresholded peaks"})

    def test_preserved_reference_payload_is_exact_and_complete(self):
        payload, digest = RUNNER.extract_payload(RUNNER.REFERENCE_REPORT)
        RUNNER.validate_payload(payload)
        self.assertEqual(digest, RUNNER.REFERENCE_JSON_SHA256)
        self.assertEqual(payload["conditions"], ["K562", "HepG2"])
        self.assertEqual(len(payload["points"]), 1019)
        self.assertEqual(len(payload["aggregate"]["motifs"]), 1009)
        samples = {
            sample["name"]
            for motif in payload["aggregate"]["motifs"]
            for condition in motif["conditions"]
            for sample in condition["samples"]
        }
        self.assertEqual(
            samples,
            {"K562_rep1", "K562_rep2", "K562_rep3", "HepG2_rep1", "HepG2_rep2", "HepG2_rep3"},
        )

    def test_reference_peak_universe_checksums_are_locked(self):
        spec = json.loads(
            (ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.spec.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            spec["reference_peak_universe"],
            {
                "comparison": "K562_vs_HepG2",
                "merged_bed_md5": "bc6ab7865e946f43509b36054902e45e",
                "background_50bp_bins_md5": "059d2d6bed8a62b9382c043125b5a620",
            },
        )

    def test_reference_default_panels_match_the_standalone_report(self):
        payload, _digest = RUNNER.extract_payload(RUNNER.REFERENCE_REPORT)
        with_profiles = {motif["prefix"] for motif in payload["aggregate"]["motifs"]}
        points = [point for point in payload["points"] if point["prefix"] in with_profiles]
        positive = sorted((point for point in points if point["change"] > 0), key=lambda point: (-point["change"], point["pvalue"]))[:2]
        negative = sorted((point for point in points if point["change"] < 0), key=lambda point: (point["change"], point["pvalue"]))[:2]
        self.assertEqual(
            [point["prefix"] for point in positive + negative],
            ["GATA2_MA0036.4", "GATA5_MA0766.3", "HNF4A_MA1494.2", "ONECUT2_MA0756.3"],
        )

    def test_reference_seed_strips_only_legacy_logo_blobs(self):
        source, _digest = RUNNER.extract_payload(RUNNER.REFERENCE_REPORT)
        with tempfile.TemporaryDirectory() as tmpdir:
            output = RUNNER.seed_reference(Path(tmpdir))
            compact = json.loads(gzip.decompress(output.read_bytes()))
        self.assertEqual(compact["logos"], {})
        for key in ("title", "report_label", "conditions", "groups", "colors", "points", "aggregate", "change_label"):
            self.assertEqual(compact[key], source[key])

    def test_static_browser_uses_canonical_gzip_payloads_and_helvetica(self):
        browser = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
        source = "\n".join((browser / name).read_text(encoding="utf-8") for name in ("index.html", "styles.css", "app.js"))
        self.assertIsNone(re.search(r"(?:src|href)=[\"']https?://", source))
        self.assertIn("DecompressionStream", source)
        self.assertIn("Helvetica,Arial,sans-serif", source.replace(" ", ""))
        self.assertIn('id="condition-1"', source)
        self.assertIn('id="condition-2"', source)
        self.assertIn('value="svg"', source)
        self.assertIn('value="png"', source)
        self.assertIn('value="pdf"', source)
        self.assertNotIn("data/profiles/", source)
        self.assertNotIn("data/motif_index.json", source)

    def test_static_reference_payload_matches_preserved_scientific_digest(self):
        browser = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting"
        payload = BUILDER.read_payload(browser / "data/reports/HepG2_vs_K562.json.gz")
        self.assertEqual(BUILDER.scientific_digest(payload), BUILDER.REFERENCE_SCIENTIFIC_SHA256)
        metadata = json.loads((browser / "data/metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(sum(len(item["samples"]) for item in metadata["conditions"]), 17)

    def test_each_static_payload_is_bounded(self):
        browser = ROOT / "docs/ENCODE-Cancer-Cell-lines-Footprinting/data/reports"
        for path in browser.glob("*.json.gz"):
            self.assertLess(path.stat().st_size, 10 * 1024 * 1024, path.name)

    def test_runner_is_pair_specific_q95_and_storage_bounded(self):
        source = RUNNER_PATH.read_text(encoding="utf-8")
        self.assertIn('"background-scale", "--stat", "q95"', source)
        self.assertIn('"--normalization", "none"', source)
        self.assertIn('"--aggregate-site-set", "bound"', source)
        self.assertIn('"--plot-aggregate", "all"', source)
        self.assertIn("shutil.rmtree(work)", source)
        self.assertIn('keep = {"report_payload.json.gz", "diff_footprints_results.txt"}', source)
        self.assertNotIn("unlink(missing_ok=True)  # downloaded BAM", source)


if __name__ == "__main__":
    unittest.main()
