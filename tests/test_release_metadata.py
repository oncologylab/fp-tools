import pathlib
import tomllib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class ReleaseMetadataTest(unittest.TestCase):
    def test_project_urls_point_to_active_repository(self):
        data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        urls = data["project"]["urls"]
        self.assertEqual(urls["Homepage"], "https://github.com/oncologylab/fp-tools")
        self.assertEqual(urls["Repository"], "https://github.com/oncologylab/fp-tools")
        self.assertEqual(urls["Issues"], "https://github.com/oncologylab/fp-tools/issues")

    def test_release_checklist_documents_required_gates(self):
        checklist = (ROOT / "RELEASE_CHECKLIST.md").read_text(encoding="utf-8")
        for phrase in [
            "pip check",
            "unittest discover",
            "call-footprints --help",
            "diff-footprints --help",
            "scripts/build_release.sh",
            "twine check",
            "benchmarks/results/",
        ]:
            self.assertIn(phrase, checklist)

    def test_release_script_can_clean_missing_dist_directory(self):
        script = (ROOT / "scripts" / "build_release.sh").read_text(encoding="utf-8")
        self.assertIn("mkdir -p dist", script)
        self.assertIn("find dist -maxdepth 1 -type f -delete", script)
        self.assertIn("-m build", script)


    def test_release_metadata_files_exist(self):
        for relative in ["LICENSE", "CITATION.cff", ".zenodo.json", "environment.yml", "Dockerfile", "Makefile"]:
            self.assertTrue((ROOT / relative).exists(), f"Missing {relative}")
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertIn('repository-code: "https://github.com/oncologylab/fp-tools"', citation)
        self.assertIn("license: MIT", citation)

    def test_benchmark_manifest_schema_documentation_exists(self):
        manifest_doc = (ROOT / "benchmarks" / "manifests" / "README.md").read_text(encoding="utf-8")
        for column in ["source", "benchmark_tier", "experiment_accession", "file_accession", "checksum", "local_path"]:
            self.assertIn(f"`{column}`", manifest_doc)


if __name__ == "__main__":
    unittest.main()
