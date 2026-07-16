import pathlib
import subprocess
import sys
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

    def test_release_versions_are_consistent(self):
        data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        version = data["project"]["version"]
        self.assertEqual(version, data["tool"]["poetry"]["version"])
        namespace = {}
        exec((ROOT / "src/fp_tools/__init__.py").read_text(encoding="utf-8"), namespace)
        self.assertEqual(version, namespace["__version__"])
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertIn(f"version: {version}", citation)
        example = (ROOT / "examples/nutrient_stress_project/run_ctrl_vs_10fbs.sh").read_text(encoding="utf-8")
        self.assertIn(f'VERSION="${{FP_TOOLS_VERSION:-{version}}}"', example)

    def test_release_checklist_documents_required_gates(self):
        checklist = (ROOT / "RELEASE_CHECKLIST.md").read_text(encoding="utf-8")
        for phrase in [
            "pip check",
            "unittest discover",
            "call-footprints --help",
            "diff-footprints --help",
            "scripts/build_release.sh",
            "auditwheel",
            "manylinux",
            "twine check",
            "PyPI",
            "benchmarks/results/",
        ]:
            self.assertIn(phrase, checklist)
        self.assertIn("examples/gui_configs/call_footprints_single.yml", checklist)
        self.assertTrue((ROOT / "examples/gui_configs/call_footprints_single.yml").is_file())

    def test_release_script_can_clean_missing_dist_directory(self):
        script = (ROOT / "scripts" / "build_release.sh").read_text(encoding="utf-8")
        self.assertIn("mkdir -p dist", script)
        self.assertIn("find dist -maxdepth 1 -type f -delete", script)
        self.assertIn('export PATH="$ROOT/.venv/bin:$PATH"', script)
        self.assertIn("-m build", script)
        self.assertIn(".venv/bin/auditwheel", script)
        self.assertIn("patchelf is missing", script)
        self.assertIn('repair "$wheel"', script)
        self.assertIn("linux_x86_64", script)

    def test_build_cython_version_is_bounded(self):
        data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertIn("Cython>=3.0,<3.2", data["build-system"]["requires"])
        setup_py = (ROOT / "setup.py").read_text(encoding="utf-8")
        self.assertIn('build_dir="build/cythonized"', setup_py)


    def test_release_metadata_files_exist(self):
        for relative in ["LICENSE", "CITATION.cff", "environment.yml", "Dockerfile", "Makefile"]:
            self.assertTrue((ROOT / relative).exists(), f"Missing {relative}")
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertIn('repository-code: "https://github.com/oncologylab/fp-tools"', citation)
        self.assertIn("license: MIT", citation)

    def test_repository_guidance_is_tracked_without_redundant_gui_plan(self):
        self.assertTrue((ROOT / "AGENTS.md").is_file())
        self.assertTrue((ROOT / "DEV_PLAN.md").is_file())
        self.assertFalse((ROOT / "GUI_PLAN.md").exists())

    def test_publish_workflow_uses_trusted_publishing_and_repaired_wheels(self):
        workflow = (ROOT / ".github" / "workflows" / "publish.yml").read_text(encoding="utf-8")
        self.assertIn("auditwheel", workflow)
        self.assertIn("twine check", workflow)
        self.assertIn("pypa/gh-action-pypi-publish", workflow)
        self.assertIn("PYPI_API_TOKEN", workflow)
        self.assertNotIn("id-token: write", workflow)

    def test_console_script_smoke_helper_covers_declared_scripts(self):
        helper = ROOT / "scripts" / "smoke_console_scripts.py"
        self.assertTrue(helper.is_file())
        result = subprocess.run(
            [sys.executable, str(helper), "--bin-dir", str(pathlib.Path(sys.executable).parent)],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        for command in data["project"]["scripts"]:
            self.assertIn(f"OK\t{command}", result.stdout)

    def test_benchmark_manifest_schema_documentation_exists(self):
        manifest_doc = (ROOT / "benchmarks" / "manifests" / "README.md").read_text(encoding="utf-8")
        for column in ["source", "benchmark_tier", "experiment_accession", "file_accession", "checksum", "local_path"]:
            self.assertIn(f"`{column}`", manifest_doc)


if __name__ == "__main__":
    unittest.main()
