import importlib.util
import pathlib
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


prepare_biomed = load_module("prepare_biomedinformatics_template", ROOT / "paper" / "scripts" / "prepare_biomedinformatics_template.py")


class PaperAssetsTest(unittest.TestCase):
    def test_paper_templates_are_present(self):
        for rel_path in [
            "paper/MDPI_template_ACS.zip",
            "paper/biomedinformatics-template.dot",
            "paper/BioMedInformatics",
        ]:
            self.assertTrue((ROOT / rel_path).exists(), rel_path)

    def test_mdpi_latex_template_zip_has_required_files(self):
        names = prepare_biomed.validate_template_zip(ROOT / "paper" / "MDPI_template_ACS.zip")
        for expected in prepare_biomed.REQUIRED_TEMPLATE_FILES:
            self.assertIn(expected, names)

    def test_prepare_template_creates_manuscript_workspace(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = pathlib.Path(tmpdir) / "manuscript"
            prepare_biomed.prepare_template(ROOT / "paper" / "MDPI_template_ACS.zip", outdir)
            self.assertTrue((outdir / "main.tex").exists())
            self.assertTrue((outdir / "Definitions" / "mdpi.cls").exists())
            self.assertTrue((outdir / "figures").is_dir())
            self.assertTrue((outdir / "tables").is_dir())
            self.assertTrue((outdir / "data_availability.md").exists())


if __name__ == "__main__":
    unittest.main()
