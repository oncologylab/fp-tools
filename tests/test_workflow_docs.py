from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_normal_bulk_examples_use_automatic_cores():
    for name in ("README.md", "docs/get-started/commands/bulk-footprinting.md",
                 "docs/get-started/workflows/bulk-atac-seq.md"):
        text = (ROOT / name).read_text()
        assert "--cores 8" not in text
        assert "--cores 36" not in text
        assert "all available cores" in text


def test_packaged_bulk_example_uses_managed_reference_and_auto_cores():
    import yaml
    for name in ("examples/gui_configs/bulk_footprinting_bam.yml",
                 "src/fp_tools/resources/gui_configs/bulk_footprinting_bam.yml"):
        job = yaml.safe_load((ROOT / name).read_text())["samples"][0]
        assert job["genome"] == "hg38"
        assert not job.get("blacklist")
        assert job.get("cores") is None
