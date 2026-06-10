import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.pseudobulk import group_fragments, write_downstream_commands


class PseudobulkTest(unittest.TestCase):
    def test_group_fragments_writes_manifest_and_group_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fragments = tmp / "fragments.tsv"
            annotations = tmp / "annotations.tsv"
            outdir = tmp / "out"
            fragments.write_text(
                "chr1\t1\t5\tcellA\t2\nchr1\t10\t20\tcellB\t3\nchr1\t30\t40\tcellC\t1\n",
                encoding="utf-8",
            )
            annotations.write_text(
                "barcode\tdonor\tcell_type\ncellA\td1\tT\ncellB\td1\tT\ncellC\td2\tB\n",
                encoding="utf-8",
            )

            manifest = group_fragments(fragments, annotations, outdir, group_by=["donor", "cell_type"], min_cells=2, min_fragments=5)

            self.assertEqual(set(manifest["group"]), {"d1__T", "d2__B"})
            kept = manifest.loc[manifest["group"] == "d1__T"].iloc[0]
            filtered = manifest.loc[manifest["group"] == "d2__B"].iloc[0]
            self.assertTrue(bool(kept["passes_filters"]))
            self.assertFalse(bool(filtered["passes_filters"]))
            self.assertTrue((outdir / "d1__T.fragments.tsv").exists())
            self.assertTrue((outdir / "pseudobulk_manifest.tsv").exists())
            self.assertTrue((outdir / "fp_tools_manifest.yml").exists())

    def test_write_downstream_commands_for_kept_groups(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fragments = tmp / "fragments.tsv"
            annotations = tmp / "annotations.tsv"
            outdir = tmp / "out"
            genome_sizes = tmp / "genome.sizes"
            commands = outdir / "commands.sh"
            fragments.write_text(
                "chr1\t1\t5\tcellA\t2\nchr1\t10\t20\tcellB\t3\nchr1\t30\t40\tcellC\t1\n",
                encoding="utf-8",
            )
            annotations.write_text(
                "barcode\tdonor\tcell_type\ncellA\td1\tT\ncellB\td1\tT\ncellC\td2\tB\n",
                encoding="utf-8",
            )
            genome_sizes.write_text("chr1\t100\n", encoding="utf-8")

            manifest = group_fragments(fragments, annotations, outdir, group_by=["donor", "cell_type"], min_cells=2, min_fragments=5)
            command_path = write_downstream_commands(manifest, commands, genome_sizes=genome_sizes, cores=8)
            text = command_path.read_text(encoding="utf-8")

        self.assertIn("bedtools genomecov", text)
        self.assertIn("bedGraphToBigWig", text)
        self.assertIn("bedToBam", text)
        self.assertIn("samtools sort -@ 8", text)
        self.assertIn("d1__T", text)
        self.assertNotIn("d2__B", text)


if __name__ == "__main__":
    unittest.main()
