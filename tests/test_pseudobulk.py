import tempfile
import unittest
from pathlib import Path

import gzip
import subprocess
import sys

import pyBigWig
import pysam

from fp_tools.tools.pseudobulk import group_fragments, write_cutsite_bigwig, write_downstream_commands


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

    def test_10x_suffix_csv_chrom_filter_and_bigwig(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fragments = tmp / "fragments.tsv.gz"
            annotations = tmp / "annotations.csv"
            genome_sizes = tmp / "genome.sizes"
            outdir = tmp / "out"
            with gzip.open(fragments, "wt", encoding="utf-8") as handle:
                handle.write(
                    "chr1\t1\t5\tcellA-1\t2\n"
                    "chr1\t10\t20\tcellB-1\t3\n"
                    "chr2\t30\t40\tcellC-1\t1\n"
                    "chrM\t2\t8\tcellA-1\t1\n"
                    "chr1\t50\t60\tunmatched-1\t1\n"
                )
            annotations.write_text(
                "barcode,cell_type\ncellA,B cell\ncellB,T cell\ncellC,Mono cell\n",
                encoding="utf-8",
            )
            genome_sizes.write_text("chr1\t100\nchr2\t100\n", encoding="utf-8")

            manifest = group_fragments(
                fragments,
                annotations,
                outdir,
                group_by=["cell_type"],
                include_chroms={"chr1", "chr2"},
                exclude_chroms={"chrM"},
                compress_output=True,
                index_output=True,
                genome_sizes=genome_sizes,
                write_cutsite_bigwigs=True,
                min_cells=1,
                min_fragments=1,
            )

            self.assertEqual(set(manifest["group"]), {"B_cell", "T_cell", "Mono_cell"})
            first_fragment = Path(manifest.loc[manifest["group"] == "B_cell", "fragment_file"].iloc[0])
            self.assertTrue(first_fragment.name.endswith(".tsv.gz"))
            self.assertTrue(first_fragment.with_suffix(first_fragment.suffix + ".tbi").exists())
            with pysam.TabixFile(str(first_fragment)) as tabix:
                rows = list(tabix.fetch("chr1", 0, 10))
            self.assertEqual(len(rows), 1)
            bigwig = Path(manifest.loc[manifest["group"] == "B_cell", "cutsite_bigwig"].iloc[0])
            self.assertTrue(bigwig.exists())
            bw = pyBigWig.open(str(bigwig))
            try:
                self.assertEqual(bw.chroms()["chr1"], 100)
                self.assertGreater(sum(0 if value != value else value for value in bw.values("chr1", 0, 10)), 0)
            finally:
                bw.close()

    def test_write_cutsite_bigwig_raw_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fragments = tmp / "fragments.tsv"
            genome_sizes = tmp / "genome.sizes"
            output = tmp / "cuts.bw"
            fragments.write_text("chr1\t1\t5\tcellA\t2\n", encoding="utf-8")
            genome_sizes.write_text("chr1\t10\n", encoding="utf-8")

            write_cutsite_bigwig(fragments, output, genome_sizes, cpm=False)

            bw = pyBigWig.open(str(output))
            try:
                self.assertEqual(bw.values("chr1", 1, 2)[0], 2.0)
                self.assertEqual(bw.values("chr1", 4, 5)[0], 2.0)
            finally:
                bw.close()

    def test_pseudobulk_aggregate_plot_script_smoke(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fragments = tmp / "fragments.tsv"
            genome_sizes = tmp / "genome.sizes"
            bigwig = tmp / "Group.cutsites.cpm.bw"
            manifest = tmp / "pseudobulk_manifest.tsv"
            site_dir = tmp / "sites"
            out_prefix = tmp / "aggregate"
            site_dir.mkdir()
            fragments.write_text("chr1\t45\t55\tcellA\t2\nchr1\t50\t60\tcellB\t1\n", encoding="utf-8")
            genome_sizes.write_text("chr1\t120\n", encoding="utf-8")
            write_cutsite_bigwig(fragments, bigwig, genome_sizes, cpm=True)
            manifest.write_text(
                "group\tfragment_file\tn_cells\tn_fragments\tcutsite_bigwig\tpasses_filters\n"
                f"Group\t{fragments}\t2\t3\t{bigwig}\tTrue\n",
                encoding="utf-8",
            )
            (site_dir / "TF1.motif_peaks.bed").write_text("chr1\t49\t51\tTF1\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    "paper/scripts/plot_pseudobulk_tf_aggregates.py",
                    "--manifest",
                    str(manifest),
                    "--tf-site-dir",
                    str(site_dir),
                    "--out-prefix",
                    str(out_prefix),
                    "--groups",
                    "Group",
                    "--tfs",
                    "TF1",
                    "--flank",
                    "20",
                ],
                cwd=Path(__file__).resolve().parents[1],
                check=True,
            )

            self.assertTrue(out_prefix.with_suffix(".png").exists())
            self.assertTrue(out_prefix.with_suffix(".pdf").exists())
            self.assertTrue(out_prefix.with_suffix(".tsv").exists())

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
