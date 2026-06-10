import tempfile
import unittest
from pathlib import Path

import pandas as pd

from fp_tools.tools.tfbs_model import train_tabular_model
from fp_tools.tools.variants import read_pwm_motifs, score_variants


class VariantScoringTest(unittest.TestCase):
    def test_score_variants_checks_ref_and_candidate_overlap(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            genome = tmp / "genome.fa"
            variants = tmp / "variants.bed"
            candidates = tmp / "candidates.bed"
            motifs = tmp / "motifs.jaspar"
            out = tmp / "scored.tsv"
            genome.write_text(">chr1\nAACCGGTTAACC\n", encoding="utf-8")
            variants.write_text("chr1\t2\t3\tvar1\tC\tT\nchr1\t8\t9\tvar2\tG\tA\n", encoding="utf-8")
            candidates.write_text("chr1\t1\t5\tcand1\t4.5\n", encoding="utf-8")
            motifs.write_text(
                ">M1\tRefPattern\n"
                "A [ 10 10  0  0  0 ]\n"
                "C [  0  0 10 10  0 ]\n"
                "G [  0  0  0  0 10 ]\n"
                "T [  0  0  0  0  0 ]\n",
                encoding="utf-8",
            )

            frame = score_variants(
                variants,
                genome,
                out,
                candidate_scores=candidates,
                sequence_flank=2,
                kmer_size=2,
                motifs=[motifs],
                motif_flank=2,
            )
            self.assertTrue(out.exists())

        first = frame.loc[frame["name"] == "var1"].iloc[0]
        second = frame.loc[frame["name"] == "var2"].iloc[0]
        self.assertTrue(bool(first["ref_matches_genome"]))
        self.assertEqual(first["candidate_name"], "cand1")
        self.assertAlmostEqual(float(first["candidate_score"]), 4.5)
        self.assertEqual(first["ref_context"], "AACCG")
        self.assertEqual(first["alt_context"], "AATCG")
        self.assertAlmostEqual(float(first["ref_gc"]), 0.6)
        self.assertAlmostEqual(float(first["alt_gc"]), 0.4)
        self.assertAlmostEqual(float(first["delta_gc"]), -0.2)
        self.assertEqual(int(first["lost_kmers"]), 2)
        self.assertEqual(int(first["gained_kmers"]), 2)
        self.assertAlmostEqual(float(first["kmer_jaccard"]), 2 / 6)
        self.assertEqual(first["best_motif_id"], "M1")
        self.assertEqual(first["best_motif_name"], "RefPattern")
        self.assertLess(float(first["motif_delta_score"]), 0.0)
        self.assertEqual(first["motif_delta_direction"], "loss")
        self.assertFalse(bool(second["ref_matches_genome"]))
        self.assertFalse(bool(second["overlaps_candidate"]))

    def test_read_pwm_motifs_supports_meme_probability_matrix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            motif_path = Path(tmpdir) / "motif.meme"
            motif_path.write_text(
                "MEME version 4\n\n"
                "ALPHABET=ACGT\n\n"
                "MOTIF\tM2\tAltPattern\n"
                "letter-probability matrix: alength=4 w=3 nsites=10 E=0\n"
                "0.90 0.05 0.03 0.02\n"
                "0.01 0.92 0.04 0.03\n"
                "0.02 0.03 0.91 0.04\n",
                encoding="utf-8",
            )
            motifs = read_pwm_motifs(motif_path)

        self.assertEqual(len(motifs), 1)
        self.assertEqual(motifs[0].motif_id, "M2")
        self.assertEqual(motifs[0].name, "AltPattern")
        self.assertEqual(len(motifs[0].probabilities), 3)

    def test_score_variants_adds_model_probability_delta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            genome = tmp / "genome.fa"
            variants = tmp / "variants.bed"
            motifs = tmp / "motifs.jaspar"
            train = tmp / "train.tsv"
            model = tmp / "model.pkl"
            out = tmp / "scored.tsv"
            genome.write_text(">chr1\nAACCGGTTAACC\n", encoding="utf-8")
            variants.write_text("chr1\t2\t3\tvar1\tC\tT\n", encoding="utf-8")
            motifs.write_text(
                ">M1\tRefPattern\n"
                "A [ 10 10  0  0  0 ]\n"
                "C [  0  0 10 10  0 ]\n"
                "G [  0  0  0  0 10 ]\n"
                "T [  0  0  0  0  0 ]\n",
                encoding="utf-8",
            )
            pd.DataFrame(
                {
                    "label": [0, 0, 0, 1, 1, 1],
                    "motif_score": [1.0, 2.0, 3.0, 7.0, 8.0, 9.0],
                    "gc": [0.2, 0.3, 0.4, 0.6, 0.7, 0.8],
                }
            ).to_csv(train, sep="\t", index=False)
            train_tabular_model(train, model, feature_columns=["motif_score", "gc"], seed=11)

            frame = score_variants(
                variants,
                genome,
                out,
                sequence_flank=2,
                motifs=[motifs],
                motif_flank=2,
                tfbs_model=model,
            )

        row = frame.iloc[0]
        self.assertIn("ref_model_probability", frame.columns)
        self.assertIn("alt_model_probability", frame.columns)
        self.assertLess(float(row["model_delta_probability"]), 0.0)
        self.assertGreater(float(row["ref_model_probability"]), float(row["alt_model_probability"]))


if __name__ == "__main__":
    unittest.main()
