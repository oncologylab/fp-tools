import tempfile
import unittest
from pathlib import Path

import pyBigWig

from fp_tools.tools.tfbs_features import build_feature_table


class TfbsFeatureBuilderTest(unittest.TestCase):
    def test_build_feature_table_with_genome_signal_and_labels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            genome = tmp / "genome.fa"
            candidates = tmp / "candidates.bed"
            labels = tmp / "labels.bed"
            signal = tmp / "signal.bw"
            out = tmp / "features.tsv"
            genome.write_text(">chr1\nAACCGGTTAACC\n", encoding="utf-8")
            candidates.write_text("chr1\t2\t6\tcand1\t3.5\nchr1\t8\t10\tcand2\t1.0\n", encoding="utf-8")
            labels.write_text("chr1\t3\t5\tpositive\n", encoding="utf-8")
            bw = pyBigWig.open(str(signal), "w")
            bw.addHeader([("chr1", 12)])
            bw.addEntries("chr1", list(range(12)), values=[float(i) for i in range(12)], span=1)
            bw.close()

            frame = build_feature_table(
                candidates,
                out,
                genome=genome,
                signals=[signal],
                signal_labels=["cut"],
                labels_bed=labels,
            )
            self.assertTrue(out.exists())

        self.assertEqual(frame["site_id"].tolist(), ["cand1", "cand2"])
        self.assertEqual(frame["label"].tolist(), [1, 0])
        self.assertAlmostEqual(float(frame.loc[0, "candidate_score"]), 3.5)
        self.assertAlmostEqual(float(frame.loc[0, "gc"]), 1.0)
        self.assertAlmostEqual(float(frame.loc[0, "cut_mean"]), 3.5)
        self.assertAlmostEqual(float(frame.loc[1, "cut_max"]), 9.0)


if __name__ == "__main__":
    unittest.main()
