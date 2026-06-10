import unittest
from pathlib import Path
import tempfile

import numpy as np
import pyBigWig

from fp_tools.tools.candidates import call_candidates_from_array, generate_candidates


class CandidateGenerationTest(unittest.TestCase):
    def test_call_candidates_from_array_finds_local_maxima(self):
        values = np.array([0.0, 1.0, 0.1, 0.2, 3.0, 0.2, 2.5, 0.1])
        candidates = call_candidates_from_array(values, "chr1", 100, 108, candidate_width=4, window=1, min_score=1.0, top_n=2)
        self.assertEqual([candidate.rank for candidate in candidates], [1, 2])
        self.assertEqual([round(candidate.score, 1) for candidate in candidates], [3.0, 2.5])
        self.assertEqual((candidates[0].chrom, candidates[0].start, candidates[0].end), ("chr1", 102, 106))

    def test_call_candidates_from_array_collapses_plateaus_and_thresholds(self):
        values = np.array([0.0, 2.0, 2.0, 0.0, 0.5])
        candidates = call_candidates_from_array(values, "chr2", 10, 15, candidate_width=2, window=1, min_score=1.0, top_n=None)
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0].start, 10)
        self.assertEqual(candidates[0].source_region, "chr2:10-15")

    def test_generate_candidates_merges_relaxed_motif_sites(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            signal = tmp / "signal.bw"
            peaks = tmp / "peaks.bed"
            motifs = tmp / "relaxed_motifs.bed"
            out = tmp / "candidates.bed"

            bw = pyBigWig.open(str(signal), "w")
            bw.addHeader([("chr1", 30)])
            bw.addEntries(["chr1"] * 30, list(range(30)), ends=list(range(1, 31)), values=[float(i) for i in range(30)])
            bw.close()
            peaks.write_text("chr1\t0\t20\n", encoding="utf-8")
            motifs.write_text("chr1\t5\t8\tCTCF_weak\t4.2\t+\nchr1\t22\t25\toutside\t9\t+\n", encoding="utf-8")

            candidates = generate_candidates(
                signal,
                peaks,
                out,
                min_score=100.0,
                top_n_per_region=1,
                motif_sites=[motifs],
            )
            rows = out.read_text(encoding="utf-8").splitlines()

        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0].generator, "motif-relaxed")
        self.assertEqual(candidates[0].motif_id, "CTCF_weak")
        self.assertAlmostEqual(candidates[0].candidate_score, 7.0)
        self.assertIn("candidate_score\tmotif_id\tmotif_score\tmotif_source", rows[0])
        self.assertIn("CTCF_weak\t4.200000", rows[1])


if __name__ == "__main__":
    unittest.main()
