import argparse
import tempfile
import unittest
from pathlib import Path

import pyBigWig

from fp_tools.tools.score_bigwig import _write_candidate_bed
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.regions import RegionList


class CallFootprintsCandidateBedTest(unittest.TestCase):
    def test_output_bed_uses_ranked_local_maxima_and_filters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bw_path = tmp / "scores.bw"
            bed_path = tmp / "calls.bed"
            bw = pyBigWig.open(str(bw_path), "w")
            bw.addHeader([("chr1", 100)])
            bw.addEntries("chr1", 0, values=[0.0, 1.0, 5.0, 1.0, 0.0, 2.0, 9.0, 2.0, 0.0, 1.0, 4.0, 1.0, 0.0, 0.0, 8.0, 1.0, 0.0, 3.0, 1.0, 0.0], span=1, step=1)
            bw.close()
            regions_path = tmp / "regions.bed"
            regions_path.write_text("chr1\t0\t20\n")
            regions = RegionList().from_bed(str(regions_path))
            args = argparse.Namespace(score="footprint", min_score=4.5, call_width=6, min_distance=5, top_n=2)
            _write_candidate_bed(str(bw_path), regions, str(bed_path), args, {"chr1": 100}, FpToolsLogger("", 0))
            lines = [line for line in bed_path.read_text().splitlines() if not line.startswith("#")]

        self.assertEqual(len(lines), 2)
        self.assertIn("footprint_1", lines[0])
        self.assertTrue(lines[0].startswith("chr1\t3\t9"), lines[0])
        self.assertTrue(lines[1].startswith("chr1\t11\t17"), lines[1])


if __name__ == "__main__":
    unittest.main()
