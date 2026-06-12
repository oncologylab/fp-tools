import tempfile
import unittest
from pathlib import Path

from fp_tools.tools.plot_aggregate_batch import _discover_motifs, write_html


class PlotAggregateBatchTest(unittest.TestCase):
    def test_discovers_motifs_from_match_dir_and_writes_html(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            bed_dir = root / "CTCF_MA0139.1" / "beds"
            bed_dir.mkdir(parents=True)
            (bed_dir / "CTCF_MA0139.1_all.bed").write_text("chr1\t10\t20\tsite1\n", encoding="utf-8")
            motifs = _discover_motifs(root)
            out = root / "report.html"
            write_html({"x": [-1, 1], "samples": [{"sample": "S1", "label": "S1", "condition": "S1", "motifs": [{"prefix": "CTCF_MA0139.1", "name": "CTCF", "score": 1, "sites": 1, "profile": [0.1, 0.2]}]}]}, out, "Demo")
            html = out.read_text(encoding="utf-8")

        self.assertEqual(motifs[0]["prefix"], "CTCF_MA0139.1")
        self.assertIn("Demo", html)
        self.assertIn("Top motifs", html)


if __name__ == "__main__":
    unittest.main()
