import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from fp_tools.tools.bindetect_functions import plot_interactive_bindetect


class InteractiveBindetectHtmlTest(unittest.TestCase):
    def test_aggregate_payload_is_serialized_as_json(self):
        motif = SimpleNamespace(name="TF1", group="Bcell_up", change=1.2, pvalue=0.001, base="")
        aggregate_data = {
            "x": [-1, 0, 1],
            "motifs": [
                {
                    "prefix": "TF1_MA0001.1",
                    "name": "TF1",
                    "conditions": [
                        {"name": "Bcell", "profile": [0.2, 0.1, 0.2]},
                        {"name": "Tcell", "profile": [0.1, 0.3, 0.1]},
                    ],
                }
            ],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report.html"
            plot_interactive_bindetect([motif], ["Bcell", "Tcell"], str(out), aggregate_data=aggregate_data)
            html = out.read_text()
        self.assertIn("const aggregateData = {", html)
        self.assertIn("\"motifs\"", html)
        self.assertIn("TF1_MA0001.1", html)
        self.assertNotIn("json.dumps", html)


if __name__ == "__main__":
    unittest.main()
