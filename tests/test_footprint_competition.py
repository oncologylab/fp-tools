import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from fp_tools.tools.footprint_competition import (
    build_competition_report,
    decompose_region,
    select_band,
)
from fp_tools.utils.multiscale import write_multiscale_npz


SCALES = [8, 16, 24, 32, 64, 100, 147]


def _write_npz(tmp: Path) -> Path:
    """Two regions: one TF-scale dominated, one nucleosome-scale dominated."""

    width = 21

    def feature(active_scales: set[int], peak: float) -> dict:
        feats = {}
        for scale in SCALES:
            arr = np.zeros(width, dtype=np.float32)
            if scale in active_scales:
                arr[width // 2] = peak
            feats[scale] = arr
        return feats

    records = [
        (("chr1", 1000, 1021), feature({8, 16}, 5.0)),
        (("chr1", 2000, 2021), feature({100, 147}, 5.0)),
    ]
    path = tmp / "multiscale.npz"
    write_multiscale_npz(str(path), records, SCALES, summary_method="max")
    return path


class SelectBandTest(unittest.TestCase):
    def test_band_masking(self):
        scales = np.array(SCALES, dtype=float)
        self.assertEqual(scales[select_band(scales, 3, 30)].tolist(), [8.0, 16.0, 24.0])
        self.assertEqual(scales[select_band(scales, 120, 200)].tolist(), [147.0])


class DecomposeRegionTest(unittest.TestCase):
    def test_tf_only_when_no_nucleosome_signal(self):
        scales = np.array(SCALES, dtype=float)
        tensor = np.zeros((len(SCALES), 11), dtype=float)
        tensor[0, 5] = 4.0  # scale 8 -> TF band
        stats = decompose_region(tensor, scales, (3, 30), (120, 200))
        self.assertGreater(stats["tf_only_auc"], 0)
        self.assertEqual(stats["shared_auc"], 0.0)
        self.assertEqual(stats["dominant_component"], "tf")
        self.assertEqual(stats["competition_index"], 0.0)

    def test_shared_when_both_bands_active(self):
        scales = np.array(SCALES, dtype=float)
        tensor = np.zeros((len(SCALES), 11), dtype=float)
        tensor[0, 5] = 3.0  # TF band
        tensor[-1, 5] = 3.0  # nucleosome band (147)
        stats = decompose_region(tensor, scales, (3, 30), (120, 200))
        self.assertGreater(stats["shared_auc"], 0)
        self.assertEqual(stats["dominant_component"], "competing")
        self.assertAlmostEqual(stats["competition_index"], 1.0, places=6)


class BuildReportTest(unittest.TestCase):
    def test_report_summary_and_figure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            npz = _write_npz(tmp)
            out = tmp / "competition.tsv"
            summary_out = tmp / "summary.tsv"
            figure_out = tmp / "competition.png"

            report, summary = build_competition_report(
                npz,
                out,
                summary_output=summary_out,
                figure_output=figure_out,
            )

            self.assertTrue(out.exists() and summary_out.exists() and figure_out.exists())
            self.assertEqual(len(report), 2)
            by_region = dict(zip(report["start"], report["dominant_component"]))
            self.assertEqual(by_region[1000], "tf")
            self.assertEqual(by_region[2000], "nucleosome")

            on_disk = pd.read_csv(out, sep="\t")
            self.assertEqual(len(on_disk), 2)
            self.assertEqual(int(summary.iloc[0]["n_regions"]), 2)
            self.assertEqual(int(summary.iloc[0]["tf_dominant"]), 1)
            self.assertEqual(int(summary.iloc[0]["nucleosome_dominant"]), 1)

    def test_overlapping_bands_raise(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            npz = _write_npz(tmp)
            with self.assertRaises(ValueError):
                build_competition_report(npz, tmp / "out.tsv", tf_band=(3, 130), nuc_band=(120, 200))

    def test_band_with_no_scales_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            npz = _write_npz(tmp)
            with self.assertRaises(ValueError):
                build_competition_report(npz, tmp / "out.tsv", tf_band=(3, 5), nuc_band=(120, 200))


if __name__ == "__main__":
    unittest.main()
