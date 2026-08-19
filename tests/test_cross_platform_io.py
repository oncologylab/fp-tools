from __future__ import annotations

import tempfile
import threading
import unittest
import shutil
import subprocess
import os
import sys
from pathlib import Path

import numpy as np

from fp_tools.utils import bigwig
from fp_tools.utils.alignment import FragmentAlignment, _BamnosticRecord
from fp_tools.utils.intervals import IntervalIndex, intersect_bed
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.motifs import MotifList, _NumpyMotifScanner, _threshold_from_p


ROOT = Path(__file__).resolve().parents[1]


class CrossPlatformIoTests(unittest.TestCase):
    def test_logger_listener_does_not_require_pickling_logger(self):
        logger = FpToolsLogger("portable-test", level=0)
        logger.start_logger_queue()
        self.assertIsInstance(logger.listener, threading.Thread)
        logger.stop_logger_queue()
        self.assertFalse(logger.listener.is_alive())

    def test_bamnostic_record_normalizes_soft_clips_and_template_length(self):
        class Raw:
            cigartuples = [(4, 5), (0, 45)]
            query_alignment_start = 0
            query_alignment_end = 45
            tlen = 78

        record = _BamnosticRecord(Raw())
        self.assertEqual(record.query_alignment_start, 5)
        self.assertEqual(record.query_alignment_end, 50)
        self.assertEqual(record.query_alignment_length, 45)
        self.assertEqual(record.template_length, 78)

    def test_pybigtools_reader_matches_native_backend(self):
        source = ROOT / "test_data" / "Tcell_uncorrected.bw"
        if bigwig._pybigwig is None:
            self.skipTest("native pyBigWig comparison backend unavailable")
        with bigwig.open(source) as handle:
            chrom = next(iter(handle.chroms()))
            expected = np.asarray(handle.values(chrom, 10_000, 11_000, numpy=True))
        native = bigwig._pybigwig
        try:
            bigwig._pybigwig = None
            with bigwig.open(source) as handle:
                observed = np.asarray(handle.values(chrom, 10_000, 11_000, numpy=True))
        finally:
            bigwig._pybigwig = native
        self.assertTrue(np.allclose(expected, observed, equal_nan=True))

    def test_pybigtools_writer_supports_fp_tools_surface(self):
        native = bigwig._pybigwig
        try:
            bigwig._pybigwig = None
            with tempfile.TemporaryDirectory() as tmp:
                output = Path(tmp) / "portable.bw"
                writer = bigwig.open(output, "w")
                writer.addHeader([("chr1", 20)])
                writer.addEntries("chr1", [2, 8], values=[1.5, -2.0], span=1)
                writer.close()
                with bigwig.open(output) as reader:
                    self.assertEqual(reader.chroms(), {"chr1": 20})
                    self.assertEqual(reader.intervals("chr1"), [(2, 3, 1.5), (8, 9, -2.0)])
        finally:
            bigwig._pybigwig = native

    def test_interval_index_and_bed_intersection(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            query = tmp / "query.bed"
            regions = tmp / "regions.bed"
            output = tmp / "out.bed"
            query.write_text("chr1\t0\t5\ta\nchr1\t10\t20\tb\nchr2\t0\t2\tc\n", encoding="utf-8")
            regions.write_text("chr1\t4\t11\n", encoding="utf-8")
            index = IntervalIndex.from_bed(regions)
            self.assertTrue(index.overlaps("chr1", 0, 5))
            self.assertFalse(index.overlaps("chr1", 11, 12))
            intersect_bed(query, regions, output)
            self.assertEqual(output.read_text(encoding="utf-8"), "chr1\t0\t5\ta\nchr1\t10\t20\tb\n")

    def test_numpy_scanner_uses_frozen_moods_threshold(self):
        motifs = MotifList().from_file(
            str(ROOT / "src" / "fp_tools" / "resources" / "motifs" / "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
        )
        motif = motifs[0]
        motif.get_pssm()
        threshold = _threshold_from_p(motif.pssm, motif.bg, 1e-4)
        self.assertAlmostEqual(threshold, 7.310867122206648, places=12)
        scanner = _NumpyMotifScanner()
        scanner.set_motifs([motif.pssm], motif.bg, [threshold])
        hits = scanner.scan("ACGT" * 50)[0]
        self.assertTrue(all(hit.score >= threshold for hit in hits))

    def test_fragment_alignment_emits_both_cutsite_anchors(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.fragments.tsv"
            path.write_text("chr1\t10\t20\tcell\t2\n", encoding="utf-8")
            alignment = FragmentAlignment(str(path))
            records = list(alignment.fetch("chr1", 0, 30))
            self.assertEqual(len(records), 4)
            self.assertEqual([record.reference_start for record in records], [10, 18, 10, 18])
            self.assertEqual([record.is_reverse for record in records], [False, True, False, True])

    def test_fragment_native_atac_correct_smoke(self):
        script_name = "atac-correct.exe" if os.name == "nt" else "atac-correct"
        executable_path = Path(sys.executable).parent / script_name
        executable = str(executable_path) if executable_path.exists() else shutil.which("atac-correct")
        if executable is None:
            self.skipTest("installed atac-correct console script is unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome = tmp / "genome.fa"
            genome.write_text(">chr1\n" + "ACGT" * 500 + "\n", encoding="utf-8")
            peaks = tmp / "peaks.bed"
            peak_rows = [(200, 400), (800, 1000), (1400, 1600)]
            peaks.write_text(
                "".join(f"chr1\t{start}\t{end}\n" for start, end in peak_rows),
                encoding="utf-8",
            )
            fragments = tmp / "sample.fragments.tsv"
            fragments.write_text(
                "".join(
                    f"chr1\t{start + 20 + index * 3}\t{start + 100 + index * 3}\tcell\t1\n"
                    for start, _end in peak_rows
                    for index in range(20)
                ),
                encoding="utf-8",
            )
            outdir = tmp / "out"
            result = subprocess.run(
                [
                    executable,
                    "--fragments", str(fragments),
                    "--genome", str(genome),
                    "--peaks", str(peaks),
                    "--regions-in", str(peaks),
                    "--regions-out", str(peaks),
                    "--outdir", str(outdir),
                    "--prefix", "fragment_test",
                    "--write-tracks", "corrected",
                    "--skip-qc",
                    "--scale-corrected", "none",
                    "--cores", "1",
                    "--verbosity", "1",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            output = outdir / "fragment_test_corrected.bw"
            self.assertTrue(output.is_file())
            with bigwig.open(output) as handle:
                self.assertEqual(handle.chroms(), {"chr1": 2000})


if __name__ == "__main__":
    unittest.main()
