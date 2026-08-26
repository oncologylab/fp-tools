from __future__ import annotations

import tempfile
import threading
import time
import errno
import unittest
import shutil
import subprocess
import os
import sys
from pathlib import Path
from unittest import mock

import numpy as np

from fp_tools.utils import bigwig
from fp_tools.utils import fasta as fasta_module
from fp_tools.utils.alignment import FragmentAlignment, _BamnosticRecord
from fp_tools.utils.fasta import open_fasta
from fp_tools.utils.intervals import IntervalIndex, intersect_bed
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.motifs import MotifList, _NumpyMotifScanner, _threshold_from_p


ROOT = Path(__file__).resolve().parents[1]


class CrossPlatformIoTests(unittest.TestCase):
    def test_pybigtools_windows_write_paths_are_not_url_like(self):
        with mock.patch.object(bigwig.os, "name", "nt"):
            self.assertEqual(
                bigwig._pybigtools_write_path(r"C:\Users\runner\result.bw"),
                r"\\?\C:\Users\runner\result.bw",
            )
            self.assertEqual(
                bigwig._pybigtools_write_path(r"\\server\share\result.bw"),
                r"\\?\UNC\server\share\result.bw",
            )
            self.assertEqual(
                bigwig._pybigtools_write_path(r"\\?\C:\data\result.bw"),
                r"\\?\C:\data\result.bw",
            )

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
                writer.addEntries("chr1", 2, values=[1.5, -2.0], span=1, step=6)
                writer.close()
                with bigwig.open(output) as reader:
                    self.assertEqual(reader.chroms(), {"chr1": 20})
                    self.assertEqual(reader.intervals("chr1"), [(2, 3, 1.5), (8, 9, -2.0)])
                    self.assertEqual(reader.header()["nBasesCovered"], 2)
        finally:
            bigwig._pybigwig = native

    def test_pyfastx_index_creation_is_serialized(self):
        active = 0
        maximum_active = 0
        state_lock = threading.Lock()

        class FakePyfastx:
            @staticmethod
            def Fasta(_filename, **_kwargs):
                nonlocal active, maximum_active
                with state_lock:
                    active += 1
                    maximum_active = max(maximum_active, active)
                time.sleep(0.1)
                with state_lock:
                    active -= 1
                return object()

        with tempfile.TemporaryDirectory() as tmp:
            index = Path(tmp) / "shared.fxi"
            with mock.patch.object(fasta_module, "pyfastx", FakePyfastx):
                threads = [
                    threading.Thread(target=fasta_module._open_pyfastx, args=("genome.fa.gz", index))
                    for _ in range(4)
                ]
                for thread in threads:
                    thread.start()
                for thread in threads:
                    thread.join()
            self.assertEqual(maximum_active, 1)
            self.assertFalse(index.with_suffix(".fxi.lock").exists())

    def test_pyfastx_lock_treats_windows_access_denied_as_contention(self):
        real_open = fasta_module.os.open
        attempts = 0

        class FakePyfastx:
            @staticmethod
            def Fasta(_filename, **_kwargs):
                return object()

        def windows_open(path, flags):
            nonlocal attempts
            attempts += 1
            if attempts == 1:
                raise PermissionError(errno.EACCES, "lock file exists", str(path))
            return real_open(path, flags)

        with tempfile.TemporaryDirectory() as tmp:
            index = Path(tmp) / "shared.fxi"
            with (
                mock.patch.object(fasta_module, "pyfastx", FakePyfastx),
                mock.patch.object(fasta_module.os, "open", side_effect=windows_open),
                mock.patch.object(fasta_module.time, "sleep"),
            ):
                fasta_module._open_pyfastx("genome.fa.gz", index)
            self.assertEqual(attempts, 2)

    def test_fasta_adapter_reports_reference_length(self):
        with open_fasta(ROOT / "test_data" / "genome.fa.gz") as fasta:
            reference = fasta.references[0]
            self.assertEqual(fasta.get_reference_length(reference), fasta.lengths[0])

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
