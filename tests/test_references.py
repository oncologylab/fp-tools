import gzip
import hashlib
import io
import tempfile
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest import mock

from fp_tools.utils import references


def _compressed(payload: bytes) -> bytes:
    return gzip.compress(payload, mtime=0)


class ManagedReferenceTest(unittest.TestCase):
    def setUp(self):
        self.fasta = _compressed(b">chr1\nACGTACGT\n")
        self.blacklist = _compressed(b"chr1\t1\t2\n")
        self.manifest = {
            "test": {
                "fasta_url": "https://example.invalid/test.fa.gz",
                "fasta_md5": hashlib.md5(self.fasta).hexdigest(),
                "blacklist_url": "https://example.invalid/test.blacklist.bed.gz",
                "blacklist_md5": hashlib.md5(self.blacklist).hexdigest(),
                "tss_gtf_url": "https://example.invalid/test.gtf.gz",
                "tss_gtf_md5": "unused",
                "macs_genome_size": "1000",
            }
        }

    def _response(self, url, **_kwargs):
        if str(url).endswith("blacklist.bed.gz"):
            return io.BytesIO(self.blacklist)
        return io.BytesIO(self.fasta)

    def test_downloads_checksums_indexes_and_reuses_cache(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(
            references.REFERENCE_MANIFEST, self.manifest, clear=True
        ), mock.patch.object(
            references.urllib.request, "urlopen", side_effect=self._response
        ) as urlopen:
            result = references.resolve_analysis_reference(
                "test", reference_dir=tmp
            )
            self.assertEqual(result.fasta.read_bytes(), b">chr1\nACGTACGT\n")
            self.assertEqual(result.blacklist.read_text(), "chr1\t1\t2\n")
            self.assertTrue(Path(str(result.fasta) + ".fai").is_file())
            self.assertEqual(urlopen.call_count, 2)

            cached = references.resolve_analysis_reference(
                "test", reference_dir=tmp
            )
            self.assertEqual(cached, result)
            self.assertEqual(urlopen.call_count, 2)

    def test_corrupt_cached_asset_is_replaced(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(
            references.REFERENCE_MANIFEST, self.manifest, clear=True
        ), mock.patch.object(
            references.urllib.request, "urlopen", side_effect=self._response
        ) as urlopen:
            result = references.resolve_analysis_reference(
                "test", reference_dir=tmp
            )
            result.fasta.write_text("corrupt", encoding="utf-8")
            repaired = references.resolve_analysis_reference(
                "test", reference_dir=tmp
            )
            self.assertEqual(repaired.fasta.read_bytes(), b">chr1\nACGTACGT\n")
            self.assertEqual(urlopen.call_count, 3)

    def test_concurrent_resolution_downloads_each_asset_once(self):
        def delayed_response(url, **kwargs):
            time.sleep(0.02)
            return self._response(url, **kwargs)

        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(
            references.REFERENCE_MANIFEST, self.manifest, clear=True
        ), mock.patch.object(
            references.urllib.request, "urlopen", side_effect=delayed_response
        ) as urlopen:
            with ThreadPoolExecutor(max_workers=4) as pool:
                results = list(
                    pool.map(
                        lambda _: references.resolve_analysis_reference(
                            "test", reference_dir=tmp
                        ),
                        range(4),
                    )
                )
            self.assertEqual(len({result.fasta for result in results}), 1)
            self.assertEqual(urlopen.call_count, 2)

    def test_dry_run_resolves_paths_without_creating_cache(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(
            references.REFERENCE_MANIFEST, self.manifest, clear=True
        ), mock.patch.object(
            references.urllib.request, "urlopen"
        ) as urlopen:
            cache = Path(tmp) / "cache"
            result = references.resolve_analysis_reference(
                "test", reference_dir=cache, dry_run=True
            )
            self.assertEqual(result.fasta, cache / "test" / "test.fa")
            self.assertEqual(
                result.blacklist, cache / "test" / "test.blacklist.bed"
            )
            self.assertFalse(cache.exists())
            urlopen.assert_not_called()

    def test_blacklist_override_skips_managed_blacklist(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(
            references.REFERENCE_MANIFEST, self.manifest, clear=True
        ), mock.patch.object(
            references.urllib.request, "urlopen", side_effect=self._response
        ) as urlopen:
            custom = Path(tmp) / "custom.bed"
            custom.write_text("chr1\t2\t3\n", encoding="utf-8")
            result = references.resolve_analysis_reference(
                "test", reference_dir=Path(tmp) / "cache", blacklist=custom
            )
            self.assertEqual(result.blacklist, custom.resolve())
            self.assertEqual(urlopen.call_count, 1)

    def test_custom_fasta_never_infers_a_blacklist(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "custom.fa"
            fasta.write_text(">contig\nACGT\n", encoding="utf-8")
            result = references.resolve_analysis_reference(fasta)
            self.assertFalse(result.managed)
            self.assertIsNone(result.assembly)
            self.assertIsNone(result.blacklist)

    def test_blacklist_and_disable_flag_are_mutually_exclusive(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "custom.fa"
            blacklist = Path(tmp) / "custom.bed"
            fasta.write_text(">contig\nACGT\n", encoding="utf-8")
            blacklist.write_text("contig\t0\t1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "mutually exclusive"):
                references.resolve_analysis_reference(
                    fasta, blacklist=blacklist, no_blacklist=True
                )


if __name__ == "__main__":
    unittest.main()
