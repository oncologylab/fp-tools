import hashlib
import io
import json
import os
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from fp_tools import runtime


class RuntimeManagerTest(unittest.TestCase):
    def test_platform_keys_are_normalized(self):
        with mock.patch("platform.system", return_value="Darwin"), mock.patch(
            "platform.machine", return_value="arm64"
        ):
            self.assertEqual(runtime.platform_key(), "macos-arm64")
        with mock.patch("platform.system", return_value="Linux"), mock.patch(
            "platform.machine", return_value="aarch64"
        ):
            self.assertEqual(runtime.platform_key(), "linux-arm64")
        with mock.patch("platform.system", return_value="Windows"), mock.patch(
            "platform.machine", return_value="AMD64"
        ):
            self.assertEqual(runtime.platform_key(), "windows-x86_64")

    def test_flag_path_translation_preserves_non_path_arguments(self):
        translated = runtime.translate_flag_paths(
            ["--samples", "reads.tsv", "--genome", "hg38", "--outdir=project", "--cores", "8"],
            {"--samples", "--outdir"},
            lambda value: "/mapped/" + value,
        )
        self.assertEqual(
            translated,
            ["--samples", "/mapped/reads.tsv", "--genome", "hg38", "--outdir=/mapped/project", "--cores", "8"],
        )

    def test_flag_path_translation_maps_every_value_of_a_list_option(self):
        translated = runtime.translate_flag_paths(
            ["--motifs", "a.meme", "b.meme", "--cores", "4"],
            {"--motifs"},
            lambda value: "/mapped/" + value,
        )
        self.assertEqual(
            translated,
            ["--motifs", "/mapped/a.meme", "/mapped/b.meme", "--cores", "4"],
        )

    def test_safe_extract_rejects_path_traversal(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            archive = Path(tmpdir) / "bad.tar.gz"
            with tarfile.open(archive, "w:gz") as handle:
                info = tarfile.TarInfo("../outside")
                info.size = 1
                handle.addfile(info, io.BytesIO(b"x"))
            with self.assertRaises(runtime.RuntimeProvisionError):
                runtime._safe_extract(archive, Path(tmpdir) / "extract")

    def test_native_runtime_download_is_verified_and_cached(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "runtime.tar.gz"
            payload = root / "payload"
            (payload / "bin").mkdir(parents=True)
            for command in ("fastp", "samtools"):
                path = payload / "bin" / command
                path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
                path.chmod(0o755)
            with tarfile.open(source, "w:gz") as handle:
                for path in sorted(payload.rglob("*")):
                    handle.add(path, arcname=path.relative_to(payload))
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            manifest = {
                "schema": 1,
                "runtime_version": "test-1",
                "repository": "example.invalid/runtime",
                "components": {"core": {"commands": ["fastp", "samtools"]}},
                "artifacts": {"linux-x86_64": {"core": {"tag": "test"}}},
            }
            with mock.patch.object(runtime, "load_runtime_manifest", return_value=manifest), mock.patch.object(
                runtime, "runtime_cache_root", return_value=root / "cache"
            ), mock.patch.object(runtime, "platform_key", return_value="linux-x86_64"), mock.patch.object(
                runtime,
                "_resolve_runtime_artifact",
                return_value=(source.as_uri(), source.stat().st_size, digest),
            ):
                first = runtime.ensure_native_runtime("core")
                second = runtime.ensure_native_runtime("core")
            self.assertEqual(first.prefix, second.prefix)
            self.assertTrue((first.prefix / "bin" / "fastp").is_file())
            marker = json.loads((first.prefix / ".fp-tools-runtime.json").read_text(encoding="utf-8"))
            self.assertEqual(marker["sha256"], digest)

    def test_release_artifact_uses_public_checksum_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            archive = root / "runtime.tar.gz"
            archive.write_bytes(b"runtime")
            digest = hashlib.sha256(archive.read_bytes()).hexdigest()
            (root / "runtime.tar.gz.sha256").write_text(
                f"{digest}  runtime.tar.gz\n", encoding="utf-8"
            )
            url, size, observed = runtime._resolve_runtime_artifact(
                {
                    "release_base_url": root.as_uri(),
                    "filename": "runtime.tar.gz",
                }
            )
        self.assertEqual(url, archive.as_uri())
        self.assertEqual(size, 0)
        self.assertEqual(observed, digest)

    def test_managed_runtime_activates_relocated_ca_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "runtime"
            (prefix / "bin").mkdir(parents=True)
            bundle = prefix / "ssl" / "cacert.pem"
            bundle.parent.mkdir(parents=True)
            bundle.write_text("test certificate bundle\n", encoding="utf-8")
            activation = runtime.RuntimeActivation("managed", "core", prefix=prefix)
            with mock.patch.object(
                runtime, "ensure_native_runtime", return_value=activation
            ), mock.patch.dict(os.environ, {"PATH": "/system/bin"}, clear=True):
                observed = runtime.activate_runtime("core", "managed")
                self.assertEqual(os.environ["CURL_CA_BUNDLE"], str(bundle.resolve()))
                self.assertEqual(os.environ["SSL_CERT_FILE"], str(bundle.resolve()))
                self.assertTrue(os.environ["PATH"].startswith(str(prefix / "bin")))
            self.assertEqual(observed, activation)

    def test_download_rejects_wrong_checksum(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "source.tar.gz"
            source.write_bytes(b"runtime")
            destination = Path(tmpdir) / "download.tar.gz"
            with self.assertRaises(runtime.RuntimeProvisionError):
                runtime._download(source.as_uri(), destination, source.stat().st_size, "0" * 64)
            self.assertFalse(destination.with_suffix(".gz.part").exists())

    def test_runtime_status_does_not_download(self):
        manifest = {
            "runtime_version": "test-1",
            "components": {"core": {"commands": []}},
            "artifacts": {"linux-x86_64": {"core": {"tag": "test"}}},
        }
        with tempfile.TemporaryDirectory() as tmpdir, mock.patch.object(
            runtime, "load_runtime_manifest", return_value=manifest
        ), mock.patch.object(runtime, "runtime_cache_root", return_value=Path(tmpdir)), mock.patch.object(
            runtime, "platform_key", return_value="linux-x86_64"
        ), mock.patch.object(runtime, "_resolve_oci_layer") as resolver:
            rows = runtime.runtime_status()
        resolver.assert_not_called()
        self.assertEqual(rows[0]["installed"], "no")

    def test_non_linux_runtime_manifest_exposes_only_meme(self):
        manifest = runtime.load_runtime_manifest()
        for target in ("macos-x86_64", "macos-arm64", "windows-x86_64"):
            self.assertEqual(set(manifest["artifacts"][target]), {"meme"})
        for target in ("linux-x86_64", "linux-arm64"):
            self.assertEqual(set(manifest["artifacts"][target]), {"core", "meme", "homer"})


if __name__ == "__main__":
    unittest.main()
