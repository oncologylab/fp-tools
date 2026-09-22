"""Exercise download failure/recovery without accessing ENCODE."""

import os
from pathlib import Path
import shutil
import subprocess

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "docs/demos/data/encode/download_encode_hepg2_k562.sh"
pytestmark = pytest.mark.skipif(os.name == "nt" or not shutil.which("bash"), reason="Linux/macOS Bash helper")


def run(tmp_path, body):
    return subprocess.run(["bash", "-c", 'source "$1"\n' + body, "test", str(SCRIPT)], cwd=tmp_path, capture_output=True, text=True)


def test_missing_prerequisite_fails_before_download_or_directory_creation(tmp_path):
    result = run(tmp_path, '''
command() { if [[ "$*" == "-v samtools" ]]; then return 1; else builtin command "$@"; fi; }
curl() { touch unwanted_download; }
main
''')
    assert result.returncode != 0
    assert "Required command missing: samtools" in result.stderr
    assert list(tmp_path.iterdir()) == []


def test_valid_cached_download_is_not_fetched_again(tmp_path):
    (tmp_path / "cached.bam").write_text("valid")
    result = run(tmp_path, '''
checksum() { echo expected; }
curl() { touch unwanted_download; return 1; }
download expected cached.bam
''')
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "unwanted_download").exists()


@pytest.mark.parametrize("code", [33, 36])
def test_unsupported_resume_restarts_and_validates(tmp_path, code):
    result = run(tmp_path, f'''
checksum() {{ echo expected; }}
curl() {{
  if [[ "$*" == *"-C -"* ]]; then return {code}; fi
  touch sample.bam.part
}}
download expected sample.bam
''')
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "sample.bam").exists()
    assert not (tmp_path / "sample.bam.part").exists()


def test_corrupt_download_is_not_promoted(tmp_path):
    result = run(tmp_path, '''
checksum() { echo wrong; }
curl() { touch sample.bam.part; }
download expected sample.bam
''')
    assert result.returncode == 1
    assert "Checksum mismatch" in result.stderr
    assert not (tmp_path / "sample.bam").exists()


def test_macos_checksum_fallback(tmp_path):
    result = run(tmp_path, '''
command() { if [[ "$*" == "-v md5sum" ]]; then return 1; else builtin command "$@"; fi; }
md5() { [[ "$1" == -q ]] && echo expected; }
checksum sample.bam
''')
    assert result.returncode == 0
    assert result.stdout.strip() == "expected"
