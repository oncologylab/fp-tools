#!/usr/bin/env python3
"""Run the checksum-verified PBMC tutorial through a released desktop app."""

import argparse
import csv
import hashlib
import io
from pathlib import Path
import subprocess
import tempfile
import urllib.request
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    executable = str(args.executable.resolve())
    base = f"https://github.com/oncologylab/fp-tools/releases/download/v{args.version.removeprefix('v')}/"
    name = "fp-tools-pbmc-chr22-demo-v1.zip"
    with urllib.request.urlopen(base + name + ".sha256", timeout=60) as response:
        expected = response.read().decode().split()[0]
    with urllib.request.urlopen(base + name, timeout=120) as response:
        data = response.read()
    if hashlib.sha256(data).hexdigest() != expected:
        raise RuntimeError("PBMC demo archive checksum mismatch")
    with tempfile.TemporaryDirectory(prefix="fp tools PBMC ") as temporary:
        root = Path(temporary)
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            for member in archive.infolist():
                if not (root / member.filename).resolve().is_relative_to(root.resolve()):
                    raise RuntimeError("Invalid archive path")
            archive.extractall(root)
        demo = root / "fp-tools-pbmc-chr22-demo-v1"
        for line in (demo / "SHA256SUMS.txt").read_text().splitlines():
            digest, filename = line.split(maxsplit=1)
            if hashlib.sha256((demo / filename).read_bytes()).hexdigest() != digest:
                raise RuntimeError(f"Input checksum mismatch: {filename}")
        for command in ("plot-aggregate", "find-signature-fp", "sc-footprinting"):
            result = subprocess.run([executable, "--fp-tools-internal-command", command, "--help"],
                                    check=True, capture_output=True, text=True, encoding="utf-8")
            help_text = " ".join(result.stdout.split())
            if command == "plot-aggregate":
                assert "--output_aggregated_signals" in help_text and "--share_y" not in help_text
            else:
                assert "genomic-bin counts" in help_text and "selected" in help_text
        subprocess.run([executable, "--fp-tools-internal-command", "run-yaml-workflow",
                        "--config", "workflow.yml"], cwd=demo, check=True, timeout=3600)
        with (demo / "results/pseudobulk_footprint_manifest.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 3 and all(row["status"] == "succeeded" for row in rows), rows
        plots = demo / "results/plots/single_cell_footprinting"
        for filename in ("per_cell_footprint_signature_heatmap.svg", "knn_footprint_signature_scores.tsv",
                         "chromvar_like_motif_activity_scores.tsv", "single_cell_footprinting_summary.svg"):
            assert (plots / filename).stat().st_size > 0, filename
        print("Verified real PBMC tutorial: checksums, corrected help, three groups and per-cell outputs")


if __name__ == "__main__":
    main()
