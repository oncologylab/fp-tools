#!/usr/bin/env python3
"""Test the public Windows EXE with a fresh WSL import and cached discovery."""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile

from smoke_frozen_meme_release import _write_inputs


def distributions():
    result = subprocess.run(["wsl.exe", "--list", "--quiet"], capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(f"Could not list WSL distributions: {result.stderr!r}")
    return result.stdout.decode("utf-16-le").replace("\x00", "").split()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--case", choices=["default", "extra"], required=True)
    args = parser.parse_args()
    executable = args.executable.resolve()
    distro = "fp-tools-" + args.version.removeprefix("v").replace(".", "-")
    if distro in distributions():
        raise SystemExit("Fresh-import test requires an isolated runner without " + distro)
    with tempfile.TemporaryDirectory(prefix="fp tools Windows release ") as directory:
        root = Path(directory)
        fasta, _ = _write_inputs(root)
        environment = os.environ.copy()
        environment["FP_TOOLS_RUNTIME_CACHE"] = str(root / "runtime cache")
        try:
            for run in ("first", "cached"):
                outdir = root / (run + " discovery")
                runtime = ["--runtime", "managed"] if run == "first" else ["--runtime=managed"]
                command = [str(executable), "--fp-tools-internal-command", "discover-motifs",
                           "--fasta", str(fasta), "--outdir", str(outdir), "--method", "streme",
                           "--known-motif-db", "jaspar2026_vertebrates", "--execute", *runtime]
                if args.case == "extra":
                    command.extend(["--extra-args", "--dna", "--nmotifs", "1", "--minw", "4", "--maxw", "8", "--seed", "1"])
                result = subprocess.run(command, cwd=root, env=environment, capture_output=True,
                                        text=True, encoding="utf-8", errors="replace", timeout=1200)
                print(result.stdout, result.stderr, flush=True)
                if result.returncode:
                    raise RuntimeError(f"{run} discovery failed with exit {result.returncode}")
                required = [outdir / "streme/streme.txt", outdir / "tomtom/tomtom.tsv",
                            outdir / "motif_summary.tsv", outdir / "motif_summary.html"]
                matrices = list((outdir / "known_motifs").glob("*.meme"))
                if not matrices or any(not p.is_file() or not p.stat().st_size for p in required + matrices):
                    raise RuntimeError(f"{run} discovery omitted required outputs")
                if distro not in distributions():
                    raise RuntimeError("Managed distribution was not imported")
                disk = list((root / "runtime cache").rglob("*.vhdx"))
                if len(disk) != 1:
                    raise RuntimeError(f"Expected one private WSL virtual disk: {disk}")
                identity = (str(disk[0]), disk[0].stat().st_ctime_ns)
                if run == "first":
                    first_identity = identity
                elif identity != first_identity:
                    raise RuntimeError("Cached discovery replaced the WSL installation")
                print(f"{args.case}: {run} discovery passed all five output checks", flush=True)
        finally:
            # This invocation verified absence before creating this test-owned distro.
            if distro in distributions():
                subprocess.run(["wsl.exe", "--unregister", distro], check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
