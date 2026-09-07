#!/usr/bin/env python3
"""Exercise first-use managed MEME discovery from a released macOS app."""

from __future__ import annotations

import argparse
import os
import random
import subprocess
import tempfile
from pathlib import Path


INTERNAL_COMMAND = "--fp-tools-internal-command"


def _write_inputs(root: Path) -> tuple[Path, Path]:
    randomizer = random.Random(20260907)
    fasta = root / "candidate_sequences.fa"
    alphabet = "ACGT"
    records = []
    for index in range(80):
        sequence = "".join(randomizer.choice(alphabet) for _ in range(36))
        sequence = sequence[:15] + "ACGTAC" + sequence[21:]
        records.append(f">site_{index + 1}\n{sequence}")
    fasta.write_text("\n".join(records) + "\n", encoding="utf-8")

    known = root / "known.jaspar"
    known.write_text(
        ">TEST1 ACGTAC\n"
        "A [ 20 0 0 0 20 0 ]\n"
        "C [ 0 20 0 0 0 20 ]\n"
        "G [ 0 0 20 0 0 0 ]\n"
        "T [ 0 0 0 20 0 0 ]\n",
        encoding="utf-8",
    )
    return fasta, known


def _run_discovery(
    executable: Path,
    fasta: Path,
    known: Path,
    outdir: Path,
    environment: dict[str, str],
    timeout: float,
) -> None:
    command = [
        str(executable),
        INTERNAL_COMMAND,
        "discover-motifs",
        "--fasta",
        str(fasta),
        "--outdir",
        str(outdir),
        "--method",
        "streme",
        "--known-motifs",
        str(known),
        "--execute",
        "--runtime",
        "managed",
        "--extra-args",
        "--dna",
        "--nmotifs",
        "1",
        "--minw",
        "4",
        "--maxw",
        "8",
        "--seed",
        "1",
    ]
    result = subprocess.run(
        command,
        cwd=outdir.parent,
        env=environment,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    if result.returncode:
        raise SystemExit(
            "Frozen managed motif discovery failed:\n"
            f"{result.stdout}\n{result.stderr}"
        )

    required = [
        outdir / "streme" / "streme.txt",
        outdir / "known_motifs" / "known.meme",
        outdir / "tomtom" / "tomtom.tsv",
        outdir / "motif_summary.tsv",
        outdir / "motif_summary.html",
    ]
    missing = [
        str(path)
        for path in required
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise SystemExit(
            "Frozen managed motif discovery omitted required outputs: "
            + ", ".join(missing)
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("--timeout", type=float, default=900.0)
    args = parser.parse_args()

    executable = args.executable.expanduser().resolve()
    if not executable.is_file():
        raise SystemExit(f"Released desktop executable does not exist: {executable}")

    with tempfile.TemporaryDirectory(prefix="fp-tools-release-meme-") as directory:
        root = Path(directory)
        fasta, known = _write_inputs(root)
        cache = root / "runtime-cache"
        environment = os.environ.copy()
        for variable in (
            "FP_TOOLS_CA_BUNDLE",
            "SSL_CERT_FILE",
            "CURL_CA_BUNDLE",
            "REQUESTS_CA_BUNDLE",
        ):
            environment.pop(variable, None)
        environment["FP_TOOLS_RUNTIME_CACHE"] = str(cache)

        _run_discovery(
            executable,
            fasta,
            known,
            root / "first_run",
            environment,
            args.timeout,
        )
        markers = list(cache.rglob(".fp-tools-runtime.json"))
        if len(markers) != 1:
            raise SystemExit(
                f"Expected one clean-cache MEME runtime installation; found {len(markers)}"
            )
        marker_mtime = markers[0].stat().st_mtime_ns

        _run_discovery(
            executable,
            fasta,
            known,
            root / "second_run",
            environment,
            args.timeout,
        )
        if markers[0].stat().st_mtime_ns != marker_mtime:
            raise SystemExit("The second app invocation unexpectedly reinstalled the runtime")

    print("Released macOS app passed clean-cache and cached managed MEME discovery")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
