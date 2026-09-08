#!/usr/bin/env python3
"""Create a deterministic nested BAM depth/seed matrix in one source scan.

The ordinary single-output downsampler is deliberately simple, but rerunning it
for every depth and seed repeatedly scans the same large public BAM. This helper
opens every missing output together, hashes each fragment name once per seed,
and writes the read to all qualifying nested depth subsets.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from downsample_bam_fragments import (  # noqa: E402
    fragment_representative,
    fragment_uniform,
    usable_alignment,
)


@dataclass(frozen=True)
class MatrixTarget:
    target_fragments: int
    seed: int
    output_bam: Path

    @property
    def metadata_path(self) -> Path:
        return self.output_bam.with_suffix(".downsampling.json")


def depth_label(depth: int) -> str:
    if depth <= 0:
        raise ValueError("depth must be positive")
    if depth % 1_000_000 == 0:
        return f"{depth // 1_000_000}m"
    return str(depth)


def build_targets(
    outdir: str | Path,
    sample: str,
    depths: list[int],
    seeds: list[int],
) -> list[MatrixTarget]:
    if not sample.strip():
        raise ValueError("sample must not be empty")
    targets = []
    root = Path(outdir)
    for depth in sorted(set(int(value) for value in depths)):
        label = depth_label(depth)
        for seed in sorted(set(int(value) for value in seeds)):
            subset_id = f"{sample}.{label}.s{seed}"
            output = root / "signals" / sample / label / f"seed_{seed}" / f"{subset_id}.bam"
            targets.append(MatrixTarget(depth, seed, output))
    return targets


def validate_targets(targets: list[MatrixTarget], available_fragments: int) -> None:
    if available_fragments <= 0:
        raise ValueError("available_fragments must be positive")
    if not targets:
        raise ValueError("at least one depth/seed target is required")
    outputs = [target.output_bam.resolve() for target in targets]
    if len(outputs) != len(set(outputs)):
        raise ValueError("matrix contains duplicate output BAM paths")
    if any(target.target_fragments <= 0 for target in targets):
        raise ValueError("target fragment counts must be positive")
    excessive = [target.target_fragments for target in targets if target.target_fragments > available_fragments]
    if excessive:
        raise ValueError(
            f"target depth {max(excessive):,} exceeds {available_fragments:,} available fragments"
        )


def selected_depths(
    query_name: str,
    seed: int,
    target_fragments: list[int],
    available_fragments: int,
) -> list[int]:
    """Return nested selected depths for a fragment; exposed for tests."""

    value = fragment_uniform(query_name, seed)
    return [
        depth
        for depth in sorted(set(target_fragments))
        if value < min(1.0, depth / available_fragments)
    ]


def _index_bam(path: Path, threads: int) -> None:
    import pysam

    if threads > 1:
        pysam.index("-@", str(threads), str(path))
    else:
        pysam.index(str(path))


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata_identity(document: dict) -> str:
    payload = {key: value for key, value in document.items() if key != "metadata_id"}
    return sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def verified_existing_target(target: MatrixTarget, source: Path, source_digest: str,
                             available_fragments: int) -> dict | None:
    index = Path(str(target.output_bam) + ".bai")
    paths = (target.output_bam, index, target.metadata_path)
    partials = [path.with_name(path.name + ".partial") for path in (target.output_bam, target.metadata_path)]
    if any(path.exists() for path in partials):
        raise ValueError(f"incomplete matrix output must be reviewed before resume: {target.output_bam}")
    if not any(path.exists() for path in paths):
        return None
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise ValueError(f"incomplete matrix output must be reviewed before resume: {target.output_bam}")
    metadata = json.loads(target.metadata_path.read_text())
    expected = {"schema": "fp-tools-nested-depth-matrix-v1", "source_bam": str(source),
                "source_bam_sha256": source_digest, "available_fragments": available_fragments,
                "target_fragments": target.target_fragments, "seed": target.seed,
                "selection_rule": "blake2b(seed,query_name)<sampling_fraction",
                "sampling_fraction": target.target_fragments / available_fragments,
                "nested_with_same_seed": True}
    if (metadata.get("metadata_id") != metadata_identity(metadata)
            or any(metadata.get(key) != value for key, value in expected.items())
            or Path(metadata.get("output_bam", "")).resolve() != target.output_bam.resolve()):
        raise ValueError(f"matrix metadata or source differs from the requested sampling: {target.metadata_path}")
    for path, field in ((target.output_bam, "output_bam_sha256"), (index, "output_bai_sha256")):
        if file_sha256(path) != metadata.get(field):
            raise ValueError(f"matrix output checksum mismatch: {path}")
    return {**metadata, "state": "skipped_existing"}


def downsample_matrix(
    source: str | Path,
    targets: list[MatrixTarget],
    *,
    available_fragments: int,
    expected_source_sha256: str | None = None,
    index_threads: int = 1,
    progress_every_reads: int = 10_000_000,
) -> list[dict[str, object]]:
    """Write all missing matrix targets while traversing ``source`` once."""

    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - environment-dependent message
        raise RuntimeError("pysam is required to downsample BAM files") from exc
    validate_targets(targets, available_fragments)
    source_path = Path(source).resolve()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    if any(target.output_bam.resolve() == source_path for target in targets):
        raise ValueError("matrix output cannot replace the source BAM")
    source_digest = file_sha256(source_path)
    if expected_source_sha256 is not None and source_digest != expected_source_sha256:
        raise ValueError("source BAM checksum differs from the frozen input")
    pending = []
    records: list[dict[str, object]] = []
    for target in targets:
        existing = verified_existing_target(target, source_path, source_digest, available_fragments)
        if existing is None:
            pending.append(target)
        else:
            records.append(existing)
    if not pending:
        return records

    by_seed: dict[int, list[MatrixTarget]] = {}
    for target in pending:
        target.output_bam.parent.mkdir(parents=True, exist_ok=True)
        by_seed.setdefault(target.seed, []).append(target)
    for seed_targets in by_seed.values():
        seed_targets.sort(key=lambda target: target.target_fragments)

    temporary = {
        target: target.output_bam.with_name(target.output_bam.name + ".partial")
        for target in pending
    }
    selected = {target: 0 for target in pending}
    started = perf_counter()
    reads_seen = 0
    try:
        with pysam.AlignmentFile(str(source_path), "rb") as source_bam, ExitStack() as stack:
            outputs = {
                target: stack.enter_context(
                    pysam.AlignmentFile(str(temporary[target]), "wb", template=source_bam)
                )
                for target in pending
            }
            for read in source_bam.fetch(until_eof=True):
                reads_seen += 1
                if not usable_alignment(read):
                    continue
                name = str(read.query_name)
                for seed, seed_targets in by_seed.items():
                    value = fragment_uniform(name, seed)
                    for target in seed_targets:
                        if value < target.target_fragments / available_fragments:
                            outputs[target].write(read)
                            selected[target] += int(fragment_representative(read))
                if progress_every_reads and reads_seen % progress_every_reads == 0:
                    elapsed = max(perf_counter() - started, 1e-6)
                    print(
                        f"processed {reads_seen:,} alignments "
                        f"({reads_seen / elapsed:,.0f}/s)",
                        flush=True,
                    )
        for target in pending:
            temporary[target].replace(target.output_bam)
            _index_bam(target.output_bam, max(1, int(index_threads)))
            metadata = {
                "schema": "fp-tools-nested-depth-matrix-v1",
                "source_bam": str(source_path),
                "source_bam_sha256": source_digest,
                "output_bam": str(target.output_bam),
                "output_bam_sha256": file_sha256(target.output_bam),
                "output_bai_sha256": file_sha256(Path(str(target.output_bam) + ".bai")),
                "seed": target.seed,
                "available_fragments": available_fragments,
                "available_fragments_source": "provided",
                "target_fragments": target.target_fragments,
                "sampling_fraction": target.target_fragments / available_fragments,
                "selected_fragments": selected[target],
                "selection_rule": "blake2b(seed,query_name)<sampling_fraction",
                "nested_with_same_seed": True,
                "matrix_source_scans": 1,
                "completed_utc": datetime.now(timezone.utc).isoformat(),
            }
            metadata["metadata_id"] = metadata_identity(metadata)
            partial_metadata = target.metadata_path.with_name(target.metadata_path.name + ".partial")
            with partial_metadata.open("x", encoding="utf-8") as handle:
                handle.write(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
            partial_metadata.replace(target.metadata_path)
            records.append({**metadata, "state": "completed"})
    except BaseException:
        for path in temporary.values():
            if path.exists():
                path.unlink()
        raise
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--available-fragments", type=int, required=True)
    parser.add_argument("--depth", type=int, action="append", required=True)
    parser.add_argument("--seed", type=int, action="append", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path)
    parser.add_argument("--index-threads", type=int, default=2)
    parser.add_argument("--progress-every-reads", type=int, default=10_000_000)
    args = parser.parse_args(argv)
    targets = build_targets(args.outdir, args.sample, args.depth, args.seed)
    records = downsample_matrix(
        args.bam,
        targets,
        available_fragments=args.available_fragments,
        index_threads=args.index_threads,
        progress_every_reads=args.progress_every_reads,
    )
    manifest = args.manifest_out or (
        args.outdir / "signals" / args.sample / "depth_matrix_manifest.json"
    )
    manifest.parent.mkdir(parents=True, exist_ok=True)
    document = {
        "schema": "fp-tools-nested-depth-matrix-manifest-v1",
        "source_bam": str(args.bam.resolve()),
        "sample": args.sample,
        "available_fragments": args.available_fragments,
        "depths": sorted(set(args.depth)),
        "seeds": sorted(set(args.seed)),
        "records": sorted(records, key=lambda row: (int(row["target_fragments"]), int(row["seed"]))),
    }
    manifest.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for row in document["records"]:
        print(
            f"{row['state']}\t{row['target_fragments']}\t{row['seed']}\t{row['output_bam']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
