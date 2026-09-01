#!/usr/bin/env python3
"""Run the checksum-locked frozen parametric factorization experiment.

Every stage receives an immutable pre-execution manifest.  A completed stage
is skipped only after all declared outputs pass format and checksum checks.
The runner is benchmark infrastructure and does not change any public command.
"""

from __future__ import annotations

import argparse
import csv
import glob
from hashlib import md5, sha256
import json
from pathlib import Path
import shlex
import subprocess
import sys
from time import perf_counter
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from sklearn.metrics import average_precision_score, roc_auc_score

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))

from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenBiasStrengthCalibrator,
    FrozenParametricFactorization,
)


RUN_SCHEMA = "fp-tools-frozen-parametric-run-v1"
STAGE_FREEZE_SCHEMA = "fp-tools-frozen-parametric-stage-freeze-v1"
STAGE_COMPLETE_SCHEMA = "fp-tools-frozen-parametric-stage-complete-v1"


def digest_file(path: str | Path, algorithm: str = "sha256") -> str:
    algorithm = algorithm.lower()
    if algorithm not in {"sha256", "md5"}:
        raise ValueError(f"unsupported checksum algorithm: {algorithm}")
    digest = sha256() if algorithm == "sha256" else md5()  # noqa: S324 - public-data integrity
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return sha256(encoded.encode()).hexdigest()


def resolve_path(path: str | Path, root: Path = REPOSITORY) -> Path:
    value = Path(path)
    return value if value.is_absolute() else root / value


def atomic_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def validate_internal_npz(path: Path) -> None:
    with np.load(path, allow_pickle=False) as arrays:
        for key in arrays.files:
            value = np.asarray(arrays[key])
            if value.dtype == object:
                raise ValueError(f"unsafe object array in {path}: {key}")


def validate_bigwig(path: Path) -> dict[str, Any]:
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - package runtime normally supplies it
        raise RuntimeError("pyBigWig is required to validate bigWig artifacts") from exc
    handle = pyBigWig.open(str(path))
    if handle is None:
        raise ValueError(f"cannot open bigWig: {path}")
    try:
        chromosomes = handle.chroms()
        header = handle.header()
        if not chromosomes or int(header.get("nBasesCovered", 0)) <= 0:
            raise ValueError(f"bigWig has no covered bases: {path}")
        for chromosome, length in list(chromosomes.items())[:3]:
            if int(length) <= 0:
                raise ValueError(f"bigWig chromosome has invalid length: {chromosome}")
            intervals = handle.intervals(chromosome, 0, min(int(length), 1_000_000))
            if intervals:
                for start, end, value in intervals[:100]:
                    if not (0 <= int(start) < int(end) <= int(length)) or not np.isfinite(value):
                        raise ValueError(f"invalid bigWig interval in {path}")
        return {"chromosomes": len(chromosomes), "covered_bases": int(header["nBasesCovered"])}
    finally:
        handle.close()


def validate_bam(path: Path) -> dict[str, Any]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pysam is required to validate BAM artifacts") from exc
    with pysam.AlignmentFile(str(path), "rb") as handle:
        if not handle.references or any(int(length) <= 0 for length in handle.lengths):
            raise ValueError(f"BAM header has no valid references: {path}")
        reference_count = len(handle.references)
        coordinate_sorted = handle.header.to_dict().get("HD", {}).get("SO") == "coordinate"
    index = Path(str(path) + ".bai")
    if not index.is_file():
        alternate = path.with_suffix(".bai")
        if not alternate.is_file():
            raise ValueError(f"BAM index is missing: {path}")
    return {"references": reference_count, "coordinate_sorted": coordinate_sorted}


def validate_artifact(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"artifact is absent or empty: {path}")
    suffixes = path.suffixes
    metadata: dict[str, Any] = {"bytes": int(path.stat().st_size)}
    if path.suffix == ".npz":
        validate_internal_npz(path)
        sidecar = path.with_suffix(".json")
        if sidecar.is_file():
            document = json.loads(sidecar.read_text(encoding="utf-8"))
            expected = document.get("npz_sha256") or document.get("profiles_sha256")
            if expected and expected != digest_file(path):
                raise ValueError(f"NPZ sidecar checksum mismatch: {path}")
    elif path.suffix == ".json":
        json.loads(path.read_text(encoding="utf-8"))
    elif path.suffix in {".tsv", ".csv"} or suffixes[-2:] == [".tsv", ".gz"]:
        separator = "\t" if ".tsv" in suffixes or path.suffix == ".tsv" else ","
        frame = pd.read_csv(path, sep=separator, nrows=2)
        if len(frame.columns) < 1:
            raise ValueError(f"tabular artifact has no columns: {path}")
        metadata["columns"] = int(len(frame.columns))
    elif path.suffix in {".bw", ".bigWig"}:
        metadata.update(validate_bigwig(path))
    elif path.suffix == ".bam":
        metadata.update(validate_bam(path))
    elif path.suffix == ".pdf":
        if not path.read_bytes()[:5] == b"%PDF-":
            raise ValueError(f"invalid PDF header: {path}")
    return metadata


def verify_declared_input(entry: dict[str, Any]) -> dict[str, Any]:
    path = resolve_path(entry["path"])
    if not path.is_file():
        if entry.get("optional", False):
            return {"path": str(path), "status": "optional_missing"}
        raise FileNotFoundError(path)
    stat = path.stat()
    expected_bytes = entry.get("bytes")
    if expected_bytes is not None and int(expected_bytes) != int(stat.st_size):
        raise ValueError(f"input size mismatch for {path}: {stat.st_size} != {expected_bytes}")
    checksum = str(entry.get("checksum") or "")
    algorithm = str(entry.get("checksum_algorithm") or ("md5" if len(checksum) == 32 else "sha256"))
    actual = digest_file(path, algorithm) if checksum else digest_file(path)
    if checksum and actual != checksum:
        raise ValueError(f"input checksum mismatch for {path}: {actual} != {checksum}")
    validate_artifact(path)
    return {
        "path": str(path),
        "bytes": int(stat.st_size),
        "checksum_algorithm": algorithm if checksum else "sha256",
        "checksum": actual,
        "status": "verified",
    }


def verify_holdout_freeze(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != "fp-tools-parametric-holdout-freeze-v1":
        raise ValueError("unsupported holdout freeze")
    if document.get("chipped_peak_contents_read") is not False or document.get("holdout_labels_scored") is not False:
        raise ValueError("holdout freeze does not certify unopened labels")
    for key in ("study", "selector", "manifest"):
        record = document[key]
        candidate = resolve_path(record["path"])
        if digest_file(candidate) != record["sha256"]:
            raise ValueError(f"holdout freeze input changed: {candidate}")
    canonical = {key: value for key, value in document.items() if key != "freeze_id"}
    if canonical_hash(canonical) != document.get("freeze_id"):
        raise ValueError("holdout freeze ID is invalid")
    return document


def expand_command(
    command: Iterable[str],
    *,
    stage_dir: Path,
    outdir: Path,
    study: Path,
    registry: Path,
) -> list[str]:
    replacements = {
        "{python}": str(REPOSITORY / ".venv" / "bin" / "python"),
        "{root}": str(REPOSITORY),
        "{stage_dir}": str(stage_dir),
        "{outdir}": str(outdir),
        "{study}": str(study),
        "{registry}": str(registry),
    }
    output = []
    for value in command:
        text = str(value)
        for token, replacement in replacements.items():
            text = text.replace(token, replacement)
        output.append(text)
    return output


def _synthetic_profiles(
    *,
    sites_per_tf: int,
    footprint: bool,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    positions = np.arange(-100, 101, dtype=float)
    count = 2 * sites_per_tf
    sample = np.repeat("synthetic", count)
    tfs = np.repeat(["TF_narrow", "TF_broad"], sites_per_tf)
    families = np.repeat(["family_A", "family_B"], sites_per_tf)
    labels = (np.arange(count) % 2 == 0).astype(int)
    log_bias = rng.normal(0.0, 0.45, size=(count, len(positions)))
    background = rng.normal(0.0, 0.25, size=count)[:, None] * positions[None, :] / 100.0
    effect = np.zeros_like(log_bias)
    if footprint:
        narrow = -1.8 * np.exp(-0.5 * np.square(positions / 6.0))
        broad = -1.3 * np.exp(-0.5 * np.square(positions / 13.0))
        shoulders = 0.45 * (
            np.exp(-0.5 * np.square((positions - 21.0) / 5.0))
            + np.exp(-0.5 * np.square((positions + 21.0) / 5.0))
        )
        curves = np.where((tfs == "TF_narrow")[:, None], narrow, broad) + shoulders
        effect = labels[:, None] * curves
    logits = 1.25 * log_bias + background + effect
    probabilities = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
    counts = np.vstack(
        [rng.multinomial(800, probability) for probability in probabilities]
    ).astype(float)
    return positions, counts, log_bias, sample, tfs, families, labels


def run_synthetic(stage_dir: Path, seed: int) -> list[Path]:
    rng = np.random.default_rng(seed)
    control_bias = rng.normal(0.0, 0.8, size=(400, 61))
    true_strength = 1.25
    probabilities = np.exp(
        true_strength * control_bias
        - logsumexp(true_strength * control_bias, axis=1, keepdims=True)
    )
    control_counts = np.vstack([rng.multinomial(80, row) for row in probabilities]).astype(float)
    calibration = FrozenBiasStrengthCalibrator().fit(control_counts, control_bias, ["synthetic"] * 400)

    positions, counts, log_bias, samples, tfs, families, labels = _synthetic_profiles(
        sites_per_tf=80,
        footprint=True,
        seed=seed + 1,
    )
    model = FrozenParametricFactorization(
        positions,
        family_shrinkage=5.0,
        tf_shrinkage=5.0,
        seed=seed,
    ).fit(
        counts,
        log_bias,
        samples,
        tfs,
        families,
        calibration,
        max_iter=12,
    )
    result = model.predict(counts, log_bias, samples, tfs)
    auroc = float(roc_auc_score(labels, result.posterior_bound))
    auprc = float(average_precision_score(labels, result.posterior_bound))
    total_error = float(np.max(np.abs(result.expected_unbound.sum(axis=1) - counts.sum(axis=1))))

    null_positions, null_counts, null_bias, null_samples, null_tfs, null_families, _ = _synthetic_profiles(
        sites_per_tf=40,
        footprint=False,
        seed=seed + 2,
    )
    null_logits = true_strength * null_bias
    null_probabilities = np.exp(null_logits - logsumexp(null_logits, axis=1, keepdims=True))
    null_counts = 5000.0 * null_probabilities
    null_model = FrozenParametricFactorization(
        null_positions,
        family_shrinkage=5.0,
        tf_shrinkage=5.0,
        use_total_component=False,
        seed=seed,
    ).fit(
        null_counts,
        null_bias,
        null_samples,
        null_tfs,
        null_families,
        {"synthetic": true_strength},
        max_iter=6,
    )
    null_result = null_model.predict(null_counts, null_bias, null_samples, null_tfs)
    center = np.abs(null_positions) <= 15
    maximum_null_effect = float(np.max(np.abs(null_result.footprint_log_effect[:, center])))
    recovered_strength = calibration.strength("synthetic")
    rows = [
        {"metric": "true_bias_strength", "value": true_strength, "gate": np.nan, "passed": True},
        {"metric": "recovered_bias_strength", "value": recovered_strength, "gate": 0.10, "passed": abs(recovered_strength - true_strength) <= 0.10},
        {"metric": "footprint_auroc", "value": auroc, "gate": 0.80, "passed": auroc >= 0.80},
        {"metric": "footprint_auprc", "value": auprc, "gate": 0.80, "passed": auprc >= 0.80},
        {"metric": "maximum_bias_only_center_effect", "value": maximum_null_effect, "gate": 0.10, "passed": maximum_null_effect <= 0.10},
        {"metric": "maximum_conditional_total_error", "value": total_error, "gate": 1e-8, "passed": total_error <= 1e-8},
    ]
    metrics_path = stage_dir / "synthetic_metrics.tsv"
    pd.DataFrame(rows).to_csv(metrics_path, sep="\t", index=False)
    profiles_path = stage_dir / "synthetic_profiles.npz"
    np.savez_compressed(
        profiles_path,
        positions=positions,
        labels=labels,
        posterior=result.posterior_bound,
        expected_unbound=result.expected_unbound,
        expected_bound=result.expected_bound,
        footprint_log_effect=result.footprint_log_effect,
        null_footprint_log_effect=null_result.footprint_log_effect,
    )
    calibration.save(stage_dir / "synthetic_calibration", {"labels_used": False})
    model.save(stage_dir / "synthetic_factorization", {"labels_used": False})

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))
    axes[0].plot(positions, np.mean(counts[labels == 0], axis=0), label="Unbound", lw=1.8)
    axes[0].plot(positions, np.mean(counts[labels == 1], axis=0), label="Bound", lw=1.8)
    axes[0].set(title="Observed synthetic cuts", xlabel="Position (bp)", ylabel="Mean cuts")
    axes[0].legend(frameon=False)
    for tf, color in (("TF_narrow", "#2166ac"), ("TF_broad", "#b2182b")):
        selected = tfs == tf
        axes[1].plot(positions, np.mean(result.footprint_log_effect[selected], axis=0), label=tf, color=color)
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].set(title="Recovered TF effects", xlabel="Position (bp)", ylabel="Log effect")
    axes[1].legend(frameon=False)
    axes[2].hist(result.posterior_bound[labels == 0], bins=20, alpha=0.65, label="Unbound")
    axes[2].hist(result.posterior_bound[labels == 1], bins=20, alpha=0.65, label="Bound")
    axes[2].set(title=f"Posterior (AUROC {auroc:.3f})", xlabel="P(bound)", ylabel="Sites")
    axes[2].legend(frameon=False)
    figure.tight_layout()
    pdf_path = stage_dir / "synthetic_validation.pdf"
    figure.savefig(pdf_path)
    plt.close(figure)
    if not all(bool(row["passed"]) for row in rows):
        raise RuntimeError("synthetic factorization gates failed; see synthetic_metrics.tsv")
    return [
        metrics_path,
        profiles_path,
        stage_dir / "synthetic_calibration.npz",
        stage_dir / "synthetic_calibration.json",
        stage_dir / "synthetic_factorization.npz",
        stage_dir / "synthetic_factorization.json",
        pdf_path,
    ]


class StageRunner:
    def __init__(
        self,
        *,
        study_path: Path,
        holdout_freeze_path: Path,
        registry_path: Path,
        outdir: Path,
        dry_run: bool = False,
    ):
        self.study_path = study_path.resolve()
        self.holdout_freeze_path = holdout_freeze_path.resolve()
        self.registry_path = registry_path.resolve()
        self.outdir = outdir.resolve()
        self.dry_run = bool(dry_run)
        self.study = json.loads(self.study_path.read_text(encoding="utf-8"))
        self.registry = json.loads(self.registry_path.read_text(encoding="utf-8"))
        self.holdout = verify_holdout_freeze(self.holdout_freeze_path)
        if self.registry.get("schema") != "fp-tools-frozen-parametric-input-registry-v1":
            raise ValueError("unsupported input registry")
        self.stages = list(self.study["stage_order"])
        branch = subprocess.run(
            ["git", "branch", "--show-current"],
            cwd=REPOSITORY,
            text=True,
            capture_output=True,
            check=True,
        ).stdout.strip()
        if branch != "research/footprint-improvement-20260830":
            raise RuntimeError(
                "the frozen experiment may run only on research/footprint-improvement-20260830"
            )
        self.outdir.mkdir(parents=True, exist_ok=True)
        run_document = {
            "schema": RUN_SCHEMA,
            "study": {"path": str(self.study_path), "sha256": digest_file(self.study_path)},
            "holdout_freeze": {"path": str(self.holdout_freeze_path), "sha256": digest_file(self.holdout_freeze_path)},
            "registry": {"path": str(self.registry_path), "sha256": digest_file(self.registry_path)},
            "runner": {"path": str(Path(__file__).resolve()), "sha256": digest_file(Path(__file__).resolve())},
            "branch_required": "research/footprint-improvement-20260830",
            "main_merge_allowed": False,
        }
        run_document["run_id"] = canonical_hash(run_document)
        run_path = self.outdir / "run.freeze.json"
        if run_path.exists():
            existing = json.loads(run_path.read_text(encoding="utf-8"))
            if existing != run_document:
                raise ValueError("run.freeze.json differs; use a new output directory")
        elif not self.dry_run:
            atomic_json(run_path, run_document)

    def _stage_config(self, stage: str) -> dict[str, Any]:
        configuration = self.registry.get("stages", {}).get(stage)
        if not isinstance(configuration, dict):
            raise ValueError(f"stage {stage!r} is not configured in {self.registry_path}")
        return configuration

    def _dependencies(self, stage: str, configuration: dict[str, Any]) -> list[str]:
        return list(configuration.get("depends_on", []))

    def _complete_path(self, stage: str) -> Path:
        return self.outdir / stage / "stage.complete.json"

    def _verify_completion(self, stage: str) -> bool:
        path = self._complete_path(stage)
        if not path.is_file():
            return False
        document = json.loads(path.read_text(encoding="utf-8"))
        if document.get("schema") != STAGE_COMPLETE_SCHEMA:
            raise ValueError(f"invalid stage completion document: {path}")
        for record in document.get("outputs", []):
            output = Path(record["path"])
            validate_artifact(output)
            if digest_file(output) != record["sha256"]:
                raise ValueError(f"completed-stage artifact changed: {output}")
        return True

    def status(self) -> pd.DataFrame:
        rows = []
        for stage in self.stages:
            configured = stage in self.registry.get("stages", {})
            complete = self._verify_completion(stage) if configured else False
            rows.append({"stage": stage, "configured": configured, "complete": complete})
        return pd.DataFrame(rows)

    def run_stage(self, stage: str) -> str:
        if stage not in self.stages:
            raise ValueError(f"unknown stage: {stage}")
        configuration = self._stage_config(stage)
        for dependency in self._dependencies(stage, configuration):
            if not self._verify_completion(dependency):
                raise RuntimeError(f"stage {stage} requires completed stage {dependency}")
        if self._verify_completion(stage):
            return "resumed_verified"
        stage_dir = self.outdir / stage
        stage_dir.mkdir(parents=True, exist_ok=True)
        verified_inputs = [verify_declared_input(entry) for entry in configuration.get("inputs", [])]
        commands = [
            expand_command(
                command,
                stage_dir=stage_dir,
                outdir=self.outdir,
                study=self.study_path,
                registry=self.registry_path,
            )
            for command in configuration.get("commands", [])
        ]
        freeze = {
            "schema": STAGE_FREEZE_SCHEMA,
            "stage": stage,
            "study_sha256": digest_file(self.study_path),
            "holdout_freeze_id": self.holdout["freeze_id"],
            "registry_sha256": digest_file(self.registry_path),
            "runner_sha256": digest_file(Path(__file__).resolve()),
            "configuration": configuration,
            "verified_inputs": verified_inputs,
            "commands": commands,
        }
        freeze["stage_id"] = canonical_hash(freeze)
        freeze_path = stage_dir / "stage.freeze.json"
        if freeze_path.exists():
            existing = json.loads(freeze_path.read_text(encoding="utf-8"))
            if existing != freeze:
                raise ValueError(f"immutable stage freeze differs for {stage}; use a new output directory")
        elif not self.dry_run:
            atomic_json(freeze_path, freeze)
        if self.dry_run:
            for command in commands:
                print(shlex.join(command))
            return "dry_run"

        started = perf_counter()
        generated: list[Path] = []
        builtin = configuration.get("builtin")
        if builtin == "synthetic":
            generated.extend(run_synthetic(stage_dir, int(configuration.get("seed", 2026))))
        elif builtin:
            raise ValueError(f"unknown built-in stage: {builtin}")
        for index, command in enumerate(commands, start=1):
            log_path = stage_dir / f"command_{index:02d}.log"
            with log_path.open("w", encoding="utf-8") as log:
                log.write("$ " + shlex.join(command) + "\n")
                log.flush()
                completed = subprocess.run(
                    command,
                    cwd=REPOSITORY,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False,
                )
            if completed.returncode != 0:
                raise RuntimeError(f"stage {stage} command {index} failed; see {log_path}")
            generated.append(log_path)

        for pattern in configuration.get("outputs", []):
            expanded = str(pattern).replace("{stage_dir}", str(stage_dir)).replace("{outdir}", str(self.outdir))
            matches = [Path(value) for value in sorted(glob.glob(expanded, recursive=True))]
            if not matches:
                raise FileNotFoundError(f"stage {stage} output pattern matched nothing: {expanded}")
            generated.extend(matches)
        unique = sorted({path.resolve() for path in generated if path.is_file()})
        output_records = []
        for path in unique:
            details = validate_artifact(path)
            output_records.append(
                {
                    "path": str(path),
                    "sha256": digest_file(path),
                    **details,
                }
            )
        integrity_path = stage_dir / "artifact_integrity.tsv"
        with integrity_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"], delimiter="\t")
            writer.writeheader()
            for record in output_records:
                writer.writerow({key: record[key] for key in ("path", "bytes", "sha256")})
        output_records.append(
            {
                "path": str(integrity_path.resolve()),
                "sha256": digest_file(integrity_path),
                "bytes": integrity_path.stat().st_size,
            }
        )
        completion = {
            "schema": STAGE_COMPLETE_SCHEMA,
            "stage": stage,
            "stage_id": freeze["stage_id"],
            "elapsed_seconds": perf_counter() - started,
            "outputs": output_records,
        }
        atomic_json(self._complete_path(stage), completion)
        return "completed"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study",
        type=Path,
        default=REPOSITORY / "benchmarks/manifests/frozen_parametric_factorization_v1.spec.json",
    )
    parser.add_argument(
        "--holdout-freeze",
        type=Path,
        default=REPOSITORY / "benchmarks/manifests/frozen_parametric_factorization_v1.freeze.json",
    )
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--stage", default="status")
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    runner = StageRunner(
        study_path=args.study,
        holdout_freeze_path=args.holdout_freeze,
        registry_path=args.registry,
        outdir=args.outdir,
        dry_run=args.dry_run,
    )
    if args.stage == "status":
        print(runner.status().to_string(index=False))
        return 0
    stages = runner.stages if args.stage == "all" else [args.stage]
    for stage in stages:
        print(f"{stage}: {runner.run_stage(stage)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
