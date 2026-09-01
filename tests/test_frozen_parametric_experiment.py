from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import freeze_parametric_holdouts  # noqa: E402
import evaluate_parametric_factorization  # noqa: E402
import run_frozen_parametric_experiment  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def _record(
    accession: str,
    *,
    step: str,
    preferred: bool = False,
    major: int = 1,
    minor: int = 0,
) -> dict:
    return {
        "accession": accession,
        "status": "released",
        "assembly": "GRCh38",
        "file_format": "bed",
        "output_type": "conservative IDR thresholded peaks",
        "biological_replicates": [1, 2],
        "preferred_default": preferred,
        "analysis_step_version": {
            "minor_version": minor,
            "name": f"{step}-v-{major}-{minor}",
            "analysis_step": {
                "step_label": step,
                "major_version": major,
                "pipelines": [],
            },
        },
    }


def test_replicate_aware_selector_excludes_pseudoreplicates() -> None:
    replicated = _record("REPLICATED", step="tf-chip-seq-replicated-idr-step")
    pseudo = _record("PSEUDO", step="tf-chip-seq-pooled-pseudoreplicated-idr-step")
    assert freeze_parametric_holdouts.is_replicate_aware_idr(replicated)
    assert not freeze_parametric_holdouts.is_replicate_aware_idr(pseudo)
    replicated["biological_replicates"] = [1]
    assert not freeze_parametric_holdouts.is_replicate_aware_idr(replicated)


def test_locked_candidate_tie_break_is_preferred_version_accession() -> None:
    preferred = _record("A", step="tf-idr-step", preferred=True, major=1)
    newer = _record("B", step="tf-idr-step", preferred=False, major=2)
    assert (
        max([preferred, newer], key=freeze_parametric_holdouts.candidate_rank)
        is preferred
    )
    newer_preferred = _record("B", step="tf-idr-step", preferred=True, major=2)
    assert (
        max([preferred, newer_preferred], key=freeze_parametric_holdouts.candidate_rank)
        is newer_preferred
    )
    same = _record("Z", step="tf-idr-step", preferred=True, major=2)
    assert (
        max([newer_preferred, same], key=freeze_parametric_holdouts.candidate_rank)
        is same
    )


def test_committed_holdout_freeze_locks_metadata_before_labels() -> None:
    freeze_path = (
        ROOT
        / "benchmarks"
        / "manifests"
        / "frozen_parametric_factorization_v1.freeze.json"
    )
    document = json.loads(freeze_path.read_text(encoding="utf-8"))
    assert document["chipped_peak_contents_read"] is False
    assert document["holdout_labels_scored"] is False
    assert document["manifest"]["rows"] == 17
    for key in ("study", "selector", "manifest"):
        path = Path(document[key]["path"])
        if not path.is_absolute():
            path = ROOT / path
        digest = sha256(path.read_bytes()).hexdigest()
        assert digest == document[key]["sha256"]
    decisions = {
        (row["cell"], row["tf"]): row["status"] for row in document["chip_decisions"]
    }
    assert decisions[("GM23338", "REST")].startswith("ineligible")
    assert decisions[("GM23338", "NANOG")] == "ineligible_no_jaspar_motif"
    assert decisions[("SK-N-SH", "MEF2A")] == "selected_pending_power_check"


def test_factorization_geometry_score_detects_center_protection() -> None:
    positions = np.arange(-100, 101, dtype=float)
    profiles = np.zeros((2, len(positions)), dtype=float)
    profiles[0, np.abs(positions) <= 7] = -2.0
    profiles[1, np.abs(positions) <= 7] = 2.0
    scores = evaluate_parametric_factorization.geometry_score(profiles, positions)
    assert scores[0] > 0
    assert scores[1] < 0


def test_factorization_dwm_loader_accepts_verified_cache_and_rejects_parametric_label(
    tmp_path: Path,
) -> None:
    cache = tmp_path / "DWM.expected.npz"
    np.savez_compressed(
        cache,
        profiles=np.ones((2, 5)),
        valid=np.ones(2, dtype=bool),
        site_hash=np.asarray([1, 2], dtype=np.uint64),
        signal_identity=np.asarray("/run/fp_tools_dwm/sample_expected.bw:1:2"),
    )
    baseline, inputs = evaluate_parametric_factorization.load_dwm_baseline(cache)
    assert baseline["expected"].shape == (2, 5)
    assert inputs == [cache]

    invalid = tmp_path / "LOG81.expected.npz"
    np.savez_compressed(
        invalid,
        profiles=np.ones((2, 5)),
        valid=np.ones(2, dtype=bool),
        site_hash=np.asarray([1, 2], dtype=np.uint64),
        signal_identity=np.asarray("/run/loglinear81/sample_expected.bw:1:2"),
    )
    try:
        evaluate_parametric_factorization.load_dwm_baseline(invalid)
    except ValueError as error:
        assert "conventional DWM" in str(error)
    else:
        raise AssertionError("parametric baseline was mislabeled as DWM")


def test_residual_selector_uses_difficult_tasks_and_ctcf_gate() -> None:
    rows = []
    for residual, difficult_ap, ctcf_auc in (
        ("deviance", 0.70, 0.79),
        ("pearson", 0.71, 0.75),
        ("difference", 0.69, 0.80),
        ("log-ratio", 0.68, 0.80),
        ("nb-center-flank", 0.67, 0.80),
    ):
        rows.extend(
            [
                {
                    "cell": "K562",
                    "tf": "MEF2A",
                    "role": "difficult",
                    "method": f"factorized_residual_{residual}",
                    "auroc": 0.7,
                    "auprc": difficult_ap,
                },
                {
                    "cell": "K562",
                    "tf": "CTCF",
                    "role": "positive_control",
                    "method": f"factorized_residual_{residual}",
                    "auroc": ctcf_auc,
                    "auprc": 0.8,
                },
            ]
        )
    rows.extend(
        [
            {
                "cell": "K562",
                "tf": "MEF2A",
                "role": "difficult",
                "method": "DWM",
                "auroc": 0.6,
                "auprc": 0.5,
            },
            {
                "cell": "K562",
                "tf": "CTCF",
                "role": "positive_control",
                "method": "DWM",
                "auroc": 0.8,
                "auprc": 0.8,
            },
        ]
    )
    selected, summary = evaluate_parametric_factorization.select_residual(
        pd.DataFrame(rows)
    )
    # Pearson has the best difficult-task AP but violates the -0.02 CTCF gate.
    assert selected == "deviance"
    assert (
        summary.loc[summary["residual"] == "pearson", "passes_ctcf_gate"].item()
        is False
    )


def test_bigwig_integrity_rejects_empty_and_accepts_covered(tmp_path: Path) -> None:
    pybigwig = __import__("pyBigWig")
    empty = tmp_path / "empty.bw"
    handle = pybigwig.open(str(empty), "w")
    handle.addHeader([("chr1", 100)])
    handle.close()
    try:
        run_frozen_parametric_experiment.validate_bigwig(empty)
    except ValueError as error:
        assert "no covered bases" in str(error)
    else:  # pragma: no cover - documents the fail-closed contract
        raise AssertionError("empty bigWig unexpectedly passed validation")

    covered = tmp_path / "covered.bw"
    handle = pybigwig.open(str(covered), "w")
    handle.addHeader([("chr1", 100)])
    handle.addEntries(["chr1"], [10], ends=[20], values=[1.5])
    handle.close()
    result = run_frozen_parametric_experiment.validate_bigwig(covered)
    assert result["covered_bases"] == 10


def test_stage_runner_completes_and_resumes_synthetic(tmp_path: Path) -> None:
    runner = run_frozen_parametric_experiment.StageRunner(
        study_path=ROOT
        / "benchmarks/manifests/frozen_parametric_factorization_v1.spec.json",
        holdout_freeze_path=ROOT
        / "benchmarks/manifests/frozen_parametric_factorization_v1.freeze.json",
        registry_path=ROOT
        / "benchmarks/manifests/frozen_parametric_local_inputs_v1.json",
        outdir=tmp_path / "run",
    )
    assert runner.run_stage("synthetic") == "completed"
    assert runner.run_stage("synthetic") == "resumed_verified"
    completion = json.loads(
        (tmp_path / "run/synthetic/stage.complete.json").read_text()
    )
    assert (
        completion["schema"] == run_frozen_parametric_experiment.STAGE_COMPLETE_SCHEMA
    )
    assert all(Path(row["path"]).is_file() for row in completion["outputs"])
