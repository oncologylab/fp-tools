from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import freeze_parametric_holdouts  # noqa: E402


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
    assert max([preferred, newer], key=freeze_parametric_holdouts.candidate_rank) is preferred
    newer_preferred = _record("B", step="tf-idr-step", preferred=True, major=2)
    assert max([preferred, newer_preferred], key=freeze_parametric_holdouts.candidate_rank) is newer_preferred
    same = _record("Z", step="tf-idr-step", preferred=True, major=2)
    assert max([newer_preferred, same], key=freeze_parametric_holdouts.candidate_rank) is same


def test_committed_holdout_freeze_locks_metadata_before_labels() -> None:
    freeze_path = ROOT / "benchmarks" / "manifests" / "frozen_parametric_factorization_v1.freeze.json"
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
    decisions = {(row["cell"], row["tf"]): row["status"] for row in document["chip_decisions"]}
    assert decisions[("GM23338", "REST")].startswith("ineligible")
    assert decisions[("GM23338", "NANOG")] == "ineligible_no_jaspar_motif"
    assert decisions[("SK-N-SH", "MEF2A")] == "selected_pending_power_check"
