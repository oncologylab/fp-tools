from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "extract_frozen_control_windows.py"
spec = importlib.util.spec_from_file_location("extract_frozen_control_windows", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_manifest_certifies_no_model_or_label_use(tmp_path) -> None:
    study = tmp_path / "study.json"
    study.write_text(json.dumps({"status": "development_locked_holdout_unscored"}))
    bam = tmp_path / "control.bam"
    bam.write_bytes(b"fixture")
    artifact = tmp_path / "windows.npz"
    artifact.write_bytes(b"npz")
    artifact.with_suffix(".json").write_text("{}")
    table = tmp_path / "control_windows.tsv"
    frame = pd.DataFrame({"cache_npz": [str(artifact)]})
    frame.to_csv(table, sep="\t", index=False)
    manifest = module.build_manifest(
        study=study,
        source="naked_dna",
        samples=[("rep2", bam)],
        window_manifest=frame,
        output_table=table,
    )
    assert manifest["models_fitted"] is False
    assert manifest["chipped_labels_used"] is False
    assert manifest["window_artifacts"][0]["sha256"] == module.file_sha256(artifact)
