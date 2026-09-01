import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "evaluate_tf_geometry_naked_dna.py"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("evaluate_tf_geometry_naked_dna", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_frozen_upper_tail_threshold_respects_call_budget_with_ties():
    values = np.array([9.0, 9.0, 8.0, 7.0, 6.0] + list(range(15)))
    threshold = module.frozen_upper_tail_threshold(values, 0.10)
    assert np.sum(values >= threshold) <= 2


def test_frozen_upper_tail_threshold_rejects_empty_and_bad_alpha():
    with pytest.raises(ValueError, match="zero finite"):
        module.frozen_upper_tail_threshold(np.array([np.nan]), 0.05)
    with pytest.raises(ValueError, match="between zero and one"):
        module.frozen_upper_tail_threshold(np.array([1.0]), 1.0)


def test_evaluate_reports_independent_and_paired_signal_support(monkeypatch, tmp_path):
    candidate = SimpleNamespace(identifier="DWM.test", correction="DWM")
    monkeypatch.setattr(module, "candidate_from_row", lambda row: candidate)
    monkeypatch.setattr(
        module.np,
        "load",
        lambda path: {"profiles": np.zeros((4, 3), dtype=float)},
    )
    candidate_calls = iter(
        [np.array([0.0, 1.0, 2.0, 3.0]), np.array([4.0, 3.0, 2.0])]
    )
    monkeypatch.setattr(
        module, "score_candidate", lambda profiles, model: next(candidate_calls)
    )
    legacy_calls = iter(
        [np.array([0.0, 1.0, 2.0, 3.0]), np.array([4.0, np.nan, 2.0])]
    )
    monkeypatch.setattr(
        module, "score_centers", lambda sites, signal: next(legacy_calls)
    )
    monkeypatch.setattr(
        module,
        "extract_profiles",
        lambda sites, signal, flank: (
            np.zeros((3, 3), dtype=float),
            np.array([True, True, False]),
        ),
    )
    development_sites = pd.DataFrame(
        {
            "cell": ["K562"] * 4,
            "tf": ["MEF2A"] * 4,
            "chromosome_split": ["validation"] * 4,
            "chip_label": [0, 0, 1, 1],
        }
    )
    naked_sites = pd.DataFrame(
        {
            "tf": ["MEF2A"] * 3,
            "TFBS_chr": ["chr1"] * 3,
            "TFBS_start": [1, 2, 3],
            "TFBS_end": [2, 3, 4],
        }
    )
    summary, scores = module.evaluate(
        development_sites=development_sites,
        winners=pd.DataFrame(
            {"cell": ["K562"], "tf": ["MEF2A"], "candidate": ["DWM.test"]}
        ),
        development_baselines=pd.DataFrame(
            {"cell": ["K562"], "signal": [str(tmp_path / "legacy.bw")]}
        ),
        profile_cache=tmp_path,
        naked_sites={"K562": naked_sites},
        naked_corrected={"DWM": tmp_path / "corrected.bw"},
        naked_legacy=tmp_path / "legacy.bw",
        tf="MEF2A",
        split="validation",
        alpha=0.5,
        flank=1,
    )
    row = summary.iloc[0]
    assert row.candidate_valid == 2
    assert row.legacy_valid == 2
    assert row.paired_candidate_valid == 1
    assert row.paired_legacy_valid == 1
    assert scores.candidate_valid.tolist() == [True, True, False]
    assert scores.legacy_valid.tolist() == [True, False, True]
    assert scores.paired_valid.tolist() == [True, False, False]
