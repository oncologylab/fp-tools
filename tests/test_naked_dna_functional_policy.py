import importlib.util
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "evaluate_naked_dna_functional_policy.py"
spec = importlib.util.spec_from_file_location("evaluate_naked_dna_functional_policy", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_wilson_interval_contains_observed_rate():
    low, high = module.wilson_interval(5, 100)
    assert low < 0.05 < high
    empty_low, empty_high = module.wilson_interval(0, 0)
    assert math.isnan(empty_low) and math.isnan(empty_high)


def test_false_positive_summary_requires_signal_for_a_call():
    probabilities = np.array([0.99, 0.8, 0.2, np.nan])
    signal = np.array([0.0, 3.0, 4.0, 2.0])
    valid = np.array([True, True, True, True])
    summary, informative, calls = module.summarize_false_positives(
        probabilities, signal, valid, threshold=0.5
    )
    assert summary["sites_valid"] == 3
    assert summary["sites_informative"] == 2
    assert summary["false_positive_calls"] == 1
    assert summary["false_positive_rate"] == 1 / 3
    assert summary["informative_false_positive_rate"] == 1 / 2
    assert informative.tolist() == [False, True, True, False]
    assert calls.tolist() == [False, True, False, False]


def test_frozen_candidate_ids_resolve_to_implemented_models():
    strand, dwm = module._candidate_lookup()
    policy_path = (
        ROOT / "benchmarks" / "manifests" / "compact" / "functional_detector_policy_v1.tsv"
    )
    policy = pd.read_csv(policy_path, sep="\t")
    promoted = policy[policy["passes_development_gates"].astype(bool)]
    assert set(promoted["candidate_id"]).issubset(strand)
    assert set(promoted["reference_candidate_id"]).issubset(dwm)


def test_promotion_false_positive_table_maps_exact_paired_methods():
    rows = []
    for method, rate in (
        ("frozen_policy_candidate", 0.03),
        ("frozen_dwm_reference", 0.02),
    ):
        rows.append(
            {
                "cell": "K562",
                "tf": "MEF2A",
                "motif_family": "MEF2",
                "replicate": "rep1",
                "method": method,
                "false_positive_rate": rate,
                "informative_false_positive_rate": rate + 0.01,
                "sites_valid": 200,
                "sites_informative": 150,
            }
        )
    output = module.promotion_false_positive_table(pd.DataFrame(rows))
    assert set(output["candidate_id"]) == {"frozen_policy", "DWM_reference"}
    assert output["false_positive_rate"].sum() == pytest.approx(0.05)
    with pytest.raises(ValueError, match="methods differ"):
        module.promotion_false_positive_table(pd.DataFrame(rows[:1]))
