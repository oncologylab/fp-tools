from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from simulate_gamma_poisson_background import case


def test_zero_rate_is_explicitly_prior_supported_and_counts_cover_zero():
    row = case(0, .1, 'poisson_working_model', replicates=1000, rng=np.random.default_rng(1))
    assert row['coverage'] == row['prior_only_fraction'] == 1
    assert row['true_target_mean'] == 0 and row['mean_predicted_count'] == pytest.approx(.05)
    assert not row['functional_uncertainty_gate_eligible']


def test_fixed_seed_reproduces_results_and_stress_is_not_implicitly_qualified():
    kwargs = dict(background_mean=100, ratio=1, generator='overdispersed_center_stress_phi1', replicates=2000)
    first = case(**kwargs, rng=np.random.default_rng(2))
    assert first == case(**kwargs, rng=np.random.default_rng(2))
    assert first['covered'] + first['below_interval'] + first['above_interval'] == 2000
    assert first['coverage_upper_95'] < .95


def test_invalid_generator_or_scenario_rejects():
    with pytest.raises(ValueError):
        case(1, 1, 'unknown', replicates=100, rng=np.random.default_rng(1))
    with pytest.raises(ValueError):
        case(1, 1, 'poisson_working_model', replicates=0, rng=np.random.default_rng(1))
