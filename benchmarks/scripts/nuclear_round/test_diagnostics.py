import numpy as np
import pytest
from diagnostics import simulation_available


def test_missing_perturbation_is_not_zero_power():
    assert simulation_available(np.zeros((2000,2,3),int))
    assert not simulation_available(np.full((2000,2,3),-1,int))
    with pytest.raises(ValueError):simulation_available(np.zeros((2000,3,2),int))
