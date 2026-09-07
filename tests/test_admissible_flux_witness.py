import numpy as np
import pytest
from degali.addons.admissible_flux_witness import minimum_amplitude_witness


def test_already_feasible_flow_is_not_artificially_changed():
    out=minimum_amplitude_witness(np.eye(2),np.ones(2),np.array([2.,3.]))
    assert out['feasible']
    assert out['amplitudes']==pytest.approx([0.,0.])
    assert out['maximum_normalized_amplitude']==0.


def test_minimum_scaled_amplitude_satisfies_strict_constraints():
    out=minimum_amplitude_witness(np.array([[-1.,0.],[0.,1.]]),np.array([-.4,-.2]),np.array([2.,1.]))
    assert out['feasible']
    assert out['amplitudes']==pytest.approx([.4,-.2])
    assert out['maximum_normalized_amplitude']==pytest.approx(.2)
    assert out['maximum_inequality_violation']<1e-12


def test_infeasible_flux_constraints_are_not_relaxed():
    out=minimum_amplitude_witness(np.array([[1.],[-1.]]),np.array([0.,-1.]),np.ones(1))
    assert not out['feasible']
    assert 'amplitudes' not in out
