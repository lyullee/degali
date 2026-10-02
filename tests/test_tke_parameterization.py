import numpy as np
import pytest

from degali.addons.tke_parameterization import (
    dissipation_time_from_integral_scale,
    k_epsilon_dissipation_coefficient,
    tke_realizability_margin,
    tke_two_normal_rms_realizability_margin,
    tke_diffusivity_from_eddy_viscosity,
)


def test_dissipation_time_is_the_explicit_k_epsilon_identity():
    out = dissipation_time_from_integral_scale([.02, .04], [4., 9.], .2)
    np.testing.assert_allclose(out['dissipation_time'], [.05, 1./15.])
    assert out['closure_inputs_explicit']
    assert not out['adopted']
    assert not out['physical_closure_validated']


def test_tke_diffusivity_requires_its_own_transport_ratio():
    out = tke_diffusivity_from_eddy_viscosity(np.array([.6, .3]), 1.5)
    np.testing.assert_allclose(out['tke_diffusivity'], [.4, .2])
    assert out['closure_inputs_explicit']
    assert not out['adopted']


@pytest.mark.parametrize('function,args', [
    (dissipation_time_from_integral_scale, (.1, 0., .2)),
    (dissipation_time_from_integral_scale, (.1, 1., np.nan)),
    (tke_diffusivity_from_eddy_viscosity, (.1, -1.)),
    (k_epsilon_dissipation_coefficient, (0.,)),
])
def test_nonphysical_parameterizations_are_rejected(function, args):
    with pytest.raises(ValueError, match='strictly positive'):
        function(*args)


def test_standard_symbol_conversion_is_not_a_hidden_default():
    assert k_epsilon_dissipation_coefficient(.09) == pytest.approx(.09**.75)


def test_tke_screen_reports_a_psd_deficit_without_repairing_it():
    out = tke_realizability_margin(
        np.array([2.0, 5.0]), np.array([2.0, 2.0]),
        np.array([[0.0, 0.0], [3.0, 4.0]]),
    )
    np.testing.assert_allclose(out["minimum_tke"], [2.0, 5.125])
    np.testing.assert_allclose(out["realizability_margin"], [0.0, -0.125])
    assert not out["realizable"]
    assert not out["uses_isotropy"]
    assert not out["physical_closure_validated"]


def test_tke_screen_rejects_invalid_or_mismatched_fields():
    with pytest.raises(ValueError, match="one-dimensional"):
        tke_realizability_margin([[1.0]], [1.0], [[0.0, 0.0]])
    with pytest.raises(ValueError, match="matching shape"):
        tke_realizability_margin([1.0, 2.0], [1.0], [[0.0, 0.0]])


def test_two_normal_rms_screen_does_not_complete_or_repair_missing_stresses():
    out = tke_two_normal_rms_realizability_margin([6.0, 8.1], [3., 4.], [1., 0.])
    np.testing.assert_allclose(out['minimum_tke'], [5., 8.])
    np.testing.assert_allclose(out['realizability_margin'], [1., .1])
    assert out['realizable']
    assert not out['uses_isotropy']
    assert not out['physical_closure_validated']
