import math
import numpy as np
import pytest
from test_enriched_transport import phase, candidate, section, mixing, reservoir, projection, transport
from degali.addons.prescribed_normal_source import prescribed_normal_fluxes, retract_source_with_normal_stress
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments


def test_supplied_equal_normal_variance_adds_both_fluxes(projection):
    q = np.r_[math.log(2.), np.zeros(projection.basis.size)]
    result = prescribed_normal_fluxes(projection, q, axial_variance_fraction=2/3)
    assert result['axial_normal_work'] == pytest.approx(2*result['tke_flux']/3)
    assert result['axial_normal_momentum'] > 0.
    assert result['shear_covariance_not_checked']
    with pytest.raises(ValueError, match='in'):
        prescribed_normal_fluxes(projection, q, axial_variance_fraction=2.1)


def test_known_Q_and_normal_source_total_is_recovered(projection):
    par = np.zeros(projection.count)
    q = np.r_[math.log(2.), np.zeros(projection.basis.size)]
    target = FaceSplitSquareMoments(projection).moments(par)
    flux = prescribed_normal_fluxes(projection, q, axial_variance_fraction=2/3)
    target[2] += flux['axial_normal_momentum']
    target[3] += flux['tke_flux']+flux['axial_normal_work']
    old_target = projection.target.copy()
    result = retract_source_with_normal_stress(projection, par, q, total_moment_target=target,
        axial_variance_fraction=2/3, pressure_assumption='ambient_pressure_no_compensation')
    assert result['numerical_passed']
    assert result['parameters'] == pytest.approx(par, abs=1e-12)
    assert not result['physical_initialization_passed'] and not result['normal_transport_closed']
    np.testing.assert_array_equal(old_target, projection.target)


def test_Q_and_normal_fluxes_are_paid_from_unchanged_source_budget(projection):
    par = np.zeros(projection.count)
    # Smaller manufactured Q stresses the coupled adjustment without asserting
    # existence for arbitrary turbulent energy at fixed centers/geometry.
    q = np.r_[math.log(.02), np.zeros(projection.basis.size)]
    target = FaceSplitSquareMoments(projection).moments(par)
    result = retract_source_with_normal_stress(projection, par, q, total_moment_target=target,
        axial_variance_fraction=2/3, pressure_assumption='ambient_pressure_no_compensation')
    assert result['numerical_passed']
    assert max(result['original_moment_errors']) < 1e-8
    assert max(abs(result['parameters'])) > 1e-9
    inner = result['q_only_inner_retraction']['actual_mean_thermal_moments']
    flux = result['source_fluxes']
    assert inner[2] == pytest.approx(target[2]-flux['axial_normal_momentum'], rel=1e-9)
    assert inner[3] == pytest.approx(target[3]-flux['tke_flux']-flux['axial_normal_work'], rel=1e-9)


def test_pressure_compensation_is_not_silently_invented(projection):
    with pytest.raises(ValueError, match='pressure'):
        retract_source_with_normal_stress(projection, np.zeros(projection.count), np.zeros(projection.basis.size+1),
            total_moment_target=projection.target, axial_variance_fraction=2/3, pressure_assumption='unselected')
