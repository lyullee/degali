"""Source bookkeeping tests; no physical initial turbulence is selected."""

import math
import numpy as np
import pytest
from test_enriched_transport import phase, candidate, section, mixing, reservoir, projection, transport
from degali.addons.prescribed_tke_source import prescribed_q_flux, retract_source_with_prescribed_tke
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments


def test_rectangular_q_flux_matches_independent_phase_grid(projection):
    qp = np.r_[math.log(2.), np.zeros(projection.basis.size)]
    actual = prescribed_q_flux(projection, qp)
    expected = 0.
    for p in FaceSplitSquareMoments(projection).rule(np.zeros(projection.count), order=4, angular_order=8):
        expected += float(p['weights'] @ (2*np.exp(-projection.section.velocity_shape_exponent*p['q'])*p['u']))
    assert actual == pytest.approx(expected, rel=1e-10, abs=1e-10)


def test_known_total_target_is_not_double_counted(projection):
    par = np.zeros(projection.count)
    qp = np.r_[math.log(2.), np.zeros(projection.basis.size)]
    target = FaceSplitSquareMoments(projection).moments(par)
    target[3] += prescribed_q_flux(projection, qp)
    old_target, old_scales = projection.target.copy(), projection.scales.copy()
    out = retract_source_with_prescribed_tke(projection, par, qp, total_moment_target=target)
    assert out['numerical_passed']
    assert out['parameters'] == pytest.approx(par, abs=1e-12)
    assert out['actual_total_moments'] == pytest.approx(target, rel=1e-9, abs=1e-8)
    np.testing.assert_array_equal(projection.target, old_target)
    np.testing.assert_array_equal(projection.scales, old_scales)
    assert not out['physical_initialization_passed'] and not out['q_was_fitted']


def test_positive_q_reallocates_energy_without_changing_total_target(projection):
    par = np.zeros(projection.count)
    qp = np.r_[math.log(2.), np.zeros(projection.basis.size)]
    target = FaceSplitSquareMoments(projection).moments(par)
    out = retract_source_with_prescribed_tke(projection, par, qp, total_moment_target=target)
    assert out['numerical_passed']
    assert max(abs(out['parameters'])) > 1e-8
    assert max(out['moment_errors']) < 1e-8
    assert out['actual_mean_thermal_moments'][3] == pytest.approx(target[3]-out['tke_flux'], rel=1e-9, abs=1e-7)
    center = projection.prepare(np.array([0.]), np.array([0.]))
    old, new = [projection.fields(v, center) for v in (par, out['parameters'])]
    for key in ('c', 'h', 'rho'):
        assert new[key] == pytest.approx(old[key], rel=1e-12)


def test_no_implicit_tke_or_mismatched_target_accepted(projection):
    with pytest.raises(ValueError, match='explicit finite'):
        prescribed_q_flux(projection, np.zeros(3))
    with pytest.raises(TypeError):
        retract_source_with_prescribed_tke(projection, np.zeros(projection.count), np.zeros(projection.basis.size+1))
    with pytest.raises(ValueError, match='six finite'):
        retract_source_with_prescribed_tke(projection, np.zeros(projection.count), np.zeros(projection.basis.size+1),
                                          total_moment_target=np.zeros(5))
