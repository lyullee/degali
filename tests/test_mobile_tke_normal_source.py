import math
import numpy as np
import pytest
from test_enriched_transport import phase, candidate, section, mixing, reservoir, projection, transport
from degali.addons.mobile_tke_normal_source import MobileTkeNormalSource
from degali.addons.exact_transverse_geometry import exact_geometry_field_view
from degali.addons.enriched_segments import encode_enriched
from degali.addons.enriched_transport import EnrichedModalTransport
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments
from degali.addons.prescribed_normal_source import prescribed_normal_fluxes


@pytest.fixture
def source(transport):
    jp = transport.projection.section.jetplume
    for name, value in dict(ustar=.2, zr=.1, rml=0., spread_floor=False).items():
        setattr(jp, name, value)
    p, par, _ = exact_geometry_field_view(transport.projection,
                                         encode_enriched(transport.projection, transport.parameters))
    base = EnrichedModalTransport(p, par, scalar_mixing=transport.mixing, thermal_species_ratio=1.,
                                 mechanical_work='reduced_buoyancy_work_immediate_shear_heat')
    return MobileTkeNormalSource(base, np.r_[math.log(2.), np.zeros(base.size)],
        axial_variance_fraction=2/3, pressure_assumption='ambient_pressure_no_compensation')


def test_six_mobile_directions_keep_location_direction_and_Q_fixed(source):
    change = np.array([.001, -.001, .002, -.002, .0001, -.0002])
    model, encoded, _ = source.field_view(change)
    assert encoded[[3, 5, 6]] == pytest.approx(source.origin[[3, 5, 6]], abs=0.)
    prep = model.projection.prepare(np.array([0., .2, .6, 1.]), np.array([0., .4, .1, .7]))
    assert prep['psi'] @ (model.parameters[:model.size]-source.base.parameters[:model.size]) == pytest.approx(change[4]*prep['q'], abs=1e-14)
    assert not np.array_equal(model.projection.mixing.state, source.base.projection.mixing.state)


def test_actual_moments_recover_independent_fixed_source_and_normal_budget(source):
    result = source.evaluate(np.zeros(6), order=8, angular_order=8)
    p, par = source.base.projection, source.base.parameters
    expected = FaceSplitSquareMoments(p).moments(par, order=8, angular_order=8)
    flux = prescribed_normal_fluxes(p, source.q_parameters, axial_variance_fraction=source.ratio)
    expected[2] += flux['axial_normal_momentum']
    expected[3] += flux['tke_flux']+flux['axial_normal_work']
    np.testing.assert_allclose(result['moments'], expected, rtol=1e-12, atol=1e-9)
    assert sum(result['energy_terms'].values()) == pytest.approx(result['moments'][3], rel=1e-13)


def test_moving_source_jacobian_matches_all_six_actual_differences(source):
    change = np.array([.001, -.001, .002, -.002, .0001, -.0002])
    result = source.evaluate(change, order=4, angular_order=8, jacobian=True)
    scales = np.maximum(abs(result['moments']), [1e-12, 1e-12, 1., 1., 1e-3, 1.])
    differences = []
    for direction in np.eye(6):
        step = 2e-6
        plus = source.evaluate(change+step*direction, order=4, angular_order=8)['moments']
        minus = source.evaluate(change-step*direction, order=4, angular_order=8)['moments']
        differences.append((plus-minus)/(2*step))
    actual = np.column_stack(differences)/scales[:, None]
    analytic = result['jacobian']/scales[:, None]
    np.testing.assert_allclose(analytic, actual, rtol=2e-5, atol=2e-6)


def test_known_Q_normal_target_does_not_move_source(source):
    target = source.evaluate(np.zeros(6), order=8, angular_order=16)['moments']
    old_target = source.base.projection.target.copy()
    result = source.solve(total_moment_target=target)
    assert result['numerical_passed'] and max(abs(result['changes'])) < 1e-10
    assert not result['old_rates_reused'] and not result['physical_initialization_passed']
    assert not result['normal_transport_closed'] and not result['scalar_mixing_recomputed']
    np.testing.assert_array_equal(old_target, source.base.projection.target)


def test_mobile_centers_pay_for_large_Q_normal_source_without_new_energy(source):
    p, par = source.base.projection, source.base.parameters
    target = FaceSplitSquareMoments(p).moments(par, order=8, angular_order=8)
    result = source.solve(total_moment_target=target)
    assert result['numerical_passed'] and max(result['moment_errors']) < 1e-8
    assert max(abs(result['changes'][:4])) > 1e-7
    assert result['energy_terms']['tke'] > 0.
    assert sum(result['energy_terms'].values()) == pytest.approx(target[3], rel=1e-9)
    assert not result['q_was_fitted'] and not result['variance_ratio_was_fitted']


def test_invalid_pressure_Q_ratio_target_and_trust_region_rejected(source):
    base, q = source.base, source.q_parameters
    for ratio, pressure, message in ((2.1, 'ambient_pressure_no_compensation', 'Rss'),
                                     (2/3, 'unselected', 'pressure')):
        with pytest.raises(ValueError, match=message):
            MobileTkeNormalSource(base, q, axial_variance_fraction=ratio, pressure_assumption=pressure)
    with pytest.raises(ValueError, match='trust region'):
        source.field_view(np.array([.1001, 0., 0., 0., 0., 0.]))
    with pytest.raises(ValueError, match='trust region'):
        source.field_view(np.array([0., 0., 0., 0., source.upper[4]+1e-6, 0.]))
    with pytest.raises(ValueError, match='six finite targets'):
        source.solve(total_moment_target=np.zeros(6))
    with pytest.raises(ValueError, match='log-Q'):
        MobileTkeNormalSource(base, [np.nan], axial_variance_fraction=2/3,
                             pressure_assumption='ambient_pressure_no_compensation')
