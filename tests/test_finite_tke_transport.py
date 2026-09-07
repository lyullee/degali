"""Manufactured finite-TKE operator tests, not chosen real-trial inputs."""

import math
import numpy as np
import pytest
from test_enriched_transport import phase, candidate, section, mixing, reservoir, projection, transport
from degali.addons.finite_tke_transport import FiniteTkeModalTransport
from degali.addons.enriched_transport import EnrichedModalTransport
from degali.addons.enriched_segments import encode_enriched
from degali.addons.exact_transverse_geometry import exact_geometry_field_view


def make_model(transport, **changes):
    supplied = dict(ambient_tke=.07, tke_diffusivity=.3, dissipation_time=.2,
                    circulation_amplitudes=np.zeros(4))
    supplied.update(changes)
    parameters = np.r_[math.log(2.), np.zeros(transport.size)]
    return FiniteTkeModalTransport(transport, parameters, **supplied)


@pytest.fixture
def finite(transport):
    return make_model(transport)


@pytest.mark.parametrize('key,value', [('ambient_tke', -.1), ('ambient_tke', np.nan),
    ('tke_diffusivity', 0.), ('tke_diffusivity', np.inf), ('dissipation_time', -.2),
    ('dissipation_time', None), ('circulation_amplitudes', [0., 0.])])
def test_invalid_or_missing_inputs_rejected(transport, key, value):
    with pytest.raises(ValueError):
        make_model(transport, **{key: value})
    with pytest.raises(TypeError):
        FiniteTkeModalTransport(transport, np.zeros(transport.size+1))


def test_callable_fields_and_pointwise_positivity(transport):
    model = make_model(transport, tke_diffusivity=lambda d: .2+.01*d['q'],
                       dissipation_time=lambda d: .3+.02*d['tke'])
    a, b = np.array([.1, .4]), np.array([.2, .7])
    d = model.local(a, b)
    assert d['chi_k'] == pytest.approx(.2+.01*d['q'])
    assert d['dissipation'] == pytest.approx(transport.area*d['q_density']/(.3+.02*d['tke']))
    model.dissipation_time = lambda d: np.ones(len(d['q'])+1)
    with pytest.raises(ValueError, match='matching field'):
        model.local(a, b)
    model.dissipation_time = lambda d: 1.-100*d['q']
    with pytest.raises(ValueError, match='positive'):
        model.local(a, b)


def test_tke_shape_or_underflow_is_rejected_not_clipped(finite):
    finite.parameters[1:] = 1.
    with pytest.raises(ValueError, match='trust region'):
        finite.local(np.array([1.]), np.array([1.]))
    finite.parameters[1:] = 0.
    finite.parameters[0] = -1000.
    with pytest.raises(ValueError, match='strictly positive'):
        finite.local(np.array([.3]), np.array([.4]))


def test_q_is_separate_from_thermal_eos_and_grad_k_matches_values(finite):
    a, b = np.array([.17, .38, .62]), np.array([.23, .41, .79])
    finite.parameters[2] = 2e-4
    d = finite.local(a, b)
    base = finite.base.local(a, b)
    for key in ('rho', 'c', 'h', 'y'):
        np.testing.assert_array_equal(d[key], base[key])
    step = 1e-6
    for axis in range(2):
        plus = finite.local(a+step*(axis == 0), b+step*(axis == 1))
        minus = finite.local(a-step*(axis == 0), b-step*(axis == 1))
        for key, gradient in (('q_density', 'grad_q'), ('tke', 'grad_k')):
            assert d[gradient][:, axis] == pytest.approx((plus[key]-minus[key])/(2*step), rel=3e-6, abs=1e-8)


def test_q_advective_jacobian_matches_actual_all_state_difference(transport, monkeypatch):
    p = transport.projection
    for name, value in dict(ustar=.2, zr=.1, rml=0., spread_floor=False).items():
        monkeypatch.setattr(p.section.jetplume, name, value, raising=False)
    encoded = encode_enriched(p, transport.parameters)
    view, params, _ = exact_geometry_field_view(p, encoded)
    mean = EnrichedModalTransport(view, params, scalar_mixing=transport.mixing,
        thermal_species_ratio=1., mechanical_work='reduced_buoyancy_work_immediate_shear_heat')
    model = make_model(mean)
    a, b = np.array([.17, .38, .62]), np.array([.23, .41, .79])
    direction = np.random.default_rng(691).normal(0., .03, model.count)
    step, actual = 1e-5, []
    for sign in (1, -1):
        shifted, sp, _ = exact_geometry_field_view(p, encoded+sign*step*direction[:mean.count])
        mb = EnrichedModalTransport(shifted, sp, scalar_mixing=transport.mixing,
            thermal_species_ratio=1., mechanical_work='reduced_buoyancy_work_immediate_shear_heat')
        trial = make_model(mb)
        trial.parameters += sign*step*direction[mean.count:]
        local = trial.local(a, b)
        actual.append(mb.area*local['q_density']*local['u'])
    expected = model.local(a, b)['bq'] @ direction
    assert expected == pytest.approx((actual[0]-actual[1])/(2*step), rel=2e-8, abs=2e-8)


def test_all_rows_retained_and_natural_face_energy_identity(finite):
    out = finite.assemble(order=4, angular_order=8, split_angles=False)
    assert finite.count == 32 and len(finite.active) == 30
    assert out['matrix'].shape == (30, 32)
    assert out['linear_scaled_error'] < 1e-8
    assert out['matrix_condition'] < 1e12
    led = out['ledger_coefficients']
    actual = led['heat_boundary']+led['mean_kinetic_boundary']+led['tke_boundary']
    expected = (.5*finite.base.projection.mixing.wind**2+finite.ambient_tke)*led['mass_boundary']
    assert actual == pytest.approx(expected, abs=1e-9)
    assert max(out['ledgers']['scaled_errors'].values()) < 2e-5
    assert not out['adopted'] and not out['field_scored'] and not out['closure_inputs_validated']


def test_dissipation_exchanges_equal_and_opposite_energy(transport):
    a = make_model(transport, dissipation_time=.2)
    b = make_model(transport, dissipation_time=.4)
    ao, bo = [v.assemble(order=4, angular_order=4, split_angles=False, solve=False) for v in (a, b)]
    np.testing.assert_array_equal(ao['matrix'], bo['matrix'])
    np.testing.assert_array_equal(ao['source'], bo['source'])
    assert ao['dissipation'] == pytest.approx(2*bo['dissipation'])
    delta = ao['right']-bo['right']
    assert delta[5:5+a.size] == pytest.approx(np.zeros(a.size), abs=1e-12)
    assert delta[5+a.size:5+2*a.size] == pytest.approx(-delta[6+2*a.size:], abs=1e-8)
    assert delta[5+2*a.size] == pytest.approx(-bo['dissipation'], abs=1e-8)


def test_ambient_tke_enters_source_and_q_boundary_once(transport):
    a, b = [make_model(transport, ambient_tke=k) for k in (0., .2)]
    ao, bo = [v.assemble(order=4, angular_order=4, split_angles=False, solve=False) for v in (a, b)]
    assert bo['source'][4]-ao['source'][4] == pytest.approx(.2*ao['source'][0], abs=1e-8)
    np.testing.assert_array_equal(ao['advective_jacobian'], bo['advective_jacobian'])
    np.testing.assert_array_equal(ao['matrix'][:5+2*a.size], bo['matrix'][:5+2*a.size])
    np.testing.assert_array_equal(ao['ledger_coefficients']['heat_boundary'], bo['ledger_coefficients']['heat_boundary'])
    assert bo['ledger_coefficients']['tke_boundary'] == pytest.approx(.2*bo['ledger_coefficients']['mass_boundary'])


def test_batch_invariance_and_independent_diagnostics(transport):
    model = make_model(transport, circulation_amplitudes=[1e-4, -1e-4, 1e-4, 0.])
    one = model.assemble(order=4, angular_order=8, split_angles=False, batch_size=1)
    many = model.assemble(order=4, angular_order=8, split_angles=False, batch_size=8)
    assert np.max(abs(one['matrix']-many['matrix'])/np.maximum(abs(many['matrix']), 1.)) < 1e-10
    assert np.max(abs(one['rates']-many['rates'])/np.maximum(abs(many['rates']), 1.)) < 1e-7
    diag = model.independent_diagnostics(many['rates'], angles=9, order=8)
    assert diag['minimum_tke'] > 0.
    assert diag['tke_edge_scaled_error'] > 0.
    assert diag['maximum_nonradial_stress_fraction'] < 1e-9
    assert not diag['physical_closure_passed'] and not diag['field_scored']


def test_manufactured_steady_q_recovers_immediate_heat_identity(finite):
    # Manufactured ledger, not an asserted equilibrium of a real profile.
    n = finite.count+1
    names = ('heat_axial', 'mean_kinetic_axial', 'tke_axial', 'production',
             'heat_boundary', 'mean_kinetic_boundary', 'tke_boundary', 'mass_boundary')
    coeff = {key: np.zeros(n) for key in names}
    for name, value in dict(heat_axial=9., heat_boundary=-2., production=7.,
                            mean_kinetic_axial=-7.).items():
        coeff[name][0] = value
    out = dict(ledger_coefficients=coeff, dissipation=7., work=0., source=np.zeros(5))
    result = finite.energy_ledgers(out, np.zeros(finite.count))
    assert result['residuals'] == dict(heat=0., tke=0., mean_kinetic=0.)
    assert result['total_budget'] == 0.


def test_normal_stress_ratio_enters_normal_stress_rows_and_ledger(transport):
    base = make_model(transport, normal_stress_ratio=0.)
    altered = make_model(transport, normal_stress_ratio=0.2)
    a = base.assemble(order=4, angular_order=4, split_angles=False, solve=False)
    b = altered.assemble(order=4, angular_order=4, split_angles=False, solve=False)
    assert np.all(a['ledger_coefficients']['normal_stress_boundary'] == 0.)
    assert np.any(np.abs(b['ledger_coefficients']['normal_stress_boundary']) > 0.)
    assert not np.allclose(a['advective_jacobian'], b['advective_jacobian'])
    assert np.max(abs(a['ledger_coefficients']['mean_kinetic_boundary']
                     - b['ledger_coefficients']['mean_kinetic_boundary'])) > 0.


def test_normal_stress_ratio_field_is_validated(transport):
    model = make_model(transport, normal_stress_ratio=.5)
    model.normal_stress_ratio = lambda d: np.full_like(d['q'], 3.0)
    with pytest.raises(ValueError, match='normal stress ratio'):
        model.local(np.array([.1]), np.array([.2]))


def test_source_terms_flags_forward_to_independent_source(transport):
    model = make_model(transport)
    base = model.assemble(order=4, angular_order=4, split_angles=False,
                          solve=False)
    expanded = model.assemble(order=4, angular_order=4, split_angles=False, solve=False,
                              include_work=True, include_thermal_relaxation=True,
                              include_phase_transition=True, include_turbulent_heat_exchange=True)
    assert expanded['source'][4] != pytest.approx(base['source'][4])
    assert expanded['source'][0] == pytest.approx(base['source'][0])
    assert expanded['thermal_source'] == ('dissipation_only+buoyancy_work+thermal_relaxation'
                                         '+phase_transition+turbulent_heat_exchange')
