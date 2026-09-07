"""State movement, dimensional scaling and strict integration-domain guards."""

import numpy as np
import pytest
from test_enriched_transport import phase, candidate, section, mixing, reservoir, projection, transport
from degali.addons.enriched_transport import shifted_advective_moments, advective_moments
from degali.addons.enriched_segments import (encode_enriched, decode_enriched, shifted_field_view,
    scaled_mixing, forward_shape_horizon, assert_transport_domain, EnrichedShortSegment)


def test_enriched_state_roundtrip_and_field_move(transport):
    p, par = transport.projection, transport.parameters
    initial = encode_enriched(p, par)
    state, modes = decode_enriched(initial, p.count)
    assert state == pytest.approx(p.mixing.state)
    assert np.array_equal(modes, par)
    direction = np.random.default_rng(211).normal(0., .01, len(initial))
    ds = 1e-4
    view, moved = shifted_field_view(p, initial+ds*direction)
    actual = advective_moments(view, moved, order=4, angular_order=4)
    expected = shifted_advective_moments(p, par, direction, ds, order=4, angular_order=4)
    assert actual == pytest.approx(expected, rel=1e-12, abs=1e-9)
    assert np.array_equal(encode_enriched(p, par), initial)
    assert view._edges is not p._edges


@pytest.mark.parametrize("policy,expected",[("equilibrium_velocity_width",1.5),("constant_geometric_diffusivity",.25)])
def test_scaling_units_and_reference_mass_are_not_rewritten(transport, policy, expected):
    old = transport.projection.mixing.state
    state = old.copy()
    state[2] *= 4.
    state[4] *= 3.
    supplied, factor = scaled_mixing(transport.mixing, old, state, policy)
    q = np.array([.1, .4, 1.])
    before = transport.mixing.evaluate(q)
    assert factor == pytest.approx(expected)
    assert supplied.evaluate(q)[:, 0] == pytest.approx(before[:, 0]*expected)
    assert supplied.evaluate(q)[:, 1] == pytest.approx(before[:, 1])
    assert np.array_equal(before, transport.mixing.evaluate(q))


def test_forward_horizon_distinguishes_inward_from_outward_motion(transport):
    p = transport.projection
    par = np.zeros(p.count)
    par[0] = 1e-3
    rates = np.zeros(transport.count)
    rates[7] = 1.
    outward = forward_shape_horizon(p, par, rates)
    inward = forward_shape_horizon(p, par, -rates)
    assert outward > 0. and inward > outward
    assert forward_shape_horizon(p, par, np.zeros_like(rates)) == np.inf


def valid_guard_record():
    return dict(edge_defects=dict(hydrogen=.01, heat=.02, momentum=.03), minimum_chi_momentum=1.,
        maximum_outward_mass=-1., curvature_half_width=.001, linear_scaled_error=0.,
        weak_heat_scaled_error=0., mass_boundary_error=0., source=np.zeros(5))


@pytest.mark.parametrize("key,value",[("minimum_chi_momentum",-1.),("maximum_outward_mass",0.),
    ("curvature_half_width",.1),("mass_boundary_error",1e-4),("linear_scaled_error",np.nan)])
def test_domain_rejection_is_not_hidden(key, value):
    out = valid_guard_record()
    assert_transport_domain(out)
    out[key] = value
    with pytest.raises(ValueError):
        assert_transport_domain(out)


def test_midpoint_integrator_checks_endpoint_before_accepting_it():
    driver = object.__new__(EnrichedShortSegment)
    driver.initial, driver.count = np.array([0.]), 1
    def evaluate(value):
        if value[0] > .4:
            raise ValueError("test endpoint outside the domain")
        return dict(**valid_guard_record(), rates=np.array([1.]), maximum_log_shape=0., mixing_amplitude_ratio=1.)
    driver.evaluate = evaluate
    out = driver.integrate(1., 2)
    assert not out["completed"]
    assert out["reached_m"] == 0.
    assert out["parameters"] == pytest.approx([0.])
    assert out["failure_stage"] == "endpoint"


def test_midpoint_rates_and_sources_share_the_same_integration():
    driver = object.__new__(EnrichedShortSegment)
    driver.initial, driver.count = np.array([1.]), 1
    def evaluate(value):
        out = valid_guard_record()
        out["source"] = np.full(5, value[0])
        return dict(**out, rates=value.copy(), maximum_log_shape=0., mixing_amplitude_ratio=1.)
    driver.evaluate = evaluate
    out = driver.integrate(.1, 4)
    assert out["completed"]
    assert out["cumulative_sources"] == pytest.approx(np.full(5, out["parameters"][0]-1.), abs=1e-14)
    assert out["rhs_calls"] == 9


def test_invalid_encoded_state_and_policy_rejected(transport):
    initial = encode_enriched(transport.projection, transport.parameters)
    with pytest.raises(ValueError):
        decode_enriched(initial[:-1], transport.projection.count)
    initial[0] = 1000.
    with pytest.raises(ValueError):
        decode_enriched(initial, transport.projection.count)
    with pytest.raises(ValueError, match="explicit"):
        scaled_mixing(transport.mixing, transport.projection.mixing.state, transport.projection.mixing.state, "guess")
