"""Square symmetry, phase consistency and exact-moment retraction tests."""

import numpy as np
import pytest
from scipy.special import roots_legendre

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from test_transverse_mixing import mixing
from test_reservoir_thermal import reservoir
from degali.addons.edge_enrichment import SquareEvenBasis, FixedTransportEdgeProjection


@pytest.fixture
def projection(reservoir):
    out = reservoir.evaluate()
    return FixedTransportEdgeProjection(reservoir, out, degree=4)


@pytest.mark.parametrize("degree,count", [(4, 8), (6, 15)])
def test_basis_is_square_symmetric_center_zero_and_orthonormal(degree, count):
    basis = SquareEvenBasis(degree)
    assert basis.size == count
    v, da, db = basis.values(np.array([0.]), np.array([0.]))
    assert v == pytest.approx(np.zeros_like(v), abs=1e-14)
    assert da == pytest.approx(np.zeros_like(da), abs=1e-14)
    assert db == pytest.approx(np.zeros_like(db), abs=1e-14)
    x, w = roots_legendre(16)
    a, b = np.meshgrid(.5*(x+1), .5*(x+1), indexing="ij")
    p = basis.values(a.ravel(), b.ravel())[0]
    weight = .25*np.outer(w, w).ravel()
    assert p.T@(weight[:, None]*p) == pytest.approx(np.eye(count), abs=1e-12)
    assert basis.values(-b.ravel(), a.ravel())[0] == pytest.approx(p, abs=1e-12)


def test_basis_normal_derivatives_match_coordinate_differences():
    basis = SquareEvenBasis(6)
    a, b, step = np.array([.13, .48, .89]), np.array([.25, .51, .77]), 1e-6
    v, da, db = basis.values(a, b)
    assert da == pytest.approx((basis.values(a+step, b)[0]-basis.values(a-step, b)[0])/(2*step), rel=1e-7, abs=1e-8)
    assert db == pytest.approx((basis.values(a, b+step)[0]-basis.values(a, b-step)[0])/(2*step), rel=1e-7, abs=1e-8)


def test_zero_correction_reproduces_actual_phase_fields(projection):
    a, b = np.array([0., .13, .48, .89, 1.]), np.array([0., .25, .51, .77, 1.])
    p = projection.prepare(a, b)
    d = projection.fields(np.zeros(projection.count), p)
    old = projection.mixing.local(p["q"])
    for key in ("rho", "y", "c", "h"):
        assert d[key] == pytest.approx(old[key], rel=1e-9, abs=1e-9)


def test_zero_edge_correction_reproduces_frozen_flux_defects(projection):
    p = projection.edge(257)
    old = projection.reservoir.edge_fields(p["q"], projection.baseline["family"])
    at = np.array([1., projection.baseline["gamma"]])
    expected = np.r_[old["species_gradient_defect"]@at/p["scale_c"], old["enthalpy_gradient_defect"]@at/p["scale_h"]]
    assert projection.edge_residual(np.zeros(projection.count), 257) == pytest.approx(expected, abs=1e-9)


def test_analytic_six_moment_jacobian_matches_independent_values(projection):
    rng = np.random.default_rng(20260905)
    parameters = rng.normal(0., 2e-5, projection.count)
    values, analytic = projection.moments(parameters, 48, jacobian=True)
    step = 1e-6
    numerical = np.column_stack([(projection.moments(parameters+step*e, 48)-projection.moments(parameters-step*e, 48))/(2*step) for e in np.eye(projection.count)])
    assert np.max(abs(analytic-numerical)/np.maximum(abs(analytic), 1.)) < 2e-5


def test_retraction_preserves_all_six_moments_and_center(projection):
    parameters = np.random.default_rng(17).normal(0., 1e-4, projection.count)
    result = projection.retract(parameters, 48)
    values = projection.moments(result, 48)
    assert np.max(abs(values-projection.target)/projection.scales) <= 1e-9
    center = projection.fields(result, projection.prepare(np.array([0.]), np.array([0.])))
    assert center["c"][0] == pytest.approx(projection.mixing.state[0]*projection.mixing.state[1], rel=1e-12)
    assert center["h"][0] == pytest.approx(projection.hc, rel=1e-12)


def test_correction_keeps_frozen_transport_and_signs(projection):
    original = {k: projection.edge(257)[k].copy() for k in ("fm", "chi", "target_h")}
    parameters = np.random.default_rng(21).normal(0., 1e-4, projection.count)
    fields = projection.fields(parameters, projection.grid(48))
    assert np.all(fields["rho"] > 0.) and np.all((fields["y"] > 0.) & (fields["y"] < 1.))
    assert np.all(fields["h"]*projection.hc > 0.)
    for k, v in original.items():
        assert np.array_equal(v, projection.edge(257)[k])
    assert not hasattr(projection, "solve")
    assert not hasattr(projection, "derivatives")


def test_excessive_or_invalid_correction_is_not_clipped(projection):
    with pytest.raises(ValueError, match="trust region"):
        projection.fields(np.ones(projection.count), projection.grid(16))
    with pytest.raises(ValueError):
        projection.fields(np.zeros(projection.count+1), projection.grid(16))
    with pytest.raises(ValueError):
        FixedTransportEdgeProjection(projection.reservoir, {"valid": False}, degree=4)


@pytest.mark.parametrize("degree", [True, 1, 2.5])
def test_invalid_basis_degree_is_rejected(degree):
    with pytest.raises(ValueError):
        SquareEvenBasis(degree)
