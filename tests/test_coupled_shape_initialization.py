"""Accelerated initialization is checked against the separate moving-grid law."""

import numpy as np
import pytest
from test_enriched_transport import phase, candidate, section, mixing, reservoir, projection, transport
from degali.addons.coupled_shape_initialization import FixedMeshCoupledTransport, SimultaneousShapeInitializer


def mesh_for(transport, **kwargs):
    return FixedMeshCoupledTransport(transport.projection, transport.parameters,
        scalar_mixing=transport.mixing, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work_immediate_shear_heat", **kwargs)


def test_vectorized_operator_matches_independent_phase_split_operator(transport):
    mesh = mesh_for(transport, order=4, angular_order=4, face_samples=5, split_angles=True)
    got = mesh.evaluate(transport.parameters)
    expected = transport.assemble(order=4, angular_order=4)
    for key in ("rates", "matrix", "right", "source"):
        assert np.max(abs(got[key]-expected[key])/np.maximum(abs(expected[key]), 1.)) < 1e-9
    assert got["linear_scaled_error"] < 1e-8
    assert got["weak_heat_scaled_error"] < 1e-8
    assert abs(got["mass_boundary_error"]) < 1e-8


def test_diagnostic_face_rays_do_not_change_rates(transport):
    one = mesh_for(transport, order=4, angular_order=6, face_samples=3).evaluate(transport.parameters)
    two = mesh_for(transport, order=4, angular_order=6, face_samples=9).evaluate(transport.parameters)
    assert one["rates"] == pytest.approx(two["rates"], rel=1e-9, abs=1e-8)
    assert one["moments"] == pytest.approx(two["moments"], rel=1e-12, abs=1e-9)


def test_six_moment_derivatives_use_actual_current_phase(transport):
    mesh = mesh_for(transport, order=4, angular_order=6, face_samples=3)
    par = np.random.default_rng(41).normal(0., 1e-5, transport.projection.count)
    direction = np.random.default_rng(42).normal(0., .01, len(par))
    out = mesh.evaluate(par)
    ds = 1e-6
    plus, minus = mesh.evaluate(par+ds*direction), mesh.evaluate(par-ds*direction)
    numerical = (plus["moments"]-minus["moments"])/(2*ds)
    assert max(abs(numerical-out["moment_jacobian"]@direction)/transport.projection.scales) < 1e-6
    assert max(abs(plus["rates"]-minus["rates"])) > 1e-7


def test_objective_and_diagnostics_do_not_clip_negative_stress(transport):
    mesh = mesh_for(transport, order=4, angular_order=6, face_samples=3)
    out = mesh.evaluate(transport.parameters)
    assert np.all(np.isfinite(out["objective"]))
    assert len(out["objective"]) > 3*3
    assert out["edge_residual"].shape == (9, 3)
    assert set(out["edge_defects"]) == {"hydrogen", "heat", "momentum"}
    assert out["minimum_chi_momentum"] < 0.  # Fixture's rejection remains visible.
    assert np.min(out["objective"][9:]) < 0.


def test_coupled_retraction_enforces_all_six_mesh_moments(transport):
    mesh = mesh_for(transport, order=4, angular_order=8, face_samples=3)
    fit = SimultaneousShapeInitializer(mesh)
    initial = np.random.default_rng(12).normal(0., 1e-5, transport.projection.count)
    parameters, out = fit.retract(initial)
    assert max(abs(out["moments"]-transport.projection.target)/transport.projection.scales) <= 1e-9
    assert max(abs(fit.constraint.bounds_matrix@parameters)) <= .0999+1e-9


def test_invalid_shapes_and_unselected_physics_are_rejected(transport):
    mesh = mesh_for(transport, order=4, angular_order=4, face_samples=3)
    with pytest.raises(ValueError, match="trust region"):
        mesh.evaluate(np.ones(transport.projection.count))
    with pytest.raises(ValueError, match="assumptions"):
        FixedMeshCoupledTransport(transport.projection, transport.parameters,
            scalar_mixing=transport.mixing, thermal_species_ratio=1., mechanical_work="unspecified")
