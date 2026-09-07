"""Degree-six/37-state coverage, without modifying the frozen pilot tests."""

import numpy as np
import pytest
from test_enriched_transport import phase, candidate, section, mixing, reservoir
from degali.addons.edge_enrichment import FixedTransportEdgeProjection
from degali.addons.enriched_transport import PrescribedRadialMixing, EnrichedModalTransport
from degali.addons.coupled_shape_initialization import FixedMeshCoupledTransport


def test_degree6_vectorized_operator_matches_independent_operator(reservoir):
    baseline = reservoir.evaluate()
    p = FixedTransportEdgeProjection(reservoir, baseline, degree=6)
    par = np.zeros(p.count)
    supplied = PrescribedRadialMixing.from_weak_baseline(p, order=8)
    kwargs = dict(scalar_mixing=supplied, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work_immediate_shear_heat")
    model = EnrichedModalTransport(p, par, **kwargs)
    mesh = FixedMeshCoupledTransport(p, par, order=4, angular_order=4,
        face_samples=3, split_angles=True, **kwargs)
    actual, expected = mesh.evaluate(par), model.assemble(order=4, angular_order=4)
    assert model.count == 37
    assert len(model.active) == 35
    assert actual["matrix"].shape == (35, 37)
    for k in ("rates", "matrix", "right", "source"):
        assert np.max(abs(actual[k]-expected[k])/np.maximum(abs(expected[k]), 1.)) < 1e-8
    assert actual["weak_heat_scaled_error"] < 1e-8
    assert abs(actual["mass_boundary_error"]) < 1e-8
    assert np.linalg.matrix_rank(actual["moment_jacobian"]) == 6
