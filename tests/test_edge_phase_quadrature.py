"""Independent square integration with rebuilt non-Gaussian phase rays."""

import numpy as np
import pytest
from test_edge_enrichment import phase, candidate, section, mixing, reservoir, projection
from degali.addons.edge_phase_quadrature import PhaseSplitSquareMoments


def test_zero_enrichment_reproduces_phase_split_radial_target(projection):
    quad = PhaseSplitSquareMoments(projection)
    actual = quad.moments(np.zeros(projection.count), order=16, angular_order=128)
    assert max(abs(actual-projection.target)/projection.scales) < 1e-8


def test_nonzero_phase_split_refinement_and_tensor_values(projection):
    quad = PhaseSplitSquareMoments(projection)
    parameters = np.random.default_rng(34).normal(0., 2e-5, projection.count)
    coarse = quad.moments(parameters, order=16, angular_order=64)
    fine = quad.moments(parameters, order=16, angular_order=128)
    tensor = projection.moments(parameters, 192)
    assert max(abs(coarse-fine)/projection.scales) < 1e-7
    assert max(abs(fine-tensor)/projection.scales) < 2e-5
    assert max(abs(fine-projection.target)/projection.scales) > 1e-8


def test_each_new_profile_rebuilds_its_phase_partitions(projection):
    quad = PhaseSplitSquareMoments(projection)
    angles = np.array([.12, .58])
    zero = np.zeros(projection.count)
    base = quad.partitions(zero, angles)
    altered = zero.copy()
    altered[2] = 1e-4
    changed = quad.partitions(altered, angles)
    assert any(not np.array_equal(a, b) for a, b in zip(base, changed))
    for angle, knots in zip(angles, changed):
        assert knots[0] == 0.
        assert knots[-1] == pytest.approx(projection.q0/np.cos(angle)**2)
        assert np.all(np.diff(knots) > 0.)


def test_unsupported_nonmonotone_ray_is_rejected(projection, monkeypatch):
    quad = PhaseSplitSquareMoments(projection)
    monkeypatch.setattr(quad, "coordinates", lambda parameters, q, phi: (300.-q, .1+.01*q))
    with pytest.raises(ValueError, match="monotonic"):
        quad.partitions(np.zeros(projection.count), np.array([.2]))
