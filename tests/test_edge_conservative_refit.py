"""Face event splitting and actual-moment constrained step checks."""

import numpy as np
import pytest
from test_edge_enrichment import phase, candidate, section, mixing, reservoir, projection
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments, ConservativeEdgeRefit


def test_face_split_zero_profile_reproduces_reference(projection):
    q = FaceSplitSquareMoments(projection)
    z = np.zeros(projection.count)
    low = q.moments(z, order=4, angular_order=4)
    high = q.moments(z, order=8, angular_order=8, probes=2049)
    assert max(abs(high-projection.target)/projection.scales) < 1e-8
    assert max(abs(high-low)/projection.scales) < 1e-8


def test_face_events_replay_table_knots_and_probe_refinement(projection):
    q = FaceSplitSquareMoments(projection)
    z = np.random.default_rng(12).normal(0., 1e-5, projection.count)
    knots = q.angular_knots(z)
    refined = q.angular_knots(z, probes=2049)
    assert knots == pytest.approx(refined, abs=1e-10)
    ti, y = q.face_coordinates(z, knots[1:-1])
    inv = projection.section.phase_inverse
    error = np.minimum(np.min(abs(ti[:, None]-inv.t_grid), axis=1),
                       np.min(abs(y[:, None]-inv.y_grid), axis=1))
    assert max(error) < 1e-8


def test_integrated_shape_direction_matches_independent_difference(projection):
    q = FaceSplitSquareMoments(projection)
    z = np.random.default_rng(31).normal(0., 2e-5, projection.count)
    direction = np.random.default_rng(9).normal(0., .1, projection.count)
    actual, jac = q.moments(z, order=4, angular_order=4, jacobian=True)
    h = 1e-5
    fd = (q.moments(z+h*direction, order=4, angular_order=4)-q.moments(z-h*direction, order=4, angular_order=4))/(2*h)
    assert max(abs(fd-jac@direction)/projection.scales) < 2e-5


def test_step_respects_linear_moments_and_interior_shape_limit(projection):
    refit = ConservativeEdgeRefit(projection)
    rng = np.random.default_rng(7)
    j = rng.normal(size=(6, projection.count))
    residual = rng.normal(0., 1e-4, 6)
    edge_j = np.eye(projection.count)
    edge = rng.normal(0., .1, projection.count)
    par = np.zeros(projection.count)
    dx, info = refit.step(par, residual, j, edge, edge_j)
    assert max(abs(j@dx+residual)) <= 1e-8
    assert max(abs(refit.bounds_matrix@dx)) <= .1
    assert max(abs(dx)) <= .020000001


def test_infeasible_constrained_step_is_rejected(projection):
    refit = ConservativeEdgeRefit(projection)
    with pytest.raises(ValueError, match="infeasible"):
        refit.step(np.zeros(projection.count), np.ones(6), np.zeros((6, projection.count)),
                   np.ones(projection.count), np.eye(projection.count))
