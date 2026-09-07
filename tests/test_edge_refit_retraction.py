"""Independent end-to-end conservation and anisotropic moment checks."""

import numpy as np
from test_edge_enrichment import phase, candidate, section, mixing, reservoir, projection
from degali.addons.edge_conservative_refit import FaceSplitSquareMoments, ConservativeEdgeRefit


def test_actual_nonlinear_retraction_reduces_conservation_error(projection):
    refit = ConservativeEdgeRefit(projection)
    initial = np.random.default_rng(94).normal(0., 3e-5, projection.count)
    result, history = refit.retract(initial)
    independent = refit.quadrature.moments(result, order=16, angular_order=16, probes=2049)
    assert max(abs(independent-projection.target)/projection.scales) < 1e-8
    assert history[-1]["maximum_error"] < history[0]["maximum_error"]
    assert np.max(abs(refit.bounds_matrix@result)) <= .1


def test_anisotropic_enthalpy_second_moment_matches_cartesian_integral(projection):
    # Keep area fixed while making the two geometric widths different.
    projection.mixing.sy *= 1.7
    projection.mixing.sn /= 1.7
    parameters = np.random.default_rng(51).normal(0., 3e-5, projection.count)
    polar = FaceSplitSquareMoments(projection).moments(parameters)
    cartesian = projection.moments(parameters, 192)
    assert abs(polar[5]-cartesian[5])/max(abs(polar[5]), 1.) < 1e-9
