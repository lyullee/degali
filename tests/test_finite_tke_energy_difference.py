import numpy as np
import pytest
pytest.importorskip('mpmath')
from test_enriched_transport import phase, candidate, section, mixing, reservoir, projection, transport
from test_paired_exact_kinetic_difference import exact_transport
from test_finite_tke_transport import make_model
from degali.addons.finite_tke_energy_difference import actual_tke_moment_difference
from degali.addons.phase_radial_quadrature import gauss_rule


def test_actual_q_zero_direction_is_exact(exact_transport):
    model = make_model(exact_transport)
    result = actual_tke_moment_difference(model, np.zeros(model.count), 1e-6, order=32)
    np.testing.assert_array_equal(result['derivatives'], np.zeros(model.tke_count))


def test_actual_q_all_moments_match_independent_jacobian_and_refine(exact_transport):
    model = make_model(exact_transport)
    model.parameters[2] = 1e-4
    rates = np.random.default_rng(706).normal(0., .02, model.count)
    x, w = gauss_rule(48)
    a, b = np.meshgrid(x, x, indexing='ij')
    local = model.local(a.ravel(), b.ravel())
    weight = (8*model.base.projection.q0*np.outer(w, w)).ravel()
    analytic = (weight*(local['bq'] @ rates)) @ local['theta_q']
    outputs = [actual_tke_moment_difference(model, rates, step, order=order)['derivatives']
               for step, order in ((1e-5, 32), (5e-6, 64))]
    for result in outputs:
        assert result == pytest.approx(analytic, rel=2e-10, abs=2e-10)


def test_actual_q_rejects_bad_direction_or_step(exact_transport):
    model = make_model(exact_transport)
    with pytest.raises(ValueError, match='matching finite'):
        actual_tke_moment_difference(model, np.zeros(model.count-1), 1e-6)
    with pytest.raises(ValueError, match='positive'):
        actual_tke_moment_difference(model, np.zeros(model.count), 0.)
