"""Exact quadrature symmetry must not change the integration rule."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.energy_crosswind import IndependentEnergyCrosswind


@pytest.mark.parametrize("points", [8, 16, 32, 64, 128, 256, 512])
def test_symmetric_quadrature_equals_original_tensor_rule(points):
    model = IndependentEnergyCrosswind.__new__(IndependentEnergyCrosswind)
    model.k = SimpleNamespace(delta=2.15)
    model._quadrature_cache = {}
    exponent, weight = model._quadrature(points)
    nodes, weights = np.polynomial.legendre.leggauss(points)
    a, b = np.meshgrid(model.k.delta * nodes, model.k.delta * nodes, indexing="ij")
    wa, wb = np.meshgrid(weights, weights, indexing="ij")
    original_exponent = math.pi / 8.0 * (a*a + b*b)
    original_weight = math.pi / 4.0 * model.k.delta**2 * wa * wb
    assert len(exponent) <= points**2 // 4
    assert np.array_equal(exponent, np.unique(original_exponent))
    for function in (
        lambda x: np.ones_like(x),
        lambda x: np.exp(-x),
        lambda x: np.exp(-1.16**2*x),
        lambda x: np.maximum(0.4 - np.exp(-x), 0.0),
        lambda x: np.exp(-x) / (1.2 - 0.5*np.exp(-x)),
    ):
        old = np.sum(original_weight * function(original_exponent))
        new = np.sum(weight * function(exponent))
        assert new == pytest.approx(old, rel=2e-14, abs=1e-14)
