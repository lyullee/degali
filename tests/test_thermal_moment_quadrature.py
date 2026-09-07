"""Adaptive response convergence, independent faces, and failure reporting."""

import math

import numpy as np
import pytest

from test_enthalpy_profile import phase, candidate
from test_buoyancy_profile import section
from test_thermal_moments import STATE
from degali.addons.thermal_moments import EnthalpyMomentOperators
from degali.addons.thermal_moment_quadrature import adaptive_unit_diffusion_response


@pytest.mark.parametrize("tolerance", [1e-6, 1e-7])
def test_adaptive_response_matches_independent_smooth_face_integrals(section, tolerance):
    op = EnthalpyMomentOperators(section)
    op.radial_diffusion_coefficient = lambda state, q: -3.*np.exp(-.8*q)
    section.jetplume._split = lambda a, *_: (2.*math.sqrt(a), .5*math.sqrt(a))
    sy, sn = section.section_widths(STATE)
    def fields(y, n):
        q = .5*((y/sy)**2+(n/sn)**2)
        coefficient = op.radial_diffusion_coefficient(STATE, q)
        return 0., coefficient*y/sy**2, coefficient*n/sn**2, 0.
    independent = op.geometry(STATE).budget(fields, points=64)
    result = adaptive_unit_diffusion_response(op, STATE, tolerance=tolerance)
    assert result.converged
    assert max(result.estimated_scaled_errors) <= tolerance
    assert result.response.transverse_volume == pytest.approx(independent.transverse_volume, rel=1e-7)
    assert result.response.outward_boundary_flux == pytest.approx(independent.outward_boundary_flux, rel=1e-7)
    assert result.response.derivative == pytest.approx(independent.derivative, rel=1e-7)


def test_subdivision_limit_is_reported_not_silently_accepted(section):
    op = EnthalpyMomentOperators(section)
    op.radial_diffusion_coefficient = lambda state, q: np.sin(777.*q)
    result = adaptive_unit_diffusion_response(op, STATE, max_subdivisions=2)
    assert not result.converged
    assert any(status != "converged" for status in result.component_status)


@pytest.mark.parametrize("tolerance", [0., -1., math.nan, math.inf, 1.])
def test_invalid_accuracy_request_rejected(section, tolerance):
    with pytest.raises(ValueError, match="tolerance"):
        adaptive_unit_diffusion_response(EnthalpyMomentOperators(section), STATE, tolerance=tolerance)
