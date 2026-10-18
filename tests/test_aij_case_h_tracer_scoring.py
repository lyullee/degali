from __future__ import annotations

import numpy as np
import pytest

from tools.score_aij_case_h_resolved_tracer import (
    _source_exit_velocity_m_s,
    _source_sigma_for_run,
    _ethylene_ppmv_from_mass_concentration,
    _sample_fluid,
)


def test_ethylene_mass_to_idealgas_ppmv_is_zero_at_zero_and_finite():
    assert _ethylene_ppmv_from_mass_concentration(
        0.0, temperature_k=293.0, pressure_pa=101_325.0,
    ) == 0.0
    value = _ethylene_ppmv_from_mass_concentration(
        1.0e-4, temperature_k=293.0, pressure_pa=101_325.0,
    )
    assert 0.0 < value < 1_000.0
    with pytest.raises(ValueError, match="nonnegative"):
        _ethylene_ppmv_from_mass_concentration(
            -1.0, temperature_k=293.0, pressure_pa=101_325.0,
        )


def test_reported_case_h_flow_and_four_mm_port_imply_exit_velocity():
    assert _source_exit_velocity_m_s(5.83e-6) == pytest.approx(0.4639366591, rel=1.0e-10)


def test_source_sigma_can_be_varied_independently_of_grid_spacing():
    assert _source_sigma_for_run(0.025, None) == pytest.approx(0.025)
    assert _source_sigma_for_run(0.025, 0.001) == pytest.approx(0.001)
    with pytest.raises(ValueError, match="positive and finite"):
        _source_sigma_for_run(0.025, 0.0)
    with pytest.raises(ValueError, match="positive and finite"):
        _source_sigma_for_run(float("nan"), None)


def test_trilinear_sampling_interpolates_and_fails_closed_on_solid_stencil():
    x = y = z = np.asarray([0.0, 1.0])
    field = np.arange(8.0).reshape((2, 2, 2))
    solid = np.zeros_like(field, dtype=bool)
    assert _sample_fluid(x, y, z, field, solid, (0.5, 0.5, 0.5)) == pytest.approx(3.5)
    solid[1, 1, 1] = True
    assert _sample_fluid(x, y, z, field, solid, (0.5, 0.5, 0.5)) is None
    assert _sample_fluid(x, y, z, field, solid, (2.0, 0.5, 0.5)) is None
