from __future__ import annotations

import pytest

from tools.run_ffi_control_transient_smoke import build_config


def test_ffi_control_smoke_uses_physical_axis_mapping_and_keeps_validation_withheld():
    config = build_config(duration_s=0.01, cell_size_m=0.025)
    assert config.domain_m == (5.8, 0.9, 0.8)
    assert config.cells == (232, 36, 32)
    assert config.obstacles == ()
    assert config.wind_speed_m_s == 0.0
    assert config.nonstreamwise_wall_boundary == "all_no_slip_tunnel"
    assert config.active_momentum_solver == "finite_volume_projection"
    assert config.binary_mixture_closure.species_name == "hydrogen"
    assert config.source_schedule.source_id == "FFI-control-A-published-LES-source-rate-no-fit:low-mach-handoff-M0.3"
    assert config.source_schedule.rate_kg_s[0] == 0.086e-3
    assert config.source_momentum_velocity_m_s[2] < 0.0
    assert config.source_position_m[2] == pytest.approx(0.342, abs=0.001)
    assert config.source_sigma_m == pytest.approx(0.03848, rel=1.0e-3)
    assert 4.0 * config.source_sigma_m / 0.025 > 4.0
    assert config.scalar_time_step_s == pytest.approx(0.005)
