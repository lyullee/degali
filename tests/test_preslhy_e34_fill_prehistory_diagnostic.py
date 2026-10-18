from pathlib import Path
import runpy

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
BUILD_TOOL = runpy.run_path(
    ROOT / "tools/build_preslhy_e34_fill_prehistory_diagnostic.py"
)
RUN_TOOL = runpy.run_path(
    ROOT / "tools/run_preslhy_e34_fill_prehistory_diagnostic.py"
)


def test_fixed_primary_frame_block_history_uses_no_concentration_input():
    time = np.arange(0.0, 20.1, 0.1)
    u = np.full_like(time, 1.0)
    v = np.full_like(time, 0.2)
    w = np.full_like(time, 0.1)

    history, diagnostics = BUILD_TOOL["_fixed_frame_block_history"](
        time,
        u,
        v,
        w,
        start_s=0.0,
        end_s=20.0,
        stream_axis_uv=(1.0, 0.0),
        cross_axis_uv=(0.0, 1.0),
        block_duration_s=10.0,
    )

    assert history.time_s == (0.0, 10.0, 20.0)
    assert history.u_m_s == (1.0, 1.0, 1.0)
    assert history.v_m_s == pytest.approx((0.2, 0.2, 0.2))
    assert history.w_m_s == pytest.approx((0.1, 0.1, 0.1))
    assert diagnostics["block_streamwise_nonpositive_count"] == 0


def test_diagnostic_numerical_gates_fail_closed():
    diagnostics = {
        "maximum_relative_mass_residual": 1.0e-7,
        "maximum_relative_energy_residual": 1.0e-7,
        "maximum_variable_density_continuity_residual_kg_m3_s": 1.0e-8,
        "maximum_lattice_mach_number": 0.10,
        "maximum_allowed_lattice_mach_number": 0.15,
        "active_scalar_applicability_passed": True,
        "thermodynamic_corrector_nonconverged_steps": 0,
    }
    gates = RUN_TOOL["_numerical_gates"](diagnostics)
    assert all(gates.values())

    diagnostics["maximum_relative_mass_residual"] = 1.0e-4
    gates = RUN_TOOL["_numerical_gates"](diagnostics)
    assert gates["relative_hydrogen_mass_residual"] is False
    assert not all(gates.values())
