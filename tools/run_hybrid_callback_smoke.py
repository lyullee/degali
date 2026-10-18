"""Run a small concurrent fast/resolved callback-coupling smoke.

This is deliberately a hook/ledger test on two overlapping toy domains.  It
does not claim a facility prediction or an accuracy score.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.hybrid_flux_packets import CartesianInterfacePlane, boundary_cell_mask
from degali.addons.hybrid_solver_bridge import SynchronizedSolverBridge
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


def run(output: Path) -> dict[str, object]:
    domain = (4.0, 4.0, 4.0)
    cells = (8, 8, 8)
    duration = 0.002
    dt = 0.001
    obstacles = (AxisAlignedCuboid(1.5, 2.0, 1.5, 2.5, 0.0, 1.0, "smoke block"),)
    zero_source = SourceRateSchedule((0.0, duration), (0.0, 0.0))
    fast_config = FastObstacleTransport3DConfig(
        domain_m=domain, cells=cells, duration_s=duration, time_step_s=dt,
        wind_speed_m_s=2.0, diffusivity_m2_s=0.03,
        source_schedule=zero_source, source_position_m=(0.5, 2.0, 0.4),
        source_sigma_m=0.15, obstacles=obstacles,
    )
    resolved_config = ResolvedObstacleFlow3DConfig(
        domain_m=domain, cells=cells, wind_speed_m_s=2.0,
        kinematic_viscosity_m2_s=0.5, flow_steps=12, minimum_flow_steps=12,
        lattice_inlet_speed=0.03, duration_s=duration,
        scalar_time_step_s=dt, scalar_diffusivity_m2_s=0.03,
        source_schedule=zero_source, source_position_m=(0.5, 2.0, 0.4),
        source_sigma_m=0.15, obstacles=obstacles, smagorinsky_constant=0.10,
    )
    bridge = SynchronizedSolverBridge(time_step_s=dt, timeout_s=60.0)
    bridge.register_side(
        name="outer", plane=CartesianInterfacePlane(0, -1, 1.0, 0.25),
        target_mask=boundary_cell_mask((8, 8, 8), axis=2, side="min"),
        cell_volume_m3=(domain[0] / cells[0]) * (domain[1] / cells[1]) * (domain[2] / cells[2]),
    )
    bridge.register_side(
        name="patch", plane=CartesianInterfacePlane(0, 0, -1.0, 0.25),
        target_mask=boundary_cell_mask((8, 8, 8), axis=2, side="max"),
        cell_volume_m3=(domain[0] / cells[0]) * (domain[1] / cells[1]) * (domain[2] / cells[2]),
    )
    outer = bridge.callbacks_for("outer")
    patch = bridge.callbacks_for("patch")
    started = time.perf_counter()

    def run_fast():
        return solve_fast_obstacle_transport_3d(
            fast_config,
            face_flux_callback=outer["face_flux_callback"],
            interface_state_callback=outer["interface_state_callback"],
            interface_state_mass_callback=outer["interface_state_mass_callback"],
            interface_momentum_callback=outer["interface_momentum_callback"],
        )

    def run_resolved():
        return solve_resolved_obstacle_flow_3d(
            resolved_config,
            face_flux_callback=patch["face_flux_callback"],
            interface_conserved_state_callback=patch["interface_conserved_state_callback"],
            interface_conserved_mass_callback=patch["interface_conserved_mass_callback"],
            interface_momentum_callback=patch["interface_momentum_callback"],
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        fast_future = executor.submit(run_fast)
        resolved_future = executor.submit(run_resolved)
        fast_error = resolved_error = None
        try:
            fast_result = fast_future.result()
        except Exception as exc:  # pragma: no cover - diagnostic CLI path
            fast_error = repr(exc)
            fast_result = None
        try:
            resolved_result = resolved_future.result()
        except Exception as exc:  # pragma: no cover - diagnostic CLI path
            resolved_error = repr(exc)
            resolved_result = None
    if fast_error or resolved_error:
        raise RuntimeError(
            f"callback smoke failed: fast={fast_error}, resolved={resolved_error}"
        )
    bridge.finalize_step(1)
    record = {
        "schema": "degali.hybrid-callback-smoke.v1",
        "status": "conditional_callback_ledger_smoke",
        "quantitative_ready": False,
        "duration_s": duration,
        "time_step_s": dt,
        "shared_clock_steps": 2,
        "exchanged_steps": bridge.exchanged_steps,
        "wall_clock_s": time.perf_counter() - started,
        "fast_diagnostics": {
            "maximum_mass_residual_kg": fast_result.diagnostics["maximum_mass_residual_kg"],
            "interface_momentum_callback_applied": fast_result.diagnostics[
                "interface_momentum_callback_applied"
            ],
            "interface_state_mass_callback_applied": fast_result.diagnostics[
                "interface_state_mass_callback_applied"
            ],
        },
        "resolved_diagnostics": {
            "maximum_relative_mass_residual": resolved_result.diagnostics[
                "maximum_relative_mass_residual"
            ],
            "interface_conserved_state_callback_applied": resolved_result.diagnostics[
                "interface_conserved_state_callback_applied"
            ],
            "interface_momentum_callback_applied": resolved_result.diagnostics[
                "interface_momentum_callback_applied"
            ],
            "interface_conserved_mass_callback_applied": resolved_result.diagnostics[
                "interface_conserved_mass_callback_applied"
            ],
        },
        "warnings": [
            "overlapping toy domains used only to verify callback ordering and ledger exchange",
            "fast outer path uses nominal density and normal momentum flux only",
            "no facility accuracy or full-event validation claim",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


if __name__ == "__main__":
    destination = Path("outputs/hybrid-callback-smoke-2026-10-14.json")
    print(json.dumps(run(destination), indent=2))
