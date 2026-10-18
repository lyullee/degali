"""Isolate carrier-flow time-step sensitivity for the active audit case."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.fast_obstacle_transport_3d import (  # noqa: E402
    FastObstacleReceptor3D,
    _receptor_stencil,
)
from degali.addons.resolved_obstacle_flow_3d import _solve_flow  # noqa: E402
from audit_active_thermal_corrector_time import (  # noqa: E402
    LATTICE_SPEEDS,
    OUTPUT,
    _config,
)


RECEPTORS = (
    FastObstacleReceptor3D("windward", (1.20, 1.20, 0.50)),
    FastObstacleReceptor3D("roof", (1.70, 1.20, 1.20)),
    FastObstacleReceptor3D("side", (1.70, 0.60, 0.50)),
    FastObstacleReceptor3D("lee", (2.40, 1.20, 0.50)),
)


def _sample(field: np.ndarray, stencil) -> float:
    return float(sum(weight * field[tuple(index)] for *index, weight in stencil))


def _relative_l2(left: np.ndarray, right: np.ndarray, fluid: np.ndarray) -> float:
    left_values = left[fluid]
    right_values = right[fluid]
    return float(np.linalg.norm(left_values - right_values) / max(
        np.linalg.norm(left_values), np.linalg.norm(right_values), 1.0e-15,
    ))


def main() -> int:
    fields = []
    cases = []
    for lattice_speed in LATTICE_SPEEDS:
        config = _config(lattice_speed)
        started = time.perf_counter()
        flow = _solve_flow(config)
        runtime = time.perf_counter() - started
        (
            x, y, z, solid, u, v, w, _u_face, _v_face, _w_face,
            _pressure, _pressure_coefficient, _eddy_viscosity,
            residual_history, _density, flow_dt, base_tau, maximum_tau,
            actual_steps, _raw_divergence, _divergence_before,
            divergence_after, _projection_iterations, _distribution,
            statistical_residual, _average_count,
        ) = flow
        receptor_values = {}
        for receptor in RECEPTORS:
            stencil = _receptor_stencil(x, y, z, solid, receptor)
            vector = np.asarray([
                _sample(component, stencil) for component in (u, v, w)
            ])
            receptor_values[receptor.label] = {
                "u_m_s": float(vector[0]),
                "v_m_s": float(vector[1]),
                "w_m_s": float(vector[2]),
                "speed_m_s": float(np.linalg.norm(vector)),
            }
        cases.append({
            "lattice_inlet_speed": lattice_speed,
            "physical_time_step_s": flow_dt,
            "flow_steps": actual_steps,
            "physical_flow_horizon_s": actual_steps * flow_dt,
            "base_relaxation_time": base_tau,
            "maximum_relaxation_time": maximum_tau,
            "final_instantaneous_residual": float(residual_history[-1, 1]),
            "statistical_residual": statistical_residual,
            "maximum_projected_divergence_s_1": divergence_after,
            "runtime_s": runtime,
            "receptors": receptor_values,
        })
        fields.append((u, v, w, solid))

    fluid = ~fields[-1][3]
    metrics = {}
    for name, index in (("u", 0), ("v", 1), ("w", 2)):
        coarse_medium = _relative_l2(
            fields[0][index], fields[1][index], fluid,
        )
        medium_fine = _relative_l2(
            fields[1][index], fields[2][index], fluid,
        )
        metrics[name] = {
            "coarse_medium_relative_l2": coarse_medium,
            "medium_fine_relative_l2": medium_fine,
            "observed_order": (
                None if min(coarse_medium, medium_fine) <= 0.0
                else math.log(coarse_medium / medium_fine) / math.log(2.0)
            ),
            "successive_difference_decreased": medium_fine < coarse_medium,
        }
    speeds = [np.sqrt(u*u + v*v + w*w) for u, v, w, _ in fields]
    coarse_medium = _relative_l2(speeds[0], speeds[1], fluid)
    medium_fine = _relative_l2(speeds[1], speeds[2], fluid)
    metrics["speed"] = {
        "coarse_medium_relative_l2": coarse_medium,
        "medium_fine_relative_l2": medium_fine,
        "observed_order": (
            None if min(coarse_medium, medium_fine) <= 0.0
            else math.log(coarse_medium / medium_fine) / math.log(2.0)
        ),
        "successive_difference_decreased": medium_fine < coarse_medium,
    }
    gate_reasons = [
        f"{name} field difference is non-monotonic"
        for name, metric in metrics.items()
        if not metric["successive_difference_decreased"]
    ]
    record = {
        "schema": "degali.active-time-carrier-flow-diagnostic.v1",
        "status": (
            "accepted_component_time" if not gate_reasons
            else "withheld_component_time"
        ),
        "scope": "carrier flow only on a fixed grid; no scalar or buoyancy",
        "cases_coarse_to_fine": cases,
        "field_metrics": metrics,
        "gate_reasons": gate_reasons,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "active_time_carrier_flow_diagnostic.json"
    figure_path = OUTPUT / "active_time_carrier_flow_diagnostic.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    time_steps = [case["physical_time_step_s"] for case in cases]
    for receptor in RECEPTORS:
        axes[0].plot(
            time_steps,
            [case["receptors"][receptor.label]["speed_m_s"] for case in cases],
            "o-", label=receptor.label,
        )
    axes[0].invert_xaxis()
    axes[0].set_xlabel("physical time step (s), finer →")
    axes[0].set_ylabel("carrier speed (m/s)")
    axes[0].set_title("Fixed receptor carrier speed")
    axes[0].legend()
    names = list(metrics)
    axes[1].bar(
        np.arange(len(names)) - 0.18,
        [metrics[name]["coarse_medium_relative_l2"] for name in names],
        width=0.36, label="coarse-medium",
    )
    axes[1].bar(
        np.arange(len(names)) + 0.18,
        [metrics[name]["medium_fine_relative_l2"] for name in names],
        width=0.36, label="medium-fine",
    )
    axes[1].set_xticks(np.arange(len(names)), names)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("relative L2 difference")
    axes[1].set_title("Carrier-flow time sensitivity")
    axes[1].legend()
    figure.suptitle("Active-case carrier-flow time isolation — " + record["status"])
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "field_metrics": metrics,
        "gate_reasons": gate_reasons,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
