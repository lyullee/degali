"""Time-refinement isolation for the finite-volume momentum/projection path."""

from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.resolved_obstacle_flow_3d import (  # noqa: E402
    solve_resolved_obstacle_flow_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule  # noqa: E402
from audit_active_thermal_corrector_time import (  # noqa: E402
    DURATION_S,
    DX_M,
    LATTICE_SPEEDS,
    OUTPUT,
    WIND_M_S,
    _config,
)


TIME_LEVEL_LATTICE_SPEEDS = (*LATTICE_SPEEDS, 0.00375)


def _relative_l2(left: np.ndarray, right: np.ndarray, fluid: np.ndarray) -> float:
    return float(np.linalg.norm(left[fluid] - right[fluid]) / max(
        np.linalg.norm(left[fluid]), np.linalg.norm(right[fluid]), 1.0e-15,
    ))


def main() -> int:
    base = _config(0.015)
    results = []
    cases = []
    for speed in TIME_LEVEL_LATTICE_SPEEDS:
        dt = speed * DX_M / WIND_M_S
        config = replace(
            base,
            scalar_time_step_s=dt,
            source_schedule=SourceRateSchedule(
                (0.0, DURATION_S), (0.0, 0.0),
            ),
            source_temperature_k=base.ambient_temperature_k,
            gravity_m_s2=0.0,
            active_momentum_solver="finite_volume_projection",
            effective_schmidt_number=0.7,
            effective_prandtl_number=0.85,
        )
        result = solve_resolved_obstacle_flow_3d(config)
        results.append(result)
        cases.append({
            "time_step_s": dt,
            "steps": int(round(DURATION_S / dt)),
            "runtime_s": result.diagnostics["total_runtime_s"],
            "maximum_momentum_balance_residual_kg_m_s": result.diagnostics[
                "maximum_momentum_balance_residual_kg_m_s"
            ],
            "maximum_continuity_residual_kg_m3_s": result.diagnostics[
                "maximum_variable_density_continuity_residual_kg_m3_s"
            ],
        })
    fluid = ~results[-1].solid_mask
    metrics = {}
    for name, getter in (
        ("u", lambda result: result.velocity_u_m_s),
        ("v", lambda result: result.velocity_v_m_s),
        ("w", lambda result: result.velocity_w_m_s),
        ("speed", lambda result: np.sqrt(
            result.velocity_u_m_s**2
            + result.velocity_v_m_s**2
            + result.velocity_w_m_s**2
        )),
    ):
        fields = [getter(result) for result in results]
        pair_differences = [
            _relative_l2(left, right, fluid)
            for left, right in zip(fields, fields[1:])
        ]
        metrics[name] = {
            "successive_relative_l2": pair_differences,
            "successive_observed_orders": [
                (
                    None if min(left, right) <= 0.0
                    else math.log(left / right) / math.log(2.0)
                )
                for left, right in zip(pair_differences, pair_differences[1:])
            ],
            "all_successive_differences_decreased": all(
                right < left
                for left, right in zip(pair_differences, pair_differences[1:])
            ),
            "finest_pair_recovered": pair_differences[-1] < pair_differences[-2],
        }
    reasons = [
        f"{name} field difference is non-monotonic"
        for name, metric in metrics.items()
        if not metric["all_successive_differences_decreased"]
    ]
    record = {
        "schema": "degali.fv-momentum-time-diagnostic.v1",
        "status": "accepted_component_time" if not reasons else "withheld_component_time",
        "scope": (
            "finite-volume momentum plus variable-density projection with zero "
            "release, zero buoyancy and constant thermodynamic density"
        ),
        "same_initial_d3q19_carrier": True,
        "cases_coarse_to_fine": cases,
        "field_metrics": metrics,
        "gate_reasons": reasons,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "fv_momentum_time_diagnostic.json"
    figure_path = OUTPUT / "fv_momentum_time_diagnostic.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    names = list(metrics)
    figure, axis = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    positions = np.arange(len(names))
    axis.bar(
        positions - 0.18,
        [metrics[name]["successive_relative_l2"][0] for name in names],
        width=0.36, label="coarse–medium",
    )
    axis.bar(
        positions + 0.18,
        [metrics[name]["successive_relative_l2"][1] for name in names],
        width=0.36, label="medium–fine",
    )
    axis.scatter(
        positions,
        [metrics[name]["successive_relative_l2"][2] for name in names],
        marker="D", color="black", label="fine–ultrafine", zorder=3,
    )
    axis.set_xticks(positions, names)
    axis.set_yscale("log")
    axis.set_ylabel("relative L2 difference")
    axis.set_title("Finite-volume momentum/projection time isolation — " + record["status"])
    axis.legend()
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)
    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "field_metrics": metrics,
        "gate_reasons": reasons,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
