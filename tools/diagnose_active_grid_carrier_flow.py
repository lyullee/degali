"""Separate carrier-flow grid sensitivity from coupled H2 transport.

This diagnostic reuses the physical domain, obstacle, viscosity, wind and
fixed-Mach grid family from ``audit_active_thermal_corrector_grid.py`` but
stops before species/enthalpy transport.  It therefore answers whether the
roof/side concentration sensitivity is already present in the spun-up carrier
flow or is introduced later by scalar transport/observation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.fast_obstacle_transport_3d import (  # noqa: E402
    FastObstacleReceptor3D,
    _receptor_stencil,
)
from degali.addons.resolved_obstacle_flow_3d import _solve_flow  # noqa: E402

from audit_active_thermal_corrector_grid import (  # noqa: E402
    OUTPUT,
    SPACINGS_M,
    _config,
)


RECEPTORS = (
    FastObstacleReceptor3D("windward", (1.05, 1.80, 0.60)),
    FastObstacleReceptor3D("roof", (1.80, 1.80, 1.35)),
    FastObstacleReceptor3D("side", (1.80, 1.05, 0.60)),
    FastObstacleReceptor3D("lee", (2.70, 1.80, 0.60)),
)


def _sample(field: np.ndarray, stencil) -> float:
    return float(sum(weight * field[tuple(index)] for *index, weight in stencil))


def _relative_l2(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right) / max(
        np.linalg.norm(left), np.linalg.norm(right), 1.0e-15,
    ))


def main() -> int:
    records = []
    results = []
    for spacing in SPACINGS_M:
        config = _config(spacing)
        started = time.perf_counter()
        flow = _solve_flow(config)
        runtime = time.perf_counter() - started
        (
            x, y, z, solid, u, v, w, u_face, v_face, w_face,
            _pressure, _pressure_coefficient, _eddy_viscosity,
            residual_history, _lattice_density, flow_dt, base_tau, maximum_tau,
            actual_steps, _raw_divergence, divergence_before,
            divergence_after, projection_iterations, _distribution,
            statistical_residual, average_count,
        ) = flow
        receptor_values = {}
        for receptor in RECEPTORS:
            stencil = _receptor_stencil(x, y, z, solid, receptor)
            components = np.asarray((
                _sample(u, stencil), _sample(v, stencil), _sample(w, stencil),
            ))
            receptor_values[receptor.label] = {
                "u_m_s": float(components[0]),
                "v_m_s": float(components[1]),
                "w_m_s": float(components[2]),
                "speed_m_s": float(np.linalg.norm(components)),
            }
        residual = (
            None if residual_history.size == 0
            else float(residual_history[-1, 1])
        )
        records.append({
            "cell_size_m": spacing,
            "cells": list(config.cells),
            "flow_dt_s": flow_dt,
            "flow_steps": actual_steps,
            "physical_flow_horizon_s": actual_steps * flow_dt,
            "base_relaxation_time": base_tau,
            "maximum_relaxation_time": maximum_tau,
            "final_instantaneous_residual": residual,
            "statistical_residual": statistical_residual,
            "average_sample_count": average_count,
            "maximum_projected_divergence_s_1": divergence_after,
            "projection_iterations": projection_iterations,
            "runtime_s": runtime,
            "receptors": receptor_values,
        })
        results.append((x, y, z, solid, u, v, w, u_face, v_face, w_face))

    coarse = results[0]
    zz, yy, xx = np.meshgrid(coarse[2], coarse[1], coarse[0], indexing="ij")
    points = np.column_stack((zz.ravel(), yy.ravel(), xx.ravel()))
    component_samples: dict[str, list[np.ndarray]] = {
        "u": [], "v": [], "w": [], "speed": [],
    }
    valid = np.ones(points.shape[0], dtype=bool)
    for x, y, z, solid, u, v, w, *_faces in results:
        fields = {"u": u, "v": v, "w": w, "speed": np.sqrt(u*u + v*v + w*w)}
        for name, field in fields.items():
            masked = np.where(solid, np.nan, field)
            sampled = RegularGridInterpolator(
                (z, y, x), masked, bounds_error=False, fill_value=np.nan,
            )(points)
            component_samples[name].append(sampled)
            valid &= np.isfinite(sampled)

    field_metrics = {}
    for name, samples in component_samples.items():
        restricted = [sample[valid] for sample in samples]
        coarse_medium = _relative_l2(restricted[0], restricted[1])
        medium_fine = _relative_l2(restricted[1], restricted[2])
        order = (
            None if coarse_medium <= 0.0 or medium_fine <= 0.0
            else math.log(coarse_medium / medium_fine) / math.log(1.5)
        )
        field_metrics[name] = {
            "coarse_medium_relative_l2": coarse_medium,
            "medium_fine_relative_l2": medium_fine,
            "observed_order": order,
            "successive_difference_decreased": medium_fine < coarse_medium,
        }

    receptor_changes = {}
    for receptor in RECEPTORS:
        label = receptor.label
        receptor_changes[label] = {}
        for component in ("u_m_s", "v_m_s", "w_m_s", "speed_m_s"):
            medium = float(records[1]["receptors"][label][component])
            fine = float(records[2]["receptors"][label][component])
            receptor_changes[label][component] = abs(fine - medium) / max(
                abs(medium), abs(fine), 0.01,
            )

    gate_reasons = []
    for name, metric in field_metrics.items():
        if not metric["successive_difference_decreased"]:
            gate_reasons.append(f"{name} field difference is non-monotonic")
        if metric["medium_fine_relative_l2"] > 0.10:
            gate_reasons.append(f"{name} field changes by more than 10%")
    for label, changes in receptor_changes.items():
        if changes["speed_m_s"] > 0.10:
            gate_reasons.append(f"{label} speed changes by more than 10%")

    record = {
        "schema": "degali.active-grid-carrier-flow-diagnostic.v1",
        "status": (
            "accepted_component_grid"
            if not gate_reasons else "withheld_component_grid"
        ),
        "scope": (
            "carrier flow only; no H2 species, enthalpy, EOS density or buoyancy"
        ),
        "common_point_count": int(np.count_nonzero(valid)),
        "cases_coarse_to_fine": records,
        "field_metrics": field_metrics,
        "finest_pair_receptor_relative_changes": receptor_changes,
        "gate_reasons": gate_reasons,
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "active_grid_carrier_flow_diagnostic.json"
    figure_path = OUTPUT / "active_grid_carrier_flow_diagnostic.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for result, case, label in zip(results, records, ("coarse", "medium", "fine")):
        x, y, z, solid, u, v, w, *_ = result
        center = int(np.argmin(np.abs(y - 1.80)))
        speed = np.sqrt(u*u + v*v + w*w)
        speed = np.where(solid, np.nan, speed)
        axes[0].plot(x, speed[int(np.argmin(np.abs(z - 1.35))), center, :], label=label)
        axes[1].plot(
            [item.position_m[0] for item in RECEPTORS],
            [case["receptors"][item.label]["speed_m_s"] for item in RECEPTORS],
            "o-", label=label,
        )
    axes[0].set_title("Centre-plane speed near roof height")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("speed (m/s)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].set_title("Fixed receptor speed")
    axes[1].set_xticks(
        [item.position_m[0] for item in RECEPTORS],
        [item.label for item in RECEPTORS], rotation=25,
    )
    axes[1].set_ylabel("speed (m/s)")
    axes[1].legend()
    axes[1].grid(alpha=0.25)
    names = list(field_metrics)
    axes[2].bar(
        np.arange(len(names)) - 0.18,
        [field_metrics[name]["coarse_medium_relative_l2"] for name in names],
        width=0.36, label="coarse-medium",
    )
    axes[2].bar(
        np.arange(len(names)) + 0.18,
        [field_metrics[name]["medium_fine_relative_l2"] for name in names],
        width=0.36, label="medium-fine",
    )
    axes[2].axhline(0.10, color="red", linestyle="--", linewidth=1)
    axes[2].set_xticks(np.arange(len(names)), names)
    axes[2].set_ylabel("relative L2 difference")
    axes[2].set_title("Common-point field sensitivity")
    axes[2].legend()
    axes[2].grid(axis="y", alpha=0.25)
    figure.suptitle("Active-grid carrier-flow isolation — " + record["status"])
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "field_metrics": field_metrics,
        "receptor_changes": receptor_changes,
        "gate_reasons": gate_reasons,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
