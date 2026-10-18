"""Linear acoustic-pulse reflection audit for the D3Q19 pressure outlet.

The short and long domains start from the same small-amplitude, nominally
right-running density pulse.  The long-domain probe is a no-return reference
over the audit window.  Subtracting its trace from the short-domain trace
isolates the wave returning from the pressure-open outlet.  This canonical
component test does not establish multi-dimensional plume domain independence.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np

from degali.addons.resolved_obstacle_flow_3d import (
    _CS2_3D,
    _advance_d3q19,
    _equilibrium_d3q19,
    _moments_d3q19,
)


SCHEMA = "degali.pressure-open-boundary-reflection-audit.v2"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-boundary-reflection-audit-2026-10-10-v2.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-boundary-reflection-audit-2026-10-10-v2.png"
)
SHORT_NX = 180
REFERENCE_NX = 360
NY = 8
NZ = 8
SOURCE_X = 40.0
PROBE_X = 120
PULSE_SIGMA_CELLS = 8.0
PULSE_AMPLITUDE = 1.0e-4
LATTICE_VISCOSITY = 0.02
STEPS = 430
REFLECTION_AMPLITUDE_RATIO_LIMIT = 0.05


def _initial_distribution(length: int, propagation_axis: str) -> np.ndarray:
    coordinate = np.arange(length, dtype=float)
    density_line = 1.0 + PULSE_AMPLITUDE * np.exp(
        -0.5 * ((coordinate - SOURCE_X) / PULSE_SIGMA_CELLS) ** 2
    )
    sound_speed = math.sqrt(_CS2_3D)
    velocity_line = sound_speed * (density_line - 1.0) / density_line
    if propagation_axis == "x":
        density = np.broadcast_to(density_line, (NZ, NY, length)).copy()
        velocity_x = np.broadcast_to(velocity_line, density.shape).copy()
        velocity_y = np.zeros_like(density)
    elif propagation_axis == "y":
        density = np.broadcast_to(
            density_line[None, :, None], (NZ, length, NY),
        ).copy()
        velocity_x = np.zeros_like(density)
        velocity_y = np.broadcast_to(
            velocity_line[None, :, None], density.shape,
        ).copy()
    else:
        raise ValueError("propagation_axis must be 'x' or 'y'")
    zero = np.zeros_like(density)
    return _equilibrium_d3q19(density, velocity_x, velocity_y, zero)


def _trace(
    length: int, boundary_scheme: str, propagation_axis: str = "x",
) -> np.ndarray:
    distribution = _initial_distribution(length, propagation_axis)
    solid = np.zeros(distribution.shape[1:], dtype=bool)
    inlet = np.zeros(NZ)
    trace = np.empty(STEPS + 1, dtype=float)
    for step in range(STEPS + 1):
        density = np.sum(distribution, axis=0)
        if propagation_axis == "x":
            trace[step] = float(np.mean(density[:, :, PROBE_X] - 1.0))
        else:
            trace[step] = float(np.mean(density[:, PROBE_X, :] - 1.0))
        if step == STEPS:
            break
        distribution, _, _ = _advance_d3q19(
            distribution,
            solid=solid,
            inlet_u_lattice=inlet,
            lattice_viscosity=LATTICE_VISCOSITY,
            flow_boundary_scheme=boundary_scheme,
            outlet_lattice_density=1.0,
            nonstreamwise_wall_boundary=(
                "all_symmetry_free_slip"
                if propagation_axis == "x" else "all_pressure_open"
            ),
            subgrid_model="none",
            smagorinsky_constant=0.0,
            wale_constant=0.0,
            maximum_relaxation_time=1.95,
            pressure_open_nonstream_faces=("y_min", "y_max"),
            apply_streamwise_boundaries=(propagation_axis == "x"),
        )
    return trace


def classify_reflection(
    short_trace: np.ndarray,
    reference_trace: np.ndarray,
    *,
    orientation: str = "streamwise",
) -> dict[str, Any]:
    short = np.asarray(short_trace, dtype=float)
    reference = np.asarray(reference_trace, dtype=float)
    if short.shape != reference.shape or short.ndim != 1:
        raise ValueError("short and reference traces must be equal-length vectors")
    sound_speed = math.sqrt(_CS2_3D)
    direct_center = (PROBE_X - SOURCE_X) / sound_speed
    reflected_center = (
        (SHORT_NX - 1.0 - SOURCE_X)
        + (SHORT_NX - 1.0 - PROBE_X)
    ) / sound_speed
    half_window = 3.0 * PULSE_SIGMA_CELLS / sound_speed
    direct_slice = slice(
        max(0, int(math.floor(direct_center - half_window))),
        min(short.size, int(math.ceil(direct_center + half_window)) + 1),
    )
    reflected_slice = slice(
        max(0, int(math.floor(reflected_center - half_window))),
        min(short.size, int(math.ceil(reflected_center + half_window)) + 1),
    )
    incident_amplitude = float(np.max(np.abs(reference[direct_slice])))
    reflected_amplitude = float(np.max(np.abs(
        short[reflected_slice] - reference[reflected_slice]
    )))
    ratio = (
        math.inf if incident_amplitude <= 0.0
        else reflected_amplitude / incident_amplitude
    )
    passed = bool(
        math.isfinite(ratio)
        and ratio <= REFLECTION_AMPLITUDE_RATIO_LIMIT
    )
    return {
        "incident_window_steps": [direct_slice.start, direct_slice.stop - 1],
        "reflection_window_steps": [
            reflected_slice.start, reflected_slice.stop - 1,
        ],
        "incident_density_amplitude": incident_amplitude,
        "reflected_density_amplitude": reflected_amplitude,
        "reflection_amplitude_ratio": ratio,
        "reflection_amplitude_ratio_limit": (
            REFLECTION_AMPLITUDE_RATIO_LIMIT
        ),
        "passed": passed,
        "status": (
            f"canonical_{orientation}_reflection_gate_passed"
            if passed else f"canonical_{orientation}_reflection_gate_failed"
        ),
        "promotion_allowed": False,
    }


def run_audit() -> dict[str, Any]:
    started = time.perf_counter()
    short = _trace(SHORT_NX, "characteristic_nonreflecting")
    reference = _trace(REFERENCE_NX, "characteristic_nonreflecting")
    fixed_pressure_short = _trace(
        SHORT_NX, "nonequilibrium_extrapolation",
    )
    classification = classify_reflection(short, reference)
    fixed_pressure_classification = classify_reflection(
        fixed_pressure_short, reference,
    )
    lateral_short = _trace(
        SHORT_NX, "characteristic_nonreflecting", "y",
    )
    lateral_reference = _trace(
        REFERENCE_NX, "characteristic_nonreflecting", "y",
    )
    lateral_classification = classify_reflection(
        lateral_short, lateral_reference, orientation="lateral_normal",
    )
    overall_passed = bool(
        classification["passed"] and lateral_classification["passed"]
    )
    uniform = np.ones((NZ, NY, 16))
    zero = np.zeros_like(uniform)
    equilibrium = _equilibrium_d3q19(uniform, zero, zero, zero)
    solid = np.zeros_like(uniform, dtype=bool)
    equilibrium, _, _ = _advance_d3q19(
        equilibrium,
        solid=solid,
        inlet_u_lattice=np.zeros(NZ),
        lattice_viscosity=LATTICE_VISCOSITY,
        flow_boundary_scheme="characteristic_nonreflecting",
        outlet_lattice_density=1.0,
        nonstreamwise_wall_boundary="all_pressure_open",
        subgrid_model="none",
        smagorinsky_constant=0.0,
        wale_constant=0.0,
        maximum_relaxation_time=1.95,
    )
    rho, ux, uy, uz = _moments_d3q19(equilibrium, solid)
    uniform_error = max(
        float(np.max(np.abs(rho - 1.0))),
        float(np.max(np.abs(ux))),
        float(np.max(np.abs(uy))),
        float(np.max(np.abs(uz))),
    )
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **classification,
        "status": (
            "canonical_normal_reflection_gates_passed"
            if overall_passed else "canonical_normal_reflection_gate_failed"
        ),
        "passed": overall_passed,
        "streamwise_status": classification["status"],
        "lateral_normal_status": lateral_classification["status"],
        "lateral_normal_reflection_amplitude_ratio": (
            lateral_classification["reflection_amplitude_ratio"]
        ),
        "lateral_normal_incident_density_amplitude": (
            lateral_classification["incident_density_amplitude"]
        ),
        "lateral_normal_reflected_density_amplitude": (
            lateral_classification["reflected_density_amplitude"]
        ),
        "fixed_pressure_baseline_reflection_amplitude_ratio": (
            fixed_pressure_classification["reflection_amplitude_ratio"]
        ),
        "uniform_state_maximum_absolute_error": uniform_error,
        "configuration": {
            "short_domain_cells": [SHORT_NX, NY, NZ],
            "reference_domain_cells": [REFERENCE_NX, NY, NZ],
            "source_x_cell": SOURCE_X,
            "probe_x_cell": PROBE_X,
            "pulse_sigma_cells": PULSE_SIGMA_CELLS,
            "pulse_density_amplitude": PULSE_AMPLITUDE,
            "lattice_viscosity": LATTICE_VISCOSITY,
            "time_steps": STEPS,
            "streamwise_boundary": "isothermal characteristic outlet",
            "fixed_pressure_baseline_boundary": (
                "nonequilibrium extrapolation fixed-density outlet"
            ),
            "transverse_boundary": "specular symmetry for 1-D isolation",
            "lateral_test_streamwise_boundary": "periodic",
            "lateral_test_vertical_boundary": "specular symmetry",
            "lateral_test_open_faces": ["y_min", "y_max"],
        },
        "short_probe_density_perturbation": short.tolist(),
        "fixed_pressure_short_probe_density_perturbation": (
            fixed_pressure_short.tolist()
        ),
        "reference_probe_density_perturbation": reference.tolist(),
        "lateral_short_probe_density_perturbation": lateral_short.tolist(),
        "lateral_reference_probe_density_perturbation": (
            lateral_reference.tolist()
        ),
        "scope": (
            "normal_incidence_linear_D3Q19_component_tests_not_oblique_or_plume_domain_tests"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_audit(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    short = np.asarray(record["short_probe_density_perturbation"])
    fixed_pressure = np.asarray(
        record["fixed_pressure_short_probe_density_perturbation"]
    )
    reference = np.asarray(record["reference_probe_density_perturbation"])
    lateral_short = np.asarray(
        record["lateral_short_probe_density_perturbation"]
    )
    lateral_reference = np.asarray(
        record["lateral_reference_probe_density_perturbation"]
    )
    steps = np.arange(short.size)
    window = record["reflection_window_steps"]
    figure, axes = plt.subplots(
        2, 2, figsize=(12.4, 6.6), sharex=True, constrained_layout=True,
    )
    axes[0, 0].plot(
        steps, reference, label="long-domain reference", color="#2166ac",
    )
    axes[0, 0].plot(
        steps, fixed_pressure, label="short fixed-pressure outlet",
        color="#b2182b", alpha=0.75,
    )
    axes[0, 0].plot(
        steps, short, label="short characteristic outlet", color="#1b7837",
    )
    axes[0, 0].axvspan(
        window[0], window[1], alpha=0.15, color="#fdae61",
    )
    axes[0, 0].set_title("Streamwise x-max outlet")
    axes[0, 0].set_ylabel(r"Probe $\rho-1$")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(alpha=0.25)
    axes[1, 0].plot(
        steps, fixed_pressure - reference,
        color="#b2182b", alpha=0.75, label="fixed-pressure - reference",
    )
    axes[1, 0].plot(
        steps, short - reference,
        color="#1b7837", label="characteristic - reference",
    )
    axes[1, 0].axvspan(
        window[0], window[1], alpha=0.15, color="#fdae61",
    )
    axes[1, 0].set_xlabel("Lattice time step")
    axes[1, 0].set_ylabel("Short - reference")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(frameon=False, fontsize=8)

    axes[0, 1].plot(
        steps, lateral_reference, label="long-domain reference",
        color="#2166ac",
    )
    axes[0, 1].plot(
        steps, lateral_short, label="short characteristic side",
        color="#1b7837",
    )
    axes[0, 1].axvspan(
        window[0], window[1], alpha=0.15, color="#fdae61",
    )
    axes[0, 1].set_title("Lateral y-max outlet")
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].grid(alpha=0.25)
    axes[1, 1].plot(
        steps, lateral_short - lateral_reference, color="#1b7837",
        label="lateral characteristic - reference",
    )
    axes[1, 1].axvspan(
        window[0], window[1], alpha=0.15, color="#fdae61",
    )
    axes[1, 1].set_xlabel("Lattice time step")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(frameon=False, fontsize=8)
    ratio = float(record["reflection_amplitude_ratio"])
    lateral_ratio = float(record["lateral_normal_reflection_amplitude_ratio"])
    figure.suptitle(
        f"D3Q19 characteristic reflection: x={ratio:.3%}, "
        f"y={lateral_ratio:.3%} "
        f"({record['status']})",
        fontsize=11,
    )
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    args = parser.parse_args(argv)
    record = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    plot_audit(record, args.plot)
    print(json.dumps({
        key: record[key]
        for key in (
            "status", "reflection_amplitude_ratio",
            "lateral_normal_reflection_amplitude_ratio",
            "fixed_pressure_baseline_reflection_amplitude_ratio",
            "reflection_amplitude_ratio_limit",
            "uniform_state_maximum_absolute_error", "runtime_s",
        )
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
