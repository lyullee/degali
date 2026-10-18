"""Corner-reflection audit for intersecting D3Q19 pressure-open faces.

A compact linear acoustic packet travels at 45 degrees toward the x-max/y-max
corner.  A larger square domain supplies the no-return reference.  The audit
compares the original one-dimensional impedance and the angle-aware candidate
without tuning either boundary to the measured trace.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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
)


SCHEMA = "degali.pressure-open-corner-reflection-audit.v1"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-corner-reflection-audit-2026-10-10.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-corner-reflection-audit-2026-10-10.png"
)
SHORT_SIZE = 96
REFERENCE_SIZE = 192
NZ = 4
SOURCE_XY = 28.0
PROBE_XY = 60
PULSE_SIGMA_CELLS = 10.0
PULSE_WAVELENGTH_CELLS = 16.0
PULSE_AMPLITUDE = 1.0e-4
LATTICE_VISCOSITY = 0.02
STEPS = 315
REFLECTION_AMPLITUDE_RATIO_LIMIT = 0.05


def _initial_distribution(size: int) -> np.ndarray:
    coordinates = np.arange(size, dtype=float)
    yy, xx = np.meshgrid(coordinates, coordinates, indexing="ij")
    delta_x = xx - SOURCE_XY
    delta_y = yy - SOURCE_XY
    envelope = np.exp(
        -0.5 * (delta_x**2 + delta_y**2) / PULSE_SIGMA_CELLS**2
    )
    wave_number = 2.0 * math.pi / PULSE_WAVELENGTH_CELLS
    component_wave_number = wave_number / math.sqrt(2.0)
    phase = component_wave_number * (delta_x + delta_y)
    density_perturbation = PULSE_AMPLITUDE * envelope * np.cos(phase)
    density_plane = 1.0 + density_perturbation
    density = np.broadcast_to(
        density_plane[None, :, :], (NZ, size, size),
    ).copy()
    sound_speed = math.sqrt(_CS2_3D)
    component_velocity = (
        sound_speed / math.sqrt(2.0)
        * density_perturbation / density_plane
    )
    velocity_x = np.broadcast_to(
        component_velocity[None, :, :], density.shape,
    ).copy()
    velocity_y = velocity_x.copy()
    return _equilibrium_d3q19(
        density, velocity_x, velocity_y, np.zeros_like(density),
    )


def _probe_value(distribution: np.ndarray) -> float:
    density = np.sum(distribution, axis=0)
    # A 3x3 average is less sensitive to lattice phase than one point and is
    # fully inside the fluid/reference region.
    patch = density[
        :, PROBE_XY - 1:PROBE_XY + 2, PROBE_XY - 1:PROBE_XY + 2,
    ]
    return float(np.mean(patch) - 1.0)


def _trace(size: int, boundary_scheme: str) -> np.ndarray:
    distribution = _initial_distribution(size)
    solid = np.zeros((NZ, size, size), dtype=bool)
    trace = np.empty(STEPS + 1, dtype=float)
    for step in range(STEPS + 1):
        trace[step] = _probe_value(distribution)
        if step == STEPS:
            break
        distribution, _, _ = _advance_d3q19(
            distribution,
            solid=solid,
            inlet_u_lattice=np.zeros(NZ),
            lattice_viscosity=LATTICE_VISCOSITY,
            flow_boundary_scheme=boundary_scheme,
            outlet_lattice_density=1.0,
            nonstreamwise_wall_boundary="all_pressure_open",
            subgrid_model="none",
            smagorinsky_constant=0.0,
            wale_constant=0.0,
            maximum_relaxation_time=1.95,
            pressure_open_nonstream_faces=("y_max",),
            apply_streamwise_boundaries=True,
        )
    return trace


def _run_trace(specification: tuple[int, str]) -> tuple[str, np.ndarray, float]:
    size, scheme = specification
    started = time.perf_counter()
    return scheme, _trace(size, scheme), time.perf_counter() - started


def classify_corner_reflection(
    short_trace: np.ndarray,
    reference_trace: np.ndarray,
) -> dict[str, Any]:
    short = np.asarray(short_trace, dtype=float)
    reference = np.asarray(reference_trace, dtype=float)
    if short.shape != reference.shape or short.ndim != 1:
        raise ValueError("short and reference traces must be equal-length vectors")
    if np.any(~np.isfinite(short)) or np.any(~np.isfinite(reference)):
        raise ValueError("corner-reflection traces must be finite")
    sound_speed = math.sqrt(_CS2_3D)
    direct_distance = math.sqrt(2.0) * (PROBE_XY - SOURCE_XY)
    reflected_distance = math.sqrt(2.0) * (
        (SHORT_SIZE - 1.0 - SOURCE_XY)
        + (SHORT_SIZE - 1.0 - PROBE_XY)
    )
    direct_center = direct_distance / sound_speed
    reflected_center = reflected_distance / sound_speed
    half_window = 3.0 * PULSE_SIGMA_CELLS / sound_speed
    direct_slice = slice(
        max(0, int(math.floor(direct_center - half_window))),
        min(short.size, int(math.ceil(direct_center + half_window)) + 1),
    )
    # A spatially localised packet reaches each face before its centre reaches
    # the geometric corner. Begin immediately after the direct-pulse window so
    # early single-face and shared-corner returns cannot be omitted.
    reflected_slice = slice(
        direct_slice.stop,
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
    return {
        "incident_window_steps": [direct_slice.start, direct_slice.stop - 1],
        "reflection_window_steps": [
            reflected_slice.start, reflected_slice.stop - 1,
        ],
        "incident_amplitude": incident_amplitude,
        "reflected_amplitude": reflected_amplitude,
        "reflection_amplitude_ratio": ratio,
        "reflection_amplitude_ratio_limit": REFLECTION_AMPLITUDE_RATIO_LIMIT,
        "passed": bool(
            math.isfinite(ratio)
            and ratio <= REFLECTION_AMPLITUDE_RATIO_LIMIT
        ),
    }


def run_audit(workers: int = 3) -> dict[str, Any]:
    started = time.perf_counter()
    specifications = (
        (REFERENCE_SIZE, "characteristic_multidimensional"),
        (SHORT_SIZE, "characteristic_nonreflecting"),
        (SHORT_SIZE, "characteristic_multidimensional"),
    )
    if workers == 1:
        traces = [_run_trace(specification) for specification in specifications]
    else:
        with ProcessPoolExecutor(max_workers=min(workers, 3)) as executor:
            traces = list(executor.map(_run_trace, specifications))
    trace_by_key = {
        (specifications[index][0], scheme): trace
        for index, (scheme, trace, _runtime) in enumerate(traces)
    }
    runtime_by_key = {
        (specifications[index][0], scheme): runtime
        for index, (scheme, _trace_values, runtime) in enumerate(traces)
    }
    reference = trace_by_key[
        (REFERENCE_SIZE, "characteristic_multidimensional")
    ]
    schemes: dict[str, Any] = {}
    for scheme in (
        "characteristic_nonreflecting",
        "characteristic_multidimensional",
    ):
        short = trace_by_key[(SHORT_SIZE, scheme)]
        schemes[scheme] = {
            **classify_corner_reflection(short, reference),
            "runtime_s": runtime_by_key[(SHORT_SIZE, scheme)],
            "probe_trace": short.tolist(),
        }
    candidate = schemes["characteristic_multidimensional"]
    baseline = schemes["characteristic_nonreflecting"]
    passed = bool(candidate["passed"])
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "intersecting_open_face_corner_reflection_gate_passed"
            if passed
            else "intersecting_open_face_corner_reflection_gate_failed"
        ),
        "passed": passed,
        "promotion_allowed": False,
        "candidate_to_baseline_reflection_ratio": (
            float(candidate["reflection_amplitude_ratio"])
            / float(baseline["reflection_amplitude_ratio"])
        ),
        "configuration": {
            "short_domain_cells_xyz": [SHORT_SIZE, SHORT_SIZE, NZ],
            "reference_domain_cells_xyz": [
                REFERENCE_SIZE, REFERENCE_SIZE, NZ,
            ],
            "source_xy_cell": SOURCE_XY,
            "probe_xy_cell": PROBE_XY,
            "pulse_sigma_cells": PULSE_SIGMA_CELLS,
            "pulse_wavelength_cells": PULSE_WAVELENGTH_CELLS,
            "pulse_density_amplitude": PULSE_AMPLITUDE,
            "incidence_angle_degrees_from_each_face_normal": 45.0,
            "intersecting_open_faces": ["x_max", "y_max"],
            "z_boundary": "specular symmetry",
            "x_min_boundary": "zero-velocity non-equilibrium inlet",
            "y_min_boundary": "specular symmetry",
            "lattice_viscosity": LATTICE_VISCOSITY,
            "time_steps": STEPS,
        },
        "schemes": schemes,
        "reference_probe_trace": reference.tolist(),
        "reference_runtime_s": runtime_by_key[
            (REFERENCE_SIZE, "characteristic_multidimensional")
        ],
        "scope": (
            "linear_localised_two_face_corner_component_test_not_vortical_or_plume_validation"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_audit(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    reference = np.asarray(record["reference_probe_trace"], dtype=float)
    steps = np.arange(reference.size)
    window = record["schemes"]["characteristic_multidimensional"][
        "reflection_window_steps"
    ]
    figure, axes = plt.subplots(
        2, 1, figsize=(9.3, 6.8), sharex=True, constrained_layout=True,
    )
    axes[0].plot(steps, reference, color="black", label="long-domain reference")
    colors = {
        "characteristic_nonreflecting": "#d73027",
        "characteristic_multidimensional": "#1a9850",
    }
    labels = {
        "characteristic_nonreflecting": "1-D impedance",
        "characteristic_multidimensional": "angle-aware impedance",
    }
    for scheme, item in record["schemes"].items():
        trace = np.asarray(item["probe_trace"], dtype=float)
        axes[0].plot(steps, trace, color=colors[scheme], label=labels[scheme])
        axes[1].plot(
            steps, trace - reference, color=colors[scheme],
            label=(
                f"{labels[scheme]}: "
                f"{100.0 * float(item['reflection_amplitude_ratio']):.2f}%"
            ),
        )
    for axis in axes:
        axis.axvspan(window[0], window[1], color="#fdae61", alpha=0.15)
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    axes[0].set_ylabel("Probe density perturbation")
    axes[1].set_ylabel("Short - reference")
    axes[1].set_xlabel("Lattice time step")
    figure.suptitle(
        "x-max/y-max corner reflection: localised 45° acoustic packet",
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
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    record = run_audit(args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    plot_audit(record, args.plot)
    baseline = record["schemes"]["characteristic_nonreflecting"]
    candidate = record["schemes"]["characteristic_multidimensional"]
    print(json.dumps({
        "status": record["status"],
        "baseline_reflection_amplitude_ratio": baseline[
            "reflection_amplitude_ratio"
        ],
        "candidate_reflection_amplitude_ratio": candidate[
            "reflection_amplitude_ratio"
        ],
        "reflection_amplitude_ratio_limit": candidate[
            "reflection_amplitude_ratio_limit"
        ],
        "runtime_s": record["runtime_s"],
    }, indent=2, allow_nan=False))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
