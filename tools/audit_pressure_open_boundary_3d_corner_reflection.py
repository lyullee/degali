"""Three-dimensional edge/corner reflection audit for D3Q19 open faces.

Localised linear acoustic packets are directed toward x-max/y-max (scaled
control), x-max/z-max, y-max/z-max and x-max/y-max/z-max intersections.  A
doubled domain supplies the no-return reference for each orientation.  The
frozen 5% amplitude gate is applied without fitting the boundary operator to
any trace.  The two-face packet is an exact half-scale replica of the already
qualified 96-cell x/y audit, including its 8-cell half-scale wavelength.
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
    _C3,
    _CS2_3D,
    _advance_d3q19,
    _equilibrium_d3q19,
)


SCHEMA = "degali.pressure-open-3d-corner-reflection-audit.v1"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-3d-corner-reflection-audit-2026-10-10.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/pressure-open-3d-corner-reflection-audit-2026-10-10.png"
)
SHORT_SIZE = 48
REFERENCE_SIZE = 96
INACTIVE_SIZE = 4
SOURCE_CELL = 14.0
PROBE_CELL = 30
PULSE_SIGMA_CELLS = 5.0
PULSE_WAVELENGTH_CELLS = 8.0
PULSE_AMPLITUDE = 1.0e-4
LATTICE_VISCOSITY = 0.02
STEPS = 210
REFLECTION_AMPLITUDE_RATIO_LIMIT = 0.05
ORIENTATIONS = {
    "x_y_scaled_control": ("x", "y"),
    "x_z_edge": ("x", "z"),
    "y_z_edge": ("y", "z"),
    "x_y_z_corner": ("x", "y", "z"),
}


def _shape(size: int, active_axes: tuple[str, ...]) -> tuple[int, int, int]:
    lengths = {
        axis: size if axis in active_axes else INACTIVE_SIZE
        for axis in ("x", "y", "z")
    }
    return lengths["z"], lengths["y"], lengths["x"]


def _initial_distribution(
    size: int,
    active_axes: tuple[str, ...],
    pulse_wavelength_cells: float = PULSE_WAVELENGTH_CELLS,
) -> np.ndarray:
    shape = _shape(size, active_axes)
    coordinates = {
        "x": np.arange(shape[2], dtype=float)[None, None, :],
        "y": np.arange(shape[1], dtype=float)[None, :, None],
        "z": np.arange(shape[0], dtype=float)[:, None, None],
    }
    deltas = {
        axis: coordinates[axis] - SOURCE_CELL
        for axis in active_axes
    }
    radius_squared = np.zeros(shape, dtype=float)
    phase_coordinate = np.zeros(shape, dtype=float)
    for delta in deltas.values():
        radius_squared = radius_squared + delta**2
        phase_coordinate = phase_coordinate + delta
    envelope = np.exp(-0.5 * radius_squared / PULSE_SIGMA_CELLS**2)
    if not math.isfinite(pulse_wavelength_cells) or pulse_wavelength_cells <= 0.0:
        raise ValueError("pulse_wavelength_cells must be positive and finite")
    wave_number = 2.0 * math.pi / pulse_wavelength_cells
    phase = wave_number * phase_coordinate / math.sqrt(len(active_axes))
    perturbation = PULSE_AMPLITUDE * envelope * np.cos(phase)
    density = 1.0 + perturbation
    component_velocity = (
        math.sqrt(_CS2_3D) / math.sqrt(len(active_axes))
        * perturbation / density
    )
    zero = np.zeros(shape, dtype=float)
    velocity = {
        axis: component_velocity if axis in active_axes else zero
        for axis in ("x", "y", "z")
    }
    return _equilibrium_d3q19(
        density, velocity["x"], velocity["y"], velocity["z"],
    )


def _probe_value(
    distribution: np.ndarray, active_axes: tuple[str, ...],
) -> float:
    density = np.sum(distribution, axis=0)
    selections = []
    for axis in ("z", "y", "x"):
        if axis in active_axes:
            selections.append(slice(PROBE_CELL - 1, PROBE_CELL + 2))
        else:
            selections.append(slice(None))
    return float(np.mean(density[tuple(selections)]) - 1.0)


def _trace(
    size: int,
    boundary_scheme: str,
    active_axes: tuple[str, ...],
    pulse_wavelength_cells: float,
) -> np.ndarray:
    distribution = _initial_distribution(
        size, active_axes, pulse_wavelength_cells,
    )
    shape = _shape(size, active_axes)
    solid = np.zeros(shape, dtype=bool)
    nonstream_faces = tuple(
        f"{axis}_max" for axis in active_axes if axis in {"y", "z"}
    )
    trace = np.empty(STEPS + 1, dtype=float)
    for step in range(STEPS + 1):
        trace[step] = _probe_value(distribution, active_axes)
        if step == STEPS:
            break
        distribution, _, _ = _advance_d3q19(
            distribution,
            solid=solid,
            inlet_u_lattice=np.zeros(shape[0]),
            lattice_viscosity=LATTICE_VISCOSITY,
            flow_boundary_scheme=boundary_scheme,
            outlet_lattice_density=1.0,
            nonstreamwise_wall_boundary="all_pressure_open",
            subgrid_model="none",
            smagorinsky_constant=0.0,
            wale_constant=0.0,
            maximum_relaxation_time=1.95,
            pressure_open_nonstream_faces=nonstream_faces,
            apply_streamwise_boundaries="x" in active_axes,
        )
    return trace


def _run_trace(
    specification: tuple[str, int, str, float],
) -> tuple[str, int, str, np.ndarray, float]:
    orientation, size, scheme, pulse_wavelength_cells = specification
    active_axes = ORIENTATIONS[orientation]
    started = time.perf_counter()
    trace = _trace(size, scheme, active_axes, pulse_wavelength_cells)
    return orientation, size, scheme, trace, time.perf_counter() - started


def classify_corner_reflection(
    short_trace: np.ndarray,
    reference_trace: np.ndarray,
    *,
    active_axis_count: int,
) -> dict[str, Any]:
    short = np.asarray(short_trace, dtype=float)
    reference = np.asarray(reference_trace, dtype=float)
    if short.shape != reference.shape or short.ndim != 1:
        raise ValueError("short and reference traces must be equal-length vectors")
    if np.any(~np.isfinite(short)) or np.any(~np.isfinite(reference)):
        raise ValueError("corner-reflection traces must be finite")
    if active_axis_count not in {2, 3}:
        raise ValueError("active_axis_count must be two or three")
    sound_speed = math.sqrt(_CS2_3D)
    diagonal = math.sqrt(active_axis_count)
    direct_center = diagonal * (PROBE_CELL - SOURCE_CELL) / sound_speed
    reflected_center = diagonal * (
        (SHORT_SIZE - 1.0 - SOURCE_CELL)
        + (SHORT_SIZE - 1.0 - PROBE_CELL)
    ) / sound_speed
    half_window = 3.0 * PULSE_SIGMA_CELLS / sound_speed
    direct_slice = slice(
        max(0, int(math.floor(direct_center - half_window))),
        min(short.size, int(math.ceil(direct_center + half_window)) + 1),
    )
    reflected_slice = slice(
        direct_slice.stop,
        min(short.size, int(math.ceil(reflected_center + half_window)) + 1),
    )
    incident = float(np.max(np.abs(reference[direct_slice])))
    reflected = float(np.max(np.abs(
        short[reflected_slice] - reference[reflected_slice]
    )))
    ratio = math.inf if incident <= 0.0 else reflected / incident
    return {
        "incident_window_steps": [direct_slice.start, direct_slice.stop - 1],
        "reflection_window_steps": [
            reflected_slice.start, reflected_slice.stop - 1,
        ],
        "incident_amplitude": incident,
        "reflected_amplitude": reflected,
        "reflection_amplitude_ratio": ratio,
        "reflection_amplitude_ratio_limit": REFLECTION_AMPLITUDE_RATIO_LIMIT,
        "passed": bool(
            math.isfinite(ratio)
            and ratio <= REFLECTION_AMPLITUDE_RATIO_LIMIT
        ),
    }


def run_audit(
    workers: int = 3,
    pulse_wavelength_cells: float = PULSE_WAVELENGTH_CELLS,
    candidate_boundary_scheme: str = "characteristic_multidimensional",
) -> dict[str, Any]:
    if not math.isfinite(pulse_wavelength_cells) or pulse_wavelength_cells <= 0.0:
        raise ValueError("pulse_wavelength_cells must be positive and finite")
    if candidate_boundary_scheme not in {
        "characteristic_multidimensional",
        "characteristic_multidimensional_extrapolated",
        "characteristic_multidimensional_vortical_limited",
        "characteristic_nscbc_convective",
    }:
        raise ValueError("unsupported candidate_boundary_scheme")
    started = time.perf_counter()
    specifications = tuple(
        (orientation, size, scheme, pulse_wavelength_cells)
        for orientation in ORIENTATIONS
        for size, scheme in (
            (REFERENCE_SIZE, "characteristic_multidimensional"),
            (SHORT_SIZE, "characteristic_nonreflecting"),
            (SHORT_SIZE, candidate_boundary_scheme),
        )
    )
    if workers == 1:
        traces = [_run_trace(item) for item in specifications]
    else:
        with ProcessPoolExecutor(max_workers=min(workers, 3)) as executor:
            traces = list(executor.map(_run_trace, specifications))
    by_key = {
        (orientation, size, scheme): (trace, runtime)
        for orientation, size, scheme, trace, runtime in traces
    }
    orientations: dict[str, Any] = {}
    for orientation, active_axes in ORIENTATIONS.items():
        reference, reference_runtime = by_key[
            (orientation, REFERENCE_SIZE, "characteristic_multidimensional")
        ]
        schemes = {}
        for scheme in (
            "characteristic_nonreflecting",
            candidate_boundary_scheme,
        ):
            trace, runtime = by_key[(orientation, SHORT_SIZE, scheme)]
            schemes[scheme] = {
                **classify_corner_reflection(
                    trace, reference, active_axis_count=len(active_axes),
                ),
                "runtime_s": runtime,
                "probe_trace": trace.tolist(),
            }
        orientations[orientation] = {
            "active_axes": list(active_axes),
            "intersecting_open_faces": [f"{axis}_max" for axis in active_axes],
            "passed": schemes[candidate_boundary_scheme]["passed"],
            "candidate_to_baseline_reflection_ratio": (
                schemes[candidate_boundary_scheme][
                    "reflection_amplitude_ratio"
                ] / schemes["characteristic_nonreflecting"][
                    "reflection_amplitude_ratio"
                ]
            ),
            "schemes": schemes,
            "reference_probe_trace": reference.tolist(),
            "reference_runtime_s": reference_runtime,
        }
    passed = all(item["passed"] for item in orientations.values())
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "z_involving_open_edge_corner_reflection_gate_passed"
            if passed else
            "z_involving_open_edge_corner_reflection_gate_failed"
        ),
        "passed": passed,
        "promotion_allowed": False,
        "configuration": {
            "short_active_axis_cells": SHORT_SIZE,
            "reference_active_axis_cells": REFERENCE_SIZE,
            "inactive_axis_cells": INACTIVE_SIZE,
            "source_cell": SOURCE_CELL,
            "probe_cell": PROBE_CELL,
            "pulse_sigma_cells": PULSE_SIGMA_CELLS,
            "pulse_wavelength_cells": pulse_wavelength_cells,
            "candidate_boundary_scheme": candidate_boundary_scheme,
            "pulse_density_amplitude": PULSE_AMPLITUDE,
            "lattice_viscosity": LATTICE_VISCOSITY,
            "time_steps": STEPS,
        },
        "orientations": orientations,
        "scope": (
            "linear_localised_scaled_xy_control_z_edge_and_three_face_corner_"
            "component_test_"
            "not_vortical_or_plume_validation"
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_audit(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    candidate_scheme = record.get("configuration", {}).get(
        "candidate_boundary_scheme", "characteristic_multidimensional",
    )
    row_count = len(record["orientations"])
    figure, axes = plt.subplots(
        row_count, 2,
        figsize=(11.0, 2.8 * row_count),
        sharex=True,
        constrained_layout=True,
    )
    colors = {
        "characteristic_nonreflecting": "#d73027",
        candidate_scheme: "#1a9850",
    }
    labels = {
        "characteristic_nonreflecting": "1-D impedance",
        candidate_scheme: (
            "time-dependent convective/NSCBC"
            if candidate_scheme == "characteristic_nscbc_convective"
            else (
            "limited vortical/acoustic impedance"
            if candidate_scheme
            == "characteristic_multidimensional_vortical_limited"
            else (
                "extrapolated angle-aware impedance"
                if candidate_scheme
                == "characteristic_multidimensional_extrapolated"
                else "angle-aware impedance"
            )
            )
        ),
    }
    for row, (name, item) in enumerate(record["orientations"].items()):
        reference = np.asarray(item["reference_probe_trace"], dtype=float)
        steps = np.arange(reference.size)
        axes[row, 0].plot(steps, reference, color="black", label="reference")
        for scheme, values in item["schemes"].items():
            trace = np.asarray(values["probe_trace"], dtype=float)
            axes[row, 0].plot(
                steps, trace, color=colors[scheme], label=labels[scheme],
            )
            axes[row, 1].plot(
                steps, trace - reference, color=colors[scheme],
                label=(
                    f"{labels[scheme]}: "
                    f"{100.0 * values['reflection_amplitude_ratio']:.2f}%"
                ),
            )
        window = item["schemes"][candidate_scheme][
            "reflection_window_steps"
        ]
        for axis in axes[row]:
            axis.axvspan(window[0], window[1], color="#fdae61", alpha=0.15)
            axis.grid(alpha=0.25)
            axis.legend(frameon=False, fontsize=8)
        axes[row, 0].set_ylabel(f"{name}\nprobe perturbation")
        axes[row, 1].set_ylabel("short - reference")
    axes[-1, 0].set_xlabel("Lattice time step")
    axes[-1, 1].set_xlabel("Lattice time step")
    figure.suptitle("D3Q19 z-edge and three-face corner reflection", fontsize=12)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument(
        "--pulse-wavelength-cells", type=float,
        default=PULSE_WAVELENGTH_CELLS,
    )
    parser.add_argument(
        "--candidate-boundary-scheme",
        choices=(
            "characteristic_multidimensional",
            "characteristic_multidimensional_extrapolated",
            "characteristic_multidimensional_vortical_limited",
            "characteristic_nscbc_convective",
        ),
        default="characteristic_multidimensional",
    )
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    if not math.isfinite(args.pulse_wavelength_cells) or args.pulse_wavelength_cells <= 0.0:
        parser.error("--pulse-wavelength-cells must be positive and finite")
    record = run_audit(
        args.workers,
        args.pulse_wavelength_cells,
        args.candidate_boundary_scheme,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    plot_audit(record, args.plot)
    print(json.dumps({
        "status": record["status"],
        "candidate_boundary_scheme": args.candidate_boundary_scheme,
        "candidate_reflection_ratios": {
            name: item["schemes"][args.candidate_boundary_scheme][
                "reflection_amplitude_ratio"
            ]
            for name, item in record["orientations"].items()
        },
        "runtime_s": record["runtime_s"],
    }, indent=2, allow_nan=False))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
    candidate_scheme = record["configuration"]["candidate_boundary_scheme"]
