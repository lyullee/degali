"""Isolate one-way H2 scalar grid sensitivity on the active-audit flow family.

The case retains the same carrier-flow solve, source history, diffusivity,
MUSCL scheme, grid and time-step family as the active thermal audit, but turns
off enthalpy/EOS/buoyancy feedback.  Comparison with the active audit separates
scalar/observation error from thermodynamic-coupling error.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.fast_obstacle_transport_3d import FastObstacleReceptor3D  # noqa: E402
from degali.addons.resolved_obstacle_flow_3d import (  # noqa: E402
    solve_resolved_obstacle_flow_3d,
)

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
SURFACE_OFFSETS_M = (0.15, 0.30, 0.45, 0.60)
SURFACE_SWEEP_RECEPTORS = tuple(
    FastObstacleReceptor3D(f"roof_offset_{offset:.2f}", (1.80, 1.80, 1.20 + offset))
    for offset in SURFACE_OFFSETS_M
) + tuple(
    FastObstacleReceptor3D(f"side_offset_{offset:.2f}", (1.80, 1.20 - offset, 0.60))
    for offset in SURFACE_OFFSETS_M
)


def _dose(time_s: np.ndarray, trace: np.ndarray) -> float:
    return float(np.trapezoid(trace, time_s))


def _relative_change(left: float, right: float, floor: float = 0.01) -> float:
    return abs(right - left) / max(abs(left), abs(right), floor)


def _restrict_factor_two_cell_averages(field: np.ndarray) -> np.ndarray:
    """Conservatively aggregate eight aligned fine FV cells to one coarse cell."""
    nz, ny, nx = field.shape
    if nz % 2 or ny % 2 or nx % 2:
        raise ValueError("factor-two restriction requires even cell counts")
    reshaped = field.reshape(nz // 2, 2, ny // 2, 2, nx // 2, 2)
    finite = np.isfinite(reshaped)
    counts = np.sum(finite, axis=(1, 3, 5))
    sums = np.nansum(reshaped, axis=(1, 3, 5))
    restricted = np.full(counts.shape, np.nan)
    np.divide(sums, counts, out=restricted, where=counts > 0)
    return restricted


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--geometry", choices=("obstacle", "open"), default="obstacle",
    )
    parser.add_argument(
        "--grid-family", choices=("ratio1p5", "factor2"),
        default="ratio1p5",
    )
    parser.add_argument(
        "--scheme",
        choices=("upwind", "muscl_van_leer", "muscl_ssprk2"),
        default="muscl_van_leer",
    )
    parser.add_argument(
        "--effective-schmidt", type=float, default=None,
        help=(
            "replace the baseline scalar diffusivity by nu/Sc for this "
            "diagnostic; this is a closure candidate, not calibration"
        ),
    )
    arguments = parser.parse_args(argv)
    if (
        arguments.effective_schmidt is not None
        and (
            not math.isfinite(arguments.effective_schmidt)
            or arguments.effective_schmidt <= 0.0
        )
    ):
        parser.error("--effective-schmidt must be finite and positive")
    spacings = (
        SPACINGS_M if arguments.grid_family == "ratio1p5"
        else (0.30, 0.15, 0.075)
    )
    configs = []
    for spacing in spacings:
        base = _config(spacing, arguments.scheme)
        configs.append(replace(
            base,
            thermal_buoyancy_coupling=False,
            source_temperature_k=None,
            effective_schmidt_number=arguments.effective_schmidt,
            obstacles=(base.obstacles if arguments.geometry == "obstacle" else ()),
        ))
    evaluation_receptors = (
        RECEPTORS + SURFACE_SWEEP_RECEPTORS
        if arguments.geometry == "obstacle" else RECEPTORS
    )
    results = [
        solve_resolved_obstacle_flow_3d(config, receptors=evaluation_receptors)
        for config in configs
    ]
    cases = []
    for config, result, spacing in zip(configs, results, spacings):
        maximum = result.maximum_concentration_vol_percent
        cases.append({
            "cell_size_m": spacing,
            "cells": list(config.cells),
            "time_step_s": config.scalar_time_step_s,
            "runtime_s": result.diagnostics["total_runtime_s"],
            "maximum_h2_vol_percent": float(np.nanmax(maximum)),
            "four_vol_percent_volume_m3": float(
                np.count_nonzero(maximum >= 4.0) * spacing**3
            ),
            "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
            "flow_residual": result.diagnostics["final_flow_residual"],
            "receptors": {
                receptor.label: {
                    "peak_h2_vol_percent": float(np.max(
                        result.receptor_traces_vol_percent[receptor.label]
                    )),
                    "dose_h2_vol_percent_s": _dose(
                        result.time_s,
                        result.receptor_traces_vol_percent[receptor.label],
                    ),
                }
                for receptor in RECEPTORS
            },
        })

    if arguments.grid_family == "factor2":
        maximum_fields = [
            result.maximum_concentration_vol_percent for result in results
        ]
        medium_on_coarse = _restrict_factor_two_cell_averages(maximum_fields[1])
        fine_on_medium = _restrict_factor_two_cell_averages(maximum_fields[2])
        valid_coarse_medium = (
            np.isfinite(maximum_fields[0]) & np.isfinite(medium_on_coarse)
        )
        valid_medium_fine = (
            np.isfinite(maximum_fields[1]) & np.isfinite(fine_on_medium)
        )
        differences = [
            float(np.sqrt(np.mean((
                maximum_fields[0][valid_coarse_medium]
                - medium_on_coarse[valid_coarse_medium]
            )**2))),
            float(np.sqrt(np.mean((
                maximum_fields[1][valid_medium_fine]
                - fine_on_medium[valid_medium_fine]
            )**2))),
        ]
        common_point_count = {
            "coarse_medium": int(np.count_nonzero(valid_coarse_medium)),
            "medium_fine": int(np.count_nonzero(valid_medium_fine)),
        }
        field_comparison_operator = (
            "pairwise conservative factor-two restriction of finite-volume "
            "cell-average maximum fields"
        )
    else:
        coarse = results[0]
        zz, yy, xx = np.meshgrid(
            coarse.z_m, coarse.y_m, coarse.x_m, indexing="ij",
        )
        points = np.column_stack((zz.ravel(), yy.ravel(), xx.ravel()))
        sampled = []
        valid = np.ones(points.shape[0], dtype=bool)
        for result in results:
            values = RegularGridInterpolator(
                (result.z_m, result.y_m, result.x_m),
                result.maximum_concentration_vol_percent,
                bounds_error=False, fill_value=np.nan,
            )(points)
            sampled.append(values)
            valid &= np.isfinite(values)
        sampled = [values[valid] for values in sampled]
        differences = [
            float(np.sqrt(np.mean((sampled[0] - sampled[1])**2))),
            float(np.sqrt(np.mean((sampled[1] - sampled[2])**2))),
        ]
        common_point_count = int(np.count_nonzero(valid))
        field_comparison_operator = (
            "trilinear interpolation to the fixed coarse-grid centres"
        )
    observed_order = (
        None if min(differences) <= 0.0 else
        math.log(differences[0] / differences[1]) / math.log(1.5)
    )
    field_metric = {
        "coarse_medium_rms_vol_percent": differences[0],
        "medium_fine_rms_vol_percent": differences[1],
        "observed_order": observed_order,
        "successive_difference_decreased": differences[1] < differences[0],
        "comparison_operator": field_comparison_operator,
    }

    medium, fine = cases[1], cases[2]
    changes = {
        "maximum_h2": _relative_change(
            medium["maximum_h2_vol_percent"], fine["maximum_h2_vol_percent"],
        ),
        "four_vol_percent_volume": _relative_change(
            medium["four_vol_percent_volume_m3"],
            fine["four_vol_percent_volume_m3"], spacings[-1]**3,
        ),
        "receptors": {},
    }
    for receptor in RECEPTORS:
        left = medium["receptors"][receptor.label]
        right = fine["receptors"][receptor.label]
        changes["receptors"][receptor.label] = {
            "peak": _relative_change(
                left["peak_h2_vol_percent"], right["peak_h2_vol_percent"],
            ),
            "dose": _relative_change(
                left["dose_h2_vol_percent_s"],
                right["dose_h2_vol_percent_s"],
            ),
        }

    surface_offset_sweep = None
    if arguments.geometry == "obstacle":
        surface_offset_sweep = {}
        for direction in ("roof", "side"):
            surface_offset_sweep[direction] = []
            for offset in SURFACE_OFFSETS_M:
                label = f"{direction}_offset_{offset:.2f}"
                left_trace = results[1].receptor_traces_vol_percent[label]
                right_trace = results[2].receptor_traces_vol_percent[label]
                left_peak = float(np.max(left_trace))
                right_peak = float(np.max(right_trace))
                left_dose = _dose(results[1].time_s, left_trace)
                right_dose = _dose(results[2].time_s, right_trace)
                surface_offset_sweep[direction].append({
                    "offset_m": offset,
                    "medium_peak_vol_percent": left_peak,
                    "fine_peak_vol_percent": right_peak,
                    "peak_relative_change": _relative_change(
                        left_peak, right_peak,
                    ),
                    "medium_dose_vol_percent_s": left_dose,
                    "fine_dose_vol_percent_s": right_dose,
                    "dose_relative_change": _relative_change(
                        left_dose, right_dose,
                    ),
                })

    gate_reasons = []
    if not field_metric["successive_difference_decreased"]:
        gate_reasons.append("maximum-H2 field difference is non-monotonic")
    if changes["maximum_h2"] > 0.10:
        gate_reasons.append("domain maximum changes by more than 10%")
    if changes["four_vol_percent_volume"] > 0.10:
        gate_reasons.append("4 vol% volume changes by more than 10%")
    for label, metric in changes["receptors"].items():
        if metric["peak"] > 0.10:
            gate_reasons.append(f"{label} peak changes by more than 10%")
        if metric["dose"] > 0.10:
            gate_reasons.append(f"{label} dose changes by more than 10%")

    record = {
        "schema": "degali.active-grid-scalar-isolation.v1",
        "status": (
            "accepted_component_grid"
            if not gate_reasons else "withheld_component_grid"
        ),
        "scope": (
            "one-way isothermal H2 transport; enthalpy, EOS density and "
            "buoyancy feedback disabled; geometry=" + arguments.geometry
        ),
        "grid_family": arguments.grid_family,
        "spacings_m": spacings,
        "source_spatial_discretization": (
            "cell-integrated separable Gaussian over Cartesian cells"
        ),
        "scalar_advection_scheme": arguments.scheme,
        "effective_schmidt_number": arguments.effective_schmidt,
        "scalar_diffusivity_m2_s": (
            configs[0].scalar_diffusivity_m2_s
            if configs[0].effective_schmidt_number is None
            else configs[0].kinematic_viscosity_m2_s
            / configs[0].effective_schmidt_number
        ),
        "common_point_count": common_point_count,
        "cases_coarse_to_fine": cases,
        "maximum_field_metric": field_metric,
        "finest_pair_relative_changes": changes,
        "surface_offset_sweep": surface_offset_sweep,
        "gate_reasons": gate_reasons,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    geometry_suffix = "" if arguments.geometry == "obstacle" else "_open_domain"
    scheme_suffix = {
        "muscl_van_leer": "",
        "upwind": "_upwind",
        "muscl_ssprk2": "_ssprk2",
    }[arguments.scheme]
    schmidt_suffix = (
        "" if arguments.effective_schmidt is None
        else f"_sc{arguments.effective_schmidt:g}".replace(".", "p")
    )
    grid_suffix = "" if arguments.grid_family == "ratio1p5" else "_factor2"
    suffix = geometry_suffix + scheme_suffix + schmidt_suffix + grid_suffix
    json_path = OUTPUT / f"active_grid_scalar_isolation{suffix}.json"
    figure_path = OUTPUT / f"active_grid_scalar_isolation{suffix}.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    common_maximum = max(
        float(np.nanmax(result.maximum_concentration_vol_percent))
        for result in results
    )
    for axis, result, case, label in zip(
        axes, results, cases, ("coarse", "medium", "fine"),
    ):
        center = int(np.argmin(np.abs(result.y_m - 1.80)))
        field = result.maximum_concentration_vol_percent[:, center, :]
        image = axis.pcolormesh(
            result.x_m, result.z_m, field, shading="nearest", cmap="magma",
            vmin=0.0, vmax=common_maximum,
        )
        if np.nanmax(field) >= 4.0:
            axis.contour(
                result.x_m, result.z_m, field, levels=[4.0], colors="white",
                linewidths=1.1,
            )
        axis.contourf(
            result.x_m, result.z_m,
            result.solid_mask[:, center, :].astype(float),
            levels=[0.5, 1.5], colors=["0.25"],
        )
        axis.set_title(
            f"{label}: max={case['maximum_h2_vol_percent']:.2f} vol%"
        )
        axis.set_xlabel("x (m)")
        axis.set_ylabel("z (m)")
        figure.colorbar(image, ax=axis, label="maximum H2 (vol%)")
    figure.suptitle(
        "One-way isothermal scalar grid isolation ("
        + arguments.geometry + ") — " + record["status"]
    )
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "cases": cases,
        "maximum_field_metric": field_metric,
        "finest_pair_relative_changes": changes,
        "surface_offset_sweep": surface_offset_sweep,
        "gate_reasons": gate_reasons,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
