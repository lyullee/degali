"""Manufactured checks for the coupled species/enthalpy/EOS/mass-flux step.

The record accepts only the discrete coupling operator.  It is deliberately
separate from external obstacle-flow validation and cannot promote Tier 3 to
quantitative facility use.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.hydrogen_air_thermodynamics import (
    HydrogenAirActiveScalarClosure,
)
from degali.addons.low_mach_thermal_corrector_3d import (
    solve_low_mach_thermal_corrector_3d,
)


OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"


def _faces(
    shape: tuple[int, int, int], solid: np.ndarray, speed_m_s: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nz, ny, nx = shape
    u = np.full((nz, ny, nx + 1), speed_m_s)
    v = np.zeros((nz, ny + 1, nx))
    w = np.zeros((nz + 1, ny, nx))
    u[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    v[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    w[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    return u, v, w


def _case_record(name: str, result, runtime_s: float) -> dict[str, object]:
    final_change = result.iteration_history[-1]
    projection = result.mass_projection
    return {
        "name": name,
        "converged": result.converged,
        "iterations": result.iterations,
        "linear_solver": result.mass_projection.solver_method,
        "runtime_s": runtime_s,
        "final_relative_density_change": float(final_change[0]),
        "final_relative_velocity_change": float(final_change[1]),
        "maximum_local_continuity_residual_kg_m3_s": (
            projection.maximum_absolute_residual_after_kg_m3_s
        ),
        "global_continuity_residual_kg_s": (
            projection.global_mass_balance_residual_kg_s
        ),
        "hydrogen_mass_residual_kg": result.hydrogen_mass_residual_kg,
        "sensible_enthalpy_residual_j": result.sensible_enthalpy_residual_j,
        "maximum_absolute_pressure_correction_pa": float(np.max(np.abs(
            result.pressure_correction_pa
        ))),
        "minimum_temperature_k": float(np.min(
            result.thermodynamic_state.temperature_k
        )),
        "maximum_relative_density_departure": float(np.max(np.abs(
            result.thermodynamic_state.relative_density_departure
        ))),
    }


def main() -> int:
    spacing = (0.2, 0.2, 0.2)
    dt = 0.01
    closure = HydrogenAirActiveScalarClosure(
        maximum_relative_density_departure=0.10,
    )

    ambient_shape = (6, 7, 10)
    ambient_solid = np.zeros(ambient_shape, dtype=bool)
    ambient_zero = np.zeros(ambient_shape)
    ambient_faces = _faces(ambient_shape, ambient_solid, 0.25)
    started = time.perf_counter()
    ambient = solve_low_mach_thermal_corrector_3d(
        ambient_zero, ambient_zero,
        hydrogen_source_increment_kg_m3=ambient_zero,
        sensible_enthalpy_source_increment_j_m3=ambient_zero,
        provisional_u_face_m_s=ambient_faces[0],
        provisional_v_face_m_s=ambient_faces[1],
        provisional_w_face_m_s=ambient_faces[2],
        solid_mask=ambient_solid,
        spacing_m=spacing,
        time_step_s=dt,
        species_diffusivity_m2_s=0.01,
        thermal_diffusivity_m2_s=0.01,
        closure=closure,
        advection_scheme="muscl_van_leer",
    )
    ambient_runtime = time.perf_counter() - started

    shape = (7, 8, 12)
    solid = np.zeros(shape, dtype=bool)
    solid[:5, 3:6, 6:8] = True
    zero = np.zeros(shape)
    source = np.zeros(shape)
    source[1:3, 3:5, 2:4] = 2.0e-4
    source_energy = source * closure.source_specific_sensible_enthalpy_j_kg(80.0)
    faces = _faces(shape, solid, 0.25)
    started = time.perf_counter()
    cold = solve_low_mach_thermal_corrector_3d(
        zero, zero,
        hydrogen_source_increment_kg_m3=source,
        sensible_enthalpy_source_increment_j_m3=source_energy,
        provisional_u_face_m_s=faces[0],
        provisional_v_face_m_s=faces[1],
        provisional_w_face_m_s=faces[2],
        solid_mask=solid,
        spacing_m=spacing,
        time_step_s=dt,
        species_diffusivity_m2_s=0.01,
        thermal_diffusivity_m2_s=0.01,
        closure=closure,
        advection_scheme="muscl_van_leer",
        maximum_iterations=12,
        relative_density_tolerance=1.0e-7,
        relative_velocity_tolerance=1.0e-7,
    )
    cold_runtime = time.perf_counter() - started

    cases = [
        _case_record("ambient_uniform_flow", ambient, ambient_runtime),
        _case_record("cold_hydrogen_source_with_block", cold, cold_runtime),
    ]
    thresholds = {
        "maximum_local_continuity_residual_kg_m3_s": 1.0e-8,
        "maximum_global_continuity_residual_kg_s": 1.0e-9,
        "maximum_hydrogen_mass_residual_kg": 1.0e-12,
        "maximum_sensible_enthalpy_residual_j": 1.0e-8,
    }
    accepted = bool(
        all(case["converged"] for case in cases)
        and max(abs(float(case["maximum_local_continuity_residual_kg_m3_s"]))
                for case in cases)
        <= thresholds["maximum_local_continuity_residual_kg_m3_s"]
        and max(abs(float(case["global_continuity_residual_kg_s"]))
                for case in cases)
        <= thresholds["maximum_global_continuity_residual_kg_s"]
        and max(abs(float(case["hydrogen_mass_residual_kg"])) for case in cases)
        <= thresholds["maximum_hydrogen_mass_residual_kg"]
        and max(abs(float(case["sensible_enthalpy_residual_j"])) for case in cases)
        <= thresholds["maximum_sensible_enthalpy_residual_j"]
    )
    record = {
        "schema": "degali.low-mach-thermal-corrector-verification.v1",
        "status": "accepted_numerical_operator_only" if accepted else "failed",
        "scope": (
            "discrete species/enthalpy/EOS/mass-flux fixed point; not physical "
            "validation of obstacle flow or cryogenic source thermodynamics"
        ),
        "thresholds": thresholds,
        "cases": cases,
        "solid_face_penetration_m_s": float(max(
            np.max(np.abs(cold.u_face_m_s[:, :, 1:-1][
                solid[:, :, :-1] | solid[:, :, 1:]
            ])),
            np.max(np.abs(cold.v_face_m_s[:, 1:-1, :][
                solid[:, :-1, :] | solid[:, 1:, :]
            ])),
            np.max(np.abs(cold.w_face_m_s[1:-1, :, :][
                solid[:-1, :, :] | solid[1:, :, :]
            ])),
        )),
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "low_mach_thermal_corrector_verification.json"
    figure_path = OUTPUT / "low_mach_thermal_corrector_verification.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    center = 4
    extent = (0.0, shape[2] * spacing[0], 0.0, shape[0] * spacing[2])
    state = cold.thermodynamic_state
    fields = (
        (100.0 * state.hydrogen_mole_fraction[:, center, :], "H2 (vol%)"),
        (state.temperature_k[:, center, :], "temperature (K)"),
        (state.mixture_density_kg_m3[:, center, :], "density (kg/m3)"),
        (cold.pressure_correction_pa[:, center, :], "pressure correction (Pa)"),
        (cold.w_face_m_s[:-1, center, :], "vertical face velocity (m/s)"),
        (np.abs(cold.mass_projection.continuity_residual_after_kg_m3_s[
            :, center, :
        ]), "|continuity residual| (kg/m3/s)"),
    )
    figure, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for axis, (field, label) in zip(axes.flat, fields):
        shown = np.asarray(field, dtype=float).copy()
        shown[solid[:, center, :]] = np.nan
        if "residual" in label:
            shown = np.maximum(shown, 1.0e-16)
            image = axis.imshow(
                shown, origin="lower", extent=extent, aspect="auto",
                cmap="viridis", norm="log",
            )
        else:
            image = axis.imshow(
                shown, origin="lower", extent=extent, aspect="auto",
                cmap="coolwarm",
            )
        axis.contour(
            np.linspace(spacing[0] / 2.0, extent[1] - spacing[0] / 2.0, shape[2]),
            np.linspace(spacing[2] / 2.0, extent[3] - spacing[2] / 2.0, shape[0]),
            solid[:, center, :].astype(float), levels=[0.5], colors="black",
        )
        axis.set_title(label)
        axis.set_xlabel("x (m)")
        axis.set_ylabel("z (m)")
        figure.colorbar(image, ax=axis)
    figure.suptitle(
        "Low-Mach thermal corrector manufactured case — numerical operator only"
    )
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "cases": cases,
        "solid_face_penetration_m_s": record["solid_face_penetration_m_s"],
    }, indent=2))
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
