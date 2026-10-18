"""Isolate time refinement of the scalar/EOS/mass-flux corrector.

The fully coupled D3Q19 audit changes the lattice relaxation time when the
physical time step is changed on a fixed grid.  This diagnostic instead uses
one physical carrier-flow field for every level, so that the remaining
differences come from species/enthalpy transport, the ideal-mixture closure,
and the variable-density mass-flux projection.
"""

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

from degali.addons.hydrogen_air_thermodynamics import (  # noqa: E402
    HydrogenAirActiveScalarClosure,
)
from degali.addons.low_mach_thermal_corrector_3d import (  # noqa: E402
    solve_low_mach_thermal_corrector_3d,
)
from degali.addons.resolved_obstacle_flow_3d import (  # noqa: E402
    _solve_flow,
)
from degali.addons.fast_obstacle_transport_3d import (  # noqa: E402
    _source_weights,
)
from audit_active_thermal_corrector_time import (  # noqa: E402
    DURATION_S,
    LATTICE_SPEEDS,
    OUTPUT,
    WIND_M_S,
    _config,
)


EFFECTIVE_SCHMIDT = 0.70
EFFECTIVE_PRANDTL = 0.85


def _rms(left: np.ndarray, right: np.ndarray, fluid: np.ndarray) -> float:
    return float(np.sqrt(np.mean((left[fluid] - right[fluid]) ** 2)))


def _metric(fields: list[np.ndarray], fluid: np.ndarray) -> dict[str, object]:
    coarse_medium = _rms(fields[0], fields[1], fluid)
    medium_fine = _rms(fields[1], fields[2], fluid)
    return {
        "coarse_medium_rms": coarse_medium,
        "medium_fine_rms": medium_fine,
        "observed_order": (
            None
            if min(coarse_medium, medium_fine) <= 0.0
            else math.log(coarse_medium / medium_fine) / math.log(2.0)
        ),
        "successive_difference_decreased": medium_fine < coarse_medium,
    }


def main() -> int:
    carrier_config = _config(0.015)
    (
        x, y, z, solid, _u, _v, _w, u_face, v_face, w_face,
        _pressure, _pressure_coefficient, _eddy_viscosity,
        _residual_history, _density, _flow_dt, _base_tau, _maximum_tau,
        _actual_steps, _raw_divergence, _divergence_before,
        _divergence_after, _projection_iterations, _distribution,
        _statistical_residual, _average_count,
    ) = _solve_flow(carrier_config)
    dx = float(x[1] - x[0])
    cell_volume = dx**3
    weights = _source_weights(
        x, y, z, solid, carrier_config.source_position_m,
        carrier_config.source_sigma_m,
    )
    closure = HydrogenAirActiveScalarClosure(
        ambient_temperature_k=carrier_config.ambient_temperature_k,
        ambient_pressure_pa=carrier_config.ambient_pressure_pa,
        ambient_air_density_kg_m3=carrier_config.ambient_air_density_kg_m3,
        dry_air_heat_capacity_j_kg_k=(
            carrier_config.dry_air_heat_capacity_j_kg_k
        ),
        hydrogen_heat_capacity_j_kg_k=(
            carrier_config.hydrogen_heat_capacity_j_kg_k
        ),
        gravity_m_s2=carrier_config.gravity_m_s2,
        maximum_relative_density_departure=(
            carrier_config.maximum_active_scalar_density_departure
        ),
    )
    source_specific_enthalpy = closure.source_specific_sensible_enthalpy_j_kg(
        float(carrier_config.source_temperature_k)
    )
    species_diffusivity = (
        carrier_config.kinematic_viscosity_m2_s / EFFECTIVE_SCHMIDT
    )
    thermal_diffusivity = (
        carrier_config.kinematic_viscosity_m2_s / EFFECTIVE_PRANDTL
    )
    final_h2_fields: list[np.ndarray] = []
    maximum_h2_fields: list[np.ndarray] = []
    final_temperature_fields: list[np.ndarray] = []
    cases: list[dict[str, object]] = []
    for lattice_speed in LATTICE_SPEEDS:
        dt = lattice_speed * dx / WIND_M_S
        steps = round(DURATION_S / dt)
        concentration = np.zeros(solid.shape, dtype=float)
        enthalpy = np.zeros(solid.shape, dtype=float)
        maximum_mole_fraction = np.zeros(solid.shape, dtype=float)
        maximum_mass_residual = 0.0
        maximum_energy_residual = 0.0
        nonconverged_steps = 0
        started = time.perf_counter()
        for step in range(steps):
            start = step * dt
            end = (step + 1) * dt
            source_mass = carrier_config.source_schedule.mass_between(start, end)
            source_h2 = source_mass * weights / cell_volume
            result = solve_low_mach_thermal_corrector_3d(
                concentration,
                enthalpy,
                hydrogen_source_increment_kg_m3=source_h2,
                sensible_enthalpy_source_increment_j_m3=(
                    source_h2 * source_specific_enthalpy
                ),
                provisional_u_face_m_s=u_face,
                provisional_v_face_m_s=v_face,
                provisional_w_face_m_s=w_face,
                solid_mask=solid,
                spacing_m=(dx, dx, dx),
                time_step_s=dt,
                species_diffusivity_m2_s=species_diffusivity,
                thermal_diffusivity_m2_s=thermal_diffusivity,
                closure=closure,
                advection_scheme="muscl_van_leer",
                source_time_integration="prestep",
                maximum_iterations=8,
                minimum_iterations=2,
                relative_density_tolerance=1.0e-6,
                relative_velocity_tolerance=1.0e-6,
            )
            concentration = result.hydrogen_mass_concentration_kg_m3
            enthalpy = result.sensible_enthalpy_density_j_m3
            np.maximum(
                maximum_mole_fraction,
                result.thermodynamic_state.hydrogen_mole_fraction,
                out=maximum_mole_fraction,
            )
            maximum_mass_residual = max(
                maximum_mass_residual, abs(result.hydrogen_mass_residual_kg)
            )
            maximum_energy_residual = max(
                maximum_energy_residual,
                abs(result.sensible_enthalpy_residual_j),
            )
            nonconverged_steps += int(not result.converged)
        state = closure.state_from_conserved(concentration, enthalpy)
        runtime = time.perf_counter() - started
        final_h2 = 100.0 * state.hydrogen_mole_fraction
        maximum_h2 = 100.0 * maximum_mole_fraction
        final_h2[solid] = np.nan
        maximum_h2[solid] = np.nan
        temperature = state.temperature_k.copy()
        temperature[solid] = np.nan
        final_h2_fields.append(final_h2)
        maximum_h2_fields.append(maximum_h2)
        final_temperature_fields.append(temperature)
        cases.append({
            "time_step_s": dt,
            "steps": steps,
            "runtime_s": runtime,
            "maximum_h2_vol_percent": float(np.nanmax(maximum_h2)),
            "minimum_temperature_k": float(np.nanmin(temperature)),
            "maximum_mass_residual_kg": maximum_mass_residual,
            "maximum_energy_residual_j": maximum_energy_residual,
            "corrector_nonconverged_steps": nonconverged_steps,
        })

    fluid = ~solid
    metrics = {
        "maximum_h2_vol_percent": _metric(maximum_h2_fields, fluid),
        "final_h2_vol_percent": _metric(final_h2_fields, fluid),
        "final_temperature_k": _metric(final_temperature_fields, fluid),
    }
    gate_reasons = [
        f"{name} successive difference is non-monotonic"
        for name, metric in metrics.items()
        if not metric["successive_difference_decreased"]
    ]
    if any(case["corrector_nonconverged_steps"] for case in cases):
        gate_reasons.append("one or more fixed-point corrector steps did not converge")
    record = {
        "schema": "degali.active-time-fixed-carrier-diagnostic.v1",
        "status": (
            "accepted_component_time" if not gate_reasons
            else "withheld_component_time"
        ),
        "scope": (
            "fixed physical carrier flow; species/enthalpy/EOS/variable-density "
            "mass-flux corrector only; no buoyancy-to-momentum feedback"
        ),
        "basis": {
            "carrier_lattice_inlet_speed": 0.015,
            "effective_schmidt_number": EFFECTIVE_SCHMIDT,
            "effective_prandtl_number": EFFECTIVE_PRANDTL,
            "source_time_integration": "prestep",
        },
        "cases_coarse_to_fine": cases,
        "field_metrics": metrics,
        "gate_reasons": gate_reasons,
        "interpretation": (
            "Acceptance isolates the conservative scalar, thermodynamic, and "
            "mass-flux-correction time discretisation. It does not validate "
            "the fully coupled D3Q19 buoyancy feedback."
        ),
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "active_time_fixed_carrier_diagnostic.json"
    figure_path = OUTPUT / "active_time_fixed_carrier_diagnostic.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    center = int(np.argmin(np.abs(y - 1.20)))
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for axis, field, case, label in zip(
        axes, maximum_h2_fields, cases, ("coarse", "medium", "fine")
    ):
        image = axis.pcolormesh(
            x, z, field[:, center, :], shading="nearest", cmap="magma",
            vmin=0.0, vmax=max(np.nanmax(item) for item in maximum_h2_fields),
        )
        axis.contourf(
            x, z, solid[:, center, :].astype(float),
            levels=[0.5, 1.5], colors=["0.25"],
        )
        axis.set_title(f"{label}: dt={case['time_step_s']:.3f} s")
        axis.set_xlabel("x (m)")
        axis.set_ylabel("z (m)")
        figure.colorbar(image, ax=axis, label="maximum H2 (vol%)")
    figure.suptitle("Fixed-carrier active-corrector time isolation — " + record["status"])
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
