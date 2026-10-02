"""Run the pre-registered FFI Test 6 finite wall-jet mechanism diagnostic."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from degali.addons import FiniteWallJetCrosswind, WallJetTransitionConfig
from degali.lh2 import (
    audit_lh2_independent_energy_interface,
    run_lh2_crosswind_research,
)
from degali.validation.nearfield import IndependentEnergyTrajectory, hydrogen_jet
from degali.validation.spadeadam import load


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "reference" / "spadeadam"


def _source_and_near_field():
    setup, _initial = hydrogen_jet(
        rate=0.833,
        diameter=0.0254,
        wind=2.3,
        height=0.5,
        ambient_temperature=277.15,
        relative_humidity=90.0,
        ambient_pressure=101325.0,
        storage_pressure_barg=2.53,
        roughness=0.001,
        stability="D",
        averaging=60.0,
        wind_reference_height=10.0,
        corrections=True,
        condensed_air_particle_diameter=1.0e-6,
        source_total_energy_consistency=True,
        source_pressure_thrust=False,
    )
    source = setup.axisymmetric_source
    coupled = run_lh2_crosswind_research(
        source,
        wind=2.3,
        height=0.5,
        ambient_temperature=277.15,
        ambient_pressure=101325.0,
        relative_humidity=90.0,
        roughness=0.001,
        stability="D",
        averaging=60.0,
        wind_reference_height=10.0,
        maximum_nearfield_distance=3.0,
        handoff_distance=3.0,
        minimum_mass_fraction=7.0e-4,
        radial_points=41,
        maximum_step=0.01,
        relative_tolerance=2.0e-6,
        thermodynamic_closure="phase_manifold",
        nearfield_establishment="entrained_mass",
        nearfield_energy_transport="total",
        crosswind_entrainment="source_momentum",
    )
    return setup, source, coupled


def _interface(coupled, ground_interaction):
    return audit_lh2_independent_energy_interface(
        coupled.near_field,
        coupled.handoff.model,
        streamline_distance=3.0,
        energy_transport="total",
        houf_width_mapping="velocity",
        ground_interaction=ground_interaction,
    )


def _arc_rows(trajectory, trial, wall=None):
    rows = []
    for radius in (30.0, 50.0, 100.0):
        observed = trial.arc(radius)
        state = trajectory.state_at(radius)
        if state is None:
            rows.append({
                "radius_m": radius,
                "status": "not_reached_by_conservative_integration",
                "observed_peak_max_pct": float(max(observed.values())),
                "observed_by_height_pct": {
                    str(height): float(value)
                    for height, value in observed.items()
                },
            })
            continue
        section = state[:7]
        _sy, sigma_z = trajectory.model.section_widths(state)
        predicted = {
            float(height): float(trajectory.concentration_at(radius, 0.0, height))
            for height in observed
        }
        row = {
            "radius_m": radius,
            "observed_peak_max_pct": float(max(observed.values())),
            "observed_by_height_pct": {
                str(height): float(value) for height, value in observed.items()
            },
            "predicted_max_pct": float(max(predicted.values())),
            "predicted_by_height_pct": {
                str(height): value for height, value in predicted.items()
            },
            "centre_height_m": float(section[6]),
            "sigma_z_m": float(sigma_z),
        }
        if wall is not None:
            diagnostic = wall.diagnostics(state)
            row["attachment_fraction"] = diagnostic.attachment_fraction
            row["geometric_contact_fraction"] = (
                diagnostic.geometric_contact_fraction
            )
            row["lift_richardson"] = diagnostic.lift_richardson
            row["vertical_release_fraction"] = (
                diagnostic.vertical_release_fraction
            )
        rows.append(row)
    return rows


def _run_base(interface, *, distance, step):
    if not interface.accepted:
        raise RuntimeError("independent-energy interface failed: " + "; ".join(
            interface.failure_reasons
        ))
    result = interface.model.solve(
        interface.state.copy(),
        maximum_distance=distance,
        maximum_step=step,
        relative_tolerance=2.0e-6,
        method="flux_RK4",
    )
    return result, IndependentEnergyTrajectory(interface, result)


def _run_wall(interface, *, distance, step):
    if not interface.accepted:
        raise RuntimeError("independent-energy interface failed: " + "; ".join(
            interface.failure_reasons
        ))
    wall = FiniteWallJetCrosswind(
        interface.model,
        WallJetTransitionConfig(
            critical_richardson=30.0,
            response_depths=1.0,
            wall_shear_multiplier=1.0,
        ),
    )
    initial = wall.initial_state(interface.state)
    result = wall.solve(
        initial,
        maximum_distance=distance,
        maximum_step=step,
        relative_tolerance=2.0e-6,
        method="flux_RK4",
    )
    view = SimpleNamespace(model=wall)
    return wall, result, IndependentEnergyTrajectory(view, result)


def run(*, step=0.1, refined_step=0.05, distance=210.0):
    setup, source, coupled = _source_and_near_field()
    free_interface = _interface(coupled, "free")
    ground_interface = _interface(coupled, "geometry")
    free_result, free_trajectory = _run_base(
        free_interface, distance=distance, step=step
    )
    ground_result, ground_trajectory = _run_base(
        ground_interface, distance=distance, step=step
    )
    wall, wall_result, wall_trajectory = _run_wall(
        ground_interface, distance=distance, step=step
    )
    refined_wall, refined_result, refined_trajectory = _run_wall(
        ground_interface, distance=distance, step=refined_step
    )
    trial = next(item for item in load(REFERENCE) if item.test == 6)

    interface = {
        "accepted": bool(ground_interface.accepted),
        "maximum_flux_residual": float(max(
            ground_interface.relative_residuals.values()
        )),
        "energy_quadrature_residual": float(
            ground_interface.energy_quadrature_residual
        ),
        "halfwidth_residual": float(ground_interface.halfwidth_residual),
        "temperature_residual_k": float(
            ground_interface.temperature_residual
        ),
        "failure_reasons": list(ground_interface.failure_reasons),
    }
    payload = {
        "case": "FFI Spadeadam Test 6 finite wall-jet diagnostic",
        "pre_registration": "docs/prereg-ffi-test6-finite-wall-jet.md",
        "source": {
            "temperature_k": float(source.temperature),
            "velocity_ms": float(source.velocity),
            "diameter_m": float(source.diameter),
            "density_kgm3": float(source.density),
            "hydrogen_mass_fraction": float(source.mass_fraction),
            "source_zone_x_m": float(source.x),
            "local_wind_ms": float(setup.local_source_wind),
            "velocity_to_wind_ratio": float(
                source.velocity / setup.local_source_wind
            ),
        },
        "interface": interface,
        "integration": {
            "distance_m": float(distance),
            "coarse_maximum_step_m": float(step),
            "refined_maximum_step_m": float(refined_step),
        },
        "free": {
            "maximum_balance_residual": float(
                free_result.maximum_relative_balance_residual
            ),
            "arcs": _arc_rows(free_trajectory, trial),
        },
        "permanent_ground_geometry": {
            "maximum_balance_residual": float(
                ground_result.maximum_relative_balance_residual
            ),
            "arcs": _arc_rows(ground_trajectory, trial),
        },
        "finite_wall_jet": {
            "configuration": {
                "critical_richardson": wall.config.critical_richardson,
                "response_depths": wall.config.response_depths,
                "wall_shear_multiplier": wall.config.wall_shear_multiplier,
            },
            "maximum_balance_residual": float(
                wall_result.maximum_relative_balance_residual
            ),
            "arcs": _arc_rows(wall_trajectory, trial, wall),
        },
        "finite_wall_jet_refined": {
            "maximum_balance_residual": float(
                refined_result.maximum_relative_balance_residual
            ),
            "arcs": _arc_rows(refined_trajectory, trial, refined_wall),
        },
    }
    coarse = payload["finite_wall_jet"]["arcs"]
    refined = payload["finite_wall_jet_refined"]["arcs"]
    pairs = [
        (a, b) for a, b in zip(coarse, refined)
        if "predicted_max_pct" in a and "predicted_max_pct" in b
    ]
    payload["refinement"] = {
        "maximum_relative_concentration_change": float(max(
            abs(a["predicted_max_pct"] - b["predicted_max_pct"])
            / max(abs(b["predicted_max_pct"]), 1.0e-12)
            for a, b in pairs
        )),
        "maximum_centre_height_change_m": float(max(
            abs(a["centre_height_m"] - b["centre_height_m"])
            for a, b in pairs
        )),
    }
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--step", type=float, default=0.1)
    parser.add_argument("--refined-step", type=float, default=0.05)
    parser.add_argument("--distance", type=float, default=210.0)
    args = parser.parse_args()
    payload = run(
        step=args.step,
        refined_step=args.refined_step,
        distance=args.distance,
    )
    rendered = json.dumps(payload, indent=2, ensure_ascii=False)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
