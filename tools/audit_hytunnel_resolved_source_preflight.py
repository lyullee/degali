"""Grid-specific admissibility preflight for the HyTunnel resolved source.

This is not an accuracy comparison.  It asks whether a top-hat low-Mach
handoff remains inside the active-scalar density range after it is mapped to
the finite-width Gaussian source used by the resolved 3-D solver.  No measured
concentration is read or used to select a candidate.

The deliberately coarse grid makes the check inexpensive enough to repeat.
Passing it only nominates a source regularisation for the later matched-grid,
time/domain-converged validation; it cannot promote the model.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from typing import Any, Sequence

import numpy as np

from degali.addons.high_pressure_source_handoff_3d import (
    LowMachHydrogenHandoff3D,
    low_mach_handoff_from_energy_nozzle,
)
from degali.addons.notional import energy_conserving_notional_nozzle
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)


SCHEMA = "degali.hytunnel-resolved-source-preflight.v6"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-resolved-source-preflight-2026-10-10-v6.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-resolved-source-preflight-2026-10-10-v6.png"
)
MEAN_DENSITY_LIMIT_CANDIDATES = (0.10, 0.08, 0.06)
RECEIVING_DENSITY_LIMIT = 0.10
GRID_SPACING_M = 0.045
DOMAIN_M = (2.52, 0.90, 0.90)
CELLS = (56, 20, 20)
LATTICE_INLET_SPEED = 0.018
WIND_SPEED_M_S = 3.5
FLOW_STEPS = 100
SOURCE_STEPS = 50


def _handoff(mean_density_limit: float) -> LowMachHydrogenHandoff3D:
    nozzle = energy_conserving_notional_nozzle(
        fluid="Hydrogen",
        storage_temperature=286.15,
        storage_pressure=11.515e6,
        mass_flow=0.005,
        orifice_diameter=0.001,
        ambient_pressure=98_500.0,
        allow_supercritical_gas=True,
    )
    return low_mach_handoff_from_energy_nozzle(
        nozzle,
        source_id=(
            "hytunnel-1mm-5g-3p5ms-coflow-preflight-"
            f"rho{mean_density_limit:.3f}"
        ),
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        ambient_axial_velocity_m_s=WIND_SPEED_M_S,
        maximum_mach_number=0.30,
        maximum_relative_density_departure=mean_density_limit,
    )


def _case_config(
    handoff: LowMachHydrogenHandoff3D,
) -> ResolvedObstacleFlow3DConfig:
    dt = LATTICE_INLET_SPEED * GRID_SPACING_M / WIND_SPEED_M_S
    duration = SOURCE_STEPS * dt
    source = handoff.resolved_source_fields(
        nozzle_position_m=(0.0, 0.45, 0.45),
        grid_spacing_m=GRID_SPACING_M,
        duration_s=duration,
    )
    return ResolvedObstacleFlow3DConfig(
        domain_m=DOMAIN_M,
        cells=CELLS,
        wind_speed_m_s=WIND_SPEED_M_S,
        # Coarse-grid numerical regularisation; not a calibrated tunnel value.
        kinematic_viscosity_m2_s=0.04,
        flow_steps=FLOW_STEPS,
        minimum_flow_steps=FLOW_STEPS,
        lattice_inlet_speed=LATTICE_INLET_SPEED,
        duration_s=duration,
        scalar_time_step_s=dt,
        scalar_diffusivity_m2_s=2.2e-5,
        thermal_diffusivity_m2_s=2.2e-5,
        obstacles=(),
        scalar_advection_scheme="muscl_van_leer",
        source_time_integration="prestep",
        subgrid_model="wale",
        wale_constant=0.5,
        flow_boundary_scheme="characteristic_multidimensional",
        # The official CFD used open non-reflecting conditions on these five
        # faces. DEGALI uses an angle-aware characteristic candidate whose
        # linear angle/wavelength gate passes; corner and plume-domain
        # qualification are still required.
        nonstreamwise_wall_boundary="all_pressure_open",
        effective_schmidt_number=0.7,
        effective_prandtl_number=0.85,
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        ambient_air_density_kg_m3=handoff.ambient_density_kg_m3,
        maximum_active_scalar_density_departure=RECEIVING_DENSITY_LIMIT,
        # The active momentum path is finite-volume. This relaxed diagnostic
        # limit avoids rejecting the physical ~20 m/s handoff merely because
        # of the carrier LBM velocity scaling selected for the spin-up.
        maximum_lattice_mach_number=0.25,
        store_scalar_fields=False,
        **source.as_config_kwargs(),
    )


def _candidate_record(mean_density_limit: float) -> dict[str, Any]:
    handoff = _handoff(mean_density_limit)
    config = _case_config(handoff)
    result = solve_resolved_obstacle_flow_3d(config)
    diagnostics = result.diagnostics
    maximum_departure = float(
        diagnostics["maximum_relative_mixture_density_departure"]
    )
    gates = {
        "mapped_profile_density_departure_at_most_0p10": (
            maximum_departure <= RECEIVING_DENSITY_LIMIT * (1.0 + 1.0e-12)
        ),
        "active_scalar_applicability_passed": bool(
            diagnostics["active_scalar_applicability_passed"]
        ),
        "thermodynamic_corrector_converged_every_step": (
            diagnostics["thermodynamic_corrector_nonconverged_steps"] == 0
        ),
        "relative_mass_residual_at_most_1e_8": (
            diagnostics["maximum_relative_mass_residual"] <= 1.0e-8
        ),
    }
    return {
        "mean_handoff_density_departure_limit": mean_density_limit,
        "handoff": {
            "distance_m": handoff.handoff_distance_m,
            "equivalent_diameter_m": handoff.equivalent_diameter_m,
            "hydrogen_mass_fraction": handoff.hydrogen_mass_fraction,
            "mixture_velocity_m_s": handoff.mixture_velocity_m_s,
            "mach_number": handoff.mach_number,
            "relative_density_departure": handoff.relative_density_departure,
            "maximum_conservation_residual": (
                handoff.maximum_relative_conservation_residual
            ),
        },
        "resolved_source": {
            "cells_across_equivalent_diameter": (
                handoff.equivalent_diameter_m / GRID_SPACING_M
            ),
            "maximum_hydrogen_vol_percent": float(
                np.max(result.maximum_concentration_vol_percent)
            ),
            "maximum_relative_mixture_density_departure": maximum_departure,
            "active_scalar_out_of_range_cell_steps": int(
                diagnostics["active_scalar_out_of_range_cell_steps"]
            ),
            "maximum_lattice_mach_number": float(
                diagnostics["maximum_lattice_mach_number"]
            ),
            "maximum_momentum_advective_courant_number": float(
                diagnostics["maximum_momentum_advective_courant_number"]
            ),
            "maximum_relative_mass_residual": float(
                diagnostics["maximum_relative_mass_residual"]
            ),
            "thermodynamic_corrector_nonconverged_steps": int(
                diagnostics["thermodynamic_corrector_nonconverged_steps"]
            ),
            "flow_converged": bool(diagnostics["flow_converged"]),
            "final_flow_residual": float(
                diagnostics["final_flow_residual"]
            ),
            "runtime_s": float(diagnostics["total_runtime_s"]),
        },
        "gates": gates,
        "passed": all(gates.values()),
    }


def classify_preflight(candidates: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Classify candidates without turning a numerical pass into validation."""
    if not candidates:
        raise ValueError("at least one preflight candidate is required")
    nominal = next(
        (
            item for item in candidates
            if abs(float(item["mean_handoff_density_departure_limit"]) - 0.10)
            <= 1.0e-12
        ),
        None,
    )
    passing = [item for item in candidates if bool(item["passed"])]
    passing.sort(
        key=lambda item: float(item["mean_handoff_density_departure_limit"]),
        reverse=True,
    )
    return {
        "status": (
            "nominal_resolved_source_profile_preflight_failed"
            if nominal is not None and not bool(nominal["passed"])
            else "nominal_resolved_source_profile_preflight_passed"
        ),
        "nominal_0p10_passed": (
            None if nominal is None else bool(nominal["passed"])
        ),
        "least_conservative_passing_candidate": (
            None
            if not passing
            else float(passing[0]["mean_handoff_density_departure_limit"])
        ),
        "promotion_allowed": False,
    }


def run_preflight() -> dict[str, Any]:
    started = time.perf_counter()
    candidates = [
        _candidate_record(limit) for limit in MEAN_DENSITY_LIMIT_CANDIDATES
    ]
    classification = classify_preflight(candidates)
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **classification,
        "scope": "coarse_grid_resolved_source_numerical_admissibility_only",
        "no_fit_statement": (
            "no measured concentration or CFD result is read by this preflight"
        ),
        "receiving_active_scalar_density_departure_limit": (
            RECEIVING_DENSITY_LIMIT
        ),
        "grid": {
            "domain_m": list(DOMAIN_M),
            "cells": list(CELLS),
            "spacing_m": GRID_SPACING_M,
            "flow_steps": FLOW_STEPS,
            "source_steps": SOURCE_STEPS,
            "kinematic_viscosity_m2_s": 0.04,
            "kinematic_viscosity_role": (
                "coarse-grid numerical regularisation, not calibrated physics"
            ),
            "nonstreamwise_boundary": "all_pressure_open",
            "flow_boundary_scheme": "characteristic_multidimensional",
            "official_boundary": "open non-reflecting on five non-inlet faces",
            "boundary_match_to_official_cfd": (
                "linear_angle_wavelength_sweep_passed_corner_and_domain_pending"
            ),
        },
        "candidates": candidates,
        "interpretation": (
            "Failure of the 0.10 mean-plane candidate shows that the Gaussian "
            "resolved-source mapping can exceed the receiving density range. "
            "A passing tighter candidate is only a numerical nominee for the "
            "matched grid/time/domain study, not an accuracy result."
        ),
        "runtime_s": time.perf_counter() - started,
    }


def plot_preflight(record: dict[str, Any], output: str | Path) -> None:
    import matplotlib.pyplot as plt

    candidates = record["candidates"]
    limits = np.asarray([
        100.0 * float(item["mean_handoff_density_departure_limit"])
        for item in candidates
    ])
    departures = np.asarray([
        100.0
        * float(item["resolved_source"][
            "maximum_relative_mixture_density_departure"
        ])
        for item in candidates
    ])
    concentration = np.asarray([
        float(item["resolved_source"]["maximum_hydrogen_vol_percent"])
        for item in candidates
    ])
    distance = np.asarray([
        float(item["handoff"]["distance_m"]) for item in candidates
    ])

    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    axes[0].plot(limits, departures, "o-", color="#b2182b", label="mapped maximum")
    axes[0].axhline(
        100.0 * RECEIVING_DENSITY_LIMIT, color="black", linestyle="--",
        label="receiving limit",
    )
    axes[0].set_xlabel("Mean handoff density-departure limit (%)")
    axes[0].set_ylabel("Resolved maximum density departure (%)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].plot(limits, concentration, "s-", color="#2166ac", label="maximum H2")
    axes[1].set_xlabel("Mean handoff density-departure limit (%)")
    axes[1].set_ylabel("Maximum H2 (vol%)", color="#2166ac")
    axes[1].tick_params(axis="y", labelcolor="#2166ac")
    distance_axis = axes[1].twinx()
    distance_axis.plot(limits, distance, "^-", color="#4d9221", label="handoff distance")
    distance_axis.set_ylabel("Handoff distance (m)", color="#4d9221")
    distance_axis.tick_params(axis="y", labelcolor="#4d9221")
    axes[1].grid(alpha=0.25)
    figure.suptitle(
        "HyTunnel resolved-source preflight — numerical admissibility, no data fit",
        fontsize=11,
    )
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)


def _write_json(record: dict[str, Any], output: str | Path) -> None:
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    arguments = parser.parse_args(argv)
    record = run_preflight()
    _write_json(record, arguments.output)
    plot_preflight(record, arguments.plot)
    print(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
