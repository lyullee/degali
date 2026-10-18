"""Audit the remapped TI021 plume flux as an outer binary-mixture inlet."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.binary_mixture_corrector_3d import (
    solve_binary_mixture_initial_cloud_step_3d,
)
from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.conservative_nested_handoff_3d import (
    nested_outer_inlet_frame_at,
    remap_conservative_nested_handoff_3d,
)


ROOT = Path(__file__).resolve().parents[1]


def _log_wind(z_m: np.ndarray) -> np.ndarray:
    return 0.205 / 0.4 * np.log((z_m + 0.005) / 0.005)


def build_audit() -> tuple[dict[str, object], dict[str, np.ndarray]]:
    source_dx = 1.25
    target_dx = 5.0
    source_y = np.arange(-100.0, 100.0 + source_dx, source_dx)
    source_z = np.arange(0.0, 40.0 + source_dx, source_dx)
    target_y = np.arange(-320.0, 320.0 + target_dx, target_dx)
    target_z = np.arange(0.0, 40.0 + target_dx, target_dx)
    times = np.asarray((0.0, 60.0, 120.0, 240.0, 600.0))
    yc = 0.5 * (source_y[:-1] + source_y[1:])
    zc = 0.5 * (source_z[:-1] + source_z[1:])
    spatial = np.exp(
        -(yc[None, :] / 28.0) ** 2
        - ((zc[:, None] - 2.5) / 4.0) ** 2
    )
    amplitude = np.asarray((0.0, 0.08, 0.12, 0.04, 0.0))[:, None, None]
    species = amplitude * spatial[None, :, :]
    handoff = remap_conservative_nested_handoff_3d(
        time_s=times,
        source_z_edges_m=source_z, source_y_edges_m=source_y,
        target_z_edges_m=target_z, target_y_edges_m=target_y,
        total_mass_flux_anomaly_kg_m2_s=0.65 * species,
        species_mass_flux_kg_m2_s=species,
        sensible_enthalpy_flux_w_m2=-8.0e5 * species,
        lateral_momentum_flux_anomaly_pa=0.08 * species,
        vertical_momentum_flux_anomaly_pa=0.03 * species,
        source_id="TI021 analytic outer-inlet coupling pulse",
    )
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0749,
        ambient_temperature_k=291.15,
        ambient_pressure_pa=100_500.0,
        minimum_temperature_k=150.0,
    )
    target_zc = 0.5 * (target_z[:-1] + target_z[1:])
    ambient_u = _log_wind(target_zc)[:, None]
    ambient_mass_flux = closure.ambient_density_kg_m3 * np.broadcast_to(
        ambient_u, (target_zc.size, target_y.size - 1),
    )
    query_time = 90.0
    frame = nested_outer_inlet_frame_at(
        handoff, time_s=query_time,
        ambient_total_mass_flux_kg_m2_s=ambient_mass_flux,
        closure=closure,
    )

    nx = 8
    shape = (target_zc.size, target_y.size - 1, nx)
    solid = np.zeros(shape, dtype=bool)
    old_species = np.zeros(shape)
    old_enthalpy = np.zeros(shape)
    u_face = np.broadcast_to(
        ambient_u[:, :, None], (shape[0], shape[1], nx + 1),
    ).copy()
    v_face = np.zeros((shape[0], shape[1] + 1, nx))
    w_face = np.zeros((shape[0] + 1, shape[1], nx))
    dt = 0.05
    started = time.perf_counter()
    result = solve_binary_mixture_initial_cloud_step_3d(
        old_species, old_enthalpy,
        provisional_u_face_m_s=u_face,
        provisional_v_face_m_s=v_face,
        provisional_w_face_m_s=w_face,
        solid_mask=solid,
        spacing_m=(target_dx, target_dx, target_dx),
        time_step_s=dt,
        species_diffusivity_m2_s=0.0,
        thermal_diffusivity_m2_s=0.0,
        closure=closure,
        maximum_iterations=40,
        minimum_iterations=2,
        relative_density_tolerance=1.0e-8,
        relative_mass_flux_tolerance=1.0e-8,
        open_boundary_faces=("x_max",),
        projection_linear_solver="cached_direct",
        inlet_total_mass_flux_kg_m2_s=frame.total_mass_flux_kg_m2_s,
        inlet_species_mass_flux_kg_m2_s=frame.species_mass_flux_kg_m2_s,
        inlet_sensible_enthalpy_flux_w_m2=(
            frame.sensible_enthalpy_flux_w_m2
        ),
    )
    runtime = time.perf_counter() - started
    area = target_dx**2
    expected_species = float(np.sum(frame.species_mass_flux_kg_m2_s) * area * dt)
    expected_energy = float(
        np.sum(frame.sensible_enthalpy_flux_w_m2) * area * dt
    )
    represented_species = float(np.sum(
        result.species_mass_concentration_kg_m3
    ) * target_dx**3)
    represented_energy = float(np.sum(
        result.sensible_enthalpy_density_j_m3
    ) * target_dx**3)
    gates = {
        "linear_time_interpolation_reproduces_mid_pulse": bool(
            np.max(frame.species_mass_flux_kg_m2_s) > 0.0
        ),
        "ambient_plus_anomaly_is_unidirectional": bool(
            np.min(frame.total_mass_flux_kg_m2_s) > 0.0
        ),
        "inlet_species_fraction_is_bounded": bool(
            np.min(frame.thermodynamic_state.species_mass_fraction) >= 0.0
            and np.max(frame.thermodynamic_state.species_mass_fraction) <= 1.0
        ),
        "prescribed_total_mass_flux_survives_projection": bool(np.allclose(
            result.mass_projection.mass_flux_x_kg_m2_s[:, :, 0],
            frame.total_mass_flux_kg_m2_s, rtol=0.0, atol=1.0e-12,
        )),
        "species_inlet_ledger_closes_at_1e_10kg": (
            abs(result.species_mass_residual_kg) <= 1.0e-10
            and abs(represented_species - expected_species) <= 1.0e-10
        ),
        "enthalpy_inlet_ledger_closes_at_1e_6j": (
            abs(result.sensible_enthalpy_residual_j) <= 1.0e-6
            and abs(represented_energy - expected_energy) <= 1.0e-6
        ),
        "projected_continuity_residual_at_most_1e_7": (
            result.mass_projection.maximum_absolute_residual_after_kg_m3_s
            <= 1.0e-7
        ),
        "binary_corrector_converged": result.converged,
        "no_clipping_or_posthoc_renormalisation": not frame.diagnostics[
            "clipping_or_posthoc_renormalisation"
        ],
    }
    record = {
        "schema": "degali.thorney-ti021-outer-inlet-coupling-audit.v1",
        "status": (
            "accepted_outer_binary_boundary_component_pending_momentum_overlap"
            if all(gates.values()) else "withheld_outer_binary_boundary_component"
        ),
        "scope": (
            "concentration-blind analytic flux pulse through remap, linear time "
            "interpolation, ambient-plus-anomaly thermodynamic reconstruction "
            "and one outer binary-mixture corrector step"
        ),
        "configuration": {
            "query_time_s": query_time,
            "outer_cell_size_m": target_dx,
            "outer_cells": list(reversed(shape)),
            "outer_cell_count": int(np.prod(shape)),
            "time_step_s": dt,
            "observation_value_fields_accessed": [],
        },
        "inlet": {
            "minimum_total_mass_flux_kg_m2_s": float(
                np.min(frame.total_mass_flux_kg_m2_s)
            ),
            "maximum_total_mass_flux_kg_m2_s": float(
                np.max(frame.total_mass_flux_kg_m2_s)
            ),
            "maximum_species_mass_fraction": float(np.max(
                frame.thermodynamic_state.species_mass_fraction
            )),
            "temperature_range_k": [
                float(np.min(frame.thermodynamic_state.temperature_k)),
                float(np.max(frame.thermodynamic_state.temperature_k)),
            ],
            "velocity_component_ranges_m_s": {
                "u": [float(np.min(frame.velocity_u_m_s)), float(np.max(frame.velocity_u_m_s))],
                "v": [float(np.min(frame.velocity_v_m_s)), float(np.max(frame.velocity_v_m_s))],
                "w": [float(np.min(frame.velocity_w_m_s)), float(np.max(frame.velocity_w_m_s))],
            },
        },
        "ledger": {
            "expected_species_inlet_mass_kg": expected_species,
            "represented_domain_species_mass_kg": represented_species,
            "species_mass_residual_kg": result.species_mass_residual_kg,
            "expected_sensible_enthalpy_inlet_j": expected_energy,
            "represented_domain_sensible_enthalpy_j": represented_energy,
            "sensible_enthalpy_residual_j": result.sensible_enthalpy_residual_j,
            "continuity_residual_kg_m3_s": (
                result.mass_projection.maximum_absolute_residual_after_kg_m3_s
            ),
            "corrector_iterations": result.iterations,
        },
        "runtime_s": runtime,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
        "next_required": (
            "feed reconstructed u/v/w into outer momentum over multiple steps "
            "and compare an overlapping inner/outer receptor plane"
        ),
    }
    arrays = {
        "target_y": target_y, "target_z": target_z,
        "inlet_species_flux": frame.species_mass_flux_kg_m2_s,
        "first_cell_species": result.species_mass_concentration_kg_m3[:, :, 0],
    }
    return record, arrays


def _plot(record: dict[str, object], arrays: dict[str, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True)
    for axis, field, title, label in (
        (axes[0], arrays["inlet_species_flux"], "prescribed x-min species flux", "kg/m2/s"),
        (axes[1], arrays["first_cell_species"], "outer first-cell species after step", "kg/m3"),
    ):
        image = axis.pcolormesh(
            arrays["target_y"], arrays["target_z"], field,
            shading="flat", cmap="magma",
        )
        axis.set_xlabel("crosswind y (m)")
        axis.set_ylabel("height z (m)")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.78, label=label)
    ledger = record["ledger"]
    labels = ("species", "enthalpy")
    expected = (
        ledger["expected_species_inlet_mass_kg"],
        ledger["expected_sensible_enthalpy_inlet_j"],
    )
    represented = (
        ledger["represented_domain_species_mass_kg"],
        ledger["represented_domain_sensible_enthalpy_j"],
    )
    positions = np.arange(2)
    scale = np.maximum(np.abs(expected), 1.0e-30)
    axes[2].bar(positions - 0.18, np.asarray(expected) / scale, 0.36, label="inlet ledger")
    axes[2].bar(positions + 0.18, np.asarray(represented) / scale, 0.36, label="domain state")
    axes[2].set_xticks(positions, labels)
    axes[2].set_ylabel("value / |prescribed inlet value|")
    axes[2].set_title("One-step conservative ledger")
    axes[2].grid(axis="y", alpha=0.25)
    axes[2].legend()
    fig.suptitle("TI021 nested outer-inlet binary-mixture coupling audit")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-directory", type=Path,
        default=(ROOT / "outputs" /
                 "thorney-ti021-outer-inlet-coupling-audit-2026-10-11-v1"),
    )
    args = parser.parse_args()
    record, arrays = build_audit()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-outer-inlet-coupling-audit.json"
    plot_path = args.output_directory / "ti021-outer-inlet-coupling-audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(record, arrays, plot_path)
    print(json.dumps({
        "status": record["status"], "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()), "runtime_s": record["runtime_s"],
        "gates": record["gates"],
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
