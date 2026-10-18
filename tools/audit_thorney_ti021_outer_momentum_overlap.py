"""Audit three-component nested momentum coupling and an overlap solution.

The pulse is analytic and concentration-blind: no TI021 dose, peak or timing
value is opened.  A narrow inner-extension reference and a wide nested outer
calculation use the same 5 m discretisation and receive the same conservatively
remapped inlet history.  This isolates the nesting/domain coupling from the
separate three-grid convergence question.  It closes the numerical handoff
gate only; it is not the field-experiment score.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.atmospheric_surface_layer_mixing import (
    NeutralSurfaceLayerMixing,
)
from degali.addons.binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
)
from degali.addons.conservative_nested_handoff_3d import (
    remap_conservative_nested_handoff_3d,
)
from degali.addons.nested_outer_binary_flow_3d import (
    advance_nested_outer_binary_flow_step_3d,
    initialise_nested_outer_binary_flow_3d,
)


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-outer-momentum-overlap-audit.v1"


def _log_wind(z_m: np.ndarray) -> np.ndarray:
    return 0.205 / 0.4 * np.log((z_m + 0.005) / 0.005)


def _volume_average_2x(field: np.ndarray) -> np.ndarray:
    nz, ny, nx = field.shape
    if nz % 2 or ny % 2 or nx % 2:
        raise ValueError("refined overlap shape must be divisible by two")
    return field.reshape(
        nz // 2, 2, ny // 2, 2, nx // 2, 2,
    ).mean(axis=(1, 3, 5))


def _history_at(time_s: np.ndarray, field: np.ndarray, query_s: float) -> np.ndarray:
    upper = int(np.searchsorted(time_s, query_s, side="right"))
    if upper == 0:
        return field[0]
    if upper >= time_s.size:
        return field[-1]
    lower = upper - 1
    fraction = (query_s - time_s[lower]) / (time_s[upper] - time_s[lower])
    return (1.0 - fraction) * field[lower] + fraction * field[upper]


def _make_handoff(
    *,
    target_y_edges_m: np.ndarray,
    target_z_edges_m: np.ndarray,
):
    source_dx = 1.25
    source_y = np.arange(-100.0, 100.0 + source_dx, source_dx)
    source_z = np.arange(0.0, 40.0 + source_dx, source_dx)
    time_s = np.asarray((0.0, 1.5, 3.0, 4.5, 6.0))
    amplitude = np.asarray((0.0, 0.055, 0.085, 0.035, 0.0))[:, None, None]
    yc = 0.5 * (source_y[:-1] + source_y[1:])
    zc = 0.5 * (source_z[:-1] + source_z[1:])
    spatial = np.exp(
        -(yc[None, :] / 28.0) ** 2
        - ((zc[:, None] - 5.0) / 6.0) ** 2
    )
    species = amplitude * spatial[None, :, :]
    common = dict(
        time_s=time_s,
        source_z_edges_m=source_z,
        source_y_edges_m=source_y,
        total_mass_flux_anomaly_kg_m2_s=0.6 * species,
        species_mass_flux_kg_m2_s=species,
        sensible_enthalpy_flux_w_m2=-8.0e5 * species,
        lateral_momentum_flux_anomaly_pa=4.0 * species,
        vertical_momentum_flux_anomaly_pa=2.0 * species,
    )
    handoff = remap_conservative_nested_handoff_3d(
        target_z_edges_m=target_z_edges_m,
        target_y_edges_m=target_y_edges_m,
        source_id=(
            "TI021 concentration-blind analytic pulse for outer momentum and "
            "overlap discretisation audit"
        ),
        **common,
    )
    return handoff, source_y, source_z, common


def _run_domain(
    *,
    dx: float,
    y_edges: np.ndarray,
    z_edges: np.ndarray,
    x_length_m: float,
    duration_s: float,
    dt: float,
    closure: BinaryIdealGasMixtureClosure,
    handoff,
) -> tuple[object, dict[str, np.ndarray | float | bool]]:
    nx = int(round(x_length_m / dx))
    ny = y_edges.size - 1
    nz = z_edges.size - 1
    solid = np.zeros((nz, ny, nx), dtype=bool)
    zc = 0.5 * (z_edges[:-1] + z_edges[1:])
    ambient_u = np.broadcast_to(_log_wind(zc)[:, None], (nz, ny)).copy()
    ambient_mass_flux = closure.ambient_density_kg_m3 * ambient_u
    state = initialise_nested_outer_binary_flow_3d(
        cells=(nx, ny, nz),
        closure=closure,
        ambient_u_m_s=ambient_u,
        solid_mask=solid,
    )
    mixing = NeutralSurfaceLayerMixing(
        friction_velocity_m_s=0.205,
        roughness_length_m=0.005,
        source_id="TI021 u*=0.205 m/s neutral primary overlap audit",
    )
    nu_t = mixing.eddy_viscosity_m2_s(zc, solid)
    viscosity = 1.5e-5 + nu_t
    species_diffusivity = 9.78865e-6 + nu_t / 0.7
    thermal_diffusivity = 3.09e-6 + nu_t / 0.74
    steps = int(round(duration_s / dt))
    if not np.isclose(steps * dt, duration_s):
        raise ValueError("duration_s must be an integer number of time steps")
    vol_history = []
    max_continuity = 0.0
    max_momentum_residual = 0.0
    all_converged = True
    exact_inlet_mass_flux = True
    max_lateral_flux_residual = 0.0
    max_vertical_flux_residual = 0.0
    started = time.perf_counter()
    for _ in range(steps):
        result = advance_nested_outer_binary_flow_step_3d(
            state,
            handoff=handoff,
            ambient_total_mass_flux_kg_m2_s=ambient_mass_flux,
            closure=closure,
            solid_mask=solid,
            spacing_m=(dx, dx, dx),
            time_step_s=dt,
            kinematic_viscosity_m2_s=viscosity,
            species_diffusivity_m2_s=species_diffusivity,
            thermal_diffusivity_m2_s=thermal_diffusivity,
            projection_linear_solver="cached_direct",
            maximum_corrector_iterations=40,
            minimum_corrector_iterations=2,
            relative_density_tolerance=1.0e-8,
            relative_mass_flux_tolerance=1.0e-8,
        )
        state = result.state
        vol_history.append(
            100.0 * state.thermodynamic_state.species_mole_fraction.copy()
        )
        max_continuity = max(
            max_continuity,
            result.binary_corrector.mass_projection
            .maximum_absolute_residual_after_kg_m3_s,
        )
        max_momentum_residual = max(
            max_momentum_residual,
            float(np.max(np.abs(
                result.momentum.momentum_balance_residual_kg_m_s
            ))),
        )
        all_converged = all_converged and result.binary_corrector.converged
        exact_inlet_mass_flux = exact_inlet_mass_flux and bool(np.allclose(
            result.binary_corrector.mass_projection
            .mass_flux_x_kg_m2_s[:, :, 0],
            result.inlet_frame.total_mass_flux_kg_m2_s,
            rtol=0.0,
            atol=1.0e-12,
        ))
        query = float(result.diagnostics["inlet_evaluation_time_s"])
        expected_lateral = _history_at(
            handoff.time_s,
            handoff.lateral_momentum_flux_anomaly_pa,
            query,
        )
        expected_vertical = _history_at(
            handoff.time_s,
            handoff.vertical_momentum_flux_anomaly_pa,
            query,
        )
        represented_lateral = (
            result.inlet_frame.total_mass_flux_kg_m2_s
            * result.inlet_frame.velocity_v_m_s
        )
        represented_vertical = (
            result.inlet_frame.total_mass_flux_kg_m2_s
            * result.inlet_frame.velocity_w_m_s
        )
        max_lateral_flux_residual = max(
            max_lateral_flux_residual,
            float(np.max(np.abs(represented_lateral - expected_lateral))),
        )
        max_vertical_flux_residual = max(
            max_vertical_flux_residual,
            float(np.max(np.abs(represented_vertical - expected_vertical))),
        )
    runtime = time.perf_counter() - started
    return state, {
        "vol_percent_history": np.asarray(vol_history),
        "runtime_s": runtime,
        "maximum_continuity_residual_kg_m3_s": max_continuity,
        "maximum_momentum_balance_residual_kg_m_s": max_momentum_residual,
        "all_corrector_steps_converged": all_converged,
        "prescribed_inlet_mass_flux_exact_every_step": exact_inlet_mass_flux,
        "maximum_lateral_momentum_flux_reconstruction_residual_pa": (
            max_lateral_flux_residual
        ),
        "maximum_vertical_momentum_flux_reconstruction_residual_pa": (
            max_vertical_flux_residual
        ),
    }


def build_audit() -> tuple[dict[str, object], dict[str, np.ndarray]]:
    overlap_dx = 5.0
    duration = 6.0
    dt = 0.125
    x_length = 20.0
    reference_y = np.arange(-100.0, 100.0 + overlap_dx, overlap_dx)
    reference_z = np.arange(0.0, 40.0 + overlap_dx, overlap_dx)
    outer_y = np.arange(-320.0, 320.0 + overlap_dx, overlap_dx)
    outer_z = np.arange(0.0, 40.0 + overlap_dx, overlap_dx)
    reference_handoff, _, _, _ = _make_handoff(
        target_y_edges_m=reference_y,
        target_z_edges_m=reference_z,
    )
    outer_handoff, _, _, _ = _make_handoff(
        target_y_edges_m=outer_y,
        target_z_edges_m=outer_z,
    )
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0749,
        ambient_temperature_k=291.15,
        ambient_pressure_pa=100_500.0,
        minimum_temperature_k=150.0,
    )
    reference_state, reference = _run_domain(
        dx=overlap_dx,
        y_edges=reference_y,
        z_edges=reference_z,
        x_length_m=x_length,
        duration_s=duration,
        dt=dt,
        closure=closure,
        handoff=reference_handoff,
    )
    outer_state, outer = _run_domain(
        dx=overlap_dx,
        y_edges=outer_y,
        z_edges=outer_z,
        x_length_m=x_length,
        duration_s=duration,
        dt=dt,
        closure=closure,
        handoff=outer_handoff,
    )
    refined_dx = 2.5
    refined_y = np.arange(-100.0, 100.0 + refined_dx, refined_dx)
    refined_z = np.arange(0.0, 40.0 + refined_dx, refined_dx)
    refined_handoff, _, _, _ = _make_handoff(
        target_y_edges_m=refined_y,
        target_z_edges_m=refined_z,
    )
    _, refined = _run_domain(
        dx=refined_dx,
        y_edges=refined_y,
        z_edges=refined_z,
        x_length_m=x_length,
        duration_s=duration,
        dt=dt,
        closure=closure,
        handoff=refined_handoff,
    )
    reference_vol = reference["vol_percent_history"]
    outer_y_start = int(round((reference_y[0] - outer_y[0]) / overlap_dx))
    outer_y_stop = outer_y_start + reference_vol.shape[2]
    outer_overlap = outer["vol_percent_history"][:, :, outer_y_start:outer_y_stop, :]
    if outer_overlap.shape != reference_vol.shape:
        raise RuntimeError("reference and outer overlap arrays do not align")

    coarse_zc = 0.5 * (outer_z[:-1] + outer_z[1:])
    reference_yc = 0.5 * (reference_y[:-1] + reference_y[1:])
    overlap_xc = np.arange(0.5 * overlap_dx, x_length, overlap_dx)
    receptor_indices = []
    for x_target in (2.5, 7.5, 12.5, 17.5):
        for z_target in (2.5, 7.5):
            receptor_indices.append((
                int(np.argmin(np.abs(coarse_zc - z_target))),
                int(np.argmin(np.abs(reference_yc))),
                int(np.argmin(np.abs(overlap_xc - x_target))),
            ))
    reference_traces = np.asarray([
        reference_vol[:, iz, iy, ix]
        for iz, iy, ix in receptor_indices
    ])
    outer_traces = np.asarray([
        outer_overlap[:, iz, iy, ix]
        for iz, iy, ix in receptor_indices
    ])
    reference_peak = np.max(reference_traces, axis=1)
    outer_peak = np.max(outer_traces, axis=1)
    reference_dose = np.sum(reference_traces, axis=1) * dt
    outer_dose = np.sum(outer_traces, axis=1) * dt
    active = reference_peak >= 1.0e-4
    peak_relative = (
        np.abs(outer_peak[active] - reference_peak[active])
        / reference_peak[active]
    )
    dose_relative = (
        np.abs(outer_dose[active] - reference_dose[active])
        / reference_dose[active]
    )
    refined_on_overlap = np.asarray([
        _volume_average_2x(field)
        for field in refined["vol_percent_history"]
    ])
    refined_traces = np.asarray([
        refined_on_overlap[:, iz, iy, ix]
        for iz, iy, ix in receptor_indices
    ])
    refined_peak = np.max(refined_traces, axis=1)
    refined_dose = np.sum(refined_traces, axis=1) * dt
    refined_active = refined_peak >= 1.0e-4
    grid_change_peak = (
        np.abs(outer_peak[refined_active] - refined_peak[refined_active])
        / refined_peak[refined_active]
    )
    grid_change_dose = (
        np.abs(outer_dose[refined_active] - refined_dose[refined_active])
        / refined_dose[refined_active]
    )
    maximum_peak_change = float(np.max(peak_relative))
    maximum_dose_change = float(np.max(dose_relative))
    maximum_continuity = max(
        float(reference["maximum_continuity_residual_kg_m3_s"]),
        float(outer["maximum_continuity_residual_kg_m3_s"]),
    )
    maximum_momentum_residual = max(
        float(reference["maximum_momentum_balance_residual_kg_m_s"]),
        float(outer["maximum_momentum_balance_residual_kg_m_s"]),
    )
    maximum_lateral_residual = max(
        float(reference["maximum_lateral_momentum_flux_reconstruction_residual_pa"]),
        float(outer["maximum_lateral_momentum_flux_reconstruction_residual_pa"]),
    )
    maximum_vertical_residual = max(
        float(reference["maximum_vertical_momentum_flux_reconstruction_residual_pa"]),
        float(outer["maximum_vertical_momentum_flux_reconstruction_residual_pa"]),
    )
    gates = {
        "observation_values_excluded": True,
        "three_component_momentum_flux_reconstructed_at_roundoff": (
            maximum_lateral_residual <= 1.0e-12
            and maximum_vertical_residual <= 1.0e-12
        ),
        "prescribed_inlet_mass_flux_exact_every_step": bool(
            reference["prescribed_inlet_mass_flux_exact_every_step"]
            and outer["prescribed_inlet_mass_flux_exact_every_step"]
        ),
        "all_binary_corrector_steps_converged": bool(
            reference["all_corrector_steps_converged"]
            and outer["all_corrector_steps_converged"]
        ),
        "maximum_continuity_residual_at_most_1e_7": (
            maximum_continuity <= 1.0e-7
        ),
        "maximum_momentum_balance_residual_at_most_1e_9": (
            maximum_momentum_residual <= 1.0e-9
        ),
        "cumulative_species_ledger_at_most_1e_8kg": (
            abs(reference_state.cumulative_species_mass_residual_kg) <= 1.0e-8
            and abs(outer_state.cumulative_species_mass_residual_kg) <= 1.0e-8
        ),
        "cumulative_enthalpy_ledger_at_most_1e_3j": (
            abs(reference_state.cumulative_sensible_enthalpy_residual_j) <= 1.0e-3
            and abs(outer_state.cumulative_sensible_enthalpy_residual_j) <= 1.0e-3
        ),
        "overlap_receptor_peak_change_at_most_10_percent": (
            maximum_peak_change <= 0.10
        ),
        "overlap_receptor_dose_change_at_most_10_percent": (
            maximum_dose_change <= 0.10
        ),
    }
    record = {
        "schema": SCHEMA,
        "status": (
            "accepted_outer_momentum_and_overlap_component"
            if all(gates.values())
            else "withheld_outer_momentum_and_overlap_component"
        ),
        "scope": (
            "concentration-blind analytic multi-step comparison of same-grid "
            "5 m narrow-reference and conservative wide outer domains; "
            "not a TI021 observation score"
        ),
        "configuration": {
            "reference_extension_cell_size_m": overlap_dx,
            "outer_cell_size_m": overlap_dx,
            "remap_source_cell_size_m": 1.25,
            "overlap_domain_m": [x_length, 200.0, 40.0],
            "outer_crosswind_width_m": 640.0,
            "duration_s": duration,
            "time_step_s": dt,
            "step_count_per_domain": int(round(duration / dt)),
            "primary_turbulent_schmidt_number": 0.7,
            "primary_turbulent_prandtl_number": 0.74,
            "observation_value_fields_accessed": [],
            "receptor_count": len(receptor_indices),
        },
        "momentum_and_conservation": {
            "maximum_lateral_momentum_flux_reconstruction_residual_pa": (
                maximum_lateral_residual
            ),
            "maximum_vertical_momentum_flux_reconstruction_residual_pa": (
                maximum_vertical_residual
            ),
            "maximum_momentum_balance_residual_kg_m_s": maximum_momentum_residual,
            "maximum_continuity_residual_kg_m3_s": maximum_continuity,
            "reference_cumulative_species_mass_residual_kg": (
                reference_state.cumulative_species_mass_residual_kg
            ),
            "outer_cumulative_species_mass_residual_kg": (
                outer_state.cumulative_species_mass_residual_kg
            ),
            "reference_cumulative_sensible_enthalpy_residual_j": (
                reference_state.cumulative_sensible_enthalpy_residual_j
            ),
            "outer_cumulative_sensible_enthalpy_residual_j": (
                outer_state.cumulative_sensible_enthalpy_residual_j
            ),
        },
        "overlap_receptors": {
            "coordinates_xyz_m": [
                [float(overlap_xc[ix]), float(reference_yc[iy]), float(coarse_zc[iz])]
                for iz, iy, ix in receptor_indices
            ],
            "reference_peak_vol_percent": reference_peak.tolist(),
            "outer_peak_vol_percent": outer_peak.tolist(),
            "reference_dose_vol_percent_s": reference_dose.tolist(),
            "outer_dose_vol_percent_s": outer_dose.tolist(),
            "peak_relative_change": peak_relative.tolist(),
            "dose_relative_change": dose_relative.tolist(),
            "maximum_peak_relative_change": maximum_peak_change,
            "maximum_dose_relative_change": maximum_dose_change,
            "acceptance_limit": 0.10,
        },
        "separate_grid_change_diagnostic": {
            "status": (
                "diagnostic_passes_10_percent_but_not_a_three_grid_audit"
                if max(float(np.max(grid_change_peak)), float(np.max(grid_change_dose))) <= 0.10
                else "diagnostic_rejects_2p5_to_5m_as_converged"
            ),
            "purpose": (
                "kept outside the same-grid handoff gate so nesting and grid "
                "convergence are not conflated"
            ),
            "refined_cell_size_m": refined_dx,
            "coarse_cell_size_m": overlap_dx,
            "maximum_peak_relative_change": float(np.max(grid_change_peak)),
            "maximum_dose_relative_change": float(np.max(grid_change_dose)),
            "promotion_effect": (
                "5 m is not promoted as a grid-converged production mesh by "
                "this component audit; the full three-grid TI021 audit remains"
            ),
        },
        "runtime_s": {
            "reference_extension": reference["runtime_s"],
            "outer": outer["runtime_s"],
            "refined_grid_change_diagnostic": refined["runtime_s"],
            "total": (
                float(reference["runtime_s"])
                + float(outer["runtime_s"])
                + float(refined["runtime_s"])
            ),
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
        "limitations": [
            "analytic pulse only; field-observation values were not used",
            "overlap acceptance qualifies the numerical nesting component, not R12 accuracy",
            "the separate 2.5-to-5 m diagnostic is retained and cannot promote the 5 m mesh",
            "the 600 s three-grid TI021 run and 145-sensor no-fit score remain required",
        ],
    }
    arrays = {
        "time_s": (np.arange(int(round(duration / dt))) + 1.0) * dt,
        "reference_traces": reference_traces,
        "outer_traces": outer_traces,
        "reference_final": reference_vol[-1],
        "outer_final": outer_overlap[-1],
        "overlap_xc": overlap_xc,
        "reference_yc": reference_yc,
        "coarse_zc": coarse_zc,
        "peak_relative": peak_relative,
        "dose_relative": dose_relative,
        "grid_change_peak_relative": grid_change_peak,
        "grid_change_dose_relative": grid_change_dose,
    }
    return record, arrays


def _plot(record: dict[str, object], arrays: dict[str, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), constrained_layout=True)
    peak_index = int(np.argmax(arrays["reference_traces"].max(axis=1)))
    axes[0, 0].plot(
        arrays["time_s"], arrays["reference_traces"][peak_index],
        linewidth=2, label="5 m narrow reference",
    )
    axes[0, 0].plot(
        arrays["time_s"], arrays["outer_traces"][peak_index],
        "--", linewidth=2, label="5 m nested outer",
    )
    axes[0, 0].set_xlabel("time (s)")
    axes[0, 0].set_ylabel("R12 (vol%)")
    axes[0, 0].set_title("Common overlap receptor trace")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.25)

    z_index = 0
    extent = [
        arrays["overlap_xc"][0] - 2.5, arrays["overlap_xc"][-1] + 2.5,
        arrays["reference_yc"][0] - 2.5,
        arrays["reference_yc"][-1] + 2.5,
    ]
    shared_max = max(
        float(np.max(arrays["reference_final"][z_index])),
        float(np.max(arrays["outer_final"][z_index])),
        1.0e-12,
    )
    for axis, field, title in (
        (axes[0, 1], arrays["reference_final"][z_index], "5 m narrow reference"),
        (axes[1, 0], arrays["outer_final"][z_index], "5 m nested outer"),
    ):
        image = axis.imshow(
            field,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="magma",
            vmin=0.0,
            vmax=shared_max,
        )
        axis.set_xlabel("downwind x (m)")
        axis.set_ylabel("crosswind y (m)")
        axis.set_title(title + "\nlow 5 m layer, final time")
        fig.colorbar(image, ax=axis, label="R12 (vol%)", shrink=0.8)

    positions = np.arange(2)
    peak_values = 100.0 * np.asarray((
        np.max(arrays["peak_relative"]),
        np.max(arrays["grid_change_peak_relative"]),
    ))
    dose_values = 100.0 * np.asarray((
        np.max(arrays["dose_relative"]),
        np.max(arrays["grid_change_dose_relative"]),
    ))
    axes[1, 1].bar(
        positions - 0.18, peak_values, 0.36,
        label="peak change",
    )
    axes[1, 1].bar(
        positions + 0.18, dose_values, 0.36,
        label="dose change",
    )
    axes[1, 1].axhline(10.0, color="black", linestyle="--", label="10% gate")
    axes[1, 1].set_xticks(
        positions,
        ("same-grid\nhandoff gate", "2.5→5 m\ngrid diagnostic"),
    )
    axes[1, 1].set_ylabel("relative change (%)")
    axes[1, 1].set_title("Separate nesting from grid convergence")
    axes[1, 1].grid(axis="y", alpha=0.25)
    axes[1, 1].legend(fontsize=8)
    fig.suptitle(
        "TI021 conservative outer momentum + overlap audit\n"
        "analytic concentration-blind pulse; not field-data scoring"
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-outer-momentum-overlap-audit-2026-10-11-v1"
        ),
    )
    args = parser.parse_args()
    record, arrays = build_audit()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-outer-momentum-overlap-audit.json"
    plot_path = args.output_directory / "ti021-outer-momentum-overlap-audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(record, arrays, plot_path)
    print(json.dumps({
        "status": record["status"],
        "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()),
        "runtime_s": record["runtime_s"],
        "maximum_peak_relative_change": record["overlap_receptors"][
            "maximum_peak_relative_change"
        ],
        "maximum_dose_relative_change": record["overlap_receptors"][
            "maximum_dose_relative_change"
        ],
        "gates": record["gates"],
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
