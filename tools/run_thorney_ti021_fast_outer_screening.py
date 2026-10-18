"""Run a fast, obstacle-aware outer screening pass from a durable TI021 handoff.

This is intentionally a screening path: the outer velocity is the finite-width
potential-flow closure and the scalar is a conservative finite-volume field.
The HDF5 handoff is read lazily and its species flux is imposed directly on the
x-min boundary, so the path is not a post-hoc field overlay.  It is not a
replacement for the resolved variable-density outer solver or a facility
validation result.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.conservative_handoff_hdf5_3d import ConservativeHandoffHdf5Reader3D
from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleReceptor3D,
    FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_HANDOFF = (
    ROOT / "outputs" / "thorney-ti021-inner-handoff-v3-600s-dx5m-2026-10-13-v1"
    / "ti021-inner-handoff.h5"
)


def _log_profile(z_m: float) -> float:
    return 0.205 / 0.4 * np.log((z_m + 0.005) / 0.005)


def run_screening(
    *, handoff_path: Path, output_path: Path, duration_s: float | None = None,
    cell_size_m: float = 10.0, time_step_s: float = 0.5,
    store_fields: bool = False,
) -> dict[str, object]:
    if not handoff_path.exists():
        raise FileNotFoundError(handoff_path)
    # The durable TI021 handoff is reported in the inner crosswind coordinate
    # q=[-100,100] m.  Read that plane exactly, then embed it into the larger
    # non-negative fast-domain y=[0,640] m by the documented q -> y=q+320 m
    # translation.  The remap itself remains conservative and area preserving.
    handoff_y_edges = np.linspace(-100.0, 100.0, int(round(200.0 / cell_size_m)) + 1)
    target_y_edges = handoff_y_edges
    target_z_edges = np.linspace(0.0, 40.0, int(round(40.0 / cell_size_m)) + 1)
    if not np.isclose(target_y_edges[-1], 100.0) or not np.isclose(target_z_edges[-1], 40.0):
        raise ValueError("cell_size_m must divide the 200 m handoff width and 40 m height")
    nx = int(round(440.0 / cell_size_m))
    ny = int(round(640.0 / cell_size_m))
    nz = int(round(40.0 / cell_size_m))
    closure = BinaryIdealGasMixtureClosure(
        species_name="R12",
        species_molar_mass_kg_mol=120.913e-3,
        species_heat_capacity_j_kg_k=606.0749,
        ambient_temperature_k=291.15,
        ambient_pressure_pa=100_500.0,
        minimum_temperature_k=150.0,
    )
    with ConservativeHandoffHdf5Reader3D(
        handoff_path, target_z_edges_m=target_z_edges, target_y_edges_m=target_y_edges,
    ) as reader:
        available_duration = float(reader.time_s[-1])
        run_duration = available_duration if duration_s is None else float(duration_s)
        if run_duration <= 0.0 or run_duration > available_duration + 1.0e-12:
            raise ValueError("duration_s must lie inside the handoff history")
        z_centres = 0.5 * (target_z_edges[:-1] + target_z_edges[1:])
        profile = tuple(float(_log_profile(value)) for value in z_centres)
        ambient_u = np.asarray(profile)[:, None]
        handoff_ny = target_y_edges.size - 1
        handoff_ambient_mass_flux = closure.ambient_density_kg_m3 * np.broadcast_to(
            ambient_u, (nz, handoff_ny),
        )
        config = FastObstacleTransport3DConfig(
            domain_m=(440.0, 640.0, 40.0),
            cells=(nx, ny, nz), duration_s=run_duration,
            time_step_s=time_step_s, wind_speed_m_s=3.9,
            diffusivity_m2_s=0.15,
            source_schedule=SourceRateSchedule(
                time_s=(0.0, run_duration), rate_kg_s=(0.0, 0.0),
                source_id="TI021 handoff-only fast screening",
            ),
            source_position_m=(10.0, 320.0, 2.5), source_sigma_m=10.0,
            obstacles=(AxisAlignedCuboid(
                150.0, 200.0, 260.0, 380.0, 0.0, 20.0,
                label="TI021 surrogate facility block",
            ),),
            inlet_profile_z_m=tuple(float(value) for value in z_centres),
            inlet_profile_u_m_s=profile,
            store_fields=store_fields,
        )
        receptors = (
            FastObstacleReceptor3D("upwind", (120.0, 320.0, 5.0)),
            FastObstacleReceptor3D("roof", (220.0, 320.0, 25.0)),
            FastObstacleReceptor3D("wake", (300.0, 320.0, 5.0)),
        )
        callback_count = 0

        def inlet_flux_callback(start, dt, u_face, concentration, solid):
            nonlocal callback_count
            callback_count += 1
            frame = reader.frame_at(
                time_s=min(float(start), available_duration),
                ambient_total_mass_flux_kg_m2_s=handoff_ambient_mass_flux,
                closure=closure,
            )
            # The fast scalar operator accepts incoming positive x-directed
            # species flux. Any non-positive numerical remap residue is not an
            # admissible injection and is fail-closed to zero.
            embedded = np.zeros((nz, ny), dtype=float)
            start_index = int(round((320.0 - 100.0) / cell_size_m))
            embedded[:, start_index:start_index + handoff_ny] = np.maximum(
                frame.species_mass_flux_kg_m2_s, 0.0,
            )
            return embedded

        started = time.perf_counter()
        result = solve_fast_obstacle_transport_3d(
            config, receptors=receptors,
            interface_inlet_flux_callback=inlet_flux_callback,
        )
        runtime = time.perf_counter() - started
        molar_density = closure.ambient_pressure_pa / (
            8.31446261815324 * closure.ambient_temperature_k
        )
        vol_percent = result.maximum_concentration_kg_m3 / (
            closure.species_molar_mass_kg_mol * molar_density
        ) * 100.0
        fluid = np.isfinite(vol_percent) & (~result.solid_mask)
        summary = {
            "schema": "degali.thorney-ti021-fast-outer-screening.v1",
            "status": "conditional_fast_outer_screening",
            "quantitative_prediction_allowed": False,
            "handoff_path": str(handoff_path.resolve()),
            "handoff_source_id": reader.source_id,
            "handoff_available_duration_s": available_duration,
            "duration_s": run_duration,
            "cell_size_m": [cell_size_m] * 3,
            "cells": [nx, ny, nz],
            "time_step_s": float(result.diagnostics["time_step_s"]),
            "runtime_s": runtime,
            "runtime_per_simulated_second_s": runtime / run_duration,
            "handoff_callback_count": callback_count,
            "full_field_history_stored": bool(store_fields),
            "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
            "maximum_species_proxy_vol_percent": float(np.nanmax(vol_percent[fluid])) if np.any(fluid) else 0.0,
            "maximum_instantaneous_volume_above_4_vol_percent_m3": float(
                np.count_nonzero(result.maximum_concentration_kg_m3 / (
                    closure.species_molar_mass_kg_mol * molar_density
                ) * 100.0 >= 4.0) * np.prod(result.diagnostics["cell_size_m"])
            ),
            "receptor_trace_kg_m3": {key: value.tolist() for key, value in result.receptor_traces_kg_m3.items()},
            "diagnostics": result.diagnostics,
            "warnings": [
                "R12 handoff species proxy, not an H2 field",
                "potential flow omits separation, recirculation, turbulence and buoyancy feedback",
                "fast screening is not a resolved outer validation or facility release decision",
            ],
        }
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        if config.store_fields and result.concentration_fields_kg_m3 is not None:
            np.savez_compressed(
                output_path.with_name(output_path.stem + "-field.npz"),
                x_m=result.x_m, y_m=result.y_m, z_m=result.z_m,
                time_s=result.time_s, maximum_concentration_kg_m3=result.maximum_concentration_kg_m3,
                concentration_fields_kg_m3=result.concentration_fields_kg_m3,
                solid_mask=result.solid_mask,
            )
        return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--handoff", type=Path, default=DEFAULT_HANDOFF)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--duration-s", type=float, default=None)
    parser.add_argument("--cell-size-m", type=float, default=10.0)
    parser.add_argument("--time-step-s", type=float, default=0.5)
    parser.add_argument("--store-fields", action="store_true", help="store the full time-resolved field (large)")
    args = parser.parse_args()
    print(json.dumps(run_screening(
        handoff_path=args.handoff, output_path=args.output,
        duration_s=args.duration_s, cell_size_m=args.cell_size_m,
        time_step_s=args.time_step_s,
        store_fields=args.store_fields,
    ), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
