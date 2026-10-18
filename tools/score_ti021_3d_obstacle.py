"""Diagnostic 3-D TI021 obstacle score using the conservative dense-gas solver.

This is an observation-blind reduced-order 3-D witness.  It uses the declared
TI021 source inventory and workbook wind, represents the registered
semi-circular fence as a finite-thickness solid mask, and reports the same
600 s positive-dose MAPE as the integral-layer scorer.  It is not a promoted
hydrogen or CFD result.
"""

from pathlib import Path
import argparse
import json
import math
import sys
import hashlib
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.addons.binary_integral_puff import BinaryPuffThermodynamics
from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from degali.addons.dense_prism_integral import integrate_inertial_prism
from degali.addons.finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud, initialise_regular_polygon_prism_cloud,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.transient_dense_gas_3d import (
    DenseGas3DReceptor,
    TransientDenseGas3DConfig,
    solve_transient_dense_gas_3d,
)
from degali.addons.transient_receptor import WindHistory
from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


def vol_percent_per_kg_m3(closure: BinaryIdealGasMixtureClosure) -> float:
    return (100.0 * 8.31446261815324 * closure.ambient_temperature_k
            / (closure.ambient_pressure_pa * closure.species_molar_mass_kg_mol))


def run(output: Path, *, cell_m: float = 5.0, dz_m: float = 2.0,
        duration_s: float = 600.0, obstacle: bool = True,
        source_sigma_m: float = 10.0, source_duration_s: float = 60.0,
        impermeable_ground: bool = True, source_mode: str = "scheduled",
        wind_angle_deg: float = 0.0, time_step_s: float = .5) -> dict:
    started = time.perf_counter()
    implementation_hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__), ROOT/'src/degali/addons/transient_dense_gas_3d.py')}
    if source_mode not in {"scheduled", "initial_cloud"}:
        raise ValueError("unsupported source mode")
    if duration_s < source_duration_s and source_mode == "scheduled":
        raise ValueError("scheduled run must cover the complete source duration")
    workbook = ROOT / "ti21_input_ascii.xls"
    dataset = read_thorney_ti021_workbook(workbook)
    closure = BinaryIdealGasMixtureClosure(
        "R12", .120913, 606.0748985077907,
        ambient_temperature_k=291.15, ambient_pressure_pa=100500.0,
    )
    thermo = BinaryPuffThermodynamics(closure)
    source = integrate_inertial_prism(
        RegularPolygonPrismInitialCloud((60.0, 320.0), 14.0, 0.0, 13.0, .32, 291.15),
        closure, duration_s=60.0, dt_s=.025, wind_xy_m_s=(3.4, 0.0),
        source_wind_scale=1.0, entrainment_shear=.1, entrainment_ambient=.1,
        drag_coefficient=.1, front_profile="grounded_head",
    )
    species_mass = float(source.states[-1].species_mass_kg)
    # Use the registered initial dense-cloud density for the buoyancy closure;
    # the 60 s handoff state is already nearly ambient after entrainment.
    source_density = float(source.states[0].density_kg_m3)
    # The density and species fraction must describe the SAME source state.
    source_fraction = float(source.states[0].species_mass_kg / source.states[0].total_mass_kg)

    points = np.array([(60.0 + s.local_y_m, 320.0 - s.local_x_m, s.height_m)
                       for s in dataset.dose_sensors], dtype=float)
    # The nearest sensor is downstream of the holder. A sensor-only box
    # started at x=60 m and cut the initial prism exactly in half.
    anchors = np.vstack((points[:, :2], (40., 300.), (80., 340.)))
    lower = np.min(anchors, axis=0) - 10.0
    upper = np.max(anchors, axis=0) + 10.0
    lower = np.floor(lower / cell_m) * cell_m
    upper = np.ceil(upper / cell_m) * cell_m
    lengths = upper - lower
    nx, ny = (int(round(value / cell_m)) for value in lengths)
    nz = int(round(20.0 / dz_m))
    x = lower[0] + (np.arange(nx) + .5) * cell_m
    y = lower[1] + (np.arange(ny) + .5) * cell_m
    z = (np.arange(nz) + .5) * dz_m
    solid = np.zeros((nz, ny, nx), dtype=bool)
    if obstacle:
        xx, yy = np.meshgrid(x, y, indexing="xy")
        radial = np.sqrt((xx - 60.0) ** 2 + (yy - 320.0) ** 2)
        shell = (xx >= 60.0) & (xx <= 110.0) & (radial >= 47.5) & (radial <= 52.5)
        solid = (z[:, None, None] < dataset.workbook_wall_height_m) & shell[None, :, :]

    initial = None
    if source_mode == "initial_cloud":
        cloud = RegularPolygonPrismInitialCloud(
            (60.0, 320.0), 14.0, 0.0, 13.0, .32, 291.15)
        initial = initialise_regular_polygon_prism_cloud(x, y, z, cloud=cloud, closure=closure)
        if abs(initial.relative_species_mass_error) > 1e-10:
            raise RuntimeError("initial cloud inventory is truncated by the domain")
        species_mass = float(initial.integrated_species_mass_kg)
        schedule = SourceRateSchedule((0.0, duration_s), (0.0, 0.0),
                                     source_id="TI021-initial-prism-only")
    elif source_duration_s < duration_s:
        schedule = SourceRateSchedule(
            (0.0, source_duration_s, duration_s),
            (species_mass / source_duration_s, 0.0, 0.0),
            source_id="TI021-R12-finite-prism-species-inventory")
    else:
        schedule = SourceRateSchedule((0.0, duration_s), (species_mass / duration_s, 0.0))
    # The scorer's local x-axis is the downwind axis; 270° meteorological
    # ``from`` direction transports toward positive x in this local frame.
    wind = WindHistory((0.0, duration_s), (3.4, 3.4),
                       (270.0 - wind_angle_deg, 270.0 - wind_angle_deg))
    config = TransientDenseGas3DConfig(
        domain_m=(float(lengths[0]), float(lengths[1]), 20.0),
        cells=(nx, ny, nz), duration_s=duration_s, time_step_s=time_step_s,
        source_schedule=schedule,
        wind_history=wind,
        source_position_m=(60.0 - lower[0], 320.0 - lower[1], 6.5),
        source_sigma_m=source_sigma_m,
        source_density_kg_m3=source_density,
        source_h2_mass_fraction=source_fraction,
        ambient_density_kg_m3=thermo.ambient_density,
        diffusivity_m2_s=1.0,
        buoyancy_length_m=2.0,
        maximum_buoyant_speed_m_s=1.0,
        prognostic_velocity=True,
        solid_mask_zyx=solid,
        impermeable_ground=impermeable_ground,
    )
    receptors = tuple(
        DenseGas3DReceptor(f"TI021-S{i + 1:03d}", tuple(point - np.array((lower[0], lower[1], 0.0))))
        for i, point in enumerate(points)
    )
    result = solve_transient_dense_gas_3d(config, receptors=receptors,
        initial_species_kg_m3=None if initial is None else initial.species_mass_concentration_kg_m3)
    # c_s = x_s p M_s / (R T). Air density converts to a mass mixing ratio,
    # not a mole/volume fraction; using it overstates R12 vol% by ~4.17.
    conversion = vol_percent_per_kg_m3(closure)
    traces = np.array([result.receptor_traces_kg_m3[r.label] * conversion for r in receptors])
    dose = np.trapezoid(traces, result.time_s, axis=1)
    observed = np.array([sensor.dose_percent_s for sensor in dataset.dose_sensors])
    positive = observed > 0.0
    full_horizon = bool(duration_s == 600.0)
    groups = {}
    radial_distance = np.linalg.norm(points[:,:2]-np.array((60.,320.)),axis=1)
    masks = {f'height_{height:g}m': points[:,2] == height for height in np.unique(points[:,2])}
    masks.update({'radius_le_70m': radial_distance <= 70., 'radius_gt_70m': radial_distance > 70.})
    for label, mask in masks.items():
        eligible = mask & positive
        groups[label] = {'positive_sensor_count': int(eligible.sum()),
            'dose_mape_percent': (float(np.mean(np.abs(dose[eligible]/observed[eligible]-1))*100)
                                  if full_horizon and eligible.any() else None)}
    report = {
        "schema": "degali.ti021-3d-obstacle-diagnostic.v2",
        "metric": "600s cumulative R12 dose MAPE; not instantaneous hydrogen concentration error",
        "full_600s_score": full_horizon,
        "status": "experimental_unqualified",
        "duration_s": duration_s,
        "obstacle": obstacle,
        "dose_mape_percent": (float(np.mean(np.abs(dose[positive] / observed[positive] - 1.0)) * 100.0)
                              if full_horizon else None),
        "dose_mae_vol_percent_s": (float(np.mean(np.abs(dose[positive] - observed[positive])))
                                  if full_horizon else None),
        "positive_observed_dose_sum_vol_percent_s": float(np.sum(observed[positive])),
        "positive_predicted_dose_sum_vol_percent_s": float(np.sum(dose[positive])),
        "normalized_mean_bias": (float(np.sum(dose[positive] - observed[positive]) / np.sum(observed[positive]))
                                 if full_horizon else None),
        "groups": groups,
        "sensor_count": int(len(observed)),
        "positive_observed_sensor_count": int(np.sum(positive)),
        "cell_size_m": cell_m,
        "vertical_cell_size_m": dz_m,
        "cells": [nx, ny, nz],
        "runtime_s": result.diagnostics.get("runtime_s"),
        "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        "solid_cell_count": result.diagnostics.get("solid_cell_count", 0),
        "source_species_mass_kg": species_mass,
        "source_density_kg_m3": source_density,
        "source_sigma_m": source_sigma_m,
        "source_duration_s": source_duration_s,
        "source_mode": source_mode,
        "initial_source_inventory_relative_error": (
            None if initial is None else initial.relative_species_mass_error),
        "domain_origin_xy_m": lower.tolist(),
        "source_mass_fraction": source_fraction,
        "vol_percent_per_kg_m3": conversion,
        "wind_angle_deg": wind_angle_deg,
        "impermeable_ground": impermeable_ground,
        "diagnostics": dict(result.diagnostics),
        "species_conservation_gate": bool(result.maximum_mass_residual_kg < 1e-7),
        "momentum_conservation_gate": None,
        "promotion_allowed": False,
        "sensor_rows": [
            {"label": r.label, "height_m": float(points[i,2]),
             "position_m": points[i].tolist(),
             "observed_dose_vol_percent_s": float(observed[i]),
             "predicted_dose_vol_percent_s": float(dose[i])}
            for i,r in enumerate(receptors)],
        "workbook_sha256": hashlib.sha256(workbook.read_bytes()).hexdigest(),
        "implementation_sha256": implementation_hashes,
        "reproduce_command": [sys.executable, str(Path(__file__)), '--output', str(output),
            '--cell-size-m', str(cell_m), '--vertical-cell-size-m', str(dz_m),
            '--duration-s', str(duration_s), '--source-sigma-m', str(source_sigma_m),
            '--source-duration-s', str(source_duration_s), '--source-mode', source_mode,
            '--wind-angle-deg', str(wind_angle_deg), '--time-step-s', str(time_step_s),
            *([] if obstacle else ['--no-obstacle']),
            *([] if impermeable_ground else ['--open-ground'])],
        "remaining_physics": [
            "scalar R12 surrogate, not a full binary momentum/energy mixture",
            "finite-thickness Cartesian fence mask",
            "no recoupled integral-layer handoff",
            "horizontal velocity reset to ambient; no pressure-mediated vertical bypass or separated wake",
            "no total-mixture or momentum conservation equation",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    np.savez_compressed(output.with_suffix('.npz'), time_s=result.time_s,
        concentration_vol_percent=traces, observed_dose_vol_percent_s=observed,
        predicted_dose_vol_percent_s=dose, total_mass_kg=result.total_mass_kg,
        outflow_mass_kg=result.outflow_mass_kg)
    report['runtime_including_output_s'] = time.perf_counter()-started
    output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='sensor_rows'}, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cell-size-m", type=float, default=5.0)
    parser.add_argument("--vertical-cell-size-m", type=float, default=2.0)
    parser.add_argument("--source-sigma-m", type=float, default=10.0)
    parser.add_argument("--source-duration-s", type=float, default=60.0)
    parser.add_argument("--no-obstacle", action="store_true")
    parser.add_argument("--open-ground", action="store_true")
    parser.add_argument("--source-mode", choices=("scheduled", "initial_cloud"), default="scheduled")
    parser.add_argument("--wind-angle-deg", type=float, default=0.0)
    parser.add_argument("--duration-s", type=float, default=600.0)
    parser.add_argument("--time-step-s", type=float, default=.5)
    args = parser.parse_args()
    run(args.output, cell_m=args.cell_size_m, dz_m=args.vertical_cell_size_m,
        obstacle=not args.no_obstacle, source_sigma_m=args.source_sigma_m,
        source_duration_s=args.source_duration_s, duration_s=args.duration_s,
        impermeable_ground=not args.open_ground, source_mode=args.source_mode,
        wind_angle_deg=args.wind_angle_deg, time_step_s=args.time_step_s)
