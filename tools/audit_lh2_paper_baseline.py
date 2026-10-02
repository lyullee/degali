"""Create an ignored, reproducible LH2 paper-validation audit bundle.

This program intentionally writes *derived* residuals and figures only to a
caller-selected directory.  It never copies a workbook, report, or a source
sensor table into the repository.  The frozen model path is the current
liquid-hydrogen fast path: corrected expanded source, at-sensor observation
operator, and free ground detachment.

The bundle has three distinct evidential strata:

* PRESLHY E3.5 individual concentration receptors (near field);
* the sealed E3.5 temperature-profile diagnostic (thermal residuals);
* DNV/FFI Spadeadam Tests 4 and 6 (30--100 m arc maxima).

It is deliberately not a parameter-estimation or curve-fitting program.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
E35 = ROOT / "reference" / "preslhy" / "e35_reduced.json"
THERMAL = ROOT / "reference" / "preslhy" / "model_thermal_profile_observation_2026-09-08.json"
SPADEADAM = ROOT / "reference" / "spadeadam"
HECHT = ROOT / "reference" / "lh2" / "hecht-panda-results-2026-09-04.json"

BASELINE = {
    "name": "degali-lh2-paper-baseline-2026-09-20",
    "corrections": True,
    "nearfield_observation_operator": "at-sensor, then arc maximum",
    "nearfield_source_rate": "flow_mean_gs",
    "ground_effect": False,
    "roughness_m": 0.001,
    "spadeadam_wind": "low mast as published table value",
    "spadeadam_wind_reference_height_m": 10.0,
    "fitted_parameter": False,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path.name}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _nearfield_trajectories() -> tuple[list[dict], dict[int, str]]:
    """Run the precise primary E3.5 baseline once per eligible trial."""
    from degali.validation import nearfield

    data = json.loads(E35.read_text(encoding="utf-8"))
    trajectories: list[dict] = []
    excluded: dict[int, str] = {}
    for trial in data["trials"]:
        number = int(trial["trial"])
        if trial["orientation"] != "horizontal":
            excluded[number] = "not horizontal"
            continue
        rate = float(trial["flow_mean_gs"]) / 1000.0
        if rate <= 0.005:
            excluded[number] = "no usable source flow"
            continue
        if nearfield._exit_ratio(trial, rate) < nearfield.MOMENTUM_RATIO:
            excluded[number] = "wind-steered source"
            continue
        jet, initial = nearfield.hydrogen_jet(
            rate=rate,
            diameter=float(trial["orifice_mm"]) / 1000.0,
            wind=float(trial["wind_ms"]),
            height=float(trial["release_height_m"]),
            ambient_temperature=float(trial["T_C"]) + 273.15,
            relative_humidity=float(trial["RH_pct"]),
            storage_pressure_barg=float(trial["tanker_barg"]),
            storage_temperature=nearfield._liquid_source_temperature(
                trial, "tank_saturation"
            ),
            wind_reference_height=float(trial["wind_ref_m"]),
            corrections=True,
            ground_effect=False,
        )
        trajectory = nearfield.Trajectory(
            jet.th.table,
            jet.run(initial, distmx=nearfield.STEP, smax=nearfield.REACH).rows,
        )
        if not trajectory.ok:
            excluded[number] = "plume integration failed"
            continue
        trajectories.append({"trial": trial, "trajectory": trajectory})
    return trajectories, excluded


def concentration_rows() -> tuple[list[dict], dict[int, str]]:
    rows: list[dict] = []
    trajectories, excluded = _nearfield_trajectories()
    for item in trajectories:
        trial, trajectory = item["trial"], item["trajectory"]
        number = int(trial["trial"])
        for sensor in trial["sensors"]:
            x, y, z = (float(sensor[name]) for name in ("x", "y", "z"))
            observed = float(sensor["peak"])
            predicted = trajectory.concentration_at(x, y, z)
            state = trajectory.at(x)
            supported = bool(state is not None and predicted > 0.0)
            row = {
                "trial": number,
                "serial": sensor["serial"],
                "x_m": x,
                "y_m": y,
                "z_m": z,
                "observed_peak_vol_pct": observed,
                "predicted_vol_pct": predicted,
                "plume_centre_m": state.z if state is not None else None,
                "sigma_z_m": state.sz if state is not None else None,
                "supported": supported,
                "log10_observed_over_predicted": (
                    math.log10(observed / predicted)
                    if observed > 0.0 and predicted > 0.0 else None
                ),
            }
            rows.append(row)
    return rows, excluded


def temperature_rows() -> list[dict]:
    """Flatten the sealed observed-profile comparison without reintegrating it."""
    data = json.loads(THERMAL.read_text(encoding="utf-8"))
    rows: list[dict] = []
    for variant, profiles in data["model_rows"].items():
        for profile in profiles:
            projection = profile["projection"]
            comparison = profile["comparison"]["thermal_deficit"]
            for relative, absolute, predicted in zip(
                projection["relative_height_m"],
                projection["absolute_height_m"],
                projection["thermal_deficit_K"],
            ):
                rows.append({
                    "variant": variant,
                    "trial": int(profile["trial"]),
                    "station_m": float(profile["station_m"]),
                    "relative_height_m": float(relative),
                    "absolute_height_m": float(absolute),
                    "predicted_thermal_deficit_K": float(predicted),
                    "observed_median_centre_deficit_K": float(
                        comparison["centre_amplitude"]["observed_median"]
                    ),
                    "centre_deficit_residual_K": float(
                        comparison["centre_amplitude"]["difference"]
                    ),
                    "thermal_width_residual_m": float(
                        comparison["width"]["difference"]
                    ),
                    "scope": "sealed diagnostic; not primary trajectory re-integration",
                })
    return rows


def hecht_panda_rows() -> list[dict]:
    """Flatten the public aggregate cryogenic-jet benchmark.

    This is kept as an external benchmark boundary because its published
    condition membership is unresolved; it is not merged into the E3.5 score.
    """
    data = json.loads(HECHT.read_text(encoding="utf-8"))
    metrics = (
        ("centerline_mass_slope", "centreline mass decay"),
        ("mass_half_width_slope_mm", "mass half width"),
        ("centerline_temperature_slope", "centreline temperature decay"),
        ("temperature_half_width_slope_mm", "temperature half width"),
    )
    return [{
        "metric": label,
        "observed": float(data["observed"][key]),
        "current_lh2": float(data["degali_variants"]["current_lh2"][key]),
        "conserved_energy": float(data["degali_conserved_energy"][key]),
        "current_over_observed": float(
            data["degali_variants"]["current_lh2"][key]
        ) / float(data["observed"][key]),
        "benchmark_scope": "aggregate slope; published condition membership unresolved",
    } for key, label in metrics]


def concentration_arc_statistics(rows: list[dict]) -> dict:
    """Apply the frozen at-sensor/arc-maximum operator to map rows."""
    from degali.validation.statistics import statistics

    arcs: dict[tuple[int, float], list[dict]] = {}
    for row in rows:
        arcs.setdefault((int(row["trial"]), float(row["x_m"])), []).append(row)
    paired = []
    for (trial, distance), group in sorted(arcs.items()):
        observed = max(float(row["observed_peak_vol_pct"]) for row in group)
        predicted = max(float(row["predicted_vol_pct"]) for row in group)
        if observed > 0.0 and predicted > 0.0:
            paired.append({
                "trial": trial,
                "x_m": distance,
                "observed_arc_max_vol_pct": observed,
                "predicted_arc_max_vol_pct": predicted,
            })
    result = statistics(
        [row["observed_arc_max_vol_pct"] for row in paired],
        [row["predicted_arc_max_vol_pct"] for row in paired],
    )
    return {"arc_count": len(paired), "statistics": dataclasses.asdict(result)}


def _ffi_trajectory(trial, *, wind: float, ground_effect: bool = False,
                    coefficients: dict | None = None):
    from degali.validation import nearfield, spadeadam

    jet, initial = nearfield.hydrogen_jet(
        rate=trial.rate,
        diameter=trial.orifice,
        wind=wind,
        height=spadeadam.RELEASE_HEIGHT,
        ambient_temperature=spadeadam.AMBIENT_TEMPERATURE,
        relative_humidity=spadeadam.RELATIVE_HUMIDITY,
        storage_pressure_barg=trial.line_pressure,
        roughness=BASELINE["roughness_m"],
        wind_reference_height=BASELINE["spadeadam_wind_reference_height_m"],
        corrections=True,
        ground_effect=ground_effect,
    )
    if coefficients:
        jet.k = dataclasses.replace(jet.k, **coefficients)
    return nearfield.Trajectory(
        jet.th.table,
        jet.run(initial, distmx=nearfield.STEP, smax=250.0).rows,
    )


def _source_coordinates(reading, wind_from: float) -> tuple[float, float]:
    """Project a DNV sensor to downwind/crosswind coordinates from the nozzle."""
    bearing = math.radians(reading.bearing)
    east = reading.radius * math.sin(bearing)
    north = reading.radius * math.cos(bearing)
    direction = math.radians(wind_from)
    down_east, down_north = -math.sin(direction), -math.cos(direction)
    right_east, right_north = down_north, -down_east
    return (
        east * down_east + north * down_north,
        east * right_east + north * right_north,
    )


def ffi_rows() -> tuple[list[dict], list[dict], dict]:
    from degali.validation import spadeadam
    from degali.validation.statistics import statistics

    arcs: list[dict] = []
    sensor_rows: list[dict] = []
    decomposition: dict = {"test_6_30_m": {}}
    trials = {trial.test: trial for trial in spadeadam.load(SPADEADAM)}

    for test in (4, 6):
        trial = trials[test]
        baseline = _ffi_trajectory(trial, wind=trial.wind_low)
        for radius in trial.radii():
            observed = trial.arc(radius)
            state = baseline.at(radius)
            if state is None or not observed:
                continue
            predicted = {
                height: baseline.concentration_at(radius, 0.0, height)
                for height in observed
            }
            arcs.append({
                "test": test,
                "radius_m": radius,
                "wind_low_m_s": trial.wind_low,
                "wind_high_m_s": trial.wind_high,
                "observed_arc_max_vol_pct": max(observed.values()),
                "predicted_arc_max_vol_pct": max(predicted.values()),
                "observed_over_predicted": max(observed.values()) / max(predicted.values()),
                "centre_m": state.z,
                "sigma_z_m": state.sz,
                "observed_flatness": min(observed.values()) / max(observed.values()),
                "predicted_flatness": min(predicted.values()) / max(predicted.values()),
            })
        for reading in trial.readings:
            if reading.over_range:
                continue
            x, y = _source_coordinates(reading, trial.wind_direction)
            state = baseline.at(x) if x > 0.0 else None
            prediction = (
                baseline.concentration_at(x, y, reading.height)
                if state is not None else 0.0
            )
            sensor_rows.append({
                "test": test,
                "sensor": reading.sensor,
                "radius_m": reading.radius,
                "height_m": reading.height,
                "bearing_deg": reading.bearing,
                "downwind_m": x,
                "crosswind_m": y,
                "observed_peak_vol_pct": reading.peak,
                "predicted_vol_pct": prediction,
                "downwind_supported": state is not None,
            })

    six = trials[6]
    observed = six.arc(30.0)
    observed_max = max(observed.values())
    transport_options = {
        "baseline_free_detachment": (False, None),
        "forced_ground_contact": (True, None),
        "shape_drag_2_diagnostic": (False, {"shape_drag": 2.0}),
        "vertical_shear_diagnostic": (False, {"vertical_shear": True}),
        "liftoff_ri_30_diagnostic": (False, {"liftoff_richardson": 30.0}),
    }
    for label, (ground, coefficients) in transport_options.items():
        trajectory = _ffi_trajectory(
            six, wind=six.wind_low, ground_effect=ground,
            coefficients=coefficients,
        )
        state = trajectory.at(30.0)
        maximum = max(trajectory.concentration_at(30.0, 0.0, z) for z in observed)
        decomposition["test_6_30_m"].setdefault("transport", {})[label] = {
            "centre_m": state.z,
            "predicted_arc_max_vol_pct": maximum,
            "observed_over_predicted": observed_max / maximum,
        }
    for label, wind in {
        "low_mast": six.wind_low,
        "mean_mast": 0.5 * (six.wind_low + six.wind_high),
        "high_mast": six.wind_high,
    }.items():
        trajectory = _ffi_trajectory(six, wind=wind)
        state = trajectory.at(30.0)
        maximum = max(trajectory.concentration_at(30.0, 0.0, z) for z in observed)
        decomposition["test_6_30_m"].setdefault("wind", {})[label] = {
            "wind_m_s": wind,
            "centre_m": state.z,
            "predicted_arc_max_vol_pct": maximum,
            "observed_over_predicted": observed_max / maximum,
        }
    direct = [row for row in sensor_rows if row["test"] == 6 and abs(row["radius_m"] - 30.0) < 0.1]
    decomposition["test_6_30_m"]["observation_operator"] = {
        "baseline": "arc maximum versus model centreline at the same heights",
        "reported_arc_max_vol_pct": observed_max,
        "physical_sensor_rows": len(direct),
        "downwind_supported_rows": sum(row["downwind_supported"] for row in direct),
        "maximum_model_value_at_physical_sensor_positions_vol_pct": max(
            row["predicted_vol_pct"] for row in direct
        ),
        "interpretation": (
            "A steady centreline is not a time-resolved sampling operator; "
            "the physical-sensor projection is a geometry sensitivity only."
        ),
    }
    decomposition["test_6_30_m"]["source"] = {
        "reported_mass_flow_kg_s": six.rate,
        "reported_P04_barg": six.line_pressure,
        "orifice_m": six.orifice,
        "status": "unidentified from the public table",
        "interpretation": (
            "The table has no independently measured atmospheric two-phase "
            "source state or synchronized source time history. No source "
            "sensitivity is fabricated or fitted."
        ),
    }
    stats = statistics(
        [row["observed_arc_max_vol_pct"] for row in arcs],
        [row["predicted_arc_max_vol_pct"] for row in arcs],
    )
    return arcs, sensor_rows, {
        "baseline_arc_statistics": dataclasses.asdict(stats),
        "test_6_30_m": decomposition["test_6_30_m"],
    }


def _colour(value: float, scale: float) -> str:
    """Blue for model-high, red for model-low; no plotting dependency."""
    fraction = min(1.0, abs(value) / scale)
    shade = int(235 - 140 * fraction)
    return f"rgb({235 if value >= 0 else shade},{shade},{shade if value >= 0 else 235})"


def _scatter_svg(path: Path, points: list[tuple[float, float, float]], *,
                 title: str, xlabel: str, ylabel: str, scale: float) -> None:
    """Small self-contained residual map for environments without matplotlib."""
    width, height, left, bottom = 900, 540, 85, 70
    plot_width, plot_height = width - left - 35, height - bottom - 70
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    dx = max(xmax - xmin, 1.0e-9)
    dy = max(ymax - ymin, 1.0e-9)
    circles = []
    for x, y, residual in points:
        px = left + (x - xmin) / dx * plot_width
        py = height - bottom - (y - ymin) / dy * plot_height
        circles.append(
            f'<circle cx="{px:.2f}" cy="{py:.2f}" r="5" fill="{_colour(residual, scale)}" />'
        )
    text = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;fill:#202020}}.small{{font-size:14px}}.title{{font-size:20px;font-weight:bold}}</style>
<text x="{left}" y="28" class="title">{title}</text><line x1="{left}" y1="{height-bottom}" x2="{width-35}" y2="{height-bottom}" stroke="#222"/><line x1="{left}" y1="{height-bottom}" x2="{left}" y2="70" stroke="#222"/>
<text x="{left}" y="{height-22}" class="small">{xlabel}: {xmin:.3g} to {xmax:.3g}</text><text x="{left}" y="50" class="small">{ylabel}: {ymin:.3g} to {ymax:.3g}</text>
<text x="{width-310}" y="50" class="small">red = observed higher; blue = model higher</text>{''.join(circles)}</svg>'''
    path.write_text(text, encoding="utf-8")


def figures(out: Path, concentration: list[dict], thermal: list[dict], arcs: list[dict], hecht: list[dict]) -> list[str]:
    """Render local, derived SVG maps; no source files are copied."""
    made: list[str] = []
    usable = [r for r in concentration if r["log10_observed_over_predicted"] is not None]
    path = out / "preslhy-concentration-residual-map.svg"
    _scatter_svg(path, [
        (r["x_m"], r["z_m"], r["log10_observed_over_predicted"]) for r in usable
    ], title="PRESLHY E3.5 concentration residual map", xlabel="downwind distance (m)",
       ylabel="sensor height (m)", scale=1.0)
    made.append(path.name)

    profile = [r for r in thermal if r["variant"] == "density_profile_control"]
    path = out / "preslhy-thermal-residual-map.svg"
    _scatter_svg(path, [
        (r["station_m"], r["absolute_height_m"], r["centre_deficit_residual_K"])
        for r in profile
    ], title="E3.5 sealed thermal-profile residual diagnostic", xlabel="station (m)",
       ylabel="height (m)", scale=100.0)
    made.append(path.name)

    path = out / "ffi-arc-residual-map.svg"
    _scatter_svg(path, [
        (r["radius_m"], float(r["test"]), math.log10(r["observed_over_predicted"]))
        for r in arcs
    ], title="FFI/DNV horizontal-release arc residual map", xlabel="arc radius (m)",
       ylabel="test identifier", scale=1.0)
    made.append(path.name)
    path = out / "hecht-panda-slope-residual-map.svg"
    _scatter_svg(path, [
        (float(index), math.log10(row["current_over_observed"]),
         math.log10(row["current_over_observed"]))
        for index, row in enumerate(hecht, start=1)
    ], title="Hecht–Panda cryogenic-jet slope residual boundary", xlabel="metric index",
       ylabel="log10(current / observed)", scale=1.0)
    made.append(path.name)
    return made


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        raise FileExistsError(f"output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    concentration, excluded = concentration_rows()
    thermal = temperature_rows()
    hecht = hecht_panda_rows()
    arcs, ffi_sensors, ffi = ffi_rows()
    write_csv(out / "preslhy_concentration_residuals.csv", concentration)
    write_csv(out / "preslhy_thermal_residuals.csv", thermal)
    write_csv(out / "ffi_arc_residuals.csv", arcs)
    write_csv(out / "ffi_physical_sensor_projection.csv", ffi_sensors)
    write_csv(out / "hecht-panda_slope_residuals.csv", hecht)
    images = figures(out, concentration, thermal, arcs, hecht)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "baseline": BASELINE,
        "inputs_sha256": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
            for path in (
                E35, THERMAL, SPADEADAM / "conditions.csv", SPADEADAM / "sensors.csv",
                HECHT,
                ROOT / "src/degali/validation/nearfield.py",
                ROOT / "src/degali/validation/spadeadam.py",
                Path(__file__),
            )
        },
        "concentration_residual_map": {
            "rows": len(concentration),
            "eligible_trials": sorted({row["trial"] for row in concentration}),
            "excluded_trials": excluded,
            "supported_receptors": sum(row["supported"] for row in concentration),
            "frozen_arc_operator": concentration_arc_statistics(concentration),
        },
        "thermal_residual_map": {
            "rows": len(thermal),
            "scope": "sealed diagnostic; does not re-integrate a primary trajectory",
        },
        "hecht_panda_external_benchmark": {
            "rows": len(hecht),
            "scope": "aggregate slope boundary; not merged into primary score",
            "rows_file": "hecht-panda_slope_residuals.csv",
        },
        "ffi": ffi,
        "figures": images,
        "no_parameter_fitting": True,
        "no_external_source_file_redistributed": True,
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output_dir": str(out),
        "concentration_rows": len(concentration),
        "thermal_rows": len(thermal),
        "ffi_arc_rows": len(arcs),
        "ffi_statistics": ffi["baseline_arc_statistics"],
    }, indent=2))


if __name__ == "__main__":
    main()
