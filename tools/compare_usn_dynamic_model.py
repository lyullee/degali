"""Compare the DEGALI transient 3-D model with public USN Exp20 data.

This is an uncalibrated component-level enclosure diagnostic.  It uses the
measured mass-flow and ventilation histories, published container geometry,
and published sensor coordinates.  It does not model the release table,
nozzle momentum, or a resolved ventilation inlet/open-door flow field.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.transient_dense_gas_3d import (
    DenseGas3DReceptor,
    TransientDenseGas3DConfig,
    solve_transient_dense_gas_3d,
)
from degali.addons.transient_receptor import WindHistory
from degali.validation.usn_carpark import (
    AMBIENT_HYDROGEN_DENSITY_KG_M3,
    read_usn_carpark_mat,
    threshold_metrics,
    well_mixed_concentration_history,
)


CONTAINER_M = (11.885, 2.240, 2.285)
SOURCE_POSITION_M = (5.000, 1.120, 0.250)
AMBIENT_AIR_DENSITY_KG_M3 = 1.225
SCORED_WINDOW_S = 900.0


def _event_nodes(duration_s: float) -> np.ndarray:
    """Retain early blowdown detail and use 5 s nodes thereafter."""

    early = np.arange(0.0, min(20.0, duration_s) + 1.0e-9, 1.0)
    later = np.arange(25.0, duration_s, 5.0)
    return np.unique(np.concatenate((early, later, np.array([duration_s]))))


def _model_position(raw_mm: tuple[float, float, float]) -> tuple[float, float, float]:
    """Map native USN coordinates to x-down-container, y-across, z-up."""

    native_x, native_height, native_lateral = raw_mm
    return (-native_x / 1000.0, 1.120 + native_lateral / 1000.0, native_height / 1000.0)


def _threshold(trace: np.ndarray, time_s: np.ndarray, threshold: float = 4.0) -> dict[str, float | None | bool]:
    above = np.asarray(trace) >= threshold
    hits = np.flatnonzero(above)
    duration = float(np.sum(np.diff(time_s) * above[:-1]))
    return {
        "arrival_s": None if not hits.size else float(time_s[hits[0]]),
        "final_exceedance_s": None if not hits.size else float(time_s[hits[-1]]),
        "duration_s": duration,
        "right_censored": bool(above[-1]),
    }


def _finite_correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    keep = np.isfinite(left) & np.isfinite(right)
    if np.count_nonzero(keep) < 3:
        return None
    a, b = left[keep], right[keep]
    if np.std(a) <= 0.0 or np.std(b) <= 0.0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def _score(model: np.ndarray, observed: np.ndarray) -> dict[str, float | None]:
    residual = model - observed
    return {
        "mae_vol_percent_points": float(np.mean(np.abs(residual))),
        "rmse_vol_percent_points": float(np.sqrt(np.mean(residual**2))),
        "bias_vol_percent_points": float(np.mean(residual)),
        "pearson_r": _finite_correlation(model.ravel(), observed.ravel()),
        "within_2_vol_percent_points_fraction": float(np.mean(np.abs(residual) <= 2.0)),
    }


def _jsonable(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"cannot serialize {type(value)!r}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "reference/datasets/usn-carpark-14405903/raw/Exp20.mat",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "outputs/usn-exp20-dynamic-comparison-2026-10-10",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    run = read_usn_carpark_mat(args.input)
    sensors = tuple(
        sensor for sensor in run.room_cloud_sensors
        if sensor.position_raw_mm is not None
    )
    invalid_sensors = tuple(
        sensor.sensor_id for sensor in sensors
        if not np.any(np.abs(sensor.concentration_vol_percent) > 1.0e-12)
    )
    valid_sensors = tuple(
        sensor for sensor in sensors if sensor.sensor_id not in invalid_sensors
    )

    duration_s = min(
        SCORED_WINDOW_S,
        float(run.source_time_s[-1] - run.source_t0_s),
        float(valid_sensors[0].time_s[-1] - run.source_t0_s),
    )
    nodes = _event_nodes(duration_s)
    raw_event_time = run.source_time_s - run.source_t0_s
    source_rate = np.interp(nodes, raw_event_time, np.maximum(run.mass_flow_g_s, 0.0)) / 1000.0
    source_rate[-1] = 0.0
    source_schedule = SourceRateSchedule(
        tuple(nodes), tuple(source_rate), source_id=f"{run.experiment_id}-measured-MFM",
        rate_operator="linear",
    )

    cross_section_m2 = CONTAINER_M[1] * CONTAINER_M[2]
    ventilation_flow = np.interp(nodes, raw_event_time, run.ventilation_flow_m3_h)
    ventilation_speed = np.maximum(ventilation_flow / 3600.0 / cross_section_m2, 1.0e-6)
    wind_history = WindHistory(
        tuple(nodes), tuple(ventilation_speed), tuple(np.full(nodes.size, 270.0)),
    )
    receptors = tuple(
        DenseGas3DReceptor(sensor.sensor_id, _model_position(sensor.position_raw_mm))
        for sensor in sensors
    )
    config = TransientDenseGas3DConfig(
        domain_m=CONTAINER_M,
        cells=(36, 8, 8),
        duration_s=duration_s,
        time_step_s=0.05,
        source_schedule=source_schedule,
        wind_history=wind_history,
        source_position_m=SOURCE_POSITION_M,
        source_sigma_m=0.25,
        source_density_kg_m3=AMBIENT_HYDROGEN_DENSITY_KG_M3,
        source_h2_mass_fraction=1.0,
        ambient_density_kg_m3=AMBIENT_AIR_DENSITY_KG_M3,
        diffusivity_m2_s=0.05,
        buoyancy_length_m=1.0,
        maximum_buoyant_speed_m_s=4.0,
        prognostic_velocity=False,
        store_fields=False,
        boundary_modes=("open", "closed", "closed"),
    )
    result = solve_transient_dense_gas_3d(config, receptors=receptors)
    model_time = result.time_s
    model = np.vstack([
        100.0 * result.receptor_traces_kg_m3[sensor.sensor_id]
        / AMBIENT_HYDROGEN_DENSITY_KG_M3
        for sensor in valid_sensors
    ])
    observed = np.vstack([
        np.interp(
            model_time + run.source_t0_s,
            sensor.time_s,
            sensor.concentration_vol_percent,
        )
        for sensor in valid_sensors
    ])

    observed_mean = np.mean(observed, axis=0)
    model_mean = np.mean(model, axis=0)
    observed_envelope = np.max(observed, axis=0)
    model_envelope = np.max(model, axis=0)
    zero_d_clock, zero_d = well_mixed_concentration_history(run)
    zero_d_model = np.interp(model_time + run.source_t0_s, zero_d_clock, zero_d)

    raw_keep = (raw_event_time >= 0.0) & (raw_event_time <= duration_s)
    raw_mass_kg = float(np.trapezoid(
        np.maximum(run.mass_flow_g_s[raw_keep], 0.0) / 1000.0,
        raw_event_time[raw_keep],
    ))
    published_threshold = threshold_metrics(run, sensor_scope="room_cloud")
    aggregate_score = _score(model, observed)
    mean_score = _score(model_mean, observed_mean)
    zero_d_score = _score(zero_d_model, observed_mean)

    sensor_rows = []
    for index, sensor in enumerate(valid_sensors):
        score = _score(model[index], observed[index])
        observed_peak_index = int(np.argmax(observed[index]))
        model_peak_index = int(np.argmax(model[index]))
        sensor_rows.append({
            "sensor_id": sensor.sensor_id,
            "x_m": receptors[[item.label for item in receptors].index(sensor.sensor_id)].position_m[0],
            "y_m": receptors[[item.label for item in receptors].index(sensor.sensor_id)].position_m[1],
            "z_m": receptors[[item.label for item in receptors].index(sensor.sensor_id)].position_m[2],
            **score,
            "observed_peak_vol_percent": float(observed[index, observed_peak_index]),
            "model_peak_vol_percent": float(model[index, model_peak_index]),
            "observed_peak_time_s": float(model_time[observed_peak_index]),
            "model_peak_time_s": float(model_time[model_peak_index]),
        })

    summary = {
        "schema": "degali.usn-dynamic-comparison.v1",
        "experiment_id": run.experiment_id,
        "comparison_verdict": "failed_physicality_gate_not_similar",
        "raw_file": str(run.source_path),
        "raw_sha256": run.raw_sha256,
        "comparison_scope": "GH2 confined dynamic transport component evidence; not LH2 end-to-end validation",
        "calibration": "none against concentration sensors",
        "scored_window_after_source_s": [0.0, duration_s],
        "sensor_scope": {
            "published_room_array_count": len(sensors),
            "scored_count": len(valid_sensors),
            "excluded_flatline_ids": list(invalid_sensors),
        },
        "geometry_and_mapping": {
            "container_m": CONTAINER_M,
            "source_position_m": SOURCE_POSITION_M,
            "native_to_model": "x=-native_x; y=1.120+native_lateral; z=native_height (mm converted to m)",
            "ventilation_direction": "+x from circular inlet wall toward open doors",
        },
        "model_configuration": {
            "cells": config.cells,
            "time_step_s": config.time_step_s,
            "boundary_modes_xyz": config.boundary_modes,
            "source_sigma_m": config.source_sigma_m,
            "effective_diffusivity_m2_s": config.diffusivity_m2_s,
            "buoyancy_length_m": config.buoyancy_length_m,
            "prognostic_velocity": config.prognostic_velocity,
            "ventilation_speed_definition": "measured m3/h divided by full container cross-section",
            "unresolved_features": [
                "release table and jet impingement",
                "nozzle momentum and compressible near field",
                "resolved circular ventilation inlet",
                "resolved open-door pressure/velocity field",
            ],
        },
        "source_resampling": {
            "raw_integrated_mass_kg": raw_mass_kg,
            "schedule_integrated_mass_kg": source_schedule.released_mass_kg,
            "relative_mass_difference_percent": 100.0 * (source_schedule.released_mass_kg - raw_mass_kg) / raw_mass_kg,
        },
        "dynamic_all_sensor_time_score": aggregate_score,
        "dynamic_spatial_mean_score": mean_score,
        "zero_d_spatial_mean_score_same_window": zero_d_score,
        "observed_envelope": {
            "peak_vol_percent": float(np.max(observed_envelope)),
            **_threshold(observed_envelope, model_time),
        },
        "dynamic_envelope": {
            "peak_vol_percent": float(np.max(model_envelope)),
            **_threshold(model_envelope, model_time),
        },
        "zero_d": {
            "peak_vol_percent": float(np.max(zero_d_model)),
            **_threshold(zero_d_model, model_time),
        },
        "published_room_array_threshold_check": {
            "peak_vol_percent": published_threshold.observed_peak_vol_percent,
            "arrival_s": published_threshold.first_arrival_after_source_s,
            "duration_s": published_threshold.exceedance_duration_s,
        },
        "physicality_checks": {
            "maximum_predicted_sensor_vol_percent": float(np.max(model)),
            "sensor_values_above_100_vol_percent": int(np.count_nonzero(model > 100.0)),
            "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        },
        "solver_diagnostics": dict(result.diagnostics),
        "interpretation_guardrail": (
            "A good score would support only transient GH2 transport timing/concentration in this enclosure. "
            "A poor score diagnoses missing enclosure/jet/obstacle closure and does not invalidate DEGALI's other validated paths."
        ),
        "similarity_interpretation": (
            "The enclosure run is not quantitatively similar. Predicted sensor concentration exceeds the "
            "physical 100 vol% ceiling, so the reported errors diagnose an inapplicable closure rather than "
            "constituting a validation metric."
        ),
    }

    (args.output_dir / "metrics.json").write_text(
        json.dumps(summary, indent=2, default=_jsonable), encoding="utf-8",
    )
    with (args.output_dir / "per_sensor_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sensor_rows[0]))
        writer.writeheader()
        writer.writerows(sensor_rows)

    sample_time = np.arange(0.0, duration_s + 1.0e-9, 1.0)
    with (args.output_dir / "timeseries_1s.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "time_after_source_s", "observed_spatial_mean_vol_percent",
            "dynamic_spatial_mean_vol_percent", "observed_sensor_envelope_vol_percent",
            "dynamic_sensor_envelope_vol_percent", "zero_d_well_mixed_vol_percent",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for time_value in sample_time:
            writer.writerow({
                "time_after_source_s": time_value,
                "observed_spatial_mean_vol_percent": np.interp(time_value, model_time, observed_mean),
                "dynamic_spatial_mean_vol_percent": np.interp(time_value, model_time, model_mean),
                "observed_sensor_envelope_vol_percent": np.interp(time_value, model_time, observed_envelope),
                "dynamic_sensor_envelope_vol_percent": np.interp(time_value, model_time, model_envelope),
                "zero_d_well_mixed_vol_percent": np.interp(time_value, model_time, zero_d_model),
            })
    np.savez_compressed(
        args.output_dir / "sensor_traces.npz",
        time_after_source_s=model_time,
        sensor_ids=np.asarray([sensor.sensor_id for sensor in valid_sensors]),
        observed_vol_percent=observed,
        dynamic_vol_percent=model,
        zero_d_vol_percent=zero_d_model,
    )

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(model_time, observed_mean, color="black", linewidth=2.2, label="Measured spatial mean")
    ax.plot(model_time, model_mean, color="#0072B2", linewidth=2.0, label="DEGALI dynamic 3-D")
    ax.plot(model_time, zero_d_model, color="#D55E00", linewidth=1.7, linestyle="--", label="0-D well mixed")
    ax.axhline(4.0, color="#CC0000", linestyle=":", linewidth=1.4, label="4 vol%")
    ax.set(xlabel="Time after source start (s)", ylabel="H$_2$ concentration (vol%)", title="A. Spatial mean history")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(model_time, observed_envelope, color="black", linewidth=2.2, label="Measured sensor envelope")
    ax.plot(model_time, model_envelope, color="#0072B2", linewidth=2.0, label="DEGALI dynamic envelope")
    ax.axhline(4.0, color="#CC0000", linestyle=":", linewidth=1.4)
    ax.set(xlabel="Time after source start (s)", ylabel="H$_2$ concentration (vol%)", title="B. Any-sensor 4 vol% operator")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    observed_peaks = np.asarray([row["observed_peak_vol_percent"] for row in sensor_rows])
    model_peaks = np.asarray([row["model_peak_vol_percent"] for row in sensor_rows])
    ax.scatter(observed_peaks, model_peaks, color="#0072B2", s=42, alpha=0.85)
    limit = max(10.0, float(np.max(observed_peaks)), float(np.max(model_peaks))) * 1.05
    ax.plot([0.0, limit], [0.0, limit], color="black", linestyle="--", linewidth=1.2, label="1:1")
    for row in sensor_rows:
        ax.annotate(row["sensor_id"], (row["observed_peak_vol_percent"], row["model_peak_vol_percent"]), fontsize=6, alpha=0.75)
    ax.set(xlim=(0.0, limit), ylim=(0.0, limit), xlabel="Measured sensor peak (vol%)", ylabel="Dynamic-model sensor peak (vol%)", title="C. Sensor peak agreement")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.axis("off")
    dynamic_threshold = summary["dynamic_envelope"]
    measured_threshold = summary["observed_envelope"]
    text = (
        "USN Exp20 — uncalibrated dynamic comparison\n"
        "PHYSICALITY GATE: FAIL — NOT QUANTITATIVELY SIMILAR\n\n"
        f"Scored sensors: {len(valid_sensors)} (flatline excluded: {', '.join(invalid_sensors)})\n"
        f"All sensor-time MAE / RMSE: {aggregate_score['mae_vol_percent_points']:.2f} / {aggregate_score['rmse_vol_percent_points']:.2f} vol%p\n"
        f"All sensor-time Pearson r: {aggregate_score['pearson_r']:.3f}\n"
        f"Spatial-mean MAE / RMSE: {mean_score['mae_vol_percent_points']:.2f} / {mean_score['rmse_vol_percent_points']:.2f} vol%p\n"
        f"0-D mean MAE / RMSE: {zero_d_score['mae_vol_percent_points']:.2f} / {zero_d_score['rmse_vol_percent_points']:.2f} vol%p\n\n"
        f"Measured / dynamic peak: {measured_threshold['peak_vol_percent']:.2f} / {dynamic_threshold['peak_vol_percent']:.2f} vol%\n"
        f"Measured / dynamic 4% arrival: {measured_threshold['arrival_s']:.1f} / "
        f"{dynamic_threshold['arrival_s']:.1f} s\n"
        f"Measured / dynamic 4% duration: {measured_threshold['duration_s']:.1f} / {dynamic_threshold['duration_s']:.1f} s\n\n"
        "Scope: GH2 enclosure component evidence only\n"
        "Missing: table/jet impingement and resolved inlet/door flow"
    )
    ax.text(0.02, 0.98, text, va="top", ha="left", fontsize=10.5, family="monospace")
    fig.suptitle("DEGALI transient 3-D vs public USN container experiment", fontsize=15, fontweight="bold")
    fig.savefig(args.output_dir / "usn_exp20_dynamic_comparison_fullscale.png", dpi=220)
    axes[0, 0].set_ylim(-0.5, 22.0)
    axes[0, 0].text(
        0.98, 0.95,
        f"dynamic mean max = {float(np.max(model_mean)):.1f} vol% (off scale)",
        transform=axes[0, 0].transAxes, ha="right", va="top", color="#0072B2",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#0072B2"},
    )
    axes[0, 1].set_ylim(-0.5, 20.0)
    axes[0, 1].text(
        0.98, 0.95,
        f"dynamic envelope max = {float(np.max(model_envelope)):.1f} vol% (off scale)",
        transform=axes[0, 1].transAxes, ha="right", va="top", color="#0072B2",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#0072B2"},
    )
    axes[1, 0].set_yscale("symlog", linthresh=1.0)
    axes[1, 0].set_xlim(6.0, 9.2)
    axes[1, 0].set_ylim(0.0, limit)
    fig.savefig(args.output_dir / "usn_exp20_dynamic_comparison.png", dpi=220)
    plt.close(fig)

    report = f"""# USN Exp20 dynamic comparison

This is a concentration-sensor comparison on the measured event clock. It is
not an LH2 end-to-end validation and no concentration data were used to tune
the model.

## Main result

**Verdict: failed physicality gate; the current 3-D enclosure closure is not
quantitatively similar to the experiment.** It predicts values above the
physical 100 vol% ceiling, so its error values are diagnostic rather than a
valid performance claim.

- All sensor-time MAE: {aggregate_score['mae_vol_percent_points']:.3f} vol%p
- All sensor-time RMSE: {aggregate_score['rmse_vol_percent_points']:.3f} vol%p
- All sensor-time Pearson r: {aggregate_score['pearson_r']:.4f}
- Spatial-mean MAE / RMSE: {mean_score['mae_vol_percent_points']:.3f} / {mean_score['rmse_vol_percent_points']:.3f} vol%p
- Same-window 0-D mean MAE / RMSE: {zero_d_score['mae_vol_percent_points']:.3f} / {zero_d_score['rmse_vol_percent_points']:.3f} vol%p
- Measured / dynamic envelope peak: {summary['observed_envelope']['peak_vol_percent']:.3f} / {summary['dynamic_envelope']['peak_vol_percent']:.3f} vol%
- Measured / dynamic 4 vol% arrival: {summary['observed_envelope']['arrival_s']} / {summary['dynamic_envelope']['arrival_s']} s
- Measured / dynamic 4 vol% duration: {summary['observed_envelope']['duration_s']:.3f} / {summary['dynamic_envelope']['duration_s']:.3f} s

## Fixed inputs and boundaries

- Container: {CONTAINER_M} m; source: {SOURCE_POSITION_M} m.
- Measured source mass-flow and ventilation histories; resampled source-mass difference: {summary['source_resampling']['relative_mass_difference_percent']:.4f}%.
- x is open at inlet/doors; y sidewalls and z floor/ceiling are closed.
- Flatline sensor(s) excluded from error aggregation: {', '.join(invalid_sensors)}.

## Remaining interpretation limit

The release table, impinging downward jet, compressible near field, circular
inlet jet and open-door velocity/pressure field are unresolved. The effective
diffusivity (0.05 m2/s) and Gaussian source width (0.25 m) are declared
uncalibrated closure assumptions. Therefore these numbers are a model-gap
diagnostic, not a validation pass/fail certificate.
"""
    (args.output_dir / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=2, default=_jsonable))


if __name__ == "__main__":
    main()
