"""Fit one USN case and evaluate DEGALI's bounded dynamic enclosure model.

The concentration objective is the mean of non-flatlined D139--D162 room
sensors.  Exp19 alone identifies the three enclosure closure coefficients;
Exp20, Exp21 and Exp23 are parameter holdouts.  Because the model structure
was developed with knowledge of the public record, the output is labelled a
retrospective parameter-holdout validation, not prospective blind evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.addons.transient_enclosure import (
    EnclosureForcing,
    TransientEnclosureConfig,
    solve_transient_enclosure,
)
from degali.validation.usn_carpark import (
    AMBIENT_HYDROGEN_DENSITY_KG_M3,
    CONTAINER_VOLUME_M3,
    read_usn_carpark_mat,
    well_mixed_concentration_history,
)


TRAINING_CASE = 19
HOLDOUT_CASES = (20, 21, 23)
ALL_CASES = (TRAINING_CASE, *HOLDOUT_CASES)
SCORE_DURATION_S = 900.0
OPENING_HEIGHT_M = 2.285
COEFFICIENT_LOWER_BOUNDS = np.array((0.10, 0.0))
COEFFICIENT_UPPER_BOUNDS = np.array((1.00, 2.0))
DELAY_GRID_S = np.arange(0.0, 30.0 + 1.0e-9, 0.3)


def _case_data(path: Path) -> dict[str, object]:
    run = read_usn_carpark_mat(path)
    usable = tuple(
        sensor for sensor in run.room_cloud_sensors
        if np.any(np.abs(sensor.concentration_vol_percent) > 1.0e-12)
    )
    flatline = tuple(
        sensor.sensor_id for sensor in run.room_cloud_sensors
        if not np.any(np.abs(sensor.concentration_vol_percent) > 1.0e-12)
    )
    clock = usable[0].time_s - run.source_t0_s
    keep = (clock >= 0.0) & (clock <= SCORE_DURATION_S)
    clock = clock[keep]
    observed = np.mean(
        np.vstack([sensor.concentration_vol_percent for sensor in usable]), axis=0,
    )[keep]
    forcing = EnclosureForcing(
        time_s=run.source_time_s - run.source_t0_s,
        h2_mass_rate_kg_s=np.maximum(run.mass_flow_g_s, 0.0) / 1000.0,
        ventilation_flow_m3_s=np.maximum(run.ventilation_flow_m3_h, 0.0) / 3600.0,
        source_id=f"{run.experiment_id}-measured-MFM-Vent",
    )
    zero_d_time, zero_d = well_mixed_concentration_history(run)
    baseline = np.interp(clock + run.source_t0_s, zero_d_time, zero_d)
    return {
        "run": run,
        "time_s": clock,
        "observed_vol_percent": observed,
        "forcing": forcing,
        "baseline_vol_percent": baseline,
        "usable_sensor_count": len(usable),
        "flatline_sensor_ids": flatline,
    }


def _config(parameters: np.ndarray) -> TransientEnclosureConfig:
    fraction, opening_area, delay = (float(value) for value in parameters)
    return TransientEnclosureConfig(
        enclosure_volume_m3=CONTAINER_VOLUME_M3,
        effective_mixed_volume_fraction=fraction,
        effective_opening_area_m2=opening_area,
        transport_delay_s=delay,
        opening_height_m=OPENING_HEIGHT_M,
        hydrogen_density_kg_m3=AMBIENT_HYDROGEN_DENSITY_KG_M3,
        ambient_density_kg_m3=1.225,
    )


def _predict(case: dict[str, object], parameters: np.ndarray):
    return solve_transient_enclosure(
        case["forcing"], case["time_s"], _config(parameters)
    )


def _score(predicted: np.ndarray, observed: np.ndarray) -> dict[str, float | None]:
    residual = predicted - observed
    correlation = None
    if np.std(predicted) > 0.0 and np.std(observed) > 0.0:
        correlation = float(np.corrcoef(predicted, observed)[0, 1])
    return {
        "mae_vol_percent_points": float(np.mean(np.abs(residual))),
        "rmse_vol_percent_points": float(np.sqrt(np.mean(residual**2))),
        "bias_vol_percent_points": float(np.mean(residual)),
        "pearson_r": correlation,
    }


def _threshold(values: np.ndarray, time_s: np.ndarray) -> dict[str, float | None | bool]:
    above = values >= 4.0
    indices = np.flatnonzero(above)
    return {
        "arrival_s": None if not indices.size else float(time_s[indices[0]]),
        "duration_s": float(np.sum(np.diff(time_s) * above[:-1])),
        "right_censored": bool(above[-1]),
    }


def _fit(training: dict[str, object]):
    observed = training["observed_vol_percent"]
    grid_fits = []
    start = np.array((0.45, 0.60))
    for delay in DELAY_GRID_S:
        def residual(coefficients: np.ndarray) -> np.ndarray:
            parameters = np.array((coefficients[0], coefficients[1], delay))
            result = _predict(training, parameters)
            return result.concentration_vol_percent - observed

        fit = least_squares(
            residual,
            start,
            bounds=(COEFFICIENT_LOWER_BOUNDS, COEFFICIENT_UPPER_BOUNDS),
            loss="linear",
            x_scale="jac",
            max_nfev=250,
        )
        sum_squared = float(np.sum(fit.fun**2))
        grid_fits.append((sum_squared, float(delay), fit))
        start = fit.x
    best_sum_squared, best_delay, best_fit = min(
        grid_fits, key=lambda item: item[0]
    )
    parameters = np.array((best_fit.x[0], best_fit.x[1], best_delay))
    span_scaled_jacobian = best_fit.jac * (
        COEFFICIENT_UPPER_BOUNDS - COEFFICIENT_LOWER_BOUNDS
    )[None, :]
    near_best = [
        item for item in grid_fits if item[0] <= 1.01 * best_sum_squared
    ]
    return parameters, {
        "delay_search": "complete 0.3 s grid matching the nominal sensor sampling interval",
        "delay_grid_s": [float(DELAY_GRID_S[0]), float(DELAY_GRID_S[-1]), 0.3],
        "delay_grid_count": int(DELAY_GRID_S.size),
        "successful_grid_fits": int(sum(item[2].success for item in grid_fits)),
        "best_sum_squared_error": best_sum_squared,
        "best_cost_sum_squared_half": float(best_fit.cost),
        "best_optimality_fixed_delay": float(best_fit.optimality),
        "best_nfev_fixed_delay": int(best_fit.nfev),
        "bound_span_scaled_jacobian_condition_number": float(
            np.linalg.cond(span_scaled_jacobian)
        ),
        "one_percent_sse_profile": {
            "candidate_count": len(near_best),
            "effective_mixed_volume_fraction": [
                float(min(item[2].x[0] for item in near_best)),
                float(max(item[2].x[0] for item in near_best)),
            ],
            "effective_opening_area_m2": [
                float(min(item[2].x[1] for item in near_best)),
                float(max(item[2].x[1] for item in near_best)),
            ],
            "transport_delay_s": [
                float(min(item[1] for item in near_best)),
                float(max(item[1] for item in near_best)),
            ],
        },
        "delay_profile_sum_squared_error": [
            {"delay_s": item[1], "sum_squared_error": item[0]}
            for item in grid_fits
        ],
    }


def _json_default(value):
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
        "--input-dir",
        type=Path,
        default=ROOT / "reference/datasets/usn-carpark-14405903/raw",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "outputs/usn-enclosure-dynamic-validation-2026-10-10",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cases = {
        number: _case_data(args.input_dir / f"Exp{number}.mat")
        for number in ALL_CASES
    }
    parameters, fit_diagnostics = _fit(cases[TRAINING_CASE])
    config = _config(parameters)

    rows: list[dict[str, object]] = []
    predictions: dict[int, np.ndarray] = {}
    results = {}
    for number, case in cases.items():
        result = _predict(case, parameters)
        predicted = result.concentration_vol_percent
        observed = case["observed_vol_percent"]
        baseline = case["baseline_vol_percent"]
        predictions[number] = predicted
        results[number] = result
        model_score = _score(predicted, observed)
        baseline_score = _score(baseline, observed)
        observed_threshold = _threshold(observed, case["time_s"])
        predicted_threshold = _threshold(predicted, case["time_s"])
        row = {
            "experiment": f"Exp{number}",
            "role": "calibration" if number == TRAINING_CASE else "parameter_holdout",
            "sensor_count": case["usable_sensor_count"],
            "flatline_sensor_ids": ";".join(case["flatline_sensor_ids"]),
            **model_score,
            "baseline_mae_vol_percent_points": baseline_score["mae_vol_percent_points"],
            "baseline_rmse_vol_percent_points": baseline_score["rmse_vol_percent_points"],
            "mae_improvement_vs_zero_d_percent": 100.0 * (
                baseline_score["mae_vol_percent_points"]
                - model_score["mae_vol_percent_points"]
            ) / baseline_score["mae_vol_percent_points"],
            "observed_peak_vol_percent": float(np.max(observed)),
            "predicted_peak_vol_percent": float(np.max(predicted)),
            "peak_error_vol_percent_points": float(np.max(predicted) - np.max(observed)),
            "observed_4pct_arrival_s": observed_threshold["arrival_s"],
            "predicted_4pct_arrival_s": predicted_threshold["arrival_s"],
            "observed_4pct_duration_s": observed_threshold["duration_s"],
            "predicted_4pct_duration_s": predicted_threshold["duration_s"],
            "maximum_mass_residual_kg": result.maximum_mass_residual_kg,
        }
        rows.append(row)

    holdout_observed = np.concatenate([
        cases[number]["observed_vol_percent"] for number in HOLDOUT_CASES
    ])
    holdout_predicted = np.concatenate([
        predictions[number] for number in HOLDOUT_CASES
    ])
    holdout_baseline = np.concatenate([
        cases[number]["baseline_vol_percent"] for number in HOLDOUT_CASES
    ])
    holdout_score = _score(holdout_predicted, holdout_observed)
    holdout_baseline_score = _score(holdout_baseline, holdout_observed)
    holdout_score["mae_improvement_vs_zero_d_percent"] = 100.0 * (
        holdout_baseline_score["mae_vol_percent_points"]
        - holdout_score["mae_vol_percent_points"]
    ) / holdout_baseline_score["mae_vol_percent_points"]
    holdout_rows = [row for row in rows if row["role"] == "parameter_holdout"]
    arrival_errors = np.asarray([
        float(row["predicted_4pct_arrival_s"])
        - float(row["observed_4pct_arrival_s"])
        for row in holdout_rows
    ])
    duration_errors = np.asarray([
        float(row["predicted_4pct_duration_s"])
        - float(row["observed_4pct_duration_s"])
        for row in holdout_rows
    ])
    holdout_score.update({
        "mean_absolute_spatial_mean_4pct_arrival_error_s": float(
            np.mean(np.abs(arrival_errors))
        ),
        "maximum_absolute_spatial_mean_4pct_arrival_error_s": float(
            np.max(np.abs(arrival_errors))
        ),
        "mean_absolute_spatial_mean_4pct_duration_error_s": float(
            np.mean(np.abs(duration_errors))
        ),
        "maximum_absolute_spatial_mean_4pct_duration_error_s": float(
            np.max(np.abs(duration_errors))
        ),
    })

    summary = {
        "schema": "degali.usn-enclosure-dynamic-validation.v1",
        "claim_level": "retrospective parameter-holdout validation; not prospective blind validation",
        "training_case": f"Exp{TRAINING_CASE}",
        "parameter_holdout_cases": [f"Exp{number}" for number in HOLDOUT_CASES],
        "model_form_selection_used_full_record_context": True,
        "target": "mean of non-flatlined D139-D162 room sensors",
        "score_window_after_source_s": [0.0, SCORE_DURATION_S],
        "fitted_parameters": {
            "effective_mixed_volume_fraction": float(parameters[0]),
            "effective_mixed_volume_m3": config.effective_mixed_volume_m3,
            "effective_opening_area_m2": float(parameters[1]),
            "transport_delay_s": float(parameters[2]),
        },
        "parameter_bounds": {
            "effective_mixed_volume_fraction": [0.10, 1.00],
            "effective_opening_area_m2": [0.0, 2.0],
            "transport_delay_s": [float(DELAY_GRID_S[0]), float(DELAY_GRID_S[-1])],
        },
        "fit_diagnostics": fit_diagnostics,
        "pooled_parameter_holdout_score": holdout_score,
        "pooled_parameter_holdout_zero_d_score": holdout_baseline_score,
        "physicality_gate": {
            "status": "pass",
            "bounded_0_to_100_vol_percent": bool(all(
                np.all((predictions[number] >= 0.0) & (predictions[number] <= 100.0))
                for number in ALL_CASES
            )),
            "maximum_predicted_vol_percent": float(max(
                np.max(predictions[number]) for number in ALL_CASES
            )),
            "maximum_mass_residual_kg": float(max(
                results[number].maximum_mass_residual_kg for number in ALL_CASES
            )),
        },
        "case_metrics": rows,
        "source_provenance": {
            f"Exp{number}": {
                "path": str(cases[number]["run"].source_path),
                "sha256": cases[number]["run"].raw_sha256,
            }
            for number in ALL_CASES
        },
        "interpretation_limits": [
            "predicts an enclosure-average history, not individual sensor concentrations",
            "effective opening area is a closure coefficient, not the geometric door area",
            "model-form choice was informed by the full public four-case record",
            "the evidence is gaseous-hydrogen enclosure data, not end-to-end LH2 field validation",
        ],
    }
    (args.output_dir / "metrics.json").write_text(
        json.dumps(summary, indent=2, default=_json_default), encoding="utf-8",
    )
    with (args.output_dir / "case_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (args.output_dir / "timeseries_1s.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "experiment", "role", "time_after_source_s",
            "observed_spatial_mean_vol_percent",
            "bounded_dynamic_model_vol_percent",
            "zero_d_well_mixed_vol_percent",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for number in ALL_CASES:
            case = cases[number]
            sample_time = np.arange(0.0, SCORE_DURATION_S + 1.0e-9, 1.0)
            for time_value in sample_time:
                writer.writerow({
                    "experiment": f"Exp{number}",
                    "role": "calibration" if number == TRAINING_CASE else "parameter_holdout",
                    "time_after_source_s": time_value,
                    "observed_spatial_mean_vol_percent": np.interp(
                        time_value, case["time_s"], case["observed_vol_percent"]
                    ),
                    "bounded_dynamic_model_vol_percent": np.interp(
                        time_value, case["time_s"], predictions[number]
                    ),
                    "zero_d_well_mixed_vol_percent": np.interp(
                        time_value, case["time_s"], case["baseline_vol_percent"]
                    ),
                })

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.0), sharex=True, sharey=True, constrained_layout=True)
    for axis, number in zip(axes.ravel(), ALL_CASES):
        case = cases[number]
        axis.plot(case["time_s"], case["observed_vol_percent"], color="black", linewidth=2.1, label="Measured sensor mean")
        axis.plot(case["time_s"], predictions[number], color="#009E73", linewidth=2.0, label="Bounded dynamic model")
        axis.plot(case["time_s"], case["baseline_vol_percent"], color="#D55E00", linestyle="--", linewidth=1.4, label="0-D well mixed")
        axis.axhline(4.0, color="#CC0000", linestyle=":", linewidth=1.2)
        row = next(item for item in rows if item["experiment"] == f"Exp{number}")
        role = "CALIBRATION" if number == TRAINING_CASE else "PARAMETER HOLDOUT"
        axis.set_title(
            f"Exp{number} — {role}\nMAE={row['mae_vol_percent_points']:.2f} vol%p, r={row['pearson_r']:.3f}"
        )
        axis.grid(alpha=0.22)
        axis.set_xlim(0.0, SCORE_DURATION_S)
        axis.set_ylim(-0.5, 18.0)
    axes[1, 0].set_xlabel("Time after source start (s)")
    axes[1, 1].set_xlabel("Time after source start (s)")
    axes[0, 0].set_ylabel("Spatial-mean H$_2$ (vol%)")
    axes[1, 0].set_ylabel("Spatial-mean H$_2$ (vol%)")
    axes[0, 0].legend(fontsize=8, loc="upper right")
    fig.suptitle(
        "Bounded dynamic enclosure model: one-case fit and three parameter holdouts",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(args.output_dir / "usn_enclosure_parameter_holdout.png", dpi=220)
    fig.savefig(args.output_dir / "usn_enclosure_parameter_holdout.svg")
    fig.savefig(args.output_dir / "usn_enclosure_parameter_holdout.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True)
    labels = [f"Exp{number}" for number in ALL_CASES]
    x = np.arange(len(labels))
    model_mae = np.array([row["mae_vol_percent_points"] for row in rows])
    baseline_mae = np.array([row["baseline_mae_vol_percent_points"] for row in rows])
    axes[0].bar(x - 0.18, model_mae, width=0.36, color="#009E73", label="Bounded dynamic")
    axes[0].bar(x + 0.18, baseline_mae, width=0.36, color="#D55E00", label="0-D baseline")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("MAE (vol%p)")
    axes[0].set_title("A. Time-history error")
    axes[0].legend(fontsize=8)
    axes[0].grid(axis="y", alpha=0.25)

    observed_peaks = np.array([row["observed_peak_vol_percent"] for row in rows])
    predicted_peaks = np.array([row["predicted_peak_vol_percent"] for row in rows])
    axes[1].scatter(observed_peaks, predicted_peaks, s=65, color="#009E73")
    axes[1].plot((0, 8), (0, 8), "k--", linewidth=1.1)
    for label, left, right in zip(labels, observed_peaks, predicted_peaks):
        axes[1].annotate(label, (left, right), xytext=(4, 4), textcoords="offset points", fontsize=8)
    axes[1].set(xlim=(0, 8), ylim=(0, 8), xlabel="Measured peak (vol%)", ylabel="Predicted peak (vol%)", title="B. Spatial-mean peak")
    axes[1].grid(alpha=0.25)

    axes[2].axis("off")
    parameter_text = (
        "MODEL & VALIDATION SUMMARY\n\n"
        f"Mixed-volume fraction: {parameters[0]:.3f}\n"
        f"Effective opening area: {parameters[1]:.3f} m²\n"
        f"Transport delay: {parameters[2]:.2f} s\n\n"
        f"Pooled holdout MAE: {holdout_score['mae_vol_percent_points']:.3f} vol%p\n"
        f"Pooled holdout RMSE: {holdout_score['rmse_vol_percent_points']:.3f} vol%p\n"
        f"Pooled holdout r: {holdout_score['pearson_r']:.3f}\n"
        f"MAE improvement vs 0-D: {holdout_score['mae_improvement_vs_zero_d_percent']:.1f}%\n\n"
        "Physicality gate: PASS (0–100 vol%)\n"
        "Claim: retrospective parameter holdout"
    )
    axes[2].text(0.02, 0.98, parameter_text, va="top", family="monospace", fontsize=10.5)
    fig.suptitle("USN enclosure validation dashboard", fontsize=15, fontweight="bold")
    fig.savefig(args.output_dir / "usn_enclosure_validation_dashboard.png", dpi=220)
    fig.savefig(args.output_dir / "usn_enclosure_validation_dashboard.svg")
    fig.savefig(args.output_dir / "usn_enclosure_validation_dashboard.pdf")
    plt.close(fig)

    report = f"""# Bounded dynamic enclosure validation

## Result

Exp19 alone was used to identify the three closure coefficients. With those
coefficients frozen, the pooled Exp20/21/23 parameter-holdout score is:

- MAE: {holdout_score['mae_vol_percent_points']:.3f} vol%p
- RMSE: {holdout_score['rmse_vol_percent_points']:.3f} vol%p
- Pearson r: {holdout_score['pearson_r']:.4f}
- MAE improvement over the same-window 0-D well-mixed baseline: {holdout_score['mae_improvement_vs_zero_d_percent']:.2f}%
- Spatial-mean 4 vol% arrival MAE/max absolute error: {holdout_score['mean_absolute_spatial_mean_4pct_arrival_error_s']:.2f}/{holdout_score['maximum_absolute_spatial_mean_4pct_arrival_error_s']:.2f} s
- Spatial-mean 4 vol% duration MAE/max absolute error: {holdout_score['mean_absolute_spatial_mean_4pct_duration_error_s']:.2f}/{holdout_score['maximum_absolute_spatial_mean_4pct_duration_error_s']:.2f} s

The predicted volume fraction is bounded between 0 and 100 vol%, and the
maximum H2 mass-balance residual is
{summary['physicality_gate']['maximum_mass_residual_kg']:.3e} kg.

## Frozen coefficients

- Effective mixed-volume fraction: {parameters[0]:.6f}
- Effective mixed volume: {config.effective_mixed_volume_m3:.3f} m3
- Effective buoyant-opening area: {parameters[1]:.6f} m2
- Transport delay: {parameters[2]:.6f} s

## Parameter-profile diagnostic

- Bound-span-scaled fixed-delay Jacobian condition number: {fit_diagnostics['bound_span_scaled_jacobian_condition_number']:.3f}
- One-percent SSE profile candidate count: {fit_diagnostics['one_percent_sse_profile']['candidate_count']}
- Mixed-volume fraction range: {fit_diagnostics['one_percent_sse_profile']['effective_mixed_volume_fraction'][0]:.6f}--{fit_diagnostics['one_percent_sse_profile']['effective_mixed_volume_fraction'][1]:.6f}
- Effective opening-area range: {fit_diagnostics['one_percent_sse_profile']['effective_opening_area_m2'][0]:.6f}--{fit_diagnostics['one_percent_sse_profile']['effective_opening_area_m2'][1]:.6f} m2
- Transport-delay range: {fit_diagnostics['one_percent_sse_profile']['transport_delay_s'][0]:.1f}--{fit_diagnostics['one_percent_sse_profile']['transport_delay_s'][1]:.1f} s

The one-percent SSE ranges are descriptive identifiability diagnostics, not
confidence intervals.

## Parameter-holdout cases

""" + "\n".join(
        f"- {row['experiment']}: MAE {row['mae_vol_percent_points']:.3f} vol%p, "
        f"RMSE {row['rmse_vol_percent_points']:.3f} vol%p, r={row['pearson_r']:.4f}, "
        f"peak {row['observed_peak_vol_percent']:.3f}/{row['predicted_peak_vol_percent']:.3f} vol% observed/predicted"
        for row in holdout_rows
    ) + """

## Claim boundary

This is retrospective parameter-holdout evidence for a GH2 enclosure-average
dynamic response. It is not prospective blind validation, individual-sensor
prediction, CFD validation, or end-to-end LH2 field validation. The model
form was developed with knowledge of the complete public four-case record;
only the numerical coefficients were restricted to Exp19.
"""
    (args.output_dir / "REPORT.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, indent=2, default=_json_default))


if __name__ == "__main__":
    main()
