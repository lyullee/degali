"""No-fit HyTunnel comparison for the conservative low-Mach source handoff.

This is a deliberately narrow falsification test.  It combines the declared
DEGALI high-pressure source, the conservative low-Mach handoff and the existing
axisymmetric free-jet closure, then compares only published *average* H2
concentrations at coordinates downstream of the handoff.  Tunnel confinement
is not represented, so failure rejects this shortcut rather than the resolved
3-D tunnel solver.  No observation is used to choose a source or coefficient.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np

from degali.addons.axisymmetric_jet import AxisymmetricJetSource, ConservedGaussianJet
from degali.addons.high_pressure_source_handoff_3d import (
    low_mach_handoff_from_energy_nozzle,
)
from degali.addons.notional import energy_conserving_notional_nozzle
from degali.validation.hytunnel_ps import (
    DATASET_DOI,
    read_hytunnel_concentration_workbook,
)
from degali.validation.statistics import statistics


SCHEMA = "degali.hytunnel-free-jet-handoff-comparison.v1"
WORKBOOK = "HTE2440PS002MIXED_300321.xlsx"
SHEET = "1 mm 5g H2 W = 3.5 co-flow"


def _mole_percent_from_mass_fraction(
    mass_fraction: float | np.ndarray,
    hydrogen_molar_mass: float,
    air_molar_mass: float,
) -> float | np.ndarray:
    fraction = np.asarray(mass_fraction, dtype=float)
    hydrogen = fraction / hydrogen_molar_mass
    air = (1.0 - fraction) / air_molar_mass
    result = 100.0 * hydrogen / (hydrogen + air)
    return float(result) if result.ndim == 0 else result


def _first_downward_crossing(
    x: np.ndarray, concentration: np.ndarray, threshold: float,
) -> float | None:
    for index in range(1, len(x)):
        left, right = concentration[index - 1], concentration[index]
        if left >= threshold and right < threshold:
            fraction = (left - threshold) / max(left - right, 1.0e-30)
            return float(x[index - 1] + fraction * (x[index] - x[index - 1]))
    return None


def _linear_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, object]:
    residual = predicted - observed
    correlation = (
        float(np.corrcoef(observed, predicted)[0, 1])
        if observed.size > 1 and np.std(observed) > 0.0 and np.std(predicted) > 0.0
        else None
    )
    relative = np.abs(residual) / np.maximum(np.abs(observed), 1.0e-12)
    return {
        "pair_count": int(observed.size),
        "mae_vol_pct_point": float(np.mean(np.abs(residual))),
        "rmse_vol_pct_point": float(np.sqrt(np.mean(residual**2))),
        "mean_bias_predicted_minus_observed_vol_pct_point": float(np.mean(residual)),
        "median_absolute_error_vol_pct_point": float(np.median(np.abs(residual))),
        "median_absolute_relative_error": float(np.median(relative)),
        "pearson_r": correlation,
        "observed_range_vol_pct": [float(np.min(observed)), float(np.max(observed))],
        "predicted_range_vol_pct": [float(np.min(predicted)), float(np.max(predicted))],
    }


def _threshold_classification(
    observed: np.ndarray, predicted: np.ndarray, threshold: float,
) -> dict[str, object]:
    observed_positive = observed >= threshold
    predicted_positive = predicted >= threshold
    true_positive = int(np.sum(observed_positive & predicted_positive))
    true_negative = int(np.sum(~observed_positive & ~predicted_positive))
    false_positive = int(np.sum(~observed_positive & predicted_positive))
    false_negative = int(np.sum(observed_positive & ~predicted_positive))
    return {
        "threshold_vol_pct": threshold,
        "true_positive": true_positive,
        "true_negative": true_negative,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "sensitivity": (
            true_positive / (true_positive + false_negative)
            if true_positive + false_negative else None
        ),
        "specificity": (
            true_negative / (true_negative + false_positive)
            if true_negative + false_positive else None
        ),
        "precision": (
            true_positive / (true_positive + false_positive)
            if true_positive + false_positive else None
        ),
    }


def compare_hytunnel_free_jet_handoff(
    root: str | Path,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    from CoolProp.CoolProp import PropsSI

    started = time.perf_counter()
    directory = Path(root)
    dataset = read_hytunnel_concentration_workbook(directory / WORKBOOK)
    configurations = [item for item in dataset.configurations if item.sheet == SHEET]
    if len(configurations) != 1:
        raise ValueError(f"expected one matched HyTunnel sheet, found {len(configurations)}")
    configuration = configurations[0]
    average_points = [
        point for point in dataset.points
        if point.configuration_id == configuration.configuration_id
        and point.statistic == "average"
    ]

    ambient_temperature = 286.15
    ambient_pressure = 98_500.0
    nozzle = energy_conserving_notional_nozzle(
        fluid="Hydrogen",
        storage_temperature=ambient_temperature,
        storage_pressure=11.515e6,
        mass_flow=0.005,
        orifice_diameter=0.001,
        ambient_pressure=ambient_pressure,
        allow_supercritical_gas=True,
    )
    handoff = low_mach_handoff_from_energy_nozzle(
        nozzle,
        source_id="hytunnel-1mm-5g-3p5ms-coflow",
        ambient_temperature_k=ambient_temperature,
        ambient_pressure_pa=ambient_pressure,
        ambient_axial_velocity_m_s=3.5,
        maximum_mach_number=0.30,
        maximum_relative_density_departure=0.10,
    )

    source = AxisymmetricJetSource(
        diameter=handoff.equivalent_diameter_m,
        velocity=handoff.mixture_velocity_m_s,
        density=handoff.mixture_density_kg_m3,
        temperature=handoff.mixture_temperature_k,
        mass_fraction=handoff.hydrogen_mass_fraction,
        theta=0.0,
        x=0.0,
        y=0.0,
    )
    model = ConservedGaussianJet(
        source,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure,
        ambient_density=handoff.ambient_density_kg_m3,
        fuel_molecular_weight=float(PropsSI("M", "Hydrogen")),
        ambient_molecular_weight=float(PropsSI("M", "Air")),
        fuel_heat_capacity=float(PropsSI(
            "C", "T", ambient_temperature, "P", ambient_pressure, "Hydrogen",
        )),
        ambient_heat_capacity=float(PropsSI(
            "C", "T", ambient_temperature, "P", ambient_pressure, "Air",
        )),
        ambient_coflow_velocity=3.5,
        radial_points=161,
        conservative_establishment="source_flux",
    )
    result = model.solve(
        maximum_distance=4.1,
        minimum_mass_fraction=1.0e-6,
        maximum_step=0.005,
        relative_tolerance=2.0e-7,
    )
    # The source-flux Gaussian conversion normally carries a published 6.2D
    # development coordinate.  Here the already mixed handoff is the declared
    # source plane, so translate the autonomous solution without altering any
    # state or flux.  This is an explicit virtual-origin convention, not a fit.
    coordinate_shift = handoff.handoff_distance_m - float(result.x[0])
    model_x = result.x + coordinate_shift

    selected = [
        point for point in average_points
        if model_x[0] - 1.0e-12 <= point.x_m <= model_x[-1] + 1.0e-12
    ]
    observed: list[float] = []
    predicted: list[float] = []
    selected_x: list[float] = []
    selected_y: list[float] = []
    state_fields = (
        result.velocity, result.width, result.density, result.mass_fraction,
        result.theta, model_x, result.y,
    )
    for point in selected:
        state = np.asarray([
            np.interp(point.x_m, model_x, values) for values in state_fields
        ])
        radial_mass_fraction = model.radial_state(
            state, abs(point.y_m - state[6]),
        )[0]
        predicted.append(float(_mole_percent_from_mass_fraction(
            radial_mass_fraction,
            model.fuel_molecular_weight,
            model.ambient_molecular_weight,
        )))
        observed.append(point.concentration_vol_percent)
        selected_x.append(point.x_m)
        selected_y.append(point.y_m)

    observed_array = np.asarray(observed)
    predicted_array = np.asarray(predicted)
    geometric = statistics(observed_array, predicted_array, floor=0.10)
    high = observed_array >= 4.0
    centre = np.isclose(np.asarray(selected_y), 0.0, atol=1.0e-12)
    centreline = np.asarray(_mole_percent_from_mass_fraction(
        result.mass_fraction,
        model.fuel_molecular_weight,
        model.ambient_molecular_weight,
    ))
    model_four = _first_downward_crossing(model_x, centreline, 4.0)
    model_ten = _first_downward_crossing(model_x, centreline, 10.0)
    experiment_four = 1.65
    experiment_ten = 0.60
    four_error = (
        None if model_four is None
        else (model_four - experiment_four) / experiment_four
    )
    ten_error = (
        None if model_ten is None
        else (model_ten - experiment_ten) / experiment_ten
    )
    diagnostic_gates = {
        "fac2_at_floor_0p10_at_least_0p70": geometric.fac2 >= 0.70,
        "four_vol_pct_length_absolute_relative_error_at_most_0p20": (
            four_error is not None and abs(four_error) <= 0.20
        ),
        "four_vol_pct_point_sensitivity_at_least_0p90": (
            (_threshold_classification(observed_array, predicted_array, 4.0)[
                "sensitivity"
            ] or 0.0) >= 0.90
        ),
    }
    record: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "diagnostic_component_route_rejected"
            if not all(diagnostic_gates.values())
            else "diagnostic_component_route_passed_full_3d_validation_still_required"
        ),
        "scope": "steady_axisymmetric_free_jet_downstream_of_low_mach_handoff",
        "dataset_doi": DATASET_DOI,
        "workbook": WORKBOOK,
        "sheet": SHEET,
        "observation_statistic": "published steady-window average",
        "no_fit_statement": (
            "source, entrainment, Gaussian spreading and coordinate mapping were "
            "fixed before scoring; no observation selected a coefficient"
        ),
        "exclusions": [
            "points upstream of the low-Mach handoff",
            "published maximum and minimum matrices",
            "tunnel wall and inlet-profile effects",
            "unavailable raw transient sensor series",
        ],
        "source": {
            "degali_notional_nozzle_diameter_m": nozzle.diameter,
            "official_cfd_notional_nozzle_diameter_m": 0.00725,
            "source_model_match_to_official_cfd": False,
            "handoff_distance_m": handoff.handoff_distance_m,
            "handoff_equivalent_diameter_m": handoff.equivalent_diameter_m,
            "handoff_mach_number": handoff.mach_number,
            "handoff_relative_density_departure": handoff.relative_density_departure,
            "handoff_maximum_conservation_residual": (
                handoff.maximum_relative_conservation_residual
            ),
            "gaussian_establishment": "source_flux",
            "gaussian_coordinate_shift_m": coordinate_shift,
            "first_model_coordinate_m": float(model_x[0]),
            "last_model_coordinate_m": float(model_x[-1]),
            "maximum_establishment_flux_residual": max(
                abs(value) for value in model.establishment_residuals.values()
            ),
        },
        "paired_average_concentrations": {
            "all_downstream_points": _linear_metrics(
                observed_array, predicted_array,
            ),
            "observed_at_or_above_4_vol_pct": (
                None if not np.any(high)
                else _linear_metrics(observed_array[high], predicted_array[high])
            ),
            "centreline_y_zero": (
                None if not np.any(centre)
                else _linear_metrics(observed_array[centre], predicted_array[centre])
            ),
            "logarithmic_statistics_floor_vol_pct": 0.10,
            "hanna_statistics": {
                "n": geometric.n,
                "mg_observed_over_predicted": geometric.mg,
                "vg": geometric.vg,
                "fac2": geometric.fac2,
                "fractional_bias": geometric.fb,
                "nmse": geometric.nmse,
                "acceptable": geometric.acceptable,
            },
        },
        "threshold_classification": {
            "four_vol_pct": _threshold_classification(
                observed_array, predicted_array, 4.0,
            ),
            "ten_vol_pct": _threshold_classification(
                observed_array, predicted_array, 10.0,
            ),
        },
        "centreline_threshold_lengths": {
            "experiment_4_vol_pct_m": experiment_four,
            "official_cfd_4_vol_pct_m": 2.25,
            "degali_free_jet_4_vol_pct_m": model_four,
            "degali_relative_error_vs_experiment": four_error,
            "experiment_10_vol_pct_m": experiment_ten,
            "official_cfd_10_vol_pct_m": 0.50,
            "degali_free_jet_10_vol_pct_m": model_ten,
            "degali_10_vol_pct_relative_error_vs_experiment": ten_error,
        },
        "diagnostic_gates": diagnostic_gates,
        "interpretation": (
            "A failure rejects the axisymmetric free-jet shortcut for the matched "
            "tunnel concentration field. It does not score or reject the unresolved "
            "3-D tunnel model."
        ),
        "runtime_s": time.perf_counter() - started,
    }
    arrays = {
        "selected_x_m": np.asarray(selected_x),
        "selected_y_m": np.asarray(selected_y),
        "observed_vol_pct": observed_array,
        "predicted_vol_pct": predicted_array,
        "model_x_m": model_x,
        "model_centerline_vol_pct": centreline,
        "handoff_distance_m": np.asarray([handoff.handoff_distance_m]),
    }
    return record, arrays


def plot_comparison(arrays: dict[str, np.ndarray], output: str | Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_values = np.unique(arrays["selected_x_m"])
    y_values = np.unique(arrays["selected_y_m"])
    observed = np.full((len(y_values), len(x_values)), np.nan)
    predicted = np.full_like(observed, np.nan)
    for x, y, obs, pred in zip(
        arrays["selected_x_m"], arrays["selected_y_m"],
        arrays["observed_vol_pct"], arrays["predicted_vol_pct"],
    ):
        ix = int(np.searchsorted(x_values, x))
        iy = int(np.searchsorted(y_values, y))
        observed[iy, ix] = obs
        predicted[iy, ix] = pred

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 7.7), constrained_layout=True)
    levels = [0.0, 1.0, 2.0, 4.0, 10.0, 20.0, 40.0, 60.0]
    for axis, values, title in (
        (axes[0, 0], observed, "HyTunnel average H$_2$"),
        (axes[0, 1], predicted, "DEGALI handoff + free jet"),
    ):
        image = axis.contourf(x_values, y_values, values, levels=levels, extend="max")
        axis.contour(x_values, y_values, values, levels=[4.0, 10.0], colors=("white", "black"))
        axis.axvline(arrays["handoff_distance_m"][0], color="magenta", linestyle="--")
        axis.set(title=title, xlabel="x from nozzle (m)", ylabel="radial y (m)")
        fig.colorbar(image, ax=axis, label="H$_2$ (vol%)")

    difference = predicted - observed
    limit = max(1.0, float(np.nanmax(np.abs(difference))))
    image = axes[1, 0].pcolormesh(
        x_values, y_values, difference, shading="nearest",
        cmap="coolwarm", vmin=-limit, vmax=limit,
    )
    axes[1, 0].axvline(
        arrays["handoff_distance_m"][0], color="magenta", linestyle="--",
    )
    axes[1, 0].set(
        title="Model minus experiment", xlabel="x from nozzle (m)",
        ylabel="radial y (m)",
    )
    fig.colorbar(image, ax=axes[1, 0], label="error (vol%-point)")

    centre = np.isclose(arrays["selected_y_m"], 0.0, atol=1.0e-12)
    order = np.argsort(arrays["selected_x_m"][centre])
    axes[1, 1].plot(
        arrays["selected_x_m"][centre][order],
        arrays["observed_vol_pct"][centre][order], "o", label="experiment average",
    )
    axes[1, 1].plot(
        arrays["model_x_m"], arrays["model_centerline_vol_pct"],
        "-", label="DEGALI free jet",
    )
    axes[1, 1].axhline(4.0, color="tab:orange", linestyle="--", label="4 vol%")
    axes[1, 1].axhline(10.0, color="tab:red", linestyle=":", label="10 vol%")
    axes[1, 1].axvline(
        arrays["handoff_distance_m"][0], color="magenta", linestyle="--",
        label="low-Mach handoff",
    )
    axes[1, 1].set(
        title="Centreline decay", xlabel="x from nozzle (m)", ylabel="H$_2$ (vol%)",
        xlim=(0.0, 3.6), ylim=(0.0, None),
    )
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(fontsize=8)
    fig.suptitle(
        "HyTunnel 1 mm, 5 g/s, 3.5 m/s co-flow — no-fit component test",
        fontsize=13,
    )
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_record(record: dict[str, Any], output: str | Path) -> None:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    args = parser.parse_args(argv)
    record, arrays = compare_hytunnel_free_jet_handoff(args.root)
    write_record(record, args.output)
    plot_comparison(arrays, args.plot)
    metrics = record["paired_average_concentrations"]["all_downstream_points"]
    print(
        f"wrote {args.output} and {args.plot}; status={record['status']} "
        f"n={metrics['pair_count']} MAE={metrics['mae_vol_pct_point']:.4f} vol%p"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
