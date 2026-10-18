"""Diagnose HyTunnel source-handoff and scalar-mixing sensitivity without fitting.

Four predeclared resolved cases are compared on the same measurement subset:
the frozen 6% handoff/Sc=0.7 case, two effective-Schmidt perturbations, and a
more conservative 5% density-departure handoff.  Observations are used only
after all cases have been executed and never select a coefficient or default.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from typing import Any

import numpy as np

from degali.addons.high_pressure_source_handoff_3d import (
    low_mach_handoff_from_energy_nozzle,
)
from degali.addons.notional import energy_conserving_notional_nozzle
from degali.validation.hytunnel_ps import (
    DATASET_DOI,
    read_hytunnel_concentration_workbook,
)
from degali.validation.statistics import statistics


SCHEMA = "degali.hytunnel-source-mixing-sensitivity.v1"
WORKBOOK = "HTE2440PS002MIXED_300321.xlsx"
SHEET = "1 mm 5g H2 W = 3.5 co-flow"
DEFAULT_CASES = {
    "handoff0p06_sc0p7": (
        "outputs/hytunnel-resolved-published-box-0p707-final-2026-10-11-v1.json"
    ),
    "handoff0p06_sc0p5": (
        "outputs/hytunnel-resolved-sensitivity-sc0p5-0p707-2026-10-11-v1.json"
    ),
    "handoff0p06_sc1p0": (
        "outputs/hytunnel-resolved-sensitivity-sc1p0-0p707-2026-10-11-v1.json"
    ),
    "handoff0p05_sc0p7": (
        "outputs/hytunnel-resolved-sensitivity-handoff0p05-0p707-2026-10-11-v1.json"
    ),
}
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-source-mixing-sensitivity-2026-10-11-v1.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-source-mixing-sensitivity-2026-10-11-v1.png"
)


def _handoff(limit: float):
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
        source_id=f"hytunnel-sensitivity-rho-{limit:g}",
        ambient_temperature_k=286.15,
        ambient_pressure_pa=98_500.0,
        ambient_axial_velocity_m_s=3.5,
        maximum_mach_number=0.30,
        maximum_relative_density_departure=limit,
    )


def _linear_weights(axis: np.ndarray, value: float, domain_max: float):
    if value < 0.0 or value > domain_max:
        raise ValueError("sample lies outside domain")
    if value <= axis[0]:
        left, right = 0, 1
    elif value >= axis[-1]:
        left, right = axis.size - 2, axis.size - 1
    else:
        right = int(np.searchsorted(axis, value, side="right"))
        left = right - 1
    return left, right, float((value - axis[left]) / (axis[right] - axis[left]))


def _sample(
    field: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    x_m: float,
    y_m: float,
    domain_m: tuple[float, float],
) -> float:
    ix0, ix1, wx = _linear_weights(x_axis, x_m, domain_m[0])
    iy0, iy1, wy = _linear_weights(y_axis, y_m, domain_m[1])
    lower = (1.0 - wx) * field[iy0, ix0] + wx * field[iy0, ix1]
    upper = (1.0 - wx) * field[iy1, ix0] + wx * field[iy1, ix1]
    return float((1.0 - wy) * lower + wy * upper)


def _metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, Any]:
    residual = predicted - observed
    geometric = statistics(observed, predicted, floor=0.10)
    correlation = (
        float(np.corrcoef(observed, predicted)[0, 1])
        if np.std(observed) > 0.0 and np.std(predicted) > 0.0 else None
    )
    positive_observed = observed >= 4.0
    positive_predicted = predicted >= 4.0
    tp = int(np.sum(positive_observed & positive_predicted))
    tn = int(np.sum(~positive_observed & ~positive_predicted))
    fp = int(np.sum(~positive_observed & positive_predicted))
    fn = int(np.sum(positive_observed & ~positive_predicted))
    return {
        "pair_count": int(observed.size),
        "mae_vol_pct_point": float(np.mean(np.abs(residual))),
        "rmse_vol_pct_point": float(np.sqrt(np.mean(residual**2))),
        "mean_bias_vol_pct_point": float(np.mean(residual)),
        "pearson_r": correlation,
        "fac2_floor_0p10": geometric.fac2,
        "four_vol_confusion": {
            "true_positive": tp,
            "true_negative": tn,
            "false_positive": fp,
            "false_negative": fn,
            "sensitivity": tp / (tp + fn) if tp + fn else None,
            "specificity": tn / (tn + fp) if tn + fp else None,
        },
    }


def audit_source_mixing_sensitivity(
    root: str | Path,
    case_paths: dict[str, str | Path] | None = None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    started = time.perf_counter()
    directory = Path(root)
    declared = DEFAULT_CASES if case_paths is None else case_paths
    cases: dict[str, dict[str, Any]] = {}
    for label, relative in declared.items():
        path = Path(relative)
        if not path.is_absolute():
            path = directory / path
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("schema") != "degali.hytunnel-resolved-case.v1":
            raise ValueError(f"{label} is not a resolved-case artifact")
        if not record["endpoint_stationarity"]["passed"]:
            raise ValueError(f"{label} endpoint stationarity failed")
        cases[label] = {"path": path, "record": record}

    dataset = read_hytunnel_concentration_workbook(
        directory / "reference/datasets/hytunnel-ps-4668554/raw" / WORKBOOK
    )
    matched = next(item for item in dataset.configurations if item.sheet == SHEET)
    points = [
        point for point in dataset.points
        if point.configuration_id == matched.configuration_id
        and point.statistic == "average"
    ]
    maximum_source_x = max(
        float(item["record"]["case"]["source_position_m"][0])
        for item in cases.values()
    )
    common_domain_x = min(
        float(item["record"]["configuration"]["common_region_domain_m"][0])
        for item in cases.values()
    )
    included = [
        point for point in points
        if maximum_source_x - 1.0e-12 <= point.x_m <= common_domain_x + 1.0e-12
    ]
    observed = np.asarray(
        [point.concentration_vol_percent for point in included], dtype=float,
    )
    selected_x = np.asarray([point.x_m for point in included], dtype=float)
    selected_y = np.asarray([point.y_m for point in included], dtype=float)

    results: dict[str, Any] = {}
    arrays: dict[str, np.ndarray] = {
        "selected_x_m": selected_x,
        "selected_y_m": selected_y,
        "observed_vol_pct": observed,
    }
    for label, item in cases.items():
        artifact = item["record"]
        config = artifact["configuration"]
        case = artifact["case"]
        spacing = float(config["grid_spacing_m"])
        domain = tuple(float(value) for value in config["common_region_domain_m"])
        field = np.asarray(case["common_final_centre_slice"], dtype=float)
        x_axis = (np.arange(field.shape[1]) + 0.5) * spacing
        y_axis = (np.arange(field.shape[0]) + 0.5) * spacing
        source_y = float(case["source_position_m"][1])
        predicted = np.maximum(0.0, np.asarray([
            _sample(
                field, x_axis, y_axis, point.x_m,
                source_y + point.y_m, (domain[0], domain[1]),
            )
            for point in included
        ]))
        arrays[f"predicted_{label}"] = predicted
        density_limit = float(
            config["source_candidate_mean_density_departure_limit"]
        )
        handoff = _handoff(density_limit)
        results[label] = {
            "artifact": str(item["path"]),
            "declared_inputs": {
                "handoff_density_departure_limit": density_limit,
                "effective_schmidt_number": float(
                    config.get("effective_schmidt_number", 0.7)
                ),
                "applied_scalar_diffusivity_m2_s": float(
                    config.get("applied_scalar_diffusivity_m2_s", 0.04 / 0.7)
                ),
            },
            "derived_handoff": {
                "distance_m": handoff.handoff_distance_m,
                "equivalent_diameter_m": handoff.equivalent_diameter_m,
                "radial_sigma_m": handoff.moment_matched_radial_sigma_m,
                "hydrogen_mass_fraction": handoff.hydrogen_mass_fraction,
                "entrained_air_mass_flow_kg_s": handoff.entrained_air_mass_flow_kg_s,
                "mixture_velocity_m_s": handoff.mixture_velocity_m_s,
            },
            "endpoint_stationarity": artifact["endpoint_stationarity"],
            "metrics_on_common_points": _metrics(observed, predicted),
        }

    baseline = results["handoff0p06_sc0p7"]["metrics_on_common_points"]
    mae_values = [
        item["metrics_on_common_points"]["mae_vol_pct_point"]
        for item in results.values()
    ]
    bias_values = [
        item["metrics_on_common_points"]["mean_bias_vol_pct_point"]
        for item in results.values()
    ]
    record: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "diagnostic_source_mixing_sensitivity_not_calibration",
        "promotion_allowed": False,
        "calibration_or_best_variant_selection_allowed": False,
        "dataset_doi": DATASET_DOI,
        "workbook": WORKBOOK,
        "sheet": SHEET,
        "frozen_variants": list(declared),
        "common_point_accounting": {
            "published_average_point_count": len(points),
            "common_supported_point_count": len(included),
            "maximum_source_handoff_x_m": maximum_source_x,
            "common_domain_x_max_m": common_domain_x,
            "included_x_stations_m": sorted({point.x_m for point in included}),
            "included_transverse_stations_m": sorted({point.y_m for point in included}),
        },
        "cases": results,
        "sensitivity_summary": {
            "baseline_mae_vol_pct_point": baseline["mae_vol_pct_point"],
            "mae_range_across_declared_variants_vol_pct_point": [
                min(mae_values), max(mae_values),
            ],
            "bias_range_across_declared_variants_vol_pct_point": [
                min(bias_values), max(bias_values),
            ],
            "interpretation_policy": (
                "observed error may diagnose model sensitivity but may not select "
                "Sc_eff or a handoff limit; selection requires independent evidence"
            ),
        },
        "limitations": [
            "the 45 mm grid is much coarser than the official 6.4 mm CFD mesh",
            "Sc_eff perturbations change an unvalidated numerical/turbulent mixing closure",
            "the handoff threshold is a receiving-model admissibility boundary, not a fit parameter",
            "this unobstructed GH2 component case does not validate obstacle wakes or cryogenic buoyancy",
        ],
        "runtime_s": time.perf_counter() - started,
    }
    return record, arrays


def plot_sensitivity(
    record: dict[str, Any], arrays: dict[str, np.ndarray], output: str | Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = list(record["cases"])
    short = [
        "6% / Sc 0.7", "6% / Sc 0.5", "6% / Sc 1.0", "5% / Sc 0.7"
    ]
    metrics = [record["cases"][label]["metrics_on_common_points"] for label in labels]
    figure, axes = plt.subplots(1, 3, figsize=(15.2, 4.8), constrained_layout=True)
    positions = np.arange(len(labels))
    axes[0].bar(positions - 0.18, [m["mae_vol_pct_point"] for m in metrics], 0.36, label="MAE")
    axes[0].bar(positions + 0.18, [m["rmse_vol_pct_point"] for m in metrics], 0.36, label="RMSE")
    axes[0].set_xticks(positions, short, rotation=25, ha="right")
    axes[0].set_ylabel("vol%-point")
    axes[0].set_title("Common 45-point concentration error")
    axes[0].legend(frameon=False)

    axes[1].bar(positions, [m["mean_bias_vol_pct_point"] for m in metrics], color="tab:orange")
    axes[1].axhline(0.0, color="black", linewidth=1)
    axes[1].set_xticks(positions, short, rotation=25, ha="right")
    axes[1].set_ylabel("predicted − observed (vol%-point)")
    axes[1].set_title("Systematic bias")

    x = arrays["selected_x_m"]
    y = arrays["selected_y_m"]
    centre = np.isclose(y, 0.0)
    order = np.argsort(x[centre])
    axes[2].plot(
        x[centre][order], arrays["observed_vol_pct"][centre][order],
        "ko-", label="experiment average",
    )
    for label, display in zip(labels, short):
        axes[2].plot(
            x[centre][order], arrays[f"predicted_{label}"][centre][order],
            marker="o", label=display,
        )
    axes[2].axhline(4.0, color="tab:red", linestyle=":")
    axes[2].set_xlabel("x from nozzle (m)")
    axes[2].set_ylabel("centreline H$_2$ (vol%)")
    axes[2].set_title("Common-station centreline")
    axes[2].legend(frameon=False, fontsize=8)
    figure.suptitle(
        "HyTunnel source/handoff and scalar-mixing sensitivity\n"
        "Diagnostic only — observations do not select a variant",
        fontsize=12,
    )
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    args = parser.parse_args()
    record, arrays = audit_source_mixing_sensitivity(args.root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    plot_sensitivity(record, arrays, args.plot)
    print(json.dumps({
        "status": record["status"],
        "common_points": record["common_point_accounting"]["common_supported_point_count"],
        "mae_range": record["sensitivity_summary"]["mae_range_across_declared_variants_vol_pct_point"],
        "bias_range": record["sensitivity_summary"]["bias_range_across_declared_variants_vol_pct_point"],
    }, indent=2))


if __name__ == "__main__":
    main()
