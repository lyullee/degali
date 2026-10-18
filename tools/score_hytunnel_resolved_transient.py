"""Score the frozen resolved HyTunnel transient against public measurements.

The resolved calculation starts at its low-Mach handoff plane, not at the
physical nozzle.  Consequently this tool scores only measurement stations
which are both downstream of that handoff and inside the numerical domain.
It never fills unsupported near-field points by extrapolation.  The current
frozen audit stores the maximum-over-time centre slice, whereas the public
workbook reports a steady-window average; that mismatch is explicit and keeps
the result conditional.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from typing import Any

import numpy as np

from degali.validation.hytunnel_ps import (
    DATASET_DOI,
    read_hytunnel_concentration_workbook,
)
from degali.validation.statistics import statistics


SCHEMA = "degali.hytunnel-resolved-transient-score.v3"
WORKBOOK = "HTE2440PS002MIXED_300321.xlsx"
SHEET = "1 mm 5g H2 W = 3.5 co-flow"
DEFAULT_AUDIT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-plume-domain-sensitivity-2026-10-11-official-0p707-v6.json"
)
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-resolved-transient-score-2026-10-11-v1.json"
)
DEFAULT_PLOT = (
    Path(__file__).resolve().parents[1]
    / "outputs/hytunnel-resolved-transient-score-2026-10-11-v1.png"
)


def _linear_axis_weights(
    coordinates: np.ndarray, value: float, *, domain: tuple[float, float],
) -> tuple[int, int, float]:
    """Return bracketing indices and a linear weight inside a face domain.

    Cell-centred data require a short one-sided linear extrapolation between
    the last cell centre and the physical boundary.  Values outside the
    declared face domain are rejected.
    """

    axis = np.asarray(coordinates, dtype=float)
    if axis.ndim != 1 or axis.size < 2 or np.any(np.diff(axis) <= 0.0):
        raise ValueError("coordinates must be a strictly increasing 1-D axis")
    if not np.isfinite(value) or value < domain[0] or value > domain[1]:
        raise ValueError("sample lies outside the declared physical domain")
    if value <= axis[0]:
        left, right = 0, 1
    elif value >= axis[-1]:
        left, right = axis.size - 2, axis.size - 1
    else:
        right = int(np.searchsorted(axis, value, side="right"))
        left = right - 1
    weight = (value - axis[left]) / (axis[right] - axis[left])
    return left, right, float(weight)


def bilinear_cell_centre_sample(
    field_yx: np.ndarray,
    x_centres_m: np.ndarray,
    y_centres_m: np.ndarray,
    x_m: float,
    y_m: float,
    *,
    domain_m: tuple[float, float],
) -> float:
    """Sample a cell-centred ``(y, x)`` field without clipping coordinates."""

    field = np.asarray(field_yx, dtype=float)
    x_axis = np.asarray(x_centres_m, dtype=float)
    y_axis = np.asarray(y_centres_m, dtype=float)
    if field.shape != (y_axis.size, x_axis.size):
        raise ValueError("field shape must match the y and x axes")
    if np.any(~np.isfinite(field)):
        raise ValueError("field must contain only finite values")
    ix0, ix1, wx = _linear_axis_weights(
        x_axis, float(x_m), domain=(0.0, float(domain_m[0])),
    )
    iy0, iy1, wy = _linear_axis_weights(
        y_axis, float(y_m), domain=(0.0, float(domain_m[1])),
    )
    lower = (1.0 - wx) * field[iy0, ix0] + wx * field[iy0, ix1]
    upper = (1.0 - wx) * field[iy1, ix0] + wx * field[iy1, ix1]
    return float((1.0 - wy) * lower + wy * upper)


def _linear_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, Any]:
    residual = predicted - observed
    correlation = (
        float(np.corrcoef(observed, predicted)[0, 1])
        if observed.size > 1 and np.std(observed) > 0.0 and np.std(predicted) > 0.0
        else None
    )
    relative = np.abs(residual) / np.maximum(np.abs(observed), 0.10)
    return {
        "pair_count": int(observed.size),
        "mae_vol_pct_point": float(np.mean(np.abs(residual))),
        "rmse_vol_pct_point": float(np.sqrt(np.mean(residual**2))),
        "mean_bias_predicted_minus_observed_vol_pct_point": float(np.mean(residual)),
        "median_absolute_error_vol_pct_point": float(np.median(np.abs(residual))),
        "median_absolute_relative_error_at_0p10_floor": float(np.median(relative)),
        "pearson_r": correlation,
        "observed_range_vol_pct": [float(np.min(observed)), float(np.max(observed))],
        "predicted_range_vol_pct": [float(np.min(predicted)), float(np.max(predicted))],
    }


def _threshold_classification(
    observed: np.ndarray, predicted: np.ndarray, threshold: float,
) -> dict[str, Any]:
    observed_positive = observed >= threshold
    predicted_positive = predicted >= threshold
    tp = int(np.sum(observed_positive & predicted_positive))
    tn = int(np.sum(~observed_positive & ~predicted_positive))
    fp = int(np.sum(~observed_positive & predicted_positive))
    fn = int(np.sum(observed_positive & ~predicted_positive))
    sensitivity = tp / (tp + fn) if tp + fn else None
    specificity = tn / (tn + fp) if tn + fp else None
    precision = tp / (tp + fp) if tp + fp else None
    f1 = (
        2.0 * precision * sensitivity / (precision + sensitivity)
        if precision is not None and sensitivity is not None
        and precision + sensitivity > 0.0 else None
    )
    return {
        "threshold_vol_pct": threshold,
        "true_positive": tp,
        "true_negative": tn,
        "false_positive": fp,
        "false_negative": fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "f1": f1,
        "accuracy": (tp + tn) / observed.size,
    }


def _case_slice(case: dict[str, Any]) -> tuple[np.ndarray, str]:
    if "common_final_centre_slice" in case:
        return (
            np.asarray(case["common_final_centre_slice"], dtype=float),
            "final_snapshot_horizontal_source_axis_slice",
        )
    if "common_maximum_centre_slice" in case:
        return (
            np.asarray(case["common_maximum_centre_slice"], dtype=float),
            "maximum_over_transient_horizontal_source_axis_slice",
        )
    raise ValueError("case has neither a final nor maximum centre slice")


def score_hytunnel_resolved_transient(
    root: str | Path,
    audit_path: str | Path = DEFAULT_AUDIT,
    resolved_case_path: str | Path | None = None,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    started = time.perf_counter()
    directory = Path(root)
    source = Path(audit_path)
    audit = json.loads(source.read_text(encoding="utf-8"))
    if not audit.get("passed"):
        raise ValueError("resolved audit must have passed its domain-sensitivity gate")
    configuration = audit["configuration"]
    spacing = float(configuration["grid_spacing_m"])
    common_domain = tuple(float(v) for v in configuration["common_region_domain_m"])
    if len(common_domain) != 3:
        raise ValueError("common_region_domain_m must contain three dimensions")

    workbook_path = (
        directory / "reference/datasets/hytunnel-ps-4668554/raw" / WORKBOOK
    )
    dataset = read_hytunnel_concentration_workbook(workbook_path)
    matches = [item for item in dataset.configurations if item.sheet == SHEET]
    if len(matches) != 1:
        raise ValueError(f"expected one matched HyTunnel sheet, found {len(matches)}")
    matched_configuration = matches[0]
    points = [
        point for point in dataset.points
        if point.configuration_id == matched_configuration.configuration_id
        and point.statistic == "average"
    ]

    audit_published_case = next(
        case for case in audit["cases"] if case["name"] == "published_box"
    )
    published_case = audit_published_case
    replay_evidence = None
    endpoint_stationarity = None
    replay_is_frozen_configuration = False
    if resolved_case_path is not None:
        replay_source = Path(resolved_case_path)
        replay = json.loads(replay_source.read_text(encoding="utf-8"))
        if replay.get("schema") not in {
            "degali.hytunnel-resolved-case.v1",
            "degali.hytunnel-resolved-case.v2",
        }:
            raise ValueError("resolved case artifact has an unsupported schema")
        replay_case = replay["case"]
        if replay_case.get("name") != "published_box":
            raise ValueError("resolved case artifact must contain published_box")
        for key in ("grid_spacing_m", "duration_s", "time_steps", "flow_steps"):
            if replay["configuration"][key] != configuration[key]:
                raise ValueError(f"resolved case replay changes frozen {key}")
        frozen_maximum = np.asarray(
            audit_published_case["common_maximum_centre_slice"], dtype=float,
        )
        replay_maximum = np.asarray(
            replay_case["common_maximum_centre_slice"], dtype=float,
        )
        maximum_l1 = float(
            np.sum(np.abs(replay_maximum - frozen_maximum))
            / max(np.sum(np.abs(frozen_maximum)), 1.0e-12)
        )
        replay_configuration = replay["configuration"]
        turbulent_diffusivity_model = replay_configuration.get(
            "turbulent_diffusivity_model",
            "legacy_total_kinematic_viscosity_ratio",
        )
        handoff_limit = float(
            replay_configuration.get(
                "source_candidate_mean_density_departure_limit", 0.06,
            )
        )
        fan_inlet_profile = replay_configuration.get(
            "fan_inlet_profile", "uniform",
        )
        replay_is_frozen_configuration = bool(
            turbulent_diffusivity_model
            == "legacy_total_kinematic_viscosity_ratio"
            and abs(handoff_limit - 0.06) <= 1.0e-12
            and abs(float(replay_configuration.get(
                "effective_schmidt_number", 0.7,
            )) - 0.7) <= 1.0e-12
            and fan_inlet_profile == "uniform"
        )
        replay_evidence = {
            "artifact": str(replay_source),
            "maximum_slice_normalised_l1_vs_frozen_audit": maximum_l1,
            "frozen_replay_required": replay_is_frozen_configuration,
            "maximum_slice_reproduction_limit": (
                1.0e-10 if replay_is_frozen_configuration else None
            ),
            "maximum_slice_reproduced": (
                bool(maximum_l1 <= 1.0e-10)
                if replay_is_frozen_configuration else None
            ),
            "comparison_status": (
                "same_configuration_replay"
                if replay_is_frozen_configuration
                else "candidate_configuration_not_a_frozen_replay"
            ),
            "candidate_turbulent_diffusivity_model": (
                turbulent_diffusivity_model
            ),
            "candidate_handoff_density_departure_limit": handoff_limit,
            "candidate_fan_inlet_profile": fan_inlet_profile,
        }
        if (
            replay_is_frozen_configuration
            and not replay_evidence["maximum_slice_reproduced"]
        ):
            raise ValueError("resolved case replay does not reproduce the frozen maximum slice")
        published_case = replay_case
        endpoint_stationarity = replay["endpoint_stationarity"]
    source_position = tuple(float(v) for v in published_case["source_position_m"])
    nx = int(round(common_domain[0] / spacing))
    ny = int(round(common_domain[1] / spacing))
    x_centres = (np.arange(nx, dtype=float) + 0.5) * spacing
    y_centres = (np.arange(ny, dtype=float) + 0.5) * spacing

    included = [
        point for point in points
        if point.x_m >= source_position[0] - 1.0e-12
        and point.x_m <= common_domain[0] + 1.0e-12
        and 0.0 <= source_position[1] + point.y_m <= common_domain[1]
    ]
    if not included:
        raise ValueError("no measured points overlap the resolved downstream domain")

    observed = np.asarray(
        [point.concentration_vol_percent for point in included], dtype=float,
    )
    selected_x = np.asarray([point.x_m for point in included], dtype=float)
    selected_y = np.asarray([point.y_m for point in included], dtype=float)
    cases_to_score = [published_case]
    if resolved_case_path is None or replay_is_frozen_configuration:
        cases_to_score.extend(
            case for case in audit["cases"] if case["name"] != "published_box"
        )
    case_scores: dict[str, Any] = {}
    predictions: dict[str, np.ndarray] = {}
    operators: set[str] = set()
    for case in cases_to_score:
        field, operator = _case_slice(case)
        operators.add(operator)
        if field.shape != (ny, nx):
            raise ValueError(
                f"{case['name']} centre slice shape {field.shape} is not {(ny, nx)}"
            )
        predicted = np.maximum(0.0, np.asarray([
            bilinear_cell_centre_sample(
                field,
                x_centres,
                y_centres,
                point.x_m,
                source_position[1] + point.y_m,
                domain_m=(common_domain[0], common_domain[1]),
            )
            for point in included
        ]))
        predictions[case["name"]] = predicted
        geometric = statistics(observed, predicted, floor=0.10)
        subset_masks = {
            "centreline_y_zero": np.isclose(selected_y, 0.0, atol=1.0e-12),
            "near_axis_absolute_y_at_most_0p025m": np.abs(selected_y) <= 0.025000000001,
            "outer_absolute_y_at_least_0p075m": np.abs(selected_y) >= 0.074999999999,
            "observed_at_or_above_4vol": observed >= 4.0,
        }
        case_scores[case["name"]] = {
            "observation_operator": operator,
            "all_supported_points": _linear_metrics(observed, predicted),
            "spatial_subsets": {
                label: _linear_metrics(observed[mask], predicted[mask])
                for label, mask in subset_masks.items() if np.any(mask)
            },
            "by_x_station": {
                f"{station:g}": _linear_metrics(
                    observed[np.isclose(selected_x, station)],
                    predicted[np.isclose(selected_x, station)],
                )
                for station in np.unique(selected_x)
            },
            "hanna_statistics_floor_0p10_vol_pct": {
                "n": geometric.n,
                "mg_observed_over_predicted": geometric.mg,
                "vg": geometric.vg,
                "fac2": geometric.fac2,
                "fractional_bias": geometric.fb,
                "nmse": geometric.nmse,
                "acceptable": geometric.acceptable,
            },
            "threshold_classification": {
                "four_vol_pct": _threshold_classification(
                    observed, predicted, 4.0,
                ),
                "ten_vol_pct": _threshold_classification(
                    observed, predicted, 10.0,
                ),
            },
        }

    primary = predictions["published_box"]
    domain_pair = predictions.get("expanded_box")
    prediction_domain_difference = None
    domain_prediction_comparable = bool(
        domain_pair is not None
        and case_scores["published_box"]["observation_operator"]
        == case_scores["expanded_box"]["observation_operator"]
    )
    if domain_prediction_comparable:
        prediction_domain_difference = _linear_metrics(primary, domain_pair)

    upstream_count = sum(point.x_m < source_position[0] for point in points)
    beyond_count = sum(point.x_m > common_domain[0] for point in points)
    exact_operator_match = False
    qualified_endpoint_proxy = bool(
        endpoint_stationarity is not None
        and endpoint_stationarity["passed"]
        and case_scores["published_box"]["observation_operator"]
        == "final_snapshot_horizontal_source_axis_slice"
    )
    four = case_scores["published_box"]["threshold_classification"]["four_vol_pct"]
    diagnostic_gates = {
        "qualified_steady_endpoint_proxy": qualified_endpoint_proxy,
        "four_vol_pct_sensitivity_at_least_0p90": (
            four["sensitivity"] is not None and four["sensitivity"] >= 0.90
        ),
        "four_vol_pct_specificity_at_least_0p90": (
            four["specificity"] is not None and four["specificity"] >= 0.90
        ),
        "fac2_at_least_0p70": (
            case_scores["published_box"]["hanna_statistics_floor_0p10_vol_pct"][
                "fac2"
            ] >= 0.70
        ),
    }
    record: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "resolved_transient_diagnostic_passed"
            if all(diagnostic_gates.values())
            else "resolved_transient_validation_withheld"
        ),
        "promotion_allowed": bool(all(diagnostic_gates.values())),
        "dataset_doi": DATASET_DOI,
        "workbook": WORKBOOK,
        "sheet": SHEET,
        "model_artifact": str(source),
        "resolved_endpoint_replay": replay_evidence,
        "model_duration_s": float(configuration["duration_s"]),
        "model_time_steps": int(configuration["time_steps"]),
        "source_position_m": list(source_position),
        "measurement_operator": "published_steady_window_average",
        "model_operators": sorted(operators),
        "exact_operator_match": exact_operator_match,
        "qualified_steady_endpoint_proxy": qualified_endpoint_proxy,
        "endpoint_stationarity": endpoint_stationarity,
        "no_fit_statement": (
            "source, grid, transport, boundary treatment and scoring coordinates "
            "were frozen before this comparison; no measured concentration was "
            "used to select a coefficient"
        ),
        "point_accounting": {
            "published_average_point_count": len(points),
            "supported_downstream_point_count": len(included),
            "excluded_upstream_of_handoff_count": upstream_count,
            "excluded_beyond_common_domain_count": beyond_count,
            "included_x_stations_m": sorted({point.x_m for point in included}),
            "included_transverse_stations_m": sorted({point.y_m for point in included}),
        },
        "case_scores": case_scores,
        "published_vs_expanded_prediction_difference": prediction_domain_difference,
        "published_vs_expanded_prediction_comparable": (
            domain_prediction_comparable
        ),
        "diagnostic_gates": diagnostic_gates,
        "limitations": [
            (
                "the published workbook contains steady-window averages; the "
                "resolved endpoint is only a qualified steady proxy, not the same "
                "finite measurement window"
                if qualified_endpoint_proxy else
                "the published workbook contains steady-window averages but the "
                "selected resolved artifact does not provide a qualified steady "
                "endpoint slice"
            ),
            "measurements upstream of the declared density-departure handoff are "
            "outside the resolved model scope and are not extrapolated",
            "the common comparison domain ends at 2.52 m, so 3.0 and 3.5 m "
            "measurement stations are not scored",
            "the case contains no physical obstacle and is a ventilated hydrogen "
            "transport component test, not obstacle-wake validation",
        ],
        "interpretation": (
            "This is an honest overlap-domain diagnostic. Failure with a maximum-"
            "over-time model operator cannot be interpreted as the error of a "
            "matched steady-average operator, and it does not validate obstacle wakes."
        ),
        "runtime_s": time.perf_counter() - started,
    }
    arrays = {
        "selected_x_m": selected_x,
        "selected_y_m": selected_y,
        "observed_vol_pct": observed,
        "published_box_predicted_vol_pct": primary,
        "expanded_box_predicted_vol_pct": (
            domain_pair if domain_pair is not None else np.asarray([], dtype=float)
        ),
        "source_x_m": np.asarray([source_position[0]]),
    }
    return record, arrays


def plot_score(
    arrays: dict[str, np.ndarray], output: str | Path,
    model_operator: str | None = None,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = arrays["selected_x_m"]
    y = arrays["selected_y_m"]
    observed = arrays["observed_vol_pct"]
    predicted = arrays["published_box_predicted_vol_pct"]
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    def matrix(values: np.ndarray) -> np.ndarray:
        result = np.full((y_unique.size, x_unique.size), np.nan)
        for px, py, value in zip(x, y, values):
            ix = int(np.where(np.isclose(x_unique, px))[0][0])
            iy = int(np.where(np.isclose(y_unique, py))[0][0])
            result[iy, ix] = value
        return result

    observed_grid = matrix(observed)
    predicted_grid = matrix(predicted)
    residual_grid = predicted_grid - observed_grid
    vmax = max(float(np.nanmax(observed_grid)), float(np.nanmax(predicted_grid)), 4.0)

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.8), constrained_layout=True)
    endpoint = model_operator == "final_snapshot_horizontal_source_axis_slice"
    model_title = (
        "DEGALI stationary endpoint at 0.707 s"
        if endpoint else "DEGALI maximum over 0.707 s"
    )
    model_axis_label = (
        "predicted endpoint H$_2$ (vol%)"
        if endpoint else "predicted maximum H$_2$ (vol%)"
    )
    for axis, values, title in (
        (axes[0, 0], observed_grid, "HyTunnel steady-window average"),
        (axes[0, 1], predicted_grid, model_title),
    ):
        image = axis.pcolormesh(x_unique, y_unique, values, shading="nearest", vmin=0, vmax=vmax)
        axis.contour(x_unique, y_unique, values, levels=[4.0], colors="white", linewidths=1.8)
        axis.set_title(title)
        axis.set_xlabel("x from nozzle (m)")
        axis.set_ylabel("transverse y (m)")
        fig.colorbar(image, ax=axis, label="H$_2$ (vol%)")

    limit = max(vmax, float(np.max(observed)), float(np.max(predicted)))
    axes[1, 0].scatter(observed, predicted, c=x, cmap="viridis", edgecolor="black", s=40)
    axes[1, 0].plot([0, limit], [0, limit], "k--", linewidth=1.2)
    axes[1, 0].axvline(4.0, color="tab:red", linestyle=":")
    axes[1, 0].axhline(4.0, color="tab:red", linestyle=":")
    axes[1, 0].set_xlim(0, limit)
    axes[1, 0].set_ylim(0, limit)
    axes[1, 0].set_aspect("equal", adjustable="box")
    axes[1, 0].set_xlabel("observed average H$_2$ (vol%)")
    axes[1, 0].set_ylabel(model_axis_label)
    axes[1, 0].set_title("Paired supported points (colour = x)")

    residual_image = axes[1, 1].pcolormesh(
        x_unique, y_unique, residual_grid, shading="nearest", cmap="coolwarm",
        vmin=-float(np.nanmax(np.abs(residual_grid))),
        vmax=float(np.nanmax(np.abs(residual_grid))),
    )
    axes[1, 1].set_title(
        "Predicted endpoint − observed average"
        if endpoint else "Predicted maximum − observed average"
    )
    axes[1, 1].set_xlabel("x from nozzle (m)")
    axes[1, 1].set_ylabel("transverse y (m)")
    fig.colorbar(residual_image, ax=axes[1, 1], label="error (vol%-point)")
    fig.suptitle(
        "HyTunnel 1 mm, 5 g/s, 3.5 m/s co-flow — resolved transient overlap score\n"
        + (
            "Stationary endpoint proxy; accuracy gates determine promotion"
            if endpoint else
            "Operator mismatch is intentional and keeps the result withheld"
        ),
        fontsize=13,
    )
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument(
        "--resolved-case", type=Path,
        help="optional reproduced case containing a stationary final centre slice",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plot", type=Path, default=DEFAULT_PLOT)
    args = parser.parse_args()
    record, arrays = score_hytunnel_resolved_transient(
        args.root, args.audit, args.resolved_case,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    plot_score(
        arrays,
        args.plot,
        record["case_scores"]["published_box"]["observation_operator"],
    )
    metrics = record["case_scores"]["published_box"]["all_supported_points"]
    print(
        f"{record['status']}: n={metrics['pair_count']} "
        f"MAE={metrics['mae_vol_pct_point']:.4f} vol%p "
        f"RMSE={metrics['rmse_vol_pct_point']:.4f} vol%p; "
        f"wrote {args.output} and {args.plot}"
    )


if __name__ == "__main__":
    main()
