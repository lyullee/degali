"""Pointwise scoring for published ExJet obstacle concentration profiles.

Inputs are explicit figure-digitized hydrogen mole fractions. Metrics are
reported as vol%-point differences under the ideal-gas mole/volume-fraction
equivalence, and the 4 vol% result is a sparse point classification only; it
is not a cloud-volume or continuous contour metric.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Mapping, Sequence

import numpy as np


_OBSERVATION_FIELDS = {
    "observation_id", "station_x_m", "radial_coordinate_m",
    "h2_molar_fraction", "quantity_kind",
    "concentration_uncertainty_mole_fraction",
}
_PREDICTION_FIELDS = {"observation_id", "h2_molar_fraction_predicted"}
_FOUR_VOL_PERCENT_MOLE_FRACTION = 0.04


def _finite_fraction(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a numeric mole fraction")
    result = float(value)
    if not math.isfinite(result) or not 0.0 <= result <= 1.0:
        raise ValueError(f"{name} must be finite and in [0, 1]")
    return result


def _trilinear_sample_fluid(
    x: Sequence[float],
    y: Sequence[float],
    z: Sequence[float],
    field,
    solid,
    point: tuple[float, float, float],
) -> float | None:
    """Sample a cell-centred field; refuse out-of-domain and solid stencils."""
    brackets = []
    for coordinates, value in zip((x, y, z), point):
        if value < coordinates[0] or value > coordinates[-1] or len(coordinates) < 2:
            return None
        upper = min(max(int(np.searchsorted(
            coordinates, value, side="right"
        )), 1), len(coordinates) - 1)
        lower = upper - 1
        fraction = (value - coordinates[lower]) / (
            coordinates[upper] - coordinates[lower]
        )
        brackets.append((lower, upper, fraction))
    bx, by, bz = brackets
    terms = []
    for k, wz in ((bz[0], 1.0 - bz[2]), (bz[1], bz[2])):
        for j, wy in ((by[0], 1.0 - by[2]), (by[1], by[2])):
            for i, wx in ((bx[0], 1.0 - bx[2]), (bx[1], bx[2])):
                weight = wx * wy * wz
                if weight > 1.0e-15:
                    terms.append((k, j, i, weight))
    if any(solid[k, j, i] for k, j, i, _weight in terms):
        return None
    return math.fsum(weight * float(field[k, j, i]) for k, j, i, weight in terms)


def sample_exjet_spatial_profile_predictions(
    result,
    observations: Sequence[Mapping[str, object]],
    *,
    nozzle_position_m: tuple[float, float, float],
) -> dict[str, object]:
    """Sample a resolved 3-D result at ExJet's vertical measurement plane.

    The paper's coordinate description places the radial sensors in the
    vertical plane through a horizontal jet. The explicit transform used here
    is therefore ``(x0 + station, y0, z0 + signed_radial_coordinate)``. The
    caller must provide the nozzle centre in solver coordinates; no hidden
    registration offset or nearest-cell fallback is applied. The declared
    radial-coordinate uncertainty is propagated by sampling both vertical
    endpoints and every cell-centre knot within the interval. For the
    piecewise-linear interpolation along z, these samples bound its extrema.
    The resulting min/max is a bounded sensitivity range, not a probabilistic
    confidence interval. If any sample is invalid, the central prediction is
    retained but that point's uncertainty is marked incomplete.
    """
    if not isinstance(nozzle_position_m, tuple) or len(nozzle_position_m) != 3:
        raise TypeError("nozzle_position_m must be a three-value tuple")
    origin = tuple(float(value) for value in nozzle_position_m)
    if any(not math.isfinite(value) for value in origin):
        raise ValueError("nozzle_position_m must be finite")
    field = getattr(result, "hydrogen_mole_fraction", None)
    if field is None:
        raise ValueError("resolved result does not contain an H2 mole-fraction field")
    x = tuple(float(value) for value in result.x_m)
    y = tuple(float(value) for value in result.y_m)
    z = tuple(float(value) for value in result.z_m)
    field = np.asarray(field, dtype=float)
    solid = np.asarray(result.solid_mask, dtype=bool)
    if field.shape != (len(z), len(y), len(x)) or solid.shape != field.shape:
        raise ValueError("resolved H2 field and solid mask must match (z, y, x) grid")
    predictions = []
    failed = []
    for observation in observations:
        key = str(observation.get("observation_id", ""))
        try:
            station = float(observation["station_x_m"])
            radial = float(observation["radial_coordinate_m"])
            radial_uncertainty = float(observation["coordinate_uncertainty_m"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                "each observation requires numeric station, radial coordinate, "
                "and radial-coordinate uncertainty"
            ) from exc
        if (
            not key or not math.isfinite(station) or not math.isfinite(radial)
            or not math.isfinite(radial_uncertainty) or radial_uncertainty <= 0.0
        ):
            raise ValueError("observation IDs and coordinates must be finite/non-empty")
        sample_point = (origin[0] + station, origin[1], origin[2] + radial)
        value = _trilinear_sample_fluid(
            x, y, z, field, solid, sample_point,
        )
        if value is None:
            failed.append({
                "observation_id": key,
                "sample_coordinate_solver_m": list(sample_point),
                "reason": "outside_domain_or_interpolation_stencil_intersects_solid",
            })
            continue
        radial_lower = radial - radial_uncertainty
        radial_upper = radial + radial_uncertainty
        vertical_lower = origin[2] + radial_lower
        vertical_upper = origin[2] + radial_upper
        vertical_samples = sorted({
            vertical_lower,
            vertical_upper,
            *(
                vertical_grid for vertical_grid in z
                if vertical_lower <= vertical_grid <= vertical_upper
            ),
        })
        radial_endpoint_values = [
            _trilinear_sample_fluid(
                x, y, z, field, solid,
                (origin[0] + station, origin[1], vertical_sample),
            )
            for vertical_sample in vertical_samples
        ]
        coordinate_bounds_complete = all(
            endpoint_value is not None for endpoint_value in radial_endpoint_values
        )
        if coordinate_bounds_complete:
            coordinate_values = [value, *radial_endpoint_values]
            coordinate_lower = min(coordinate_values)
            coordinate_upper = max(coordinate_values)
        else:
            coordinate_lower = coordinate_upper = None
        predictions.append({
            "observation_id": key,
            "h2_molar_fraction_predicted": value,
            "h2_molar_fraction_coordinate_lower": coordinate_lower,
            "h2_molar_fraction_coordinate_upper": coordinate_upper,
            "coordinate_uncertainty_sampling_status": (
                "complete" if coordinate_bounds_complete
                else "incomplete_interval_sample_outside_or_solid_stencil"
            ),
            "sample_coordinate_solver_m": list(sample_point),
        })
    time_values = tuple(float(value) for value in result.time_s)
    return {
        "schema": "degali.exjet-resolved-field-sampling.v1",
        "coordinate_transform": "x=x0+station; y=y0; z=z0+signed_radial_coordinate",
        "nozzle_position_solver_m": list(origin),
        "sampled_field": "hydrogen_mole_fraction",
        "sampled_time_s": time_values[-1] if time_values else None,
        "sampled_field_is_final_snapshot": True,
        "coordinate_uncertainty_propagated": bool(predictions) and all(
            item["coordinate_uncertainty_sampling_status"] == "complete"
            for item in predictions
        ),
        "interpolation": "trilinear cell-centred; any solid-stencil point withheld",
        "predictions": predictions,
        "failed_points": failed,
        "complete_pairing": not failed and len(predictions) == len(observations),
        "complete_coordinate_uncertainty_propagation": bool(predictions) and all(
            item["coordinate_uncertainty_sampling_status"] == "complete"
            for item in predictions
        ),
        "quantitative_promotion_allowed": False,
        "scope_warning": "spatial point comparison only; ExJet digitized data do not provide a transient history",
    }


def score_exjet_profile_points(
    observations: Sequence[Mapping[str, object]],
    predictions: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Score predictions paired by observation ID; never interpolate silently.

    Required observation/prediction ID sets must match exactly. Predictions
    are paired at the published station/radial coordinates represented by the
    observation IDs; callers must perform and record the field sampling step.
    """
    if not observations or not predictions:
        raise ValueError("observations and predictions must both be non-empty")

    observed_by_id: dict[str, Mapping[str, object]] = {}
    for row in observations:
        missing = sorted(_OBSERVATION_FIELDS - set(row))
        if missing:
            raise ValueError(f"observation is missing fields: {missing}")
        key = str(row["observation_id"])
        if not key or key in observed_by_id:
            raise ValueError("observation IDs must be non-empty and unique")
        observed_by_id[key] = row

    predicted_by_id: dict[str, Mapping[str, object]] = {}
    for row in predictions:
        missing = sorted(_PREDICTION_FIELDS - set(row))
        if missing:
            raise ValueError(f"prediction is missing fields: {missing}")
        key = str(row["observation_id"])
        if not key or key in predicted_by_id:
            raise ValueError("prediction IDs must be non-empty and unique")
        predicted_by_id[key] = row

    missing_predictions = sorted(set(observed_by_id) - set(predicted_by_id))
    unexpected_predictions = sorted(set(predicted_by_id) - set(observed_by_id))
    if missing_predictions or unexpected_predictions:
        raise ValueError(
            "prediction IDs must match observation IDs exactly; "
            f"missing={missing_predictions}, unexpected={unexpected_predictions}"
        )

    paired: list[dict[str, object]] = []
    for key, observation in observed_by_id.items():
        predicted = predicted_by_id[key]
        observed_fraction = _finite_fraction(
            observation["h2_molar_fraction"], "observed h2_molar_fraction",
        )
        predicted_fraction = _finite_fraction(
            predicted["h2_molar_fraction_predicted"],
            "predicted h2_molar_fraction",
        )
        uncertainty = float(observation["concentration_uncertainty_mole_fraction"])
        coordinate_uncertainty = float(observation["coordinate_uncertainty_m"])
        station = float(observation["station_x_m"])
        radial = float(observation["radial_coordinate_m"])
        if not all(math.isfinite(value) for value in (
            uncertainty, coordinate_uncertainty, station, radial,
        )):
            raise ValueError("observation coordinates and uncertainty must be finite")
        if uncertainty <= 0.0 or coordinate_uncertainty <= 0.0 or station < 0.0:
            raise ValueError("uncertainty must be positive and station non-negative")
        error_vol_pct_point = 100.0 * (predicted_fraction - observed_fraction)
        observed_vol_pct = 100.0 * observed_fraction
        predicted_vol_pct = 100.0 * predicted_fraction
        uncertainty_vol_pct_point = 100.0 * uncertainty
        raw_predicted_lower = predicted.get("h2_molar_fraction_coordinate_lower")
        raw_predicted_upper = predicted.get("h2_molar_fraction_coordinate_upper")
        coordinate_bounds_available = (
            raw_predicted_lower is not None and raw_predicted_upper is not None
        )
        if (raw_predicted_lower is None) != (raw_predicted_upper is None):
            raise ValueError("both predicted coordinate bounds must be supplied together")
        predicted_lower_raw = (
            raw_predicted_lower if coordinate_bounds_available else predicted_fraction
        )
        predicted_upper_raw = (
            raw_predicted_upper if coordinate_bounds_available else predicted_fraction
        )
        predicted_lower = _finite_fraction(
            predicted_lower_raw, "predicted coordinate lower mole fraction",
        )
        predicted_upper = _finite_fraction(
            predicted_upper_raw, "predicted coordinate upper mole fraction",
        )
        if predicted_lower > predicted_upper:
            raise ValueError("predicted coordinate bounds must be ordered")
        observed_lower = max(0.0, observed_fraction - uncertainty)
        observed_upper = min(1.0, observed_fraction + uncertainty)
        error_lower_vol_pct_point = 100.0 * (predicted_lower - observed_upper)
        error_upper_vol_pct_point = 100.0 * (predicted_upper - observed_lower)
        if observed_lower >= _FOUR_VOL_PERCENT_MOLE_FRACTION:
            observed_threshold_class = "above"
        elif observed_upper < _FOUR_VOL_PERCENT_MOLE_FRACTION:
            observed_threshold_class = "below"
        else:
            observed_threshold_class = "indeterminate_within_plot_reading_uncertainty"
        if not coordinate_bounds_available:
            predicted_threshold_class = "point_estimate_only"
        elif predicted_lower >= _FOUR_VOL_PERCENT_MOLE_FRACTION:
            predicted_threshold_class = "above"
        elif predicted_upper < _FOUR_VOL_PERCENT_MOLE_FRACTION:
            predicted_threshold_class = "below"
        else:
            predicted_threshold_class = "indeterminate_within_coordinate_uncertainty"
        paired.append({
            "observation_id": key,
            "station_x_m": station,
            "radial_coordinate_m": radial,
            "coordinate_uncertainty_m": coordinate_uncertainty,
            "quantity_kind": str(observation["quantity_kind"]),
            "observed_vol_percent": observed_vol_pct,
            "predicted_vol_percent": predicted_vol_pct,
            "error_vol_percent_point": error_vol_pct_point,
            "absolute_error_vol_percent_point": abs(error_vol_pct_point),
            "reading_uncertainty_vol_percent_point": uncertainty_vol_pct_point,
            "predicted_coordinate_lower_vol_percent": 100.0 * predicted_lower,
            "predicted_coordinate_upper_vol_percent": 100.0 * predicted_upper,
            "coordinate_uncertainty_bounds_available": coordinate_bounds_available,
            "combined_error_interval_lower_vol_percent_point": error_lower_vol_pct_point,
            "combined_error_interval_upper_vol_percent_point": error_upper_vol_pct_point,
            "compatible_with_combined_plot_and_coordinate_interval": (
                error_lower_vol_pct_point <= 0.0 <= error_upper_vol_pct_point
            ) if coordinate_bounds_available else None,
            "within_digitization_reading_uncertainty": (
                abs(error_vol_pct_point) <= uncertainty_vol_pct_point
            ),
            "observed_4_vol_percent_classification": observed_threshold_class,
            "predicted_4_vol_percent_classification": predicted_threshold_class,
        })

    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in paired:
        groups[str(row["quantity_kind"])].append(row)

    def summarize(rows: Sequence[dict[str, object]]) -> dict[str, object]:
        errors = [float(row["error_vol_percent_point"]) for row in rows]
        abs_errors = [abs(error) for error in errors]
        tp = sum(
            row["observed_4_vol_percent_classification"] == "above"
            and row["predicted_4_vol_percent_classification"] == "above"
            for row in rows
        )
        tn = sum(
            row["observed_4_vol_percent_classification"] == "below"
            and row["predicted_4_vol_percent_classification"] == "below"
            for row in rows
        )
        fp = sum(
            row["observed_4_vol_percent_classification"] == "below"
            and row["predicted_4_vol_percent_classification"] == "above"
            for row in rows
        )
        fn = sum(
            row["observed_4_vol_percent_classification"] == "above"
            and row["predicted_4_vol_percent_classification"] == "below"
            for row in rows
        )
        indeterminate = sum(
            row["observed_4_vol_percent_classification"]
            == "indeterminate_within_plot_reading_uncertainty"
            for row in rows
        )
        predicted_indeterminate = sum(
            row["predicted_4_vol_percent_classification"]
            == "indeterminate_within_coordinate_uncertainty"
            for row in rows
        )
        confusion_matrix_n = sum(
            row["observed_4_vol_percent_classification"] in {"above", "below"}
            and row["predicted_4_vol_percent_classification"] in {"above", "below"}
            for row in rows
        )
        combined_compatibility_rows = [
            row["compatible_with_combined_plot_and_coordinate_interval"]
            for row in rows
            if row["compatible_with_combined_plot_and_coordinate_interval"] is not None
        ]
        return {
            "n": len(rows),
            "mean_bias_vol_percent_point": math.fsum(errors) / len(errors),
            "mae_vol_percent_point": math.fsum(abs_errors) / len(errors),
            "rmse_vol_percent_point": math.sqrt(
                math.fsum(error * error for error in errors) / len(errors)
            ),
            "maximum_absolute_error_vol_percent_point": max(abs_errors),
            "within_plot_reading_uncertainty_fraction": sum(
                bool(row["within_digitization_reading_uncertainty"]) for row in rows
            ) / len(rows),
            "coordinate_uncertainty_propagation_fraction": sum(
                bool(row["coordinate_uncertainty_bounds_available"]) for row in rows
            ) / len(rows),
            "compatible_with_combined_plot_and_coordinate_interval_fraction": (
                sum(bool(value) for value in combined_compatibility_rows)
                / len(combined_compatibility_rows)
                if combined_compatibility_rows else None
            ),
            "combined_interval_comparison_n": len(combined_compatibility_rows),
            "pointwise_4_vol_percent_classification": {
                "threshold_basis": "mole_fraction >= 0.04; ideal-gas mole/volume equivalence",
                "true_positive": tp,
                "true_negative": tn,
                "false_positive": fp,
                "false_negative": fn,
                "indeterminate_observation_count": indeterminate,
                "indeterminate_prediction_count": predicted_indeterminate,
                "confusion_matrix_n": confusion_matrix_n,
                "scope_warning": "sparse digitized points only; not a continuous hazard contour or volume score",
            },
        }

    return {
        "schema": "degali.exjet-spatial-profile-score.v1",
        "comparison_scope": "sparse steady/published spatial concentration profiles",
        "input_concentration": "H2 mole fraction",
        "reporting_concentration": "vol% under ideal-gas mole/volume-fraction equivalence",
        "digitized_observations_are_raw_data": False,
        "dynamic_validation_supported": False,
        "coordinate_uncertainty_propagated": all(
            bool(row["coordinate_uncertainty_bounds_available"]) for row in paired
        ),
        "quantitative_promotion_allowed": False,
        "promotion_note": "Apply the separate ExJet readiness audit for source and geometry gates; this score does not validate transient response.",
        "n": len(paired),
        "by_quantity": {
            name: summarize(rows) for name, rows in sorted(groups.items())
        },
        "points": paired,
    }
