"""Diagnose which global spatial mechanism can still reduce TI021 dose MAPE.

The audit does not promote a fitted result.  It computes deliberately generous
upper bounds for two geometry-wide mechanisms while prohibiting a separate
gain or offset at each receptor:

* conservative redistribution inside exact mirrored crosswind pairs; and
* one complete precomputed wind-angle run per downwind group.

If even these bounds miss the target, a new side-split or curved-centreline
closure cannot honestly be justified as the route to ten percent.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog, minimize_scalar

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]


def _mape(prediction: np.ndarray, observed: np.ndarray) -> float:
    positive = observed > 0.0
    return float(np.mean(np.abs(prediction[positive] / observed[positive] - 1.0)) * 100.0)


def _load(path: Path, labels: list[str], observed: np.ndarray) -> tuple[dict, np.ndarray]:
    record = json.loads(path.read_text(encoding="utf-8"))
    rows = record.get("sensor_rows", [])
    if [row.get("label") for row in rows] != labels:
        raise ValueError(f"sensor order mismatch: {path}")
    candidate_observed = np.asarray(
        [row.get("observed_dose_vol_percent_s", np.nan) for row in rows], dtype=float
    )
    if not np.allclose(candidate_observed, observed, rtol=0.0, atol=1.0e-9):
        raise ValueError(f"observations mismatch: {path}")
    prediction = np.asarray(
        [row.get("predicted_dose_vol_percent_s", np.nan) for row in rows], dtype=float
    )
    if not np.all(np.isfinite(prediction)):
        raise ValueError(f"non-finite prediction: {path}")
    return record, prediction


def _pair_redistribution(
    prediction: np.ndarray, observed: np.ndarray, coordinates: np.ndarray
) -> tuple[np.ndarray, list[dict]]:
    adjusted = prediction.copy()
    used: set[int] = set()
    rows: list[dict] = []
    for left in range(len(observed)):
        if left in used or observed[left] <= 0.0:
            continue
        x, y, z = coordinates[left]
        if abs(x) <= 1.0e-9:
            continue
        candidates = np.flatnonzero(
            np.isclose(coordinates[:, 0], -x)
            & np.isclose(coordinates[:, 1], y)
            & np.isclose(coordinates[:, 2], z)
        )
        candidates = [int(index) for index in candidates if index not in used]
        if not candidates:
            continue
        right = candidates[0]
        if observed[right] <= 0.0:
            continue
        total = float(prediction[left] + prediction[right])
        breakpoints = np.clip(
            np.asarray([0.0, total, observed[left], total - observed[right]]),
            0.0,
            total,
        )
        errors = (
            np.abs(breakpoints / observed[left] - 1.0)
            + np.abs((total - breakpoints) / observed[right] - 1.0)
        )
        best = int(np.argmin(errors))
        adjusted[left] = breakpoints[best]
        adjusted[right] = total - breakpoints[best]
        used.update((left, right))
        rows.append(
            {
                "left_index": left,
                "right_index": right,
                "local_y_m": float(y),
                "height_m": float(z),
                "preserved_predicted_pair_total": total,
                "observed_pair_total": float(observed[left] + observed[right]),
            }
        )
    return adjusted, rows


def _angle_group_bound(
    models: list[dict], observed: np.ndarray, group_values: np.ndarray
) -> tuple[np.ndarray, list[dict]]:
    assembled = np.zeros_like(observed)
    choices: list[dict] = []
    for value in sorted(set(float(item) for item in group_values)):
        mask = np.isclose(group_values, value) & (observed > 0.0)
        if not np.any(mask):
            continue
        scored = []
        for model in models:
            error = float(np.mean(np.abs(model["prediction"][mask] / observed[mask] - 1.0)))
            scored.append((error, model))
        error, selected = min(scored, key=lambda item: item[0])
        assembled[np.isclose(group_values, value)] = selected["prediction"][
            np.isclose(group_values, value)
        ]
        choices.append(
            {
                "group": value,
                "positive_sensor_count": int(np.count_nonzero(mask)),
                "selected_angle_deg": selected["angle_deg"],
                "group_mape_percent": error * 100.0,
                "source": selected["path"],
            }
        )
    return assembled, choices


def _conservative_group_redistribution(
    prediction: np.ndarray,
    observed: np.ndarray,
    keys: list[tuple[float, ...]],
) -> tuple[np.ndarray, list[dict]]:
    """Minimize percentage error inside each group at fixed predicted sum."""
    adjusted = prediction.copy()
    rows: list[dict] = []
    for key in sorted(set(keys)):
        indices = np.asarray(
            [index for index, candidate in enumerate(keys)
             if candidate == key and observed[index] > 0.0], dtype=int
        )
        if indices.size < 2:
            continue
        total = float(prediction[indices].sum())
        count = int(indices.size)
        # Variables: allocated non-negative predictions then absolute relative
        # residuals.  This is the exact L1 percentage-error optimum.
        objective = np.r_[np.zeros(count), np.full(count, 1.0 / count)]
        ratio = np.diag(1.0 / observed[indices])
        upper = np.block([[ratio, -np.eye(count)], [-ratio, -np.eye(count)]])
        result = linprog(
            objective,
            A_ub=upper,
            b_ub=np.r_[np.ones(count), -np.ones(count)],
            A_eq=np.r_[np.ones(count), np.zeros(count)][None, :],
            b_eq=np.asarray([total]),
            bounds=[(0.0, None)] * (2 * count),
            method="highs",
        )
        if not result.success:
            raise RuntimeError(f"group redistribution failed for {key}: {result.message}")
        adjusted[indices] = result.x[:count]
        rows.append(
            {
                "key": list(key),
                "positive_sensor_count": count,
                "preserved_predicted_sum": total,
                "observed_sum": float(observed[indices].sum()),
            }
        )
    return adjusted, rows


def _one_scale_per_group(
    prediction: np.ndarray,
    observed: np.ndarray,
    keys: list[tuple[float, ...]],
) -> tuple[np.ndarray, list[dict]]:
    adjusted = prediction.copy()
    rows: list[dict] = []
    for key in sorted(set(keys)):
        indices = np.asarray(
            [index for index, candidate in enumerate(keys)
             if candidate == key and observed[index] > 0.0], dtype=int
        )
        if indices.size == 0:
            continue
        ratio = prediction[indices] / observed[indices]
        result = minimize_scalar(
            lambda scale: float(np.mean(np.abs(scale * ratio - 1.0))),
            bounds=(0.0, 10.0), method="bounded",
            options={"xatol": 1.0e-12},
        )
        adjusted[indices] = prediction[indices] * float(result.x)
        rows.append(
            {
                "key": list(key),
                "positive_sensor_count": int(indices.size),
                "optimal_scale": float(result.x),
                "group_mape_percent": float(result.fun * 100.0),
            }
        )
    return adjusted, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--angle-run", type=Path, action="append", default=[],
        help="same-physics run differing only in wind angle; may be repeated",
    )
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
    labels = [row["label"] for row in baseline["sensor_rows"]]
    observed = np.asarray(
        [row["observed_dose_vol_percent_s"] for row in baseline["sensor_rows"]]
    )
    candidate_record, candidate = _load(args.candidate, labels, observed)
    candidate_observation_mode = candidate_record.get(
        "observation_species_mode", "r12_component"
    )
    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
    coordinates = np.asarray(
        [(sensor.local_x_m, sensor.local_y_m, sensor.height_m)
         for sensor in dataset.dose_sensors], dtype=float
    )

    redistributed, pairs = _pair_redistribution(candidate, observed, coordinates)

    angle_by_value: dict[float, dict] = {}
    angle_paths = (
        args.angle_run
        if args.angle_run
        else sorted((ROOT / "outputs").glob("candidate-angle-*.json"))
    )
    for path in angle_paths:
        try:
            record, prediction = _load(path, labels, observed)
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        if record.get("observation_species_mode", "r12_component") != candidate_observation_mode:
            continue
        angle = float(record["wind_angle_deg"])
        item = {
            "angle_deg": angle,
            "path": str(path.resolve()),
            "prediction": prediction,
            "mape_percent": _mape(prediction, observed),
        }
        if angle not in angle_by_value or item["mape_percent"] < angle_by_value[angle]["mape_percent"]:
            angle_by_value[angle] = item
    models = list(angle_by_value.values())
    if models:
        per_arc, arc_choices = _angle_group_bound(models, observed, coordinates[:, 1])
    else:
        per_arc, arc_choices = None, []
    broad_group = np.where(
        coordinates[:, 1] < 150.0, 0.0,
        np.where(coordinates[:, 1] < 300.0, 1.0, 2.0),
    )
    if models:
        per_broad, broad_choices = _angle_group_bound(models, observed, broad_group)
    else:
        per_broad, broad_choices = None, []

    tower_keys = [
        (float(x), float(y)) for x, y in coordinates[:, :2]
    ]
    arc_keys = [(float(y),) for y in coordinates[:, 1]]
    vertical_redistributed, vertical_groups = _conservative_group_redistribution(
        candidate, observed, tower_keys
    )
    tower_scaled, tower_scale_rows = _one_scale_per_group(
        candidate, observed, tower_keys
    )
    arc_scaled, arc_scale_rows = _one_scale_per_group(
        candidate, observed, arc_keys
    )

    relative = np.abs(candidate / np.where(observed > 0.0, observed, 1.0) - 1.0)
    worst = np.argsort(np.where(observed > 0.0, relative, -1.0))[-15:][::-1]
    report = {
        "schema": "degali.ti021-residual-structure-audit.v1",
        "status": "diagnostic_bounds_not_promoted_models",
        "candidate": str(args.candidate.resolve()),
        "observation_species_mode": candidate_observation_mode,
        "candidate_mape_percent": _mape(candidate, observed),
        "conservative_mirrored_pair_redistribution": {
            "pair_count": len(pairs),
            "mape_percent": _mape(redistributed, observed),
            "target_10_percent_reached": _mape(redistributed, observed) <= 10.0,
            "contract": "preserve each model-predicted mirrored-pair total; optimize only its left/right split",
        },
        "one_angle_per_exact_downwind_arc_bound": {
            "mape_percent": None if per_arc is None else _mape(per_arc, observed),
            "target_10_percent_reached": False if per_arc is None else _mape(per_arc, observed) <= 10.0,
            "choices": arc_choices,
            "status": "unavailable_no_matching_observable_runs" if per_arc is None else "available",
        },
        "one_angle_per_three_downwind_bands_bound": {
            "mape_percent": None if per_broad is None else _mape(per_broad, observed),
            "target_10_percent_reached": False if per_broad is None else _mape(per_broad, observed) <= 10.0,
            "choices": broad_choices,
            "status": "unavailable_no_matching_observable_runs" if per_broad is None else "available",
        },
        "conservative_vertical_redistribution_within_each_tower": {
            "tower_count": len(vertical_groups),
            "mape_percent": _mape(vertical_redistributed, observed),
            "target_10_percent_reached": _mape(vertical_redistributed, observed) <= 10.0,
            "contract": "preserve each horizontal tower's predicted sum and optimize only its allocation among heights",
        },
        "one_scalar_per_exact_downwind_arc_bound": {
            "mape_percent": _mape(arc_scaled, observed),
            "target_10_percent_reached": _mape(arc_scaled, observed) <= 10.0,
            "groups": arc_scale_rows,
        },
        "one_scalar_per_horizontal_tower_bound": {
            "mape_percent": _mape(tower_scaled, observed),
            "target_10_percent_reached": _mape(tower_scaled, observed) <= 10.0,
            "group_count": len(tower_scale_rows),
            "warning": "high-flexibility diagnostic; not an admissible model",
        },
        "worst_candidate_sensors": [
            {
                "label": labels[index],
                "local_x_m": float(coordinates[index, 0]),
                "local_y_m": float(coordinates[index, 1]),
                "height_m": float(coordinates[index, 2]),
                "observed": float(observed[index]),
                "predicted": float(candidate[index]),
                "absolute_percentage_error": float(relative[index] * 100.0),
            }
            for index in worst
        ],
        "interpretation_guard": (
            "These are optimistic mechanism bounds fitted on TI021 and are not validation results."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in report.items() if key != "worst_candidate_sensors"}, indent=2))


if __name__ == "__main__":
    main()
