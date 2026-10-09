"""Audit the FFI/DNV Test 6 exact sensor operator from local public extracts.

The CSVs remain outside the package. This tool reads a user-held
``reference/spadeadam`` directory, projects the jet onto the reported sensor
coordinates, and prints aggregate values only.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from degali import project_lh2_jet_to_sensors
from degali.validation.spadeadam import load


def _coordinates(trial, readings):
    direction = math.radians(trial.wind_direction)
    down_east, down_north = -math.sin(direction), -math.cos(direction)
    right_east, right_north = down_north, -down_east
    points, observed = [], []
    for reading in readings:
        bearing = math.radians(reading.bearing)
        east = reading.radius * math.sin(bearing)
        north = reading.radius * math.cos(bearing)
        points.append((
            east * down_east + north * down_north,
            east * right_east + north * right_north,
            reading.height,
        ))
        observed.append(reading.peak)
    return points, observed


def _operator_table(readings, predicted, *, group: str) -> list[dict[str, object]]:
    """Summarise the finite sensor operator without emitting raw rows.

    ``group="radius"`` is the arc-max operator: the reported peak and the
    exact-coordinate prediction are maximised over all retained heights and
    bearings on one arc.  ``group="height"`` keeps the same finite sensor
    projection but reports one row per sensor height, which makes the
    sensor-height effect visible without copying the source table into the
    generated artifact.
    """
    if group not in {"radius", "height"}:
        raise ValueError("group must be 'radius' or 'height'")
    if len(readings) != len(predicted):
        raise ValueError("readings and predicted must have equal lengths")

    buckets: dict[float, dict[str, object]] = {}
    for reading, value in zip(readings, predicted):
        key_value = float(reading.radius if group == "radius" else reading.height)
        key_value = round(key_value, 6)
        bucket = buckets.setdefault(
            key_value,
            {"observed": [], "predicted": []},
        )
        bucket["observed"].append(float(reading.peak))
        bucket["predicted"].append(float(value))

    rows: list[dict[str, object]] = []
    for key_value in sorted(buckets):
        bucket = buckets[key_value]
        observed = max(bucket["observed"])
        projected = max(bucket["predicted"])
        rows.append({
            ("radius_m" if group == "radius" else "height_m"): key_value,
            "sensor_count": len(bucket["observed"]),
            "observed_arc_max_vol_pct": observed,
            "projected_sensor_max_vol_pct": projected,
            "observed_over_projected": (
                None if projected <= 0.0 else observed / projected
            ),
            "observed_is_lower_bound": True,
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=None,
                        help="local reference/spadeadam directory")
    parser.add_argument("--radius", type=float, default=30.0,
                        help="arc radius to audit in metres (default: 30)")
    args = parser.parse_args()
    trial = next(item for item in load(args.root) if item.test == 6)
    readings = [
        reading for reading in trial.readings
        if abs(reading.radius - args.radius) < 0.5 and not reading.over_range
    ]
    if not readings:
        raise SystemExit(f"no usable Test 6 readings near radius {args.radius:g} m")
    points, _observed = _coordinates(trial, readings)
    projection = project_lh2_jet_to_sensors(
        rate=trial.rate, wind=trial.wind_low, height=0.5,
        orifice=trial.orifice, storage_pressure=1.013 + trial.line_pressure,
        ambient_temperature=277.15, relative_humidity=90.0,
        points_m=points, max_distance=max(1.0, args.radius),
    )
    predicted = [100.0 * value for value in projection.mole_fractions]
    arc_max_table = _operator_table(readings, predicted, group="radius")
    sensor_height_table = _operator_table(readings, predicted, group="height")
    arc_row = arc_max_table[0]
    payload = {
        "test": 6,
        "radius_m": args.radius,
        "sensor_count": arc_row["sensor_count"],
        "observed_arc_max_vol_pct": arc_row["observed_arc_max_vol_pct"],
        "projected_sensor_max_vol_pct": arc_row["projected_sensor_max_vol_pct"],
        "observed_over_projected": arc_row["observed_over_projected"],
        "operator_definition": {
            "arc_max": (
                "maximum over retained sensor heights and bearings at one radius"
            ),
            "sensor_height": (
                "exact-coordinate projection grouped by reported sensor height"
            ),
            "observed_is_lower_bound": True,
        },
        "arc_max_table": arc_max_table,
        "sensor_height_table": sensor_height_table,
        "warnings": list(projection.warnings),
        "source_files_remain_external": True,
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
