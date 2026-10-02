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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=None,
                        help="local reference/spadeadam directory")
    args = parser.parse_args()
    trial = next(item for item in load(args.root) if item.test == 6)
    readings = [
        reading for reading in trial.readings
        if abs(reading.radius - 30.0) < 0.5 and not reading.over_range
    ]
    points, observed = _coordinates(trial, readings)
    projection = project_lh2_jet_to_sensors(
        rate=trial.rate, wind=trial.wind_low, height=0.5,
        orifice=trial.orifice, storage_pressure=1.013 + trial.line_pressure,
        ambient_temperature=277.15, relative_humidity=90.0,
        points_m=points, max_distance=30.0,
    )
    predicted = [100.0 * value for value in projection.mole_fractions]
    payload = {
        "test": 6,
        "sensor_count": len(readings),
        "observed_arc_max_vol_pct": max(observed),
        "projected_sensor_max_vol_pct": max(predicted),
        "observed_over_projected": max(observed) / max(predicted),
        "warnings": list(projection.warnings),
        "source_files_remain_external": True,
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
