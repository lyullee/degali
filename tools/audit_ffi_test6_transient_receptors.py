"""Replay FFI/DNV Test 6 receptors from a user-supplied wind history.

The public reports contain plotted time histories but no machine-readable raw
channel export. This tool requires an explicit CSV instead of inventing a
trace from published summary values. Expensive steady plume solutions are
created once at the requested table winds and linearly interpolated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from degali.addons import (
    FixedReceptor,
    SteadyPlumeTable,
    WindHistory,
    replay_fixed_receptors,
)
from degali.validation.nearfield import STEP, Trajectory, hydrogen_jet
from degali.validation.spadeadam import load


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "reference" / "spadeadam"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_wind(path: Path) -> WindHistory:
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append((
                float(row["time_s"]),
                float(row["wind_speed_ms"]),
                float(row["wind_direction_from_deg"]),
            ))
    if not rows:
        raise ValueError("wind CSV contains no samples")
    return WindHistory(
        time_s=[row[0] for row in rows],
        speed_m_s=[row[1] for row in rows],
        direction_from_deg=[row[2] for row in rows],
    )


def _steady_sampler(trial, wind: float, maximum_distance: float):
    setup, initial = hydrogen_jet(
        rate=trial.rate,
        diameter=trial.orifice,
        wind=wind,
        height=0.5,
        ambient_temperature=277.15,
        relative_humidity=90.0,
        ambient_pressure=101325.0,
        storage_pressure_barg=trial.line_pressure,
        corrections=True,
        ground_effect=False,
    )
    result = setup.run(
        initial,
        distmx=STEP,
        smax=max(40.0, 2.5 * maximum_distance),
    )
    trajectory = Trajectory(setup.th.table, result.rows)
    if not trajectory.ok:
        raise RuntimeError(f"steady trajectory failed at {wind:g} m/s")

    def sample(alongwind_m: float, crosswind_m: float, height_m: float) -> float:
        if trajectory.at(alongwind_m) is None:
            raise ValueError(
                f"steady trajectory at {wind:g} m/s does not reach "
                f"{alongwind_m:g} m"
            )
        return max(
            0.0,
            float(trajectory.concentration_at(
                alongwind_m, crosswind_m, height_m
            )) / 100.0,
        )

    return sample


def _receptors(trial, response_t90_s: float | None):
    output = []
    seen = set()
    for reading in trial.readings:
        if reading.sensor in seen or reading.over_range:
            continue
        seen.add(reading.sensor)
        bearing = np.radians(reading.bearing)
        output.append(FixedReceptor(
            name=reading.sensor,
            east_m=float(reading.radius * np.sin(bearing)),
            north_m=float(reading.radius * np.cos(bearing)),
            height_m=float(reading.height),
            response_t90_s=response_t90_s,
        ))
    return output


def run(
    wind_csv: Path,
    *,
    table_winds_m_s: list[float],
    response_t90_s: float | None,
    start_s: float,
    end_s: float,
    reference_root: Path = REFERENCE,
) -> dict:
    history = _read_wind(wind_csv)
    time, measured_speed, _direction = history.arrays()
    table_winds = sorted(set(float(value) for value in table_winds_m_s))
    if not table_winds:
        raise ValueError("at least one table wind is required")
    if table_winds[0] > float(np.min(measured_speed)) or (
        table_winds[-1] < float(np.max(measured_speed))
    ):
        raise ValueError("steady table must bracket every measured wind speed")
    trial = next(item for item in load(reference_root) if item.test == 6)
    maximum_distance = max(reading.radius for reading in trial.readings)
    table = SteadyPlumeTable(
        table_winds,
        [
            _steady_sampler(trial, wind, maximum_distance)
            for wind in table_winds
        ],
    )
    traces = replay_fixed_receptors(
        history, _receptors(trial, response_t90_s), table
    )
    rows = []
    for trace in traces:
        stats = trace.statistics(start_s, end_s)
        rows.append({
            "sensor": trace.receptor.name,
            "east_m": trace.receptor.east_m,
            "north_m": trace.receptor.north_m,
            "height_m": trace.receptor.height_m,
            "true_mean_vol_pct": 100.0 * stats.true_mean,
            "true_max_vol_pct": 100.0 * stats.true_maximum,
            "indicated_mean_vol_pct": 100.0 * stats.indicated_mean,
            "indicated_max_vol_pct": 100.0 * stats.indicated_maximum,
        })
    return {
        "case": "FFI/DNV outdoor Test 6 transient receptor replay",
        "status": "research_only_observation_operator",
        "wind_input": {
            "path": str(wind_csv.resolve()),
            "sha256": _sha256(wind_csv),
            "samples": int(time.size),
            "speed_range_m_s": [
                float(np.min(measured_speed)), float(np.max(measured_speed))
            ],
        },
        "steady_table_winds_m_s": table_winds,
        "response_t90_s": response_t90_s,
        "statistics_window_s": [start_s, end_s],
        "limitations": [
            "steady spatial solutions are replayed; plume storage is absent",
            "wind direction rotates receptors but does not deform the plume",
            "no source time history or travel-time delay is inferred",
            "no parameter is fitted to Test 6 concentration observations",
        ],
        "receptors": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("wind_csv", type=Path)
    parser.add_argument(
        "--table-winds",
        default="0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6",
        help="comma-separated steady-solution wind grid in m/s",
    )
    parser.add_argument(
        "--response-t90",
        type=float,
        default=6.0,
        help="declared sensor 90%% response time; use 0 to bypass",
    )
    parser.add_argument("--start", type=float, default=20.0)
    parser.add_argument("--end", type=float, default=140.0)
    parser.add_argument("--reference-root", type=Path, default=REFERENCE)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    winds = [float(value) for value in args.table_winds.split(",")]
    payload = run(
        args.wind_csv,
        table_winds_m_s=winds,
        response_t90_s=(None if args.response_t90 == 0.0 else args.response_t90),
        start_s=args.start,
        end_s=args.end,
        reference_root=args.reference_root,
    )
    rendered = json.dumps(payload, indent=2)
    if args.output is None:
        print(rendered)
    else:
        args.output.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
