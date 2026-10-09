"""Replay a supplied FFI Test 6 wind history through steady LH2 plumes.

This is an observation operator, not a transient dispersion solver.  The
steady free-field ``hydrogen_jet`` calculation remains unchanged: one plume
is built for each declared wind-table speed, then fixed Test 6 receptors are
rotated into the instantaneous wind frame and sampled.  A declared causal
first-order response may be applied after the true signal is calculated.

The public report contains plotted time histories but not a redistributable
machine-readable export.  Consequently this tool requires the caller to
provide the CSV and records its SHA-256 rather than embedding the trace in
the output artifact.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Mapping, Sequence

import numpy as np

from degali.addons.field_json import strict_json_loads
from degali.addons.transient_receptor import (
    FixedReceptor,
    SteadyPlumeTable,
    WindHistory,
    replay_fixed_receptors,
)
from degali.validation.nearfield import STEP, Trajectory, hydrogen_jet
from degali.validation.ffi_source_state import ffi_reference_provenance
from degali.validation.spadeadam import (
    AMBIENT_TEMPERATURE,
    DEFAULT_ROOT,
    RELEASE_HEIGHT,
    RELATIVE_HUMIDITY,
    load,
)


SCHEMA = "degali.ffi-test6-transient-receptor-execution.v1"
HISTORY_COLUMNS = (
    "time_s", "wind_speed_ms", "wind_direction_from_deg",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite(value: str, name: str, row_number: int) -> float:
    try:
        parsed = float(value.strip())
    except (AttributeError, ValueError) as error:
        raise ValueError(
            f"wind-history row {row_number} has a non-numeric {name}"
        ) from error
    if not math.isfinite(parsed):
        raise ValueError(f"wind-history row {row_number} has non-finite {name}")
    return parsed


def read_wind_history_csv(path: str | Path) -> tuple[WindHistory, int]:
    """Read the exact three-column wind-history contract.

    Extra, missing, duplicate, blank or malformed columns are rejected before
    the model is called.  The returned row count excludes the header.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    with source.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("wind-history CSV is empty") from error
        normalized = tuple(item.strip() for item in header)
        if normalized != HISTORY_COLUMNS or len(set(normalized)) != len(normalized):
            raise ValueError(
                "wind-history CSV header must be exactly "
                + ",".join(HISTORY_COLUMNS)
            )
        times: list[float] = []
        speeds: list[float] = []
        directions: list[float] = []
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(HISTORY_COLUMNS):
                raise ValueError(
                    f"wind-history row {row_number} has {len(row)} fields; "
                    f"expected {len(HISTORY_COLUMNS)}"
                )
            if any(not item.strip() for item in row):
                raise ValueError(f"wind-history row {row_number} contains a blank field")
            times.append(_finite(row[0], HISTORY_COLUMNS[0], row_number))
            speeds.append(_finite(row[1], HISTORY_COLUMNS[1], row_number))
            directions.append(_finite(row[2], HISTORY_COLUMNS[2], row_number))
    if not times:
        raise ValueError("wind-history CSV contains no data rows")
    history = WindHistory(times, speeds, directions)
    history.arrays()
    return history, len(times)


def parse_table_winds(text: str) -> tuple[float, ...]:
    """Parse a strictly increasing positive comma-separated wind table."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("--table-winds must contain positive comma-separated speeds")
    values: list[float] = []
    for item in text.split(","):
        try:
            value = float(item.strip())
        except ValueError as error:
            raise ValueError("--table-winds must contain only finite numbers") from error
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("--table-winds speeds must be finite and positive")
        values.append(value)
    if not values or any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("--table-winds speeds must be strictly increasing")
    return tuple(values)


def history_diagnostics(history: WindHistory) -> dict[str, object]:
    """Summarize the supplied time axis without assigning a quality grade."""
    time, speed, direction = history.arrays()
    if time.size > 1:
        intervals = time[1:] - time[:-1]
        minimum_dt = float(np.min(intervals))
        maximum_dt = float(np.max(intervals))
        uniform = bool(np.allclose(intervals, intervals[0], rtol=1.0e-6, atol=1.0e-12))
    else:
        minimum_dt = None
        maximum_dt = None
        uniform = True
    return {
        "sample_count": int(time.size),
        "time_start_s": float(time[0]),
        "time_end_s": float(time[-1]),
        "minimum_dt_s": minimum_dt,
        "maximum_dt_s": maximum_dt,
        "uniform_sampling": uniform,
        "wind_speed_min_m_s": float(np.min(speed)),
        "wind_speed_max_m_s": float(np.max(speed)),
        "direction_from_deg_first": float(direction[0]),
        "direction_from_deg_last": float(direction[-1]),
    }


def _test6_receptors(trial) -> tuple[FixedReceptor, ...]:
    receptors = []
    for reading in trial.readings:
        bearing = math.radians(reading.bearing)
        receptors.append(FixedReceptor(
            name=reading.sensor,
            east_m=float(reading.radius * math.sin(bearing)),
            north_m=float(reading.radius * math.cos(bearing)),
            height_m=float(reading.height),
            response_t90_s=None,
        ))
    if not receptors:
        raise ValueError("Test 6 has no reference receptors")
    return tuple(receptors)


def _steady_table(
    trial,
    wind_speeds: Sequence[float],
    *,
    corrections: bool,
    maximum_radius_m: float,
) -> SteadyPlumeTable:
    """Build one steady trajectory per declared wind-table speed."""
    trajectories: list[Trajectory] = []
    integration_limit = max(2.0 * maximum_radius_m, maximum_radius_m + 20.0)
    for wind in wind_speeds:
        setup, initial = hydrogen_jet(
            rate=trial.rate,
            diameter=trial.orifice,
            wind=wind,
            height=RELEASE_HEIGHT,
            ambient_temperature=AMBIENT_TEMPERATURE,
            relative_humidity=RELATIVE_HUMIDITY,
            storage_pressure_barg=trial.line_pressure,
            wind_reference_height=1.5,
            corrections=corrections,
        )
        trajectory = Trajectory(
            setup.th.table,
            setup.run(y0=initial, distmx=STEP, smax=integration_limit).rows,
        )
        if not trajectory.ok or trajectory.rows[-1, 0] < maximum_radius_m:
            raise ValueError(
                f"steady plume at wind {wind:g} m/s does not reach the "
                f"farthest Test 6 receptor ({maximum_radius_m:g} m)"
            )
        trajectories.append(trajectory)

    return SteadyPlumeTable(
        wind_speeds,
        tuple(
            lambda alongwind, crosswind, height, trajectory=trajectory: (
                trajectory.concentration_at(alongwind, crosswind, height) / 100.0
            )
            for trajectory in trajectories
        ),
    )


def _record_replay(
    trial,
    traces,
    *,
    start_s: float,
    end_s: float,
    response_t90_s: float | None,
) -> list[dict[str, object]]:
    readings = {reading.sensor: reading for reading in trial.readings}
    rows: list[dict[str, object]] = []
    for trace in traces:
        reading = readings[trace.receptor.name]
        stats = trace.statistics(start_s, end_s)
        rows.append({
            "sensor": trace.receptor.name,
            "radius_m": float(reading.radius),
            "bearing_deg": float(reading.bearing),
            "height_m": float(reading.height),
            "observed_mean_vol_pct": None if reading.over_range else float(reading.mean),
            "observed_peak_vol_pct": None if reading.over_range else float(reading.peak),
            "observed_over_range": bool(reading.over_range),
            "observed_is_lower_bound": True,
            "response_t90_s": response_t90_s,
            "window": {
                "start_s": stats.start_s,
                "end_s": stats.end_s,
                "sample_count": stats.sample_count,
                "true_mean_vol_pct": 100.0 * stats.true_mean,
                "true_maximum_vol_pct": 100.0 * stats.true_maximum,
                "indicated_mean_vol_pct": 100.0 * stats.indicated_mean,
                "indicated_maximum_vol_pct": 100.0 * stats.indicated_maximum,
            },
        })
    return rows


def run_replay(
    history: WindHistory,
    *,
    trial_number: int = 6,
    reference_root: str | Path | None = None,
    table_winds: Sequence[float],
    start_s: float,
    end_s: float,
    response_t90_s: float | None,
    corrections: bool = True,
) -> tuple[object, list[dict[str, object]]]:
    """Run the replay and return the source trial plus compact sensor rows."""
    if trial_number != 6:
        raise ValueError("this replay boundary is defined for FFI Test 6 only")
    if not math.isfinite(start_s) or not math.isfinite(end_s) or end_s < start_s:
        raise ValueError("replay window must be finite and end after its start")
    if response_t90_s is not None and (
        not math.isfinite(response_t90_s) or response_t90_s <= 0.0
    ):
        raise ValueError("response_t90_s must be positive or None")
    wind_speeds = tuple(float(value) for value in table_winds)
    if not wind_speeds or any(
        not math.isfinite(value) or value <= 0.0
        for value in wind_speeds
    ) or any(b <= a for a, b in zip(wind_speeds, wind_speeds[1:])):
        raise ValueError("table_winds must be strictly increasing positive values")
    history.arrays()
    trial = next(
        (item for item in load(reference_root or DEFAULT_ROOT) if item.test == trial_number),
        None,
    )
    if trial is None:
        raise ValueError("reference data contains no FFI Test 6")
    receptors = _test6_receptors(trial)
    maximum_radius = max(float(reading.radius) for reading in trial.readings)
    table = _steady_table(
        trial,
        wind_speeds,
        corrections=corrections,
        maximum_radius_m=maximum_radius,
    )
    configured = tuple(
        FixedReceptor(
            receptor.name,
            receptor.east_m,
            receptor.north_m,
            receptor.height_m,
            response_t90_s=response_t90_s,
        )
        for receptor in receptors
    )
    traces = replay_fixed_receptors(
        history,
        configured,
        table,
    )
    return trial, _record_replay(
        trial,
        traces,
        start_s=start_s,
        end_s=end_s,
        response_t90_s=response_t90_s,
    )


def build_payload(
    history_path: str | Path,
    *,
    table_winds: Sequence[float],
    start_s: float,
    end_s: float,
    response_t90_s: float | None,
    corrections: bool,
    reference_root: str | Path | None = None,
) -> dict[str, object]:
    source = Path(history_path).resolve()
    history, row_count = read_wind_history_csv(source)
    diagnostics = history_diagnostics(history)
    trial, receptors = run_replay(
        history,
        reference_root=reference_root,
        table_winds=table_winds,
        start_s=start_s,
        end_s=end_s,
        response_t90_s=response_t90_s,
        corrections=corrections,
    )
    indicated_peaks = [row["window"]["indicated_maximum_vol_pct"] for row in receptors]
    true_peaks = [row["window"]["true_maximum_vol_pct"] for row in receptors]
    return {
        "schema": SCHEMA,
        "input": {
            "history_path": str(source),
            "history_sha256": _sha256(source),
            "history_row_count": row_count,
            "history_diagnostics": diagnostics,
            "reference_root": str(Path(reference_root or DEFAULT_ROOT).resolve()),
            "reference_provenance": ffi_reference_provenance(reference_root or DEFAULT_ROOT),
            "test": trial.test,
            "table_winds_m_s": [float(value) for value in table_winds],
            "corrections": bool(corrections),
            "response_t90_s": response_t90_s,
            "window": {"start_s": float(start_s), "end_s": float(end_s)},
        },
        "observation_operator": {
            "kind": "instantaneous_steady_wind_frame_with_first_order_sensor_response",
            "source_history_is_time_varying": False,
            "plume_storage_included": False,
            "observed_values_are_lower_bounds": True,
            "sensor_count": len(receptors),
        },
        "replay": {
            "status": "complete",
            "sensor_rows": receptors,
            "summary": {
                "sensor_count": len(receptors),
                "indicated_peak_max_vol_pct": max(indicated_peaks),
                "true_peak_max_vol_pct": max(true_peaks),
                "over_range_sensor_count": sum(
                    1 for row in receptors if row["observed_over_range"]
                ),
            },
            "promotion_allowed": False,
            "validation_qualified": False,
        },
    }


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _require_keys(
    value: Mapping[str, object], name: str, required: set[str],
) -> None:
    missing = sorted(required - set(value))
    unknown = sorted(set(value) - required)
    if missing or unknown:
        details = []
        if missing:
            details.append("missing=" + ", ".join(missing))
        if unknown:
            details.append("unknown=" + ", ".join(unknown))
        raise ValueError(f"{name} keys are invalid: " + "; ".join(details))


def _sha_matches(path: Path, expected: object, name: str) -> None:
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError(f"{name} must be a SHA-256 digest")
    actual = _sha256(path)
    if actual != expected.lower() or expected != expected.lower():
        raise ValueError(f"{name} does not match current file bytes: {path}")


def verify_replay_artifact(path: str | Path) -> dict[str, object]:
    """Verify a saved replay's provenance and compact report invariants.

    This is intentionally integrity-only.  It does not claim that the replay
    is validation, and it does not silently rerun a potentially expensive
    CoolProp calculation.
    """
    artifact = Path(path).resolve()
    if not artifact.is_file():
        raise FileNotFoundError(artifact)
    payload = strict_json_loads(artifact.read_text(encoding="utf-8-sig"))
    root = _mapping(payload, "transient replay")
    _require_keys(
        root,
        "transient replay",
        {"schema", "input", "observation_operator", "replay"},
    )
    if root["schema"] != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA!r}")
    input_data = _mapping(root["input"], "transient replay input")
    _require_keys(
        input_data,
        "transient replay input",
        {
            "history_path", "history_sha256", "history_row_count", "history_diagnostics", "reference_root",
            "reference_provenance", "test", "table_winds_m_s", "corrections",
            "response_t90_s", "window",
        },
    )
    history_path = Path(str(input_data["history_path"])).resolve()
    if not history_path.is_file():
        raise FileNotFoundError(history_path)
    _sha_matches(history_path, input_data["history_sha256"], "input.history_sha256")
    history, row_count = read_wind_history_csv(history_path)
    if input_data["history_row_count"] != row_count:
        raise ValueError("input.history_row_count does not match the history file")
    stored_diagnostics = _mapping(
        input_data["history_diagnostics"], "input.history_diagnostics",
    )
    expected_diagnostics = history_diagnostics(history)
    _require_keys(
        stored_diagnostics,
        "input.history_diagnostics",
        set(expected_diagnostics),
    )
    for key, expected in expected_diagnostics.items():
        actual = stored_diagnostics[key]
        if isinstance(expected, bool):
            if actual is not expected:
                raise ValueError(f"input.history_diagnostics.{key} changed")
        elif expected is None:
            if actual is not None:
                raise ValueError(f"input.history_diagnostics.{key} changed")
        elif not isinstance(actual, (int, float)) or not math.isclose(
            float(actual), float(expected), rel_tol=1.0e-12, abs_tol=1.0e-12,
        ):
            raise ValueError(f"input.history_diagnostics.{key} changed")
    reference_root = Path(str(input_data["reference_root"])).resolve()
    expected_provenance = _mapping(
        input_data["reference_provenance"], "input.reference_provenance",
    )
    actual_provenance = ffi_reference_provenance(reference_root)
    if expected_provenance != actual_provenance:
        raise ValueError("reference provenance changed since replay generation")
    if input_data["test"] != 6:
        raise ValueError("replay input.test must be 6")
    table_winds = input_data["table_winds_m_s"]
    if not isinstance(table_winds, list):
        raise ValueError("input.table_winds_m_s must be an array")
    parse_table_winds(",".join(str(value) for value in table_winds))
    if not isinstance(input_data["corrections"], bool):
        raise ValueError("input.corrections must be boolean")
    response = input_data["response_t90_s"]
    if response is not None and (
        isinstance(response, bool) or not isinstance(response, (int, float))
        or not math.isfinite(float(response)) or float(response) <= 0.0
    ):
        raise ValueError("input.response_t90_s must be positive or null")
    window = _mapping(input_data["window"], "input.window")
    _require_keys(window, "input.window", {"start_s", "end_s"})
    start_s = float(window["start_s"])
    end_s = float(window["end_s"])
    if not math.isfinite(start_s) or not math.isfinite(end_s) or end_s < start_s:
        raise ValueError("input.window is invalid")
    history_time, _history_speed, _history_direction = history.arrays()
    expected_sample_count = int(
        ((history_time >= start_s) & (history_time <= end_s)).sum()
    )
    if expected_sample_count <= 0:
        raise ValueError("input.window contains no history samples")
    operator = _mapping(root["observation_operator"], "observation operator")
    _require_keys(
        operator,
        "observation operator",
        {
            "kind", "source_history_is_time_varying", "plume_storage_included",
            "observed_values_are_lower_bounds", "sensor_count",
        },
    )
    if operator["kind"] != "instantaneous_steady_wind_frame_with_first_order_sensor_response":
        raise ValueError("observation operator kind is unsupported")
    if operator["source_history_is_time_varying"] is not False:
        raise ValueError("replay must not claim a time-varying source history")
    if operator["plume_storage_included"] is not False:
        raise ValueError("replay must not claim plume storage")
    if operator["observed_values_are_lower_bounds"] is not True:
        raise ValueError("replay must retain lower-bound observation provenance")
    replay = _mapping(root["replay"], "transient replay report")
    _require_keys(
        replay,
        "transient replay report",
        {"status", "sensor_rows", "summary", "promotion_allowed", "validation_qualified"},
    )
    if replay["status"] != "complete":
        raise ValueError("replay.status must be complete")
    if replay["promotion_allowed"] is not False or replay["validation_qualified"] is not False:
        raise ValueError("replay promotion and qualification must remain false")
    rows = replay["sensor_rows"]
    if not isinstance(rows, list) or not rows:
        raise ValueError("replay.sensor_rows must be a non-empty array")
    reference_trial = next(
        (item for item in load(reference_root) if item.test == 6), None,
    )
    if reference_trial is None:
        raise ValueError("reference data contains no FFI Test 6")
    readings = {item.sensor: item for item in reference_trial.readings}
    seen: set[str] = set()
    for index, item in enumerate(rows):
        record = _mapping(item, f"replay.sensor_rows[{index}]")
        required = {
            "sensor", "radius_m", "bearing_deg", "height_m",
            "observed_mean_vol_pct", "observed_peak_vol_pct",
            "observed_over_range", "observed_is_lower_bound",
            "response_t90_s", "window",
        }
        _require_keys(record, f"replay.sensor_rows[{index}]", required)
        sensor = record["sensor"]
        if not isinstance(sensor, str) or sensor in seen or sensor not in readings:
            raise ValueError(f"replay sensor identity is invalid: {sensor!r}")
        seen.add(sensor)
        reference = readings[sensor]
        for key, expected in (
            ("radius_m", reference.radius),
            ("bearing_deg", reference.bearing),
            ("height_m", reference.height),
        ):
            if float(record[key]) != float(expected):
                raise ValueError(f"replay {sensor} {key} does not match reference")
        if not isinstance(record["observed_over_range"], bool):
            raise ValueError(f"replay {sensor} over-range flag is invalid")
        if record["observed_is_lower_bound"] is not True:
            raise ValueError(f"replay {sensor} must retain lower-bound provenance")
        if record["response_t90_s"] != response:
            raise ValueError(f"replay {sensor} response time does not match input")
        if record["observed_over_range"]:
            if record["observed_mean_vol_pct"] is not None or record["observed_peak_vol_pct"] is not None:
                raise ValueError(f"replay {sensor} over-range observation must be null")
        else:
            for key in ("observed_mean_vol_pct", "observed_peak_vol_pct"):
                value = record[key]
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)) or float(value) < 0.0:
                    raise ValueError(f"replay {sensor} {key} is invalid")
        sensor_window = _mapping(record["window"], f"replay.sensor_rows[{index}].window")
        _require_keys(
            sensor_window,
            f"replay.sensor_rows[{index}].window",
            {
                "start_s", "end_s", "sample_count", "true_mean_vol_pct",
                "true_maximum_vol_pct", "indicated_mean_vol_pct",
                "indicated_maximum_vol_pct",
            },
        )
        if sensor_window["start_s"] != start_s or sensor_window["end_s"] != end_s:
            raise ValueError(f"replay {sensor} window does not match input")
        if (
            not isinstance(sensor_window["sample_count"], int)
            or sensor_window["sample_count"] != expected_sample_count
        ):
            raise ValueError(f"replay {sensor} sample count is invalid")
        for key in (
            "true_mean_vol_pct", "true_maximum_vol_pct",
            "indicated_mean_vol_pct", "indicated_maximum_vol_pct",
        ):
            value = sensor_window[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"replay {sensor} {key} is invalid")
    if seen != set(readings):
        raise ValueError("replay sensor rows do not cover the Test 6 reference array")
    summary = _mapping(replay["summary"], "replay.summary")
    _require_keys(
        summary,
        "replay.summary",
        {"sensor_count", "indicated_peak_max_vol_pct", "true_peak_max_vol_pct", "over_range_sensor_count"},
    )
    if operator["sensor_count"] != len(rows) or summary["sensor_count"] != len(rows):
        raise ValueError("replay summary sensor count does not match rows")
    if summary["over_range_sensor_count"] != sum(
        1 for item in rows if item["observed_over_range"]
    ):
        raise ValueError("replay summary over-range count does not match rows")
    true_peaks = [item["window"]["true_maximum_vol_pct"] for item in rows]
    indicated_peaks = [item["window"]["indicated_maximum_vol_pct"] for item in rows]
    if summary["true_peak_max_vol_pct"] != max(true_peaks):
        raise ValueError("replay summary true peak does not match rows")
    if summary["indicated_peak_max_vol_pct"] != max(indicated_peaks):
        raise ValueError("replay summary indicated peak does not match rows")
    return {
        "schema": "degali.ffi-test6-transient-receptor-verification.v1",
        "artifact_path": str(artifact),
        "artifact_sha256": _sha256(artifact),
        "history_path": str(history_path),
        "history_sha256_verified": True,
        "reference_provenance_verified": True,
        "report_recomputed": False,
        "promotion_allowed": False,
    }


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "history", type=Path, nargs="?",
        help="explicit three-column wind-history CSV",
    )
    parser.add_argument(
        "--verify", type=Path,
        help="verify a saved degali.ffi-test6-transient-receptor-execution.v1 artifact",
    )
    parser.add_argument(
        "--table-winds",
        help="strictly increasing comma-separated steady plume wind speeds (m/s)",
    )
    parser.add_argument("--response-t90", type=float, default=6.0)
    parser.add_argument("--start", type=float)
    parser.add_argument("--end", type=float)
    parser.add_argument("--reference-root", type=Path)
    parser.add_argument("--no-corrections", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.verify is not None:
            if args.history is not None or args.table_winds is not None or args.start is not None or args.end is not None:
                raise ValueError("--verify cannot be combined with replay inputs")
            print(json.dumps(verify_replay_artifact(args.verify), indent=2, ensure_ascii=False) + "\n", end="")
            return 0
        if args.history is None or args.table_winds is None or args.start is None or args.end is None:
            raise ValueError("history, --table-winds, --start and --end are required unless --verify is used")
        winds = parse_table_winds(args.table_winds)
        response = None if args.response_t90 == 0.0 else args.response_t90
        payload = build_payload(
            args.history,
            table_winds=winds,
            start_s=args.start,
            end_s=args.end,
            response_t90_s=response,
            corrections=not args.no_corrections,
            reference_root=args.reference_root,
        )
        rendered = json.dumps(
            payload, indent=2, ensure_ascii=False, allow_nan=False,
        ) + "\n"
        if args.output is None:
            print(rendered, end="")
        else:
            output = args.output
            if output.exists():
                raise FileExistsError(
                    f"refusing to overwrite transient replay: {output}"
                )
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("x", encoding="utf-8", newline="\n") as stream:
                stream.write(rendered)
            print(f"FFI Test 6 transient replay: {output}")
            print(f"sensor count: {payload['replay']['summary']['sensor_count']}")
        return 0
    except (OSError, TypeError, ValueError, RuntimeError) as error:
        print(f"audit_ffi_test6_transient_receptors: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
