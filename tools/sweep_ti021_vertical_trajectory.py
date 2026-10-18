"""Screen a small physics-based TI021 overtop vertical-trajectory grid.

Every point is a complete 600 s replay.  The two parameters describe a common
decaying upward velocity after wall impact; no receptor values or per-sensor
coefficients enter the solver.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout
import inspect
import io
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook
from tools.score_ti021_dense_layer_hybrid import run


TRAJECTORIES = (
    (0.0, 1.0),
    (0.0125, 20.0),
    (0.025, 20.0),
    (0.05, 20.0),
    (0.075, 20.0),
    (0.05, 40.0),
    (0.10, 20.0),
    (0.075, 40.0),
)


def _base_parameters(record: dict) -> dict:
    accepted = inspect.signature(run).parameters
    params = {key: value for key, value in record.items() if key in accepted}
    params.update(
        h=record["horizontal_cell_size_m"],
        dt=record["max_dt_s"],
        source_dt=record["source_time_step_s"],
        entrainment=record["dense_entrainment_model"],
    )
    return params


def _execute(item: tuple[str, dict]) -> str:
    output_text, params = item
    with redirect_stdout(io.StringIO()):
        run(Path(output_text), **params)
    return output_text


def _mape(prediction: np.ndarray, observed: np.ndarray, mask: np.ndarray) -> float:
    selected = mask & (observed > 0.0)
    return float(np.mean(np.abs(prediction[selected] / observed[selected] - 1.0)) * 100.0)


def _groups(prediction: np.ndarray, observed: np.ndarray) -> dict[str, float]:
    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
    height = np.asarray([sensor.height_m for sensor in dataset.dose_sensors])
    downwind = np.asarray([sensor.local_y_m for sensor in dataset.dose_sensors])
    result = {
        f"height_{value:g}_m": _mape(prediction, observed, np.isclose(height, value))
        for value in sorted(set(height))
    }
    for name, lower, upper in (
        ("downwind_lt_150_m", -np.inf, 150.0),
        ("downwind_150_to_300_m", 150.0, 300.0),
        ("downwind_ge_300_m", 300.0, np.inf),
    ):
        result[name] = _mape(
            prediction, observed, (downwind >= lower) & (downwind < upper)
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument(
        "--branch", choices=("marker", "conserved"), default="marker"
    )
    args = parser.parse_args()

    base = json.loads(args.base.read_text(encoding="utf-8"))
    params = _base_parameters(base)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    jobs: list[tuple[str, dict]] = []
    for velocity, rise_time in TRAJECTORIES:
        token = f"v{velocity:g}-tau{rise_time:g}".replace(".", "p")
        output = args.output_dir / f"ti021-vertical-trajectory-{token}.json"
        candidate = dict(params)
        prefix = (
            "wall_impulse_elevated"
            if args.branch == "marker"
            else "wall_impulse_elevated_conserved"
        )
        candidate[f"{prefix}_rise_velocity_m_s"] = velocity
        candidate[f"{prefix}_rise_time_s"] = rise_time
        jobs.append((str(output), candidate))

    completed: list[Path] = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {executor.submit(_execute, job): job[0] for job in jobs}
        for future in as_completed(futures):
            completed.append(Path(future.result()))

    base_observed = np.asarray(
        [row["observed_dose_vol_percent_s"] for row in base["sensor_rows"]], dtype=float
    )
    base_prediction = np.asarray(
        [row["predicted_dose_vol_percent_s"] for row in base["sensor_rows"]], dtype=float
    )
    base_groups = _groups(base_prediction, base_observed)
    rows = []
    for path in completed:
        record = json.loads(path.read_text(encoding="utf-8"))
        prediction = np.asarray(
            [row["predicted_dose_vol_percent_s"] for row in record["sensor_rows"]],
            dtype=float,
        )
        groups = _groups(prediction, base_observed)
        group_delta = {name: groups[name] - base_groups[name] for name in base_groups}
        rows.append(
            {
                "path": str(path.resolve()),
                "rise_velocity_m_s": record[f"{prefix}_rise_velocity_m_s"],
                "rise_time_s": record[f"{prefix}_rise_time_s"],
                "asymptotic_rise_displacement_m": (
                    record[f"{prefix}_rise_velocity_m_s"]
                    * record[f"{prefix}_rise_time_s"]
                ),
                "mape_percent": record["dose_mape_percent"],
                "runtime_including_output_s": record["runtime_including_output_s"],
                "group_mape_percent": groups,
                "group_delta_percentage_points": group_delta,
                "all_groups_nondegraded": all(value <= 1.0e-12 for value in group_delta.values()),
            }
        )
    rows.sort(key=lambda row: row["mape_percent"])
    summary = {
        "schema": "degali.ti021-vertical-trajectory-sweep.v1",
        "status": "diagnostic_not_promoted",
        "base": str(args.base.resolve()),
        "branch": args.branch,
        "parameter_contract": (
            "common decaying post-impact upward velocity; no receptor-specific calibration"
        ),
        "rows": rows,
        "best": rows[0],
        "target_10_percent_reached": rows[0]["mape_percent"] <= 10.0,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
