"""Run a bounded one-parameter TI021 physics sweep from an exact result."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import inspect
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.score_ti021_dense_layer_hybrid import run
from tools.sweep_ti021_vertical_trajectory import _base_parameters, _execute, _groups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--parameter", required=True)
    parser.add_argument("--value", type=float, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()

    accepted = inspect.signature(run).parameters
    if args.parameter not in accepted or args.parameter == "output":
        raise ValueError(f"not a score_ti021_dense_layer_hybrid.run parameter: {args.parameter}")
    base = json.loads(args.base.read_text(encoding="utf-8"))
    parameters = _base_parameters(base)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    jobs: list[tuple[str, dict]] = []
    for value in args.value:
        token = f"{value:g}".replace("-", "m").replace(".", "p")
        output = args.output_dir / f"ti021-{args.parameter}-{token}.json"
        candidate = dict(parameters)
        candidate[args.parameter] = value
        jobs.append((str(output), candidate))

    completed: list[Path] = []
    failures: list[dict] = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {
            executor.submit(_execute, job): {
                "path": job[0],
                "parameter_value": job[1][args.parameter],
            }
            for job in jobs
        }
        for future in as_completed(futures):
            metadata = futures[future]
            try:
                completed.append(Path(future.result()))
            except Exception as error:
                failures.append({
                    **metadata,
                    "error_type": type(error).__name__,
                    "error": str(error),
                })

    observed = np.asarray(
        [row["observed_dose_vol_percent_s"] for row in base["sensor_rows"]], dtype=float
    )
    base_prediction = np.asarray(
        [row["predicted_dose_vol_percent_s"] for row in base["sensor_rows"]], dtype=float
    )
    base_groups = _groups(base_prediction, observed)
    rows = []
    for path in completed:
        record = json.loads(path.read_text(encoding="utf-8"))
        prediction = np.asarray(
            [row["predicted_dose_vol_percent_s"] for row in record["sensor_rows"]],
            dtype=float,
        )
        groups = _groups(prediction, observed)
        deltas = {name: groups[name] - base_groups[name] for name in base_groups}
        rows.append(
            {
                "path": str(path.resolve()),
                "parameter_value": record[args.parameter],
                "mape_percent": record["dose_mape_percent"],
                "runtime_including_output_s": record["runtime_including_output_s"],
                "group_mape_percent": groups,
                "group_delta_percentage_points": deltas,
                "all_groups_nondegraded": all(value <= 1.0e-12 for value in deltas.values()),
            }
        )
    rows.sort(key=lambda row: row["mape_percent"])
    if not rows:
        raise RuntimeError(f"all sweep points failed: {failures}")
    report = {
        "schema": "degali.ti021-single-parameter-sweep.v1",
        "status": "diagnostic_not_promoted",
        "base": str(args.base.resolve()),
        "parameter": args.parameter,
        "rows": rows,
        "failures": sorted(failures, key=lambda row: row["parameter_value"]),
        "best": rows[0],
        "target_10_percent_reached": rows[0]["mape_percent"] <= 10.0,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
