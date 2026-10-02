"""Read one public E3.1 D4/80-K file and audit declared blowdown end members."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.addons.cryogenic_blowdown import (
    CryogenicBlowdownConfig,
    run_cryogenic_blowdown,
)
from degali.validation.preslhy_e31 import (
    nozzle_pressure_rise,
    read_e31_pressure_run,
    valve_open_interval,
)


STATIONS_S = (0.5, 1.0, 2.0)
TEMPERATURES_K = (77.0, 80.0, 84.0)
DISCHARGE_COEFFICIENTS = (0.6, 0.7, 0.8)
WITHDRAWALS = ("vapour", "homogeneous")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def measured_pressure_at(run, absolute_time_s: float) -> float:
    pressure = np.asarray(run.pressure.vessel_pressure_bar, dtype=float)
    time = np.asarray(run.pressure.time_s, dtype=float)
    valid = np.isfinite(time) & np.isfinite(pressure) & (pressure > 0.0)
    if not valid.any() or absolute_time_s < time[valid][0] or absolute_time_s > time[valid][-1]:
        raise ValueError("declared source-time station lies outside vessel-pressure record")
    return float(np.interp(absolute_time_s, time[valid], pressure[valid]))


def model_row(*, pressure_bar: float, temperature_k: float, cd: float, withdrawal: str) -> dict:
    config = CryogenicBlowdownConfig(
        vessel_volume_m3=2.815e-3,
        nozzle_diameter_m=4.0e-3,
        initial_temperature_k=temperature_k,
        initial_pressure_pa=pressure_bar * 1.0e5,
        discharge_coefficient=cd,
        duration_s=max(STATIONS_S),
        time_step_s=0.01,
        two_phase_withdrawal=withdrawal,
        tank_internal_diameter_m=0.160 if withdrawal == "vapour" else None,
        outlet_height_from_bottom_m=0.030 if withdrawal == "vapour" else None,
    )
    result = run_cryogenic_blowdown(config)
    state_time = np.asarray([state.time_s for state in result.states], dtype=float)
    state_pressure = np.asarray([state.pressure_pa for state in result.states], dtype=float)
    state_temperature = np.asarray([state.temperature_k for state in result.states], dtype=float)
    available = tuple(float(station) <= state_time[-1] + 1e-12 for station in STATIONS_S)
    pressures = tuple(
        None if not ok else float(np.interp(station, state_time, state_pressure) / 1.0e5)
        for station, ok in zip(STATIONS_S, available)
    )
    temperatures = tuple(
        None if not ok else float(np.interp(station, state_time, state_temperature))
        for station, ok in zip(STATIONS_S, available)
    )
    return {
        "temperature_k": temperature_k,
        "discharge_coefficient": cd,
        "two_phase_withdrawal": withdrawal,
        "termination": result.termination,
        "maximum_discrete_energy_residual_j": result.maximum_discrete_energy_residual_j,
        "model_pressure_bar": pressures,
        "model_tank_temperature_k": temperatures,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.workbook = args.workbook.resolve()
    args.output = args.output.resolve()
    if args.output.exists():
        parser.error("refusing to overwrite source-validation evidence")
    if not args.workbook.is_file():
        raise FileNotFoundError(args.workbook)
    prereg = ROOT / "docs/prereg-e31-d4-80k-blowdown-pressure.md"
    dependencies = [args.workbook, prereg, Path(__file__), ROOT / "src/degali/validation/preslhy_e31.py", ROOT / "src/degali/addons/cryogenic_blowdown.py"]
    hashes_before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    run = read_e31_pressure_run(args.workbook)
    interval = valve_open_interval(run)
    response = nozzle_pressure_rise(run, interval, required_rise_bar=2.0)
    source_start = response.response_time_s
    initial_pressure = measured_pressure_at(run, source_start)
    observed = tuple(measured_pressure_at(run, source_start + station) for station in STATIONS_S)
    rows = [
        model_row(
            pressure_bar=initial_pressure,
            temperature_k=temperature,
            cd=cd,
            withdrawal=withdrawal,
        )
        for withdrawal in WITHDRAWALS
        for temperature in TEMPERATURES_K
        for cd in DISCHARGE_COEFFICIENTS
    ]
    coverage = []
    for index, station in enumerate(STATIONS_S):
        predictions = [row["model_pressure_bar"][index] for row in rows]
        finite = [value for value in predictions if value is not None]
        coverage.append({
            "time_after_source_start_s": station,
            "measured_pressure_bar": observed[index],
            "minimum_declared_model_pressure_bar": None if not finite else min(finite),
            "maximum_declared_model_pressure_bar": None if not finite else max(finite),
            "within_declared_envelope": bool(finite and min(finite) <= observed[index] <= max(finite)),
            "available_end_members": len(finite),
        })
    hashes_after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if hashes_before != hashes_after:
        raise RuntimeError("audit input changed during calculation")
    result = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": "PRESLHY E3.1 part A, D=4 mm, nominal 80 K",
        "source_doi": "10.35097/1187",
        "source_files_distributed": False,
        "model_coefficient_fitted": False,
        "default_model_changed": False,
        "hashes": hashes_after,
        "workbook": args.workbook.name,
        "initial_vessel_pressure_metadata_bar": run.initial_vessel_pressure_bar,
        "nozzle_diameter_mm": run.nozzle_diameter_mm,
        "relay_interval": asdict(interval),
        "nozzle_pressure_response": asdict(response),
        "source_start_time_s": source_start,
        "model_initial_vessel_pressure_bar": initial_pressure,
        "stations_after_source_start_s": STATIONS_S,
        "measured_vessel_pressure_bar": observed,
        "end_members": rows,
        "declared_envelope_coverage": coverage,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps({
        "workbook": result["workbook"],
        "initial_pressure_bar": initial_pressure,
        "response_delay_s": response.delay_after_relay_s,
        "coverage": coverage,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
