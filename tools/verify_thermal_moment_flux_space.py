"""Independently reconstruct the sealed flux-space endpoint balances."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_reservoir_thermal_segments import actual_source_model, serial
from degali.addons.reservoir_thermal import ReservoirShortSegment


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path,
        default=REF / "thermal_moment_flux_space_2026-09-09.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=REF / "thermal_moment_flux_space_verification_2026-09-09.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite flux-space verification evidence")
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    dependencies = [
        args.input.resolve(),
        reduced_path,
        measured_path,
        ROOT / "src/degali/addons/reservoir_thermal.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    evidence = json.loads(args.input.read_text(encoding="utf-8"))
    if not evidence["completed"]:
        raise RuntimeError("flux-space audit is incomplete")
    trials = json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
    measured_rows = {
        row["trial"]: row
        for row in json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    }
    trial = next(row for row in trials if row["trial"] == evidence["trial"])
    jp, thermodynamics, _source = actual_source_model(trial, measured_rows)
    driver = ReservoirShortSegment(
        jp, thermodynamics, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work",
    )

    initial_parameters = np.asarray(evidence["coarse"]["parameters"][0], float)
    initial_evaluation = driver.evaluate(initial_parameters, order=16, probes=1025)
    initial = driver.flux_values(initial_evaluation["mixing"], order=16)
    scale = np.maximum(abs(initial), [1e-12, 1e-12, 1., 1., 1., 1.])
    rows = {}
    for label in ("coarse", "refined"):
        record = evidence[label]
        final_parameters = np.asarray(record["parameters"][-1], float)
        final_evaluation = driver.evaluate(final_parameters, order=16, probes=1025)
        final = driver.flux_values(final_evaluation["mixing"], order=16)
        cumulative = np.asarray(record["cumulative_sources"][-1], float)
        balance = (final-initial-cumulative)/scale
        stored = np.asarray(record["independent_order16_scaled_balance"], float)
        rows[label] = dict(
            initial_fluxes=initial,
            final_fluxes=final,
            cumulative_sources=cumulative,
            scaled_balance=balance,
            maximum_scaled_balance=float(max(abs(balance))),
            stored_balance_difference=float(max(abs(balance-stored))),
            passed=bool(max(abs(balance)) <= 1e-5 and max(abs(balance-stored)) <= 1e-12),
        )
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed during independent verification")
    payload = dict(
        completed=True,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        source_audit=str(args.input.relative_to(ROOT)),
        quadrature_order=16,
        probes=1025,
        balance_limit=1e-5,
        rows=rows,
        all_passed=all(row["passed"] for row in rows.values()),
        hashes=after,
    )
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    print(json.dumps({
        "all_passed": payload["all_passed"],
        "coarse_balance": rows["coarse"]["maximum_scaled_balance"],
        "refined_balance": rows["refined"]["maximum_scaled_balance"],
        "saved": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
