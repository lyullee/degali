"""Recompute coarse/refined six-balance residuals from sealed reach states."""

from __future__ import annotations

from datetime import datetime, timezone
import argparse
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
        "--input",
        type=Path,
        default=REF / "thermal_moment_downstream_reach_2026-09-08.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REF / "thermal_moment_downstream_reach_verification_2026-09-08.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite reach verification evidence")
    reduced = REF / "e35_reduced.json"
    measured = REF / "measured_pipe_source_2026-09-05.json"
    dependencies = [
        args.input.resolve(),
        reduced,
        measured,
        ROOT / "src/degali/addons/reservoir_thermal.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    evidence = json.loads(args.input.read_text(encoding="utf-8"))
    trials = json.loads(reduced.read_text(encoding="utf-8"))["trials"]
    measured_rows = {
        row["trial"]: row
        for row in json.loads(measured.read_text(encoding="utf-8"))["trials"]
    }
    jp, thermodynamics, _source = actual_source_model(
        next(row for row in trials if row["trial"] == 10), measured_rows
    )
    driver = ReservoirShortSegment(
        jp,
        thermodynamics,
        thermal_species_ratio=1.0,
        mechanical_work="reduced_buoyancy_work",
    )
    initial = np.asarray(evidence["coarse"]["states"][0], dtype=float)
    initial_evaluation = driver.evaluate(initial[:8], order=16, probes=1025)
    initial_values = driver.flux_values(initial_evaluation["mixing"])
    scales = np.maximum(
        np.abs(initial_values), np.asarray([1e-12, 1e-12, 1.0, 1.0, 1.0, 1.0])
    )
    rows = {}
    for name in ("coarse", "refined"):
        final = np.asarray(evidence[name]["states"][-1], dtype=float)
        evaluation = driver.evaluate(final[:8], order=16, probes=1025)
        values = driver.flux_values(evaluation["mixing"])
        balance = (values - initial_values - final[8:]) / scales
        rows[name] = {
            "final_values": values,
            "cumulative_sources": final[8:],
            "scaled_balance": balance,
            "maximum_scaled_balance": float(np.max(np.abs(balance))),
            "limit": 5.0e-4,
            "passed": bool(np.max(np.abs(balance)) <= 5.0e-4),
        }
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed during reach verification")
    payload = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "input_unmodified": True,
        "rows": rows,
        "all_balances_passed": all(row["passed"] for row in rows.values()),
        "candidate_promoted": False,
        "hashes": after,
    }
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    print(json.dumps({
        "coarse": rows["coarse"]["maximum_scaled_balance"],
        "refined": rows["refined"]["maximum_scaled_balance"],
        "all_balances_passed": payload["all_balances_passed"],
        "saved": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
