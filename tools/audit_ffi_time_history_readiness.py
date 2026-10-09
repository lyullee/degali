"""Audit whether an FFI event can enter the causal time-history operator.

This tool never fabricates a source or sensor time series from a run-average
or a steady-model setting. It records available files, hashes, headers, and
the resulting conditional/withheld gate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


SCHEMA = "degali.ffi-time-history-readiness.v1"


def fingerprint(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        rows = sum(1 for _ in reader)
    return {
        "path": str(path.resolve()),
        "sha256": digest,
        "row_count": rows,
        "columns": header,
    }


def fingerprint_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "top_level_keys": sorted(payload) if isinstance(payload, dict) else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slabx-root", type=Path, required=True)
    parser.add_argument("--paired-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--event-id", default="FFI-Test4")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    inputs = args.slabx_root / "model-comparison" / "inputs"
    conditions = args.slabx_root / "slabx-lh2" / "lh2_ffi_conditions.csv"
    high = inputs / "ffi_test4_high_anemometer_digitized_5s.csv"
    low = inputs / "ffi_test4_low_anemometer_digitized_5s.csv"
    source_state = (
        args.slabx_root / "tmp" / "urt_20261007" / "q3" / "002"
        / "ffi_test4_common_gas_handoff_v1" / "source_term_mass_balance.csv"
    )
    wind_provenance = (
        args.slabx_root / "tmp" / "urt_20261007" / "q3" / "002"
        / "ffi_state_transition_review_20260922" / "test4_measured_history_mean_no_fit"
        / "provenance.json"
    )
    for path in (args.paired_csv, conditions, high, low, source_state, wind_provenance):
        if not path.is_file():
            raise FileNotFoundError(path)

    paired = fingerprint(args.paired_csv)
    condition_record = fingerprint(conditions)
    high_record = fingerprint(high)
    low_record = fingerprint(low)
    source_state_record = fingerprint(source_state)
    wind_provenance_record = fingerprint_json(wind_provenance)
    source_columns = set(condition_record["columns"])
    source_state_columns = set(source_state_record["columns"])
    receptor_columns = set(paired["columns"])
    weather_columns = set(high_record["columns"])
    source_history_columns = {"time_s", "mass_rate_kg_s"}
    receptor_history_columns = {"sensor", "timestamp", "observed_mole_fraction"}
    source_history_available = source_history_columns <= source_columns
    source_state_history_available = {"time_s", "mass_rate_kg_s"} <= source_state_columns
    receptor_history_available = receptor_history_columns <= receptor_columns
    weather_history_available = {"start_s", "end_s", "speed_m_s", "from_bearing_deg"} <= weather_columns

    gates: list[str] = []
    reasons: list[str] = []
    actions: list[str] = []
    if not source_history_available:
        gates.append("source_history_missing")
        reasons.append(
            "the FFI source file contains run-average outflow/run duration but no time_s + mass_rate_kg_s history"
        )
        actions.append("obtain a hash-pinned source flow/pressure/temperature time series on the event clock")
    if source_state_columns and not source_state_history_available:
        reasons.append(
            "the discovered source-state handoff is a one-row P04 pressure/temperature/quality record, not q(t)"
        )
    if not receptor_history_available:
        gates.append("receptor_time_history_missing")
        reasons.append(
            "the paired receptor file contains 275 s aggregate means, not timestamped sensor observations"
        )
        actions.append("obtain timestamped receptor samples and sensor response/calibration metadata")
    if weather_history_available:
        gates.append("weather_history_available")
        reasons.append("digitized high/low weather records are available for a 25-300 s window")
    else:
        gates.append("weather_history_missing")
        reasons.append("weather history lacks the required time and vector columns")
        actions.append("obtain a hash-pinned weather vector history covering the source and observation clock")
    gates.append("common_clock_partial")
    reasons.append(
        "weather timestamps and 275 s receptor means do not establish a common source/weather/receptor clock"
    )
    actions.append("declare and hash-pin the common clock origin, offsets, and observation operator")
    status = "withheld" if not (source_history_available and receptor_history_available and weather_history_available) else "conditional"
    gates.append("time_history_execution_withheld" if status == "withheld" else "time_history_execution_conditional")
    payload = {
        "schema": SCHEMA,
        "event_id": args.event_id,
        "status": status,
        "promotion_allowed": False,
        "gate_codes": gates,
        "reasons": reasons,
        "required_actions": actions,
        "available_inputs": {
            "source_boundary": condition_record,
            "source_state_handoff": source_state_record,
            "weather_high": high_record,
            "weather_low": low_record,
            "wind_history_provenance": wind_provenance_record,
            "receptor_observations": paired,
        },
        "checks": {
            "source_history_columns": sorted(source_history_columns),
            "source_history_available": source_history_available,
            "source_state_history_available": source_state_history_available,
            "weather_history_available": weather_history_available,
            "receptor_history_columns": sorted(receptor_history_columns),
            "receptor_history_available": receptor_history_available,
            "common_clock_established": False,
        },
        "operator_contract": {
            "manifest_bound_entrypoint": "run_manifest_bound_time_history",
            "packet_operator": "replay_fixed_receptors_with_source_history",
            "steady_275s_setting_is_not_used_as_source_history": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(json.dumps({"status": status, "gate_codes": gates, "output": str(args.output.resolve())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
