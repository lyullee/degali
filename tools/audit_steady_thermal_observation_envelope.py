"""Compare stored steady thermal predictions to raw PRESLHY time envelopes.

This is deliberately a read-only observation audit: it neither re-integrates
a trajectory nor changes a transport coefficient.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"
sys.path.insert(0, str(ROOT / "src"))

from degali.validation.thermal_profile_moments import (
    best_integer_lag,
    empirical_scalar_envelope,
    read_paired_profiles,
)

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REF / "steady_thermal_observation_envelope_2026-09-17.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite observation evidence")

    report = ROOT / "tmp/pdfs/PRESLHY_D3.6_Summary_of_Rainout_Experiments_V1.20.pdf"
    model_path = REF / "model_thermal_profile_observation_2026-09-08.json"
    trials = {
        10: (REF / "raw/trial_10_13-09-2019alldata.xlsx", (30, 78)),
        23: (REF / "raw/trial_23_18-09-2019alldata.xlsx", (28, 125)),
    }
    dependencies = [
        report,
        model_path,
        ROOT / "docs/prereg-steady-thermal-observation-envelope.md",
        ROOT / "src/degali/validation/thermal_profile_moments.py",
        ROOT / "tests/test_thermal_observation_envelope.py",
        Path(__file__).resolve(),
        *(path for path, _ in trials.values()),
    ]
    hashes_before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    models = json.loads(model_path.read_text(encoding="utf-8"))["model_rows"]
    predicted = {
        (variant, row["trial"], float(row["station_m"])):
        row["comparison"]["thermal_deficit"]["centre_amplitude"]["predicted"]
        for variant, rows in models.items()
        for row in rows
    }
    comparison = {}
    streamwise_decay = {}
    for trial, (workbook, rows) in trials.items():
        paired = read_paired_profiles(workbook, report, flexlogger_rows=rows)
        comparison[str(trial)] = {}
        for station in (1.78, 4.0):
            centre_trace = paired["profiles"][str(station)]["thermal_deficit"][:, 2]
            selected = centre_trace[centre_trace >= 5.0]
            if selected.size < 2:
                raise ValueError("frozen centre signal gate left too few samples")
            comparison[str(trial)][str(station)] = {
                variant: asdict(empirical_scalar_envelope(
                    selected, predicted[(variant, trial, station)]
                ))
                for variant in models
            }
        upstream = paired["profiles"]["1.78"]["thermal_deficit"][:, 2]
        downstream = paired["profiles"]["4.0"]["thermal_deficit"][:, 2]
        lag = best_integer_lag(upstream, downstream, maximum_seconds=5)
        if lag.seconds:
            paired_upstream = upstream[:-lag.seconds]
            paired_downstream = downstream[lag.seconds:]
        else:
            paired_upstream, paired_downstream = upstream, downstream
        selected = (paired_upstream >= 5.0) & (paired_downstream >= 5.0)
        observed_ratio = paired_downstream[selected] / paired_upstream[selected]
        if observed_ratio.size < 2:
            raise ValueError("frozen paired signal gates left too few samples")
        streamwise_decay[str(trial)] = {
            "lag_seconds": lag.seconds,
            "lag_correlation": lag.correlation,
            "valid_pairs": int(observed_ratio.size),
            "variants": {
                variant: asdict(empirical_scalar_envelope(
                    observed_ratio,
                    predicted[(variant, trial, 4.0)]
                    / predicted[(variant, trial, 1.78)],
                ))
                for variant in models
            },
        }

    hashes_after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if hashes_before != hashes_after:
        raise RuntimeError("an audit input changed while calculating")
    output = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": "PRESLHY E3.5",
        "source_doi": "10.35097/1481",
        "source_files_distributed": False,
        "trajectory_reintegrated": False,
        "coefficient_fitted": False,
        "default_model_changed": False,
        "hashes": hashes_after,
        "comparison": comparison,
        "streamwise_temperature_deficit_decay": streamwise_decay,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
    print(json.dumps({
        "comparison": comparison,
        "streamwise_temperature_deficit_decay": streamwise_decay,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
