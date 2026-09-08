"""Read-only paired PRESLHY thermal/species profile-moment audit."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from degali.validation.thermal_profile_moments import (
    read_paired_profiles,
    summarize_paired_profile_moments,
)


ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def arrays_to_lists(value):
    if isinstance(value, dict):
        return {key: arrays_to_lists(item) for key, item in value.items()}
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report",
        type=Path,
        default=ROOT / "tmp/pdfs/PRESLHY_D3.6_Summary_of_Rainout_Experiments_V1.20.pdf",
    )
    parser.add_argument(
        "--raw-root", type=Path, default=ROOT / "reference/preslhy/raw"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "reference/preslhy/observed_thermal_profile_moments_2026-09-08.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite observed profile-moment evidence")

    trials = {
        10: ("trial_10_13-09-2019alldata.xlsx", (30, 78)),
        23: ("trial_23_18-09-2019alldata.xlsx", (28, 125)),
    }
    prereg = ROOT / "docs/prereg-observed-thermal-profile-moments.md"
    observation_operator = (
        ROOT / "src/degali/validation/thermal_profile_moments.py"
    )
    operator_tests = ROOT / "tests/test_thermal_profile_moments.py"
    dependencies = [
        args.report,
        prereg,
        observation_operator,
        operator_tests,
        Path(__file__).resolve(),
    ]
    results = {}
    for trial, (name, rows) in trials.items():
        workbook = args.raw_root / name
        dependencies.append(workbook)
        paired = read_paired_profiles(
            workbook, args.report, flexlogger_rows=rows
        )
        summary = summarize_paired_profile_moments(paired)
        results[str(trial)] = dict(
            workbook=str(workbook.relative_to(ROOT)),
            workbook_sha256=digest(workbook),
            source_rows_inclusive=list(rows),
            clipped_values=paired["clipped_values"],
            **summary,
        )

    hashes = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    output = dict(
        completed=True,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        source_dataset="PRESLHY E3.5",
        source_doi="10.35097/1481",
        source_license="CC BY-SA 4.0",
        source_files_distributed=False,
        model_coefficient_fitted=False,
        trajectory_reintegrated=False,
        diffusivity_ratio_adopted=False,
        hashes=hashes,
        trials=results,
    )
    if hashes != {str(path.relative_to(ROOT)): digest(path) for path in dependencies}:
        raise RuntimeError("an audit input changed while reading")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(arrays_to_lists(output), stream, indent=2, allow_nan=False)
    compact = {
        trial: {
            "lag": data["lag"],
            "common_delay_passed": data["common_delay_passed"],
            "paired_profiles": data["paired_profiles"],
            "station_summary": data["station_summary"],
            "growth_summary": data["growth_summary"],
        }
        for trial, data in results.items()
    }
    print(json.dumps(compact, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
