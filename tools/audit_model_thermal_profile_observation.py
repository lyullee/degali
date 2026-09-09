"""Project stored models through the paired five-sensor observation operator."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from degali.validation.thermal_profile_moments import (
    project_model_profile,
    read_paired_profiles,
    truncated_profile_moment,
)
from run_preslhy_ambient_profile_audit import replay


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distribution(values) -> dict:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or not array.size or not np.all(np.isfinite(array)):
        raise ValueError("a nonempty finite one-dimensional distribution is required")
    return {
        "count": int(array.size),
        "minimum": float(np.min(array)),
        "p25": float(np.percentile(array, 25)),
        "median": float(np.median(array)),
        "p75": float(np.percentile(array, 75)),
        "maximum": float(np.max(array)),
    }


def observed_profile_summary(coordinate, values, *, centre_gate: float) -> dict:
    """Summarize amplitudes and direct moments without fitting outer tails."""
    z = np.asarray(coordinate, dtype=float)
    profiles = np.asarray(values, dtype=float)
    if profiles.ndim != 2 or profiles.shape[1:] != z.shape:
        raise ValueError("observed profiles must be time by coordinate")
    if not np.all(np.isfinite(profiles)) or np.any(profiles < 0.0):
        raise ValueError("observed profiles must be finite and nonnegative")
    centre = int(np.argmin(np.abs(z)))
    if abs(z[centre]) > 1.0e-12:
        raise ValueError("observed coordinates must include the release axis")
    selected = profiles[profiles[:, centre] >= centre_gate]
    if not len(selected):
        raise ValueError("no observed profile passes the frozen centre gate")
    moments = [truncated_profile_moment(z, profile) for profile in selected]
    return {
        "usable_profiles": len(moments),
        "centre_amplitude": distribution(selected[:, centre]),
        "zeroth_moment": distribution([moment.zeroth for moment in moments]),
        "centroid": distribution([moment.centre for moment in moments]),
        "variance": distribution([moment.variance for moment in moments]),
        "width": distribution([moment.width for moment in moments]),
    }


def model_summary(projection, observed: dict) -> dict:
    centre = int(np.argmin(np.abs(projection.coordinate)))
    fields = {
        "thermal_deficit": {
            "centre_amplitude": float(projection.thermal_deficit[centre]),
            "zeroth_moment": projection.thermal_moment.zeroth,
            "centroid": projection.thermal_moment.centre,
            "variance": projection.thermal_moment.variance,
            "width": projection.thermal_moment.width,
        },
        "hydrogen": {
            "centre_amplitude": float(projection.hydrogen[centre]),
            "zeroth_moment": projection.hydrogen_moment.zeroth,
            "centroid": projection.hydrogen_moment.centre,
            "variance": projection.hydrogen_moment.variance,
            "width": projection.hydrogen_moment.width,
        },
    }
    comparisons = {}
    for scalar, scalar_fields in fields.items():
        comparisons[scalar] = {}
        for name, predicted in scalar_fields.items():
            reference = observed[scalar][name]
            median = reference["median"]
            entry = {
                "predicted": predicted,
                "observed_median": median,
                "observed_p25": reference["p25"],
                "observed_p75": reference["p75"],
                "difference": predicted - median,
                "within_observed_iqr": bool(
                    reference["p25"] <= predicted <= reference["p75"]
                ),
            }
            if predicted > 0.0 and median > 0.0 and name != "centroid":
                entry["predicted_over_observed"] = predicted / median
                entry["absolute_log_ratio"] = abs(math.log(predicted / median))
            else:
                entry["predicted_over_observed"] = None
                entry["absolute_log_ratio"] = None
            comparisons[scalar][name] = entry
    return {
        "projection": {
            "station_m": projection.station,
            "release_height_m": projection.release_height,
            "relative_height_m": projection.coordinate.tolist(),
            "absolute_height_m": projection.absolute_height.tolist(),
            "ambient_temperature_K": projection.ambient_temperature,
            "temperature_K": projection.temperature.tolist(),
            "thermal_deficit_K": projection.thermal_deficit.tolist(),
            "hydrogen_vol_pct": projection.hydrogen.tolist(),
            "thermal_moment": asdict(projection.thermal_moment),
            "hydrogen_moment": asdict(projection.hydrogen_moment),
        },
        "comparison": comparisons,
    }


def aggregate(rows: list[dict]) -> dict:
    output = {}
    for scalar in ("thermal_deficit", "hydrogen"):
        output[scalar] = {}
        for field in ("centre_amplitude", "variance"):
            errors = [
                row["comparison"][scalar][field]["absolute_log_ratio"]
                for row in rows
            ]
            if any(error is None for error in errors):
                output[scalar][field] = {"complete": False}
                continue
            output[scalar][field] = {
                "complete": True,
                "profiles": len(errors),
                "median_absolute_log_ratio": float(np.median(errors)),
                "maximum_absolute_log_ratio": float(np.max(errors)),
                "within_observed_iqr": sum(
                    row["comparison"][scalar][field]["within_observed_iqr"]
                    for row in rows
                ),
            }
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REF / "model_thermal_profile_observation_2026-09-08.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite model-observation evidence")

    report = ROOT / "tmp/pdfs/PRESLHY_D3.6_Summary_of_Rainout_Experiments_V1.20.pdf"
    reduced_path = REF / "e35_reduced.json"
    observed_audit = REF / "observed_thermal_profile_moments_2026-09-08.json"
    model_paths = {
        "density_profile_control": REF / "phase_ambient_consistency_field_complete_2026-09-05.json",
        "same_width_enthalpy_profile": REF / "gaussian_enthalpy_profile_allflux_field_complete_2026-09-05.json",
    }
    trial_files = {
        10: (REF / "raw/trial_10_13-09-2019alldata.xlsx", (30, 78)),
        23: (REF / "raw/trial_23_18-09-2019alldata.xlsx", (28, 125)),
    }
    dependencies = [
        report,
        reduced_path,
        observed_audit,
        ROOT / "docs/prereg-model-thermal-profile-observation.md",
        ROOT / "docs/prereg-observed-thermal-profile-moments.md",
        ROOT / "src/degali/validation/thermal_profile_moments.py",
        ROOT / "src/degali/validation/nearfield.py",
        ROOT / "src/degali/addons/axisymmetric_jet.py",
        ROOT / "src/degali/addons/energy_crosswind.py",
        ROOT / "src/degali/addons/enthalpy_profile.py",
        ROOT / "tests/test_thermal_profile_moments.py",
        ROOT / "tests/test_model_thermal_profile_observation.py",
        ROOT / "tools/run_preslhy_ambient_profile_audit.py",
        Path(__file__).resolve(),
        *model_paths.values(),
        *(path for path, _ in trial_files.values()),
    ]
    hashes_before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    reduced = json.loads(reduced_path.read_text(encoding="utf-8"))
    trial_records = {trial["trial"]: trial for trial in reduced["trials"]}
    fields = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in model_paths.items()
    }

    observed = {}
    model_rows = {name: [] for name in fields}
    replay_checks = {}
    for number, (workbook, rows) in trial_files.items():
        trial = trial_records[number]
        paired = read_paired_profiles(workbook, report, flexlogger_rows=rows)
        observed[str(number)] = {}
        for station in (1.78, 4.0):
            profiles = paired["profiles"][str(station)]
            station_observed = {
                "thermal_deficit": observed_profile_summary(
                    paired["coordinate"], profiles["thermal_deficit"], centre_gate=5.0
                ),
                "hydrogen": observed_profile_summary(
                    paired["coordinate"], profiles["hydrogen"], centre_gate=1.0
                ),
            }
            observed[str(number)][str(station)] = station_observed
        for variant, field in fields.items():
            trajectory, check = replay(field, trial)
            replay_checks[f"{variant}_{number}"] = check
            if trajectory is None:
                raise ValueError(f"stored {variant} trajectory is missing Trial {number}")
            for station in (1.78, 4.0):
                projection = project_model_profile(
                    trajectory, station, trial["release_height_m"]
                )
                row = {
                    "variant": variant,
                    "trial": number,
                    "station_m": station,
                    **model_summary(
                        projection, observed[str(number)][str(station)]
                    ),
                }
                model_rows[variant].append(row)

    aggregates = {variant: aggregate(rows) for variant, rows in model_rows.items()}
    control = aggregates["density_profile_control"]["thermal_deficit"]
    candidate = aggregates["same_width_enthalpy_profile"]["thermal_deficit"]
    componentwise = all(
        candidate[field]["median_absolute_log_ratio"]
        <= control[field]["median_absolute_log_ratio"]
        for field in ("centre_amplitude", "variance")
    )
    strictly_better = any(
        candidate[field]["median_absolute_log_ratio"]
        < control[field]["median_absolute_log_ratio"]
        for field in ("centre_amplitude", "variance")
    )
    hashes_after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if hashes_before != hashes_after:
        raise RuntimeError("an audit input changed while calculating")
    payload = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": "PRESLHY E3.5",
        "source_doi": "10.35097/1481",
        "source_files_distributed": False,
        "trajectory_reintegrated": False,
        "coefficient_fitted": False,
        "thermal_width_transport_tested": False,
        "hashes": hashes_after,
        "replay_checks": replay_checks,
        "observed": observed,
        "model_rows": model_rows,
        "aggregates": aggregates,
        "same_width_enthalpy_dominates_control": bool(componentwise and strictly_better),
        "promotion_allowed": False,
        "limitations": [
            "A steady prediction is compared with the observed release-window distribution.",
            "Temperature deficit is an observation proxy, not transported enthalpy.",
            "The five-point variance is truncated and contains no fitted outer tail.",
            "The enthalpy-profile path keeps thermal and species widths equal.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
    print(json.dumps({
        "aggregates": aggregates,
        "same_width_enthalpy_dominates_control": payload["same_width_enthalpy_dominates_control"],
        "saved": str(args.output),
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
