"""Diagnose the first physical gate after a stopped direct-flux march."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_reservoir_thermal_segments import actual_source_model, serial
from degali.addons.reservoir_thermal import ReservoirShortSegment, decode_section
from degali.addons.thermal_moment_flux_march import (
    FluxSpaceThermalMomentMarch,
    ThermalMomentFluxInverter,
)


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trial", type=int, default=10)
    parser.add_argument(
        "--attempt", type=Path,
        default=REF / "thermal_moment_observation_validation_first_attempt_2026-09-10.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=REF / "thermal_moment_observation_stop_diagnosis_2026-09-10.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite stop diagnosis")
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    dependencies = [
        args.attempt, reduced_path, measured_path,
        ROOT / "src/degali/addons/reservoir_thermal.py",
        ROOT / "src/degali/addons/thermal_moment_flux_march.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    attempt = json.loads(args.attempt.read_text(encoding="utf-8"))
    trial = next(
        row for row in json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
        if row["trial"] == args.trial
    )
    measured = {
        row["trial"]: row
        for row in json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    }
    source_row = next(row for row in attempt["rows"] if row["trial"] == args.trial)
    coarse = source_row["coarse"]
    parameters = np.asarray(coarse["parameters"], float)[-1]
    fluxes = np.asarray(coarse["fluxes"], float)[-1]
    state, beta = decode_section(parameters)
    jp, thermodynamics, source = actual_source_model(trial, measured)
    driver = ReservoirShortSegment(
        jp, thermodynamics, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work",
    )
    accepted = driver.evaluate(parameters, order=8, probes=257)
    accepted_summary = {
        "state": state,
        "thermal_width_ratio": beta,
        "weak_budget_scaled_error": accepted["weak_budget_scaled_error"],
        "minimum_chi_species": accepted["minimum_chi_species"],
        "minimum_chi_momentum": accepted["minimum_chi_momentum"],
        "maximum_outward_mass": accepted["maximum_outward_mass"],
        "incoming": accepted["incoming"],
        "positive_diffusion": accepted["positive_diffusion"],
        "curvature_half_width": accepted["curvature_half_width"],
        "tangent_maximum_scaled_residual": accepted["family"].maximum_scaled_residual,
    }
    inverter = ThermalMomentFluxInverter(
        jp, thermodynamics, order=8, probes=257,
        quadrature_points=128, tolerance=1e-8,
    )
    marcher = FluxSpaceThermalMomentMarch(
        inverter, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work",
        stage_inverse_tolerance=2e-5,
    )
    primary = np.r_[fluxes, state[5:7]]
    probes = []
    for step in (.005, .0025, .00125, .000625, .0003125, .00015625, .000078125):
        diagnostics = dict(
            rhs_calls=0, maximum_inverse_residual=0., maximum_weak_residual=0.,
            minimum_sampled_diffusivity=math.inf, maximum_edge_heat_defect=0.,
        )
        entry = {"step_m": step}
        try:
            next_primary, next_parameters, _source_increment = marcher._rk4(
                primary, step, parameters, diagnostics,
            )
            next_state, next_beta = decode_section(next_parameters)
            entry.update(
                passed=True,
                next_downwind_m=float(next_primary[6]),
                next_state=next_state,
                next_thermal_width_ratio=next_beta,
            )
        except (ArithmeticError, FloatingPointError, RuntimeError, ValueError,
                np.linalg.LinAlgError) as error:
            entry.update(passed=False, failure=str(error))
        entry["diagnostics_before_failure"] = diagnostics
        probes.append(entry)
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed during stop diagnosis")
    payload = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "trial": args.trial,
        "source": source,
        "source_attempt": str(args.attempt.relative_to(ROOT)),
        "accepted_section": accepted_summary,
        "one_step_probes": probes,
        "step_probe_changes_acceptance": False,
        "coefficient_fitted": False,
        "hashes": after,
    }
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    print(json.dumps({
        "accepted_section": accepted_summary,
        "one_step_probes": probes,
        "saved": str(args.output),
    }, default=serial, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
