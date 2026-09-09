"""Locate the momentum-diffusivity sign loss in the stopped Trial 10 path."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_reservoir_thermal_segments import actual_source_model, serial
from degali.addons.reservoir_thermal import ReservoirThermalMoments, decode_section
from degali.addons.thermal_moment_flux_march import (
    FluxSpaceThermalMomentMarch,
    ThermalMomentFluxInverter,
)


ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference/preslhy"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def radial_summary(match, *, order=8, probes=257):
    closure = ReservoirThermalMoments(
        match.mixing, thermal_species_ratio=1.,
        mechanical_work="reduced_buoyancy_work",
    )
    evaluated = closure.evaluate(order=order, probes=probes)
    partition = match.mixing.partition
    q = np.unique(np.r_[
        np.linspace(0., partition.qmax, 2049),
        .5*(partition.knots[:-1]+partition.knots[1:]),
    ])
    fields = closure.shear.radial_fields(
        q, evaluated["family"], thermal_species_ratio=1., order=order,
    )
    at = np.array([1., evaluated["gamma"]])
    chi_species = fields["chi_species"]@at
    chi_momentum = fields["chi_momentum"]@at
    tail_fraction = np.array([.9, .99, .999, .9999, .99999, 1.])
    tail_q = partition.qmax*tail_fraction
    tail_fields = closure.shear.radial_fields(
        tail_q, evaluated["family"], thermal_species_ratio=1., order=order,
    )
    tail_species = tail_fields["chi_species"]@at
    tail_momentum = tail_fields["chi_momentum"]@at
    tail_angle = np.where(
        tail_q <= partition.q0,
        2.*math.pi,
        2.*math.pi-8.*np.arctan(np.sqrt(tail_q/partition.q0-1.)),
    )

    def minimum(values, affine):
        index = int(np.argmin(values))
        normalized = math.sqrt(2.*q[index])
        return {
            "value": float(values[index]),
            "q": float(q[index]),
            "normalized_radius": normalized,
            "lateral_distance_m": normalized*match.mixing.sy,
            "normal_distance_m": normalized*match.mixing.sn,
            "affine_origin": float(affine[index, 0]),
            "affine_response": float(affine[index, 1]),
            "density_kg_m3": float(fields["local"]["rho"][index]),
            "velocity_gradient_q": float(fields["uq"][index]),
            "radial_stress": float(fields["stress"][index]@at),
        }

    return {
        "valid": evaluated["valid"],
        "state": match.state,
        "thermal_width_ratio": match.thermal_width_ratio,
        "thermal_width_log_rate_per_m": evaluated["gamma"],
        "inverse_maximum_scaled_residual": float(
            np.max(np.abs(match.scaled_residuals))
        ),
        "weak_budget_scaled_error": evaluated["weak_budget_scaled_error"],
        "incoming": evaluated["incoming"],
        "positive_diffusion": evaluated["positive_diffusion"],
        "maximum_outward_mass": evaluated["maximum_outward_mass"],
        "curvature_half_width": evaluated["curvature_half_width"],
        "minimum_species": minimum(chi_species, fields["chi_species"]),
        "minimum_momentum": minimum(chi_momentum, fields["chi_momentum"]),
        "negative_species_samples": int(np.count_nonzero(chi_species < 0.)),
        "negative_momentum_samples": int(np.count_nonzero(chi_momentum < 0.)),
        "radial_samples": int(len(q)),
        "outer_tail": [
            {
                "q_over_qmax": float(fraction),
                "angular_measure_rad": float(angle),
                "chi_species": float(species),
                "chi_momentum": float(momentum),
            }
            for fraction, angle, species, momentum in zip(
                tail_fraction, tail_angle, tail_species, tail_momentum
            )
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--attempt", type=Path,
        default=REF / "thermal_moment_observation_validation_first_attempt_2026-09-10.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=REF / "momentum_diffusivity_collapse_2026-09-10.json",
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite momentum-diffusivity diagnosis")
    reduced_path = REF / "e35_reduced.json"
    measured_path = REF / "measured_pipe_source_2026-09-05.json"
    dependencies = [
        args.attempt, reduced_path, measured_path,
        ROOT / "docs/prereg-momentum-diffusivity-collapse.md",
        ROOT / "src/degali/addons/reservoir_thermal.py",
        ROOT / "src/degali/addons/shear_thermal.py",
        ROOT / "src/degali/addons/thermal_moment_flux_march.py",
        ROOT / "src/degali/addons/transverse_mixing.py",
        Path(__file__).resolve(),
    ]
    before = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    attempt = json.loads(args.attempt.read_text(encoding="utf-8"))
    trial = next(
        row for row in json.loads(reduced_path.read_text(encoding="utf-8"))["trials"]
        if row["trial"] == 10
    )
    measured = {
        row["trial"]: row
        for row in json.loads(measured_path.read_text(encoding="utf-8"))["trials"]
    }
    row = next(row for row in attempt["rows"] if row["trial"] == 10)
    coarse = row["coarse"]
    parameters = np.asarray(coarse["parameters"], float)[-1]
    fluxes = np.asarray(coarse["fluxes"], float)[-1]
    state, _beta = decode_section(parameters)
    jp, thermodynamics, source = actual_source_model(trial, measured)
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
    accepted_match = inverter.match(
        primary[:6], parameters, position=primary[6:8], tolerance=2e-5,
    )
    diagnostics = dict(
        rhs_calls=0, maximum_inverse_residual=0., maximum_weak_residual=0.,
        minimum_sampled_diffusivity=math.inf, maximum_edge_heat_defect=0.,
    )
    first_rhs, first_parameters = marcher._rhs(primary, parameters, diagnostics)
    half_primary = primary+.0025*first_rhs
    half_match = inverter.match(
        half_primary[:6], first_parameters, position=half_primary[6:8],
        tolerance=2e-5,
    )
    payload = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "trial": 10,
        "source": source,
        "source_attempt": str(args.attempt.relative_to(ROOT)),
        "rejected_nominal_step_m": .005,
        "accepted_section": radial_summary(accepted_match),
        "rejected_half_stage": radial_summary(half_match),
        "coefficient_fitted": False,
        "acceptance_changed": False,
        "hashes": before,
    }
    after = {str(path.relative_to(ROOT)): digest(path) for path in dependencies}
    if before != after:
        raise RuntimeError("an input changed during radial diagnosis")
    payload["hashes"] = after
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, default=serial, indent=2, allow_nan=False)
    print(json.dumps(payload, default=serial, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
