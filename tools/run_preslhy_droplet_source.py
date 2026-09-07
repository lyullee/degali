"""Freeze post-flash liquid/vapour source diagnostics for PRESLHY nozzles."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.addons.lh2_droplets import (  # noqa: E402
    critical_droplet_diameter_for_phase_delay,
    critical_diffusivity_for_phase_delay,
    flashing_hydrogen_droplet_source,
    gasflow_phase_relaxation_coefficient,
    homogeneous_equilibrium_hydrogen_source,
    minimum_heat_limited_hydrogen_evaporation_time,
)
from degali.addons.cryogenic_air import (  # noqa: E402
    particle_relaxation_time,
    particle_terminal_velocity,
)


DEFAULT_REDUCED = ROOT / "reference" / "preslhy" / "e35_reduced.json"
DEFAULT_SOURCE = (
    ROOT / "reference" / "preslhy" / "measured_pipe_source_2026-09-05.json"
)
DEFAULT_OUTPUT = (
    ROOT / "reference" / "preslhy"
    / "measured_pipe_droplet_source_2026-09-05.json"
)

# NBS Monograph 168, Table 25: parahydrogen self-diffusion in
# orthohydrogen at 20.4 K and 0.1 MPa.  It is a scale comparison, not a claim
# that self-diffusion equals the effective vapour-to-droplet coefficient.
REFERENCE_HYDROGEN_DIFFUSIVITY = 0.8016e-6


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reduced", type=Path, default=DEFAULT_REDUCED)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--hydrogen-spin-isomer",
        choices=("normal", "para", "ortho"),
        default="normal",
    )
    args = parser.parse_args()

    hydrogen_species = {
        "normal": "Hydrogen",
        "para": "ParaHydrogen",
        "ortho": "OrthoHydrogen",
    }[args.hydrogen_spin_isomer]

    reduced = json.loads(args.reduced.read_text(encoding="utf-8"))
    measured = json.loads(args.source.read_text(encoding="utf-8"))
    conditions = {int(row["trial"]): row for row in reduced["trials"]}
    from CoolProp.CoolProp import PropsSI

    records = []
    for row in measured["trials"]:
        trial = int(row["trial"])
        if row["topology"] != "nozzle":
            continue
        condition = conditions[trial]
        common = dict(
            mass_flow=row["pressure_loss_mass_flow_g_s"] / 1000.0,
            orifice_diameter=condition["orifice_mm"] / 1000.0,
            upstream_temperature=row["tc3_temperature_k"]["median"],
            upstream_pressure=101325.0 + row["pt2_barg"]["median"] * 1.0e5,
            ambient_temperature=273.15 + condition["T_C"],
            ambient_pressure=101325.0,
        )
        sensitivities = {}
        air_density = float(PropsSI(
            "D", "T", common["ambient_temperature"],
            "P", common["ambient_pressure"], "Air"
        ))
        air_viscosity = float(PropsSI(
            "V", "T", common["ambient_temperature"],
            "P", common["ambient_pressure"], "Air"
        ))
        air_conductivity = float(PropsSI(
            "L", "T", common["ambient_temperature"],
            "P", common["ambient_pressure"], "Air"
        ))
        air_cp = float(PropsSI(
            "C", "T", common["ambient_temperature"],
            "P", common["ambient_pressure"], "Air"
        ))
        air_prandtl = air_cp * air_viscosity / air_conductivity
        for coefficient in (10.0, 15.0, 20.0):
            result = flashing_hydrogen_droplet_source(
                **common,
                droplet_size_coefficient=coefficient,
                hydrogen_species=hydrogen_species,
            )
            result_record = asdict(result)
            result_record["liquid_mass_fraction"] = result.liquid_mass_fraction
            if result.droplet_diameter is not None:
                diameter = result.droplet_diameter
                maximum_relative_reynolds = (
                    air_density * result.postflash_velocity * diameter
                    / air_viscosity
                )
                stagnant_lifetime = (
                    minimum_heat_limited_hydrogen_evaporation_time(
                        diameter=diameter,
                        gas_temperature=common["ambient_temperature"],
                        droplet_temperature=result.postflash_temperature,
                        gas_thermal_conductivity=air_conductivity,
                        gas_prandtl=air_prandtl,
                        hydrogen_species=hydrogen_species,
                    )
                )
                forced_lifetime = (
                    minimum_heat_limited_hydrogen_evaporation_time(
                        diameter=diameter,
                        gas_temperature=common["ambient_temperature"],
                        droplet_temperature=result.postflash_temperature,
                        gas_thermal_conductivity=air_conductivity,
                        gas_prandtl=air_prandtl,
                        particle_reynolds=maximum_relative_reynolds,
                        hydrogen_species=hydrogen_species,
                    )
                )
                result_record["isolated_ambient_diagnostic"] = {
                    "meaning": (
                        "fastest-heating bound only; production droplets see "
                        "local cold jet gas and its finite heat capacity"
                    ),
                    "air_prandtl": air_prandtl,
                    "maximum_relative_reynolds": maximum_relative_reynolds,
                    "stagnant_evaporation_time_s": stagnant_lifetime,
                    "maximum_forced_evaporation_time_s": forced_lifetime,
                    "stagnant_evaporation_distance_m": (
                        stagnant_lifetime * result.postflash_velocity
                    ),
                    "maximum_forced_evaporation_distance_m": (
                        forced_lifetime * result.postflash_velocity
                    ),
                    "velocity_relaxation_time_s": particle_relaxation_time(
                        diameter=diameter,
                        particle_density=result.liquid_density,
                        gas_viscosity=air_viscosity,
                    ),
                    "terminal_settling_velocity_m_s": (
                        particle_terminal_velocity(
                            diameter=diameter,
                            particle_density=result.liquid_density,
                            gas_density=air_density,
                            gas_viscosity=air_viscosity,
                        )
                    ),
                }
            sensitivities[str(int(coefficient))] = result_record
        central = sensitivities["15"]
        collective_bound = homogeneous_equilibrium_hydrogen_source(
            mass_flow=central["mass_flow"],
            orifice_diameter=central["orifice_diameter"],
            upstream_temperature=central["upstream_temperature"],
            upstream_pressure=central["upstream_pressure"],
            ambient_temperature=common["ambient_temperature"],
            ambient_pressure=common["ambient_pressure"],
            hydrogen_species=hydrogen_species,
        )
        collective_endpoint = collective_bound.phase_plane
        postflash_transit_time = (
            collective_endpoint.formation_distance
            / central["postflash_velocity"]
        )
        critical_reference_diameter = (
            critical_droplet_diameter_for_phase_delay(
                diffusion_coefficient=REFERENCE_HYDROGEN_DIFFUSIVITY,
                delay_time=postflash_transit_time,
                schmidt_number=1.0,
            )
        )
        for sensitivity in sensitivities.values():
            relaxation_coefficient = gasflow_phase_relaxation_coefficient(
                droplet_diameter=sensitivity["droplet_diameter"],
                diffusion_coefficient=REFERENCE_HYDROGEN_DIFFUSIVITY,
                schmidt_number=1.0,
            )
            critical_diffusivity = critical_diffusivity_for_phase_delay(
                droplet_diameter=sensitivity["droplet_diameter"],
                delay_time=postflash_transit_time,
                schmidt_number=1.0,
            )
            sensitivity["gasflow_relaxation_applicability"] = {
                "target_delay_time_s": postflash_transit_time,
                "meaning": (
                    "stagnant-droplet effective diffusivity required for "
                    "1/c to equal the collective endpoint transit time"
                ),
                "critical_diffusivity_m2_s": critical_diffusivity,
                "reference_self_diffusivity_m2_s": (
                    REFERENCE_HYDROGEN_DIFFUSIVITY
                ),
                "critical_to_reference_diffusivity_ratio": (
                    critical_diffusivity / REFERENCE_HYDROGEN_DIFFUSIVITY
                ),
                "reference_relaxation_coefficient_1_s": (
                    relaxation_coefficient
                ),
                "reference_relaxation_time_s": 1.0 / relaxation_coefficient,
                "critical_reference_droplet_diameter_m": (
                    critical_reference_diameter
                ),
                "critical_to_correlated_diameter_ratio": (
                    critical_reference_diameter
                    / sensitivity["droplet_diameter"]
                ),
            }
        records.append(
            {
                "trial": trial,
                "source_file": row["file"],
                "ambient_temperature_k": common["ambient_temperature"],
                "collective_all_vapour_endpoint_diagnostic": {
                    "meaning": (
                        "enthalpy-limited integral endpoint; not a spatial "
                        "finite-rate solution"
                    ),
                    "entrained_air_to_hydrogen_mass_ratio": (
                        collective_endpoint.endpoint.air_ratio
                    ),
                    "hydrogen_mass_fraction": (
                        collective_endpoint.endpoint.hydrogen_mass_fraction
                    ),
                    "temperature_k": collective_endpoint.endpoint.temperature,
                    "constant_entrainment_formation_distance_m": (
                        collective_endpoint.formation_distance
                    ),
                    "postflash_transit_time_s": postflash_transit_time,
                    "velocity_m_s": collective_endpoint.velocity,
                    "diameter_m": collective_endpoint.diameter,
                    "relative_energy_residual": (
                        collective_endpoint.endpoint.relative_energy_residual
                    ),
                    "hydrogen_mass_residual": (
                        collective_bound.hydrogen_mass_residual
                    ),
                    "total_mass_residual": collective_bound.total_mass_residual,
                    "momentum_residual": collective_bound.momentum_residual,
                },
                "sensitivities": sensitivities,
            }
        )

    payload = {
        "protocol": "docs/prereg-preslhy-finite-rate-lh2-droplets.md",
        "dataset_doi": "10.35097/1481",
        "source_table": args.source.name,
        "status": "source diagnostic only; not field-scored",
        "hydrogen_spin_isomer": args.hydrogen_spin_isomer,
        "hydrogen_property_species": hydrogen_species,
        "droplet_size_coefficients": [10, 15, 20],
        "phase_relaxation_reference": {
            "source": "NBS Monograph 168, Table 25",
            "url": (
                "https://nvlpubs.nist.gov/nistpubs/Legacy/MONO/"
                "nbsmonograph168.pdf"
            ),
            "quantity": (
                "parahydrogen self-diffusion in orthohydrogen; scale "
                "comparison, not an effective spray fit"
            ),
            "temperature_k": 20.4,
            "pressure_mpa": 0.1,
            "diffusivity_m2_s": REFERENCE_HYDROGEN_DIFFUSIVITY,
        },
        "trials": records,
    }
    rendered = json.dumps(payload, indent=2, ensure_ascii=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    print(f"saved: {args.output}")
    for row in records:
        central = row["sensitivities"]["15"]
        print(
            f"trial {row['trial']}: Q={central['postflash_quality']:.6f}, "
            f"liquid={central['liquid_mass_fraction']:.6f}, "
            f"d={central['droplet_diameter'] * 1e6:.3f} um"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
