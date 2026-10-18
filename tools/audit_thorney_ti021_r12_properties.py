"""Freeze cited R12 property inputs for the TI021 full-mixture preflight."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from CoolProp.CoolProp import PropsSI, get_global_param_string


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=(
            ROOT / "outputs" / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
        ),
    )
    args = parser.parse_args()

    temperature_k = 291.15
    pressure_pa = 100_500.0
    # Matsunaga, Hori & Nagashima (1994), measured by the Taylor method.
    diffusivity_m2_s = 4.77e-10 * temperature_k**1.75
    fluid = "R12"
    coolprop = {
        "version": get_global_param_string("version"),
        "git_revision": get_global_param_string("gitrevision"),
        "fluid": fluid,
        "molar_mass_kg_mol": PropsSI("M", "T", temperature_k, "P", pressure_pa, fluid),
        "density_kg_m3": PropsSI(
            "Dmass", "T", temperature_k, "P", pressure_pa, fluid,
        ),
        "isobaric_heat_capacity_j_kg_k": PropsSI(
            "Cpmass", "T", temperature_k, "P", pressure_pa, fluid,
        ),
        "dynamic_viscosity_pa_s": PropsSI(
            "viscosity", "T", temperature_k, "P", pressure_pa, fluid,
        ),
        "thermal_conductivity_w_m_k": PropsSI(
            "conductivity", "T", temperature_k, "P", pressure_pa, fluid,
        ),
        "prandtl_number": PropsSI(
            "Prandtl", "T", temperature_k, "P", pressure_pa, fluid,
        ),
    }
    coolprop["thermal_diffusivity_m2_s"] = (
        coolprop["thermal_conductivity_w_m_k"]
        / (
            coolprop["density_kg_m3"]
            * coolprop["isobaric_heat_capacity_j_kg_k"]
        )
    )
    gates = {
        "event_temperature_inside_measured_diffusivity_range": (
            273.15 <= temperature_k <= 453.15
        ),
        "event_pressure_is_atmospheric_within_1pct": (
            abs(pressure_pa / 101_325.0 - 1.0) <= 0.01
        ),
        "coolprop_r12_molar_mass_matches_workbook_species": (
            abs(coolprop["molar_mass_kg_mol"] / 0.12091 - 1.0) <= 5.0e-4
        ),
        "positive_finite_transport_properties": all(
            value > 0.0
            for key, value in coolprop.items()
            if key not in {"version", "git_revision", "fluid"}
        ) and diffusivity_m2_s > 0.0,
    }
    record = {
        "schema": "degali.thorney-ti021-r12-property-audit.v1",
        "status": (
            "accepted_molecular_property_component"
            if all(gates.values()) else "withheld_molecular_property_component"
        ),
        "scope": (
            "molecular R12-air diffusion and pure-R12 ambient gas properties; "
            "does not qualify turbulent scalar mixing"
        ),
        "event_state": {
            "temperature_k": temperature_k,
            "pressure_pa": pressure_pa,
        },
        "r12_air_binary_diffusion": {
            "diffusivity_m2_s": diffusivity_m2_s,
            "measurement_equation": "D = 4.77e-10 * T^1.75 m2/s",
            "temperature_range_k": [273.15, 453.15],
            "pressure_basis": "atmospheric pressure",
            "reported_experimental_uncertainty_percent": 2.0,
            "citation": {
                "authors": "Naoki Matsunaga, Morio Hori, Akira Nagashima",
                "year": 1994,
                "title": (
                    "Measurements of the Mutual Diffusion Coefficients of Gases "
                    "by the Taylor Method, 4th Report"
                ),
                "doi": "10.1299/kikaib.60.504",
                "url": (
                    "https://www.jstage.jst.go.jp/article/kikaib1979/60/570/"
                    "60_570_504/_article"
                ),
            },
            "frozen_sensitivity_bounds_m2_s": [
                0.98 * diffusivity_m2_s,
                1.02 * diffusivity_m2_s,
            ],
        },
        "pure_r12_ambient_properties": coolprop,
        "property_software_provenance": {
            "model": "CoolProp R12 backend",
            "reference_database": (
                "NIST REFPROP provides R12 thermodynamic and transport-property "
                "models over 116.099-525 K; CoolProp values are recorded with "
                "exact software revision for reproducibility"
            ),
            "nist_refprop_doi": "10.18434/T4/1502528",
            "nist_refprop_url": "https://www.nist.gov/srd/refprop",
            "coolprop_citation_url": "https://coolprop.org/coolprop/Citation.html",
        },
        "solver_policy": {
            "species_molar_mass_kg_mol": coolprop["molar_mass_kg_mol"],
            "species_heat_capacity_j_kg_k": coolprop[
                "isobaric_heat_capacity_j_kg_k"
            ],
            "molecular_species_diffusivity_m2_s": diffusivity_m2_s,
            "initial_temperature_k": temperature_k,
            "initial_sensible_enthalpy_density_j_m3": 0.0,
            "initial_state_cp_sensitivity": (
                "none: initial gas and carrier are at the enthalpy reference "
                "temperature, so initial density and buoyancy do not depend on cp"
            ),
            "thermal_diffusivity_policy": (
                "energy field remains identically zero for the isothermal TI021 "
                "initial condition; pure-R12 thermal diffusivity is diagnostic, "
                "not a binary-mixture thermal closure"
            ),
            "turbulent_diffusivity_policy": (
                "not inferred from molecular diffusion; turbulent mixing closure "
                "must be frozen independently of TI021 concentration observations"
            ),
        },
        "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": record["status"],
        "output": str(args.output.resolve()),
        "r12_air_diffusivity_m2_s": diffusivity_m2_s,
        "r12_cp_j_kg_k": coolprop["isobaric_heat_capacity_j_kg_k"],
        "coolprop_version": coolprop["version"],
        "gates": gates,
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
