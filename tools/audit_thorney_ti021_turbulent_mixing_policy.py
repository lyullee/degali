"""Freeze a concentration-blind TI021 atmospheric mixing envelope.

Only the event friction velocity, roughness length and cited prior ranges enter
this policy.  The TI021 concentration, dose and arc observations are neither
loaded nor scored, so they cannot select the turbulent coefficients.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.atmospheric_surface_layer_mixing import (
    NeutralSurfaceLayerMixing,
)


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-turbulent-mixing-policy.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_record(
    *,
    friction_velocity_m_s: float = 0.205,
    roughness_length_m: float = 0.005,
    molecular_species_diffusivity_m2_s: float = 9.78865e-6,
    molecular_thermal_diffusivity_m2_s: float = 3.09e-6,
    heights_m: np.ndarray | None = None,
    workbook_sha256: str | None = None,
    report_sha256: str | None = None,
    property_audit_sha256: str | None = None,
) -> dict[str, object]:
    """Return the frozen neutral reference and mandatory sensitivity ensemble."""

    z = (
        np.arange(0.5, 16.0, 1.0, dtype=float)
        if heights_m is None else np.asarray(heights_m, dtype=float)
    )
    if z.ndim != 1 or z.size < 2 or np.any(~np.isfinite(z)) or np.any(z <= 0.0):
        raise ValueError("heights_m must be a positive finite one-dimensional grid")
    if np.any(np.diff(z) <= 0.0):
        raise ValueError("heights_m must be strictly increasing")
    for value, name in (
        (molecular_species_diffusivity_m2_s, "molecular species diffusivity"),
        (molecular_thermal_diffusivity_m2_s, "molecular thermal diffusivity"),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")

    closure = NeutralSurfaceLayerMixing(
        friction_velocity_m_s=friction_velocity_m_s,
        roughness_length_m=roughness_length_m,
        source_id="TI021 workbook meteorology: u_star=0.205 m/s, z0=0.005 m",
    )
    solid = np.zeros((z.size, 1, 1), dtype=bool)
    nu_t = closure.eddy_viscosity_m2_s(z, solid)[:, 0, 0]
    schmidt_values = (0.2, 0.7, 1.3)
    prandtl_values = (0.74, 1.0)
    species = {
        f"Sc_t={value:g}": (
            molecular_species_diffusivity_m2_s + nu_t / value
        ).tolist()
        for value in schmidt_values
    }
    thermal = {
        f"Pr_t={value:g}": (
            molecular_thermal_diffusivity_m2_s + nu_t / value
        ).tolist()
        for value in prandtl_values
    }
    ensemble = [
        {
            "member_id": f"neutral-Sc{schmidt:g}-Pr{prandtl:g}",
            "turbulent_schmidt_number": schmidt,
            "turbulent_prandtl_number": prandtl,
            "role": (
                "primary_prior" if schmidt == 0.7 and prandtl == 0.74
                else "mandatory_prior_sensitivity"
            ),
        }
        for schmidt in schmidt_values
        for prandtl in prandtl_values
    ]
    expected = closure.von_karman_constant * friction_velocity_m_s * z
    gates = {
        "no_concentration_dose_or_arc_input_loaded": True,
        "friction_velocity_and_roughness_declared_before_scoring": True,
        "neutral_similarity_formula_reproduced": bool(
            np.allclose(nu_t, expected, rtol=0.0, atol=1.0e-15)
        ),
        "all_momentum_species_and_thermal_diffusivities_positive": bool(
            np.all(nu_t > 0.0)
            and all(np.all(np.asarray(values) > 0.0) for values in species.values())
            and all(np.all(np.asarray(values) > 0.0) for values in thermal.values())
        ),
        "literature_range_endpoints_retained": bool(
            schmidt_values == (0.2, 0.7, 1.3)
            and prandtl_values == (0.74, 1.0)
        ),
        "stability_correction_not_inferred_without_monin_obukhov_length": True,
        "all_six_prior_members_required_for_no_fit_run": len(ensemble) == 6,
    }
    return {
        "schema": SCHEMA,
        "status": (
            "accepted_concentration_blind_turbulent_mixing_envelope_component"
            if all(gates.values())
            else "withheld_turbulent_mixing_policy_component"
        ),
        "scope": (
            "neutral atmospheric background mixing and a mandatory prior "
            "Sc_t/Pr_t sensitivity envelope; not obstacle-wake turbulence "
            "validation and not a concentration-selected calibration"
        ),
        "evidence_provenance": {
            "workbook_sha256": workbook_sha256,
            "report_sha256": report_sha256,
            "molecular_property_audit_sha256": property_audit_sha256,
            "concentration_observation_fields_accessed": [],
        },
        "event_inputs": {
            "friction_velocity_m_s": friction_velocity_m_s,
            "roughness_length_m": roughness_length_m,
            "reported_stability_class": "D/E",
            "monin_obukhov_length_m": None,
            "stability_policy": (
                "use the neutral phi_m=1 reference only; do not manufacture "
                "a stable correction from the categorical D/E label"
            ),
        },
        "closure": closure.record(),
        "prior_sources": [
            {
                "use": "neutral Monin-Obukhov eddy-viscosity relation",
                "citation": (
                    "van der Laan et al. (2017), A new k-epsilon model "
                    "consistent with Monin-Obukhov similarity theory"
                ),
                "doi": "10.1002/we.2017",
                "url": "https://onlinelibrary.wiley.com/doi/10.1002/we.2017",
            },
            {
                "use": "predeclared turbulent Schmidt sensitivity bounds",
                "citation": (
                    "Tominaga and Stathopoulos (2007), Turbulent Schmidt "
                    "numbers for CFD pollutant dispersion"
                ),
                "doi": "10.1016/j.atmosenv.2007.06.054",
                "url": (
                    "https://www.sciencedirect.com/science/article/abs/pii/"
                    "S1352231007006036"
                ),
                "reported_non_universal_range_used": [0.2, 1.3],
            },
            {
                "use": "neutral turbulent Prandtl prior endpoints",
                "citation": "Kansas/Businger-Dyer neutral similarity values",
                "values": [0.74, 1.0],
                "traceability_note": (
                    "values are frozen as a prior envelope, not inferred from "
                    "the TI021 observations"
                ),
            },
        ],
        "frozen_policy": {
            "primary_turbulent_schmidt_number": 0.7,
            "required_turbulent_schmidt_numbers": list(schmidt_values),
            "primary_turbulent_prandtl_number": 0.74,
            "required_turbulent_prandtl_numbers": list(prandtl_values),
            "selection_rule": (
                "run and retain all six prior combinations; the primary is "
                "reported first but may not be replaced after viewing scores"
            ),
            "promotion_rule": (
                "report every member and the complete error envelope; no member "
                "is promoted as calibrated from TI021"
            ),
            "ensemble": ensemble,
        },
        "profile": {
            "height_m": z.tolist(),
            "eddy_viscosity_m2_s": nu_t.tolist(),
            "molecular_species_diffusivity_m2_s": (
                molecular_species_diffusivity_m2_s
            ),
            "effective_species_diffusivity_m2_s": species,
            "molecular_thermal_diffusivity_m2_s": (
                molecular_thermal_diffusivity_m2_s
            ),
            "effective_thermal_diffusivity_m2_s": thermal,
        },
        "gates": gates,
        "all_gates_pass": all(gates.values()),
        "limitations": [
            "neutral reference only because Monin-Obukhov length is unavailable",
            "no obstacle-generated turbulence production closure",
            "gradient-diffusion Sc_t and Pr_t values are priors, not universal constants",
            "full transient grid/time/domain convergence and no-fit scoring remain required",
        ],
    }


def _plot(record: dict[str, object], output: Path) -> None:
    profile = record["profile"]
    z = np.asarray(profile["height_m"], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), constrained_layout=True)
    axes[0].plot(profile["eddy_viscosity_m2_s"], z, color="black", linewidth=2)
    axes[0].set_xlabel(r"$\nu_t$ (m$^2$/s)")
    axes[0].set_ylabel("height z (m)")
    axes[0].set_title("Neutral momentum mixing")
    for label, values in profile["effective_species_diffusivity_m2_s"].items():
        axes[1].plot(values, z, linewidth=2, label=label)
    axes[1].set_xlabel(r"$D_{eff}$ (m$^2$/s)")
    axes[1].set_title("Mandatory species envelope")
    axes[1].legend(fontsize=8)
    for label, values in profile["effective_thermal_diffusivity_m2_s"].items():
        axes[2].plot(values, z, linewidth=2, label=label)
    axes[2].set_xlabel(r"$\alpha_{eff}$ (m$^2$/s)")
    axes[2].set_title("Mandatory thermal envelope")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.suptitle(
        "TI021 concentration-blind turbulent-mixing policy\n"
        "prior envelope only — no dose or concentration used"
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--property-audit", type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
        ),
    )
    parser.add_argument(
        "--output-directory", type=Path,
        default=(
            ROOT / "outputs"
            / "thorney-ti021-turbulent-mixing-policy-2026-10-11-v1"
        ),
    )
    args = parser.parse_args()
    workbook = args.workbook.resolve()
    report = args.report.resolve()
    property_path = args.property_audit.resolve()
    property_record = json.loads(property_path.read_text(encoding="utf-8"))
    if property_record.get("status") != "accepted_molecular_property_component":
        raise ValueError("R12 molecular-property audit is not accepted")
    record = build_record(
        molecular_species_diffusivity_m2_s=float(
            property_record["solver_policy"]["molecular_species_diffusivity_m2_s"]
        ),
        molecular_thermal_diffusivity_m2_s=float(
            property_record["pure_r12_ambient_properties"][
                "thermal_diffusivity_m2_s"
            ]
        ),
        workbook_sha256=_sha256(workbook),
        report_sha256=_sha256(report),
        property_audit_sha256=_sha256(property_path),
    )
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-turbulent-mixing-policy.json"
    plot_path = args.output_directory / "ti021-turbulent-mixing-policy.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(record, plot_path)
    print(json.dumps({
        "status": record["status"],
        "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()),
        "ensemble_members": len(record["frozen_policy"]["ensemble"]),
        "gates": record["gates"],
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
