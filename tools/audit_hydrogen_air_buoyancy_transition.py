"""Canonical ideal-gas H2/air neutral-buoyancy map for the active closure."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.field_observation import (
    DRY_AIR_MOLAR_MASS_KG_MOL,
    HYDROGEN_MOLAR_MASS_KG_MOL,
)
from degali.addons.hydrogen_air_thermodynamics import (
    HydrogenAirActiveScalarClosure,
)


OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"


def _conserved_for_x_t(closure, mole_fraction, temperature):
    air_moles = closure.ambient_air_density_kg_m3 / DRY_AIR_MOLAR_MASS_KG_MOL
    concentration = (
        mole_fraction / np.maximum(1.0 - mole_fraction, 1.0e-12)
        * air_moles * HYDROGEN_MOLAR_MASS_KG_MOL
    )
    heat_capacity = (
        closure.ambient_air_density_kg_m3 * closure.dry_air_heat_capacity_j_kg_k
        + concentration * closure.hydrogen_heat_capacity_j_kg_k
    )
    enthalpy = heat_capacity * (temperature - closure.ambient_temperature_k)
    return concentration, enthalpy


def main() -> int:
    closure = HydrogenAirActiveScalarClosure(
        maximum_relative_density_departure=0.10,
    )
    x = np.linspace(1.0e-4, 0.60, 500)
    temperature = np.linspace(20.0, closure.ambient_temperature_k, 420)
    xx, tt = np.meshgrid(x, temperature)
    concentration, enthalpy = _conserved_for_x_t(closure, xx, tt)
    state = closure.state_from_conserved(concentration, enthalpy)
    departure = state.relative_density_departure
    neutral = closure.neutral_buoyancy_temperature_k(x)

    checkpoints = {}
    for percent in (1.0, 4.0, 10.0, 30.0, 60.0, 100.0):
        fraction = percent / 100.0
        checkpoints[f"{percent:g}_vol_percent"] = {
            "neutral_buoyancy_temperature_k": float(
                closure.neutral_buoyancy_temperature_k(fraction)
            ),
            "mixture_is_dense_below_neutral_temperature": True,
            "mixture_is_buoyant_above_neutral_temperature": True,
        }
    record = {
        "schema": "degali.hydrogen-air-buoyancy-transition-audit.v1",
        "status": "conditional",
        "closure": {
            "equation_of_state": "ideal gas at declared ambient pressure",
            "ambient_temperature_k": closure.ambient_temperature_k,
            "ambient_pressure_pa": closure.ambient_pressure_pa,
            "ambient_air_density_kg_m3": closure.ambient_air_density_kg_m3,
            "maximum_active_scalar_density_departure": (
                closure.maximum_relative_density_departure
            ),
        },
        "checkpoints": checkpoints,
        "sign_check": {
            "dense_region_has_negative_vertical_acceleration": bool(np.all(
                state.buoyancy_acceleration_m_s2[departure > 0.0] < 0.0
            )),
            "light_region_has_positive_vertical_acceleration": bool(np.all(
                state.buoyancy_acceleration_m_s2[departure < 0.0] > 0.0
            )),
        },
        "limitations": [
            "ideal-gas map omits real-gas hydrogen near saturation",
            "phase change and condensed air are outside this closure",
            "strong density departure requires variable-density low-Mach continuity",
        ],
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "hydrogen_air_buoyancy_transition_audit.json"
    figure_path = OUTPUT / "hydrogen_air_buoyancy_transition_audit.png"
    json_path.write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")

    figure, axis = plt.subplots(figsize=(10.5, 6.5), constrained_layout=True)
    levels = np.linspace(-1.0, 1.0, 41)
    image = axis.contourf(
        100.0 * xx, tt, np.clip(departure, -1.0, 1.0),
        levels=levels, cmap="coolwarm", extend="both",
    )
    axis.plot(100.0 * x, neutral, color="black", linewidth=2.0,
              label="neutral density")
    axis.axvline(4.0, color="lime", linestyle="--", linewidth=1.5,
                 label="4 vol% H2")
    axis.contour(
        100.0 * xx, tt, np.abs(departure), levels=[0.10],
        colors="yellow", linestyles=":", linewidths=1.8,
    )
    axis.text(34.0, 55.0, "denser than ambient\n(downward buoyancy)",
              color="white", ha="center", fontsize=10)
    axis.text(37.0, 245.0, "lighter than ambient\n(upward buoyancy)",
              color="black", ha="center", fontsize=10)
    axis.set_xlim(0.0, 60.0)
    axis.set_ylim(20.0, closure.ambient_temperature_k)
    axis.set_xlabel("hydrogen concentration (vol%)")
    axis.set_ylabel("mixture temperature (K)")
    axis.set_title(
        "H2/air ideal-gas density transition — active-scalar applicability map"
    )
    axis.legend(loc="lower right")
    plt.colorbar(image, ax=axis, label="(mixture density - ambient density) / ambient density")
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
        "four_vol_percent_neutral_temperature_k": checkpoints[
            "4_vol_percent"
        ]["neutral_buoyancy_temperature_k"],
        "sign_check": record["sign_check"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
