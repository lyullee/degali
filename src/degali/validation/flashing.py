"""Equivalent source terms for flashing pressurised releases.

A liquefied gas stored under its own vapour pressure and released to
atmosphere is outside DEGADIS entirely.  The liquid flashes at the orifice,
throws an aerosol, and entrains air while the droplets evaporate.  DEGADIS's
jet model handles a single-phase buoyant jet; its ground model starts from a
pool.  Neither is a flashing jet.

The accepted treatment is to hand the dispersion model the plume state at the
point where no liquid remains, and let it disperse from there.  Some trials
come with such a source term published -- the SMEDIS exercise issued them for
Desert Tortoise and two Lathen trials -- but most do not, and this module
computes one.

The calculation
---------------
Two conditions close the problem.  Writing :math:`r` for the mass of air
entrained per unit mass of contaminant, an adiabatic energy balance between
the stored liquid and the diluted plume gives

.. math::
    h_{store} + r\\, h_{air}(T_a) = h_{vap}(T) + r\\, h_{air}(T)

and the liquid has *just* disappeared when the vapour is exactly saturated at
its own partial pressure,

.. math:: y\\, p_{amb} = p_{sat}(T)

with :math:`y` the contaminant mole fraction implied by :math:`r`.  Two
equations, two unknowns.  Note the resulting temperature lies well below the
substance's normal boiling point -- ammonia comes out near 205 K against a
boiling point of 240 K -- because the plume is dilute and the partial pressure
is a fraction of ambient.

For liquid hydrogen that endpoint is near 20 K, where entrained nitrogen and
oxygen cannot remain gaseous.  The optional phase-safe boundary continues the
same adiabatic plug-flow balance until N2, O2 and Ar are all below saturation
at their component partial pressures.  It is a model-form sensitivity bound,
not the default: validation shows that replacing the whole condensed-air zone
by one warm all-gas source worsens residual variance and plume-centre height.

What it assumes
---------------
* All the flashed liquid evaporates, none rains out.  For a fine aerosol from
  a small orifice that is reasonable; for a large low-momentum release it is
  not, and rainout would leave a pool the model would also have to handle.
* Entrainment is adiabatic: no heat from the ground over the flashing length.
* The plume is at ambient pressure once expanded.

Validation
----------
Run against the two Desert Tortoise trials for which SMEDIS published a source
term, the calculation gives 204.9 K and 13.6 mole per cent against their
205 K and 13 per cent.  That agreement is what licenses applying it to trials
where no source term was issued.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq

from ..core.constants import WMA


# The three constituents that make up more than 99.9 % of dry air.  The
# fractions are the standard dry-air mole fractions, normalised here because
# the omitted CO2 and trace gases otherwise leave their sum just below one.
# Keeping the components separate is essential at cryogenic temperature: it
# is each component's *partial* pressure, not one atmosphere, that determines
# whether it can remain gaseous.
_BULK_DRY_AIR = (
    ("Nitrogen", 0.78084, 28.0134),
    ("Oxygen", 0.20946, 31.9988),
    ("Argon", 0.00934, 39.948),
)
_AIR_MOLE_SUM = sum(x for _name, x, _mw in _BULK_DRY_AIR)
_BULK_DRY_AIR = tuple(
    (name, x / _AIR_MOLE_SUM, mw) for name, x, mw in _BULK_DRY_AIR
)
_BULK_DRY_AIR_MW = sum(x * mw for _name, x, mw in _BULK_DRY_AIR)
_BULK_DRY_AIR_MASS = tuple(
    (name, x, mw, x * mw / _BULK_DRY_AIR_MW)
    for name, x, mw in _BULK_DRY_AIR
)

#: Optional extra floor on the search for the equivalent-source temperature.
#: The physical bound is the fluid's own triple point -- below it there is no
#: saturation line to sit on, and the assumption that the aerosol evaporated
#: rather than froze stops holding -- so that is what governs by default.
#:
#: A fixed floor cannot serve both: an equivalent source for ammonia lands
#: near 205 K, one for hydrogen near 20 K, and any value that excludes
#: nonsense for the first excludes the answer for the second.
MINIMUM_TEMPERATURE = 0.0


@dataclass(frozen=True)
class FlashResult:
    """Equivalent plume handoff after H2 evaporation or phase-safe heating."""

    temperature: float  #: K
    mole_fraction: float  #: contaminant
    air_ratio: float  #: kg air entrained per kg contaminant
    mass_fraction: float  #: contaminant, kg/kg
    flash_fraction: float  #: fraction vaporised at the orifice
    density: float | None = None  #: plume density at the equivalent source
    #: Density of the flashed contaminant *at the orifice plane*, before any
    #: air has been entrained, kg/m**3.  This is the state a jet model should
    #: start from when it is going to do the entrainment itself; the
    #: equivalent-source density is further downstream, where the plume is
    #: already much wider than the orifice.
    orifice_density: float | None = None
    #: True when the equivalent source was advanced beyond hydrogen
    #: evaporation to the first state where bulk dry air can all be gaseous.
    bulk_air_phase_safe: bool = False
    #: Bulk-air component on its phase boundary at the phase-safe handoff.
    limiting_air_species: str | None = None
    #: Largest component partial-pressure / saturation-pressure ratio.
    max_air_saturation_ratio: float | None = None

    @property
    def molar_percent(self) -> float:
        return self.mole_fraction * 100.0


def _vapour_enthalpy(fluid: str, temp: float, pressure: float, props) -> float:
    """Contaminant enthalpy as vapour, avoiding the saturation boundary.

    At the ambient pressure the contaminant may be exactly at saturation, and
    CoolProp declines a single-phase query there rather than picking a side.
    Above the critical temperature there is no saturation line at all and the
    plain pressure-temperature call is the right one.
    """
    # Only the saturation lookup may legitimately fail -- at or above the
    # critical point there is no saturation line -- and the fallback is the
    # plain pressure-temperature call. Catching anything wider would hide a
    # bad fluid name or a bad temperature behind a number that looks fine.
    try:
        critical = float(props("Tcrit", fluid))
    except Exception as exc:  # an unknown fluid, not a saturation edge
        raise ValueError(f"no critical temperature for {fluid!r}") from exc

    if temp < critical:
        try:
            p_sat = props("P", "T", temp, "Q", 1, fluid)
        except Exception:
            p_sat = None  # below the triple point: no saturation line either
        if p_sat is not None:
            return float(
                props("H", "T", temp, "P", min(pressure, 0.999 * p_sat), fluid)
            )
    return float(props("H", "T", temp, "P", pressure, fluid))


def _gas_enthalpy(fluid: str, temp: float, pressure: float, props) -> float:
    """Mass enthalpy on the gas branch at a component partial pressure."""
    return float(
        props("H", "T|gas", temp, "P", max(float(pressure), 1.0), fluid)
    )


def _bulk_air_phase_safe_endpoint(
    *, h_store: float, ambient_temperature: float, ambient_pressure: float,
    molecular_weight: float, props,
) -> tuple[float, float, float, str, float]:
    """First adiabatic state where N2, O2 and Ar can all be gaseous.

    This collapses the thermodynamically inaccessible part of a cryogenic
    hydrogen jet into a plug-flow control volume.  At a candidate
    temperature, the entrained dry-air mass follows from component enthalpy
    conservation.  Its composition then supplies the component partial
    pressures.  The exit is the first temperature for which none exceeds its
    saturation pressure.

    Temporary condensate inside the skipped zone needs no solid-property
    model because every bulk-air component is gaseous again at the handoff;
    enthalpy is a state function.  The corresponding momentum closure is
    applied by ``JetPlume.expanded_source_start`` rather than here.
    """
    h_air_in = sum(
        w * _gas_enthalpy(name, ambient_temperature, x * ambient_pressure, props)
        for name, x, _mw, w in _BULK_DRY_AIR_MASS
    )

    def composition(air_ratio: float):
        contaminant_moles = 1.0 / molecular_weight
        air_moles = air_ratio / _BULK_DRY_AIR_MW
        total_moles = contaminant_moles + air_moles
        y_hydrogen = contaminant_moles / total_moles
        y_air = air_moles / total_moles
        return y_hydrogen, tuple(
            (name, x * y_air) for name, x, _mw, _w in _BULK_DRY_AIR_MASS
        )

    def energy_residual(air_ratio: float, temp: float) -> float:
        y_hydrogen, y_components = composition(air_ratio)
        h_hydrogen = _gas_enthalpy(
            "Hydrogen", temp, y_hydrogen * ambient_pressure, props
        )
        partials = dict(y_components)
        h_air_out = sum(
            w * _gas_enthalpy(
                name, temp, partials[name] * ambient_pressure, props
            )
            for name, _x, _mw, w in _BULK_DRY_AIR_MASS
        )
        return h_store + air_ratio * h_air_in - (
            h_hydrogen + air_ratio * h_air_out
        )

    def state(temp: float):
        low, high = 1.0e-10, 1.0
        f_low = energy_residual(low, temp)
        f_high = energy_residual(high, temp)
        while f_low * f_high > 0.0 and high < 1.0e4:
            high *= 2.0
            f_high = energy_residual(high, temp)
        if f_low * f_high > 0.0:
            raise ValueError(
                f"no positive bulk-air energy balance at {temp:.2f} K"
            )
        ratio = brentq(
            lambda r: energy_residual(r, temp), low, high, xtol=1.0e-11
        )
        y_hydrogen, y_components = composition(ratio)
        saturation = {}
        for name, y_component in y_components:
            p_sat = float(props("P", "T", temp, "Q", 1, name))
            saturation[name] = y_component * ambient_pressure / p_sat
        limiting = max(saturation, key=saturation.get)
        return ratio, y_hydrogen, limiting, saturation[limiting]

    # The H2 is unambiguously gaseous above its critical point.  Start just
    # above nitrogen's triple point as well: below it a deliberately
    # supersaturated gas-state enthalpy query is not numerically defined by
    # every CoolProp backend.  Ordinary LH2 releases cross the phase-safe
    # boundary higher than this, near 68 K.
    lower = max(
        float(props("Tcrit", "Hydrogen")) + 0.1,
        float(props("Ttriple", "Nitrogen")) + 0.1,
    )
    upper = min(ambient_temperature - 1.0, 120.0)
    if lower >= upper:
        raise ValueError("ambient temperature is too low for a phase-safe source")

    def dew_residual(temp: float) -> float:
        return state(temp)[3] - 1.0

    try:
        temp = brentq(dew_residual, lower, upper, xtol=1.0e-9)
    except ValueError as exc:
        raise ValueError(
            "no bulk-air phase boundary between "
            f"{lower:.1f} and {upper:.1f} K"
        ) from exc
    ratio, y_hydrogen, limiting, saturation = state(temp)
    return float(temp), float(ratio), float(y_hydrogen), limiting, float(saturation)


def equivalent_source(
    fluid: str,
    *,
    storage_temperature: float,
    ambient_temperature: float,
    ambient_pressure: float,
    molecular_weight: float,
    minimum_temperature: float = MINIMUM_TEMPERATURE,
    bulk_air_phase_safe: bool = False,
    storage_pressure: float | None = None,
) -> FlashResult:
    """Plume state where a flashing jet's liquid has fully evaporated.

    Parameters
    ----------
    fluid
        CoolProp fluid name.
    storage_temperature, ambient_temperature
        K.  By default the release is assumed to be saturated liquid at the
        storage temperature, which is what a vessel under its own vapour
        pressure holds.
    storage_pressure
        Optional independently measured absolute pressure, Pa.  With this
        supplied, the stored enthalpy is evaluated at ``T,P`` instead of on
        the saturation curve.  This represents pressure-driven subcooled
        liquid without inventing a saturation temperature.
    ambient_pressure
        Pa.
    molecular_weight
        kg/kmol.
    bulk_air_phase_safe
        For hydrogen, advance the plug-flow source beyond completion of H2
        evaporation to the first adiabatic state where the N2, O2 and Ar in
        entrained dry air can all remain gaseous at their mixture partial
        pressures.  This removes the invalid 20 K ideal-air segment without
        introducing a fitted handoff temperature.  Off by default so the
        validated ammonia route and historical hydrogen results are unchanged.

    Raises
    ------
    ValueError
        If no solution exists between ``minimum_temperature`` and ambient --
        which happens when the substance is not actually a liquefied gas at
        these conditions, and the caller should not be using this route.
    """
    from CoolProp.CoolProp import PropsSI as props

    if storage_pressure is None:
        h_store = float(props("H", "T", storage_temperature, "Q", 0, fluid))
    else:
        storage_pressure = float(storage_pressure)
        if storage_pressure <= 0.0:
            raise ValueError("storage pressure must be positive")
        h_store = float(props(
            "H", "T", storage_temperature, "P", storage_pressure, fluid
        ))
    limiting_air_species = None
    max_air_saturation_ratio = None

    if bulk_air_phase_safe:
        if fluid.lower() != "hydrogen":
            raise ValueError(
                "bulk-air phase-safe source is implemented for Hydrogen only"
            )
        temp, ratio, y, limiting_air_species, max_air_saturation_ratio = (
            _bulk_air_phase_safe_endpoint(
                h_store=h_store,
                ambient_temperature=ambient_temperature,
                ambient_pressure=ambient_pressure,
                molecular_weight=molecular_weight,
                props=props,
            )
        )
    else:
        triple = float(props("Ttriple", fluid))
        critical = float(props("Tcrit", fluid))
        floor = max(minimum_temperature, triple + 0.5)
        # The saturation condition only has meaning below the critical point.
        # For ammonia the critical temperature is 405 K and ambient is nowhere
        # near it, so the whole range is available. For hydrogen it is 33 K,
        # far below ambient: above that the fluid is supercritical, there is no
        # saturation line to sit on, and any remaining liquid has necessarily
        # gone. So the search is capped there.
        ceiling = min(ambient_temperature - 1.0, critical - 0.05)
        if floor >= ceiling:
            raise ValueError(
                f"{fluid} has no sub-critical liquid range above {floor:.1f} K"
            )
        cp_air = float(
            props("C", "T", ambient_temperature, "P", ambient_pressure, "Air")
        )

        def state(temp: float) -> tuple[float, float]:
            h_c = _vapour_enthalpy(fluid, temp, ambient_pressure, props)
            ratio = (h_store - h_c) / (cp_air * (temp - ambient_temperature))
            if ratio <= 0.0:
                return ratio, -1.0
            moles = 1.0 / molecular_weight
            return ratio, moles / (moles + ratio / WMA)

        def residual(temp: float) -> float:
            ratio, y = state(temp)
            if ratio <= 0.0:
                return 1.0e6
            return y * ambient_pressure - float(
                props("P", "T", temp, "Q", 1, fluid)
            )

        try:
            temp = brentq(residual, floor, ceiling, xtol=1e-9)
        except ValueError as exc:
            raise ValueError(
                f"no equivalent source for {fluid} between {floor:.1f} and "
                f"{ceiling:.1f} K; is it a liquefied gas at these conditions?"
            ) from exc

        ratio, y = state(temp)
    wc = 1.0 / (1.0 + ratio)

    # flash fraction at the orifice, for reporting: isenthalpic expansion to
    # ambient pressure
    h_liquid = float(props("H", "P", ambient_pressure, "Q", 0, fluid))
    h_vapour = float(props("H", "P", ambient_pressure, "Q", 1, fluid))
    flash = (h_store - h_liquid) / (h_vapour - h_liquid)

    air_mw = _BULK_DRY_AIR_MW if bulk_air_phase_safe else WMA
    wm = y * molecular_weight + (1.0 - y) * air_mw
    density = ambient_pressure * wm / 8314.462618 / temp

    # homogeneous two-phase density at the orifice: the flashed vapour and the
    # liquid it left behind, at ambient pressure and no air yet
    rho_liquid = float(props("D", "P", ambient_pressure, "Q", 0, fluid))
    rho_vapour = float(props("D", "P", ambient_pressure, "Q", 1, fluid))
    x = min(max(flash, 0.0), 1.0)
    orifice_density = 1.0 / ((1.0 - x) / rho_liquid + x / rho_vapour)

    return FlashResult(
        temperature=float(temp),
        mole_fraction=float(y),
        air_ratio=float(ratio),
        mass_fraction=float(wc),
        flash_fraction=float(min(max(flash, 0.0), 1.0)),
        density=float(density),
        orifice_density=float(orifice_density),
        bulk_air_phase_safe=bulk_air_phase_safe,
        limiting_air_species=limiting_air_species,
        max_air_saturation_ratio=max_air_saturation_ratio,
    )


def plume_half_width(
    rate: float, result: FlashResult, velocity: float
) -> float:
    """Half-width of the plume at the equivalent source, m.

    From continuity: the total mass flow -- contaminant plus entrained air --
    crosses a circular section at the plume velocity.  ``velocity`` has to
    come from somewhere else; a jet decelerating into a crossflow ends up at
    something close to the wind speed, so that is a usable estimate when no
    better one exists.
    """
    total = rate * (1.0 + result.air_ratio)
    area = total / result.density / velocity
    return math.sqrt(area / math.pi)
