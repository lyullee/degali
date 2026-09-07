"""Turning a REDIPHEM trial into a DEGADIS case.

A model-evaluation study lives or dies on how its input decks were built, and
the usual failure is silent: someone transcribes conditions from a report into
a deck, makes a judgement call at each ambiguity, and the judgement calls
never appear in the results.  This module makes the mapping explicit and
executable, and records every assumption it had to make on the case it
returns.

What the trial gives directly
-----------------------------
Wind speed and its reference height, friction velocity, surface roughness,
Monin-Obukhov length, stability class, ambient temperature, pressure and
humidity, release rate and duration, and the pool diameter or nozzle
diameter.  These go straight across.

What has to be supplied
-----------------------
``SPECS.DAT`` describes an experiment, not a model run, so a few things are
not in it:

* **Contaminant properties.**  The trial names a substance; its molecular
  weight, boiling point and vapour density come from a small table here.
* **Surface temperature.**  Rarely recorded.  Defaults to ambient, which
  understates ground heating over a desert playa in August.
* **Levels of concern.**  A property of the question being asked, not of the
  experiment.  Defaults to the flammable limits for fuels and to a nominal
  pair otherwise.
* **Release temperature.**  Given for pressurised releases; for a pool spill
  it is the substance's normal boiling point.

Every one of these is recorded in :attr:`TrialCase.assumptions`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core.atmosphere import absolute_humidity, stability_defaults
from ..core.constants import ATM_TO_PA, PI, RGAS
from ..core.thermo import (
    AmbientConditions,
    GasProperties,
    LegacyBackend,
    Thermo,
)
from ..io.inp import Case, SourceTable
from .rediphem import Trial


@dataclass(frozen=True)
class Substance:
    """Properties a trial specification does not carry."""

    label: str  #: three-character DEGADIS name
    mw: float  #: kg/kmol
    boiling_point: float  #: K at one atmosphere
    ulc: float  #: upper level of concern, mole fraction
    llc: float  #: lower level of concern, mole fraction
    #: Molar heat-capacity constants for the ``CPC`` correlation,
    #: ``cp_molar = 3.33e4 + q1 p1 T**(p1-1)``.  Zero ``cpk`` means a
    #: constant heat capacity, which the deck reader expands.
    cpk: float = 0.0
    cpp: float = 1.0

    @property
    def flammable(self) -> bool:
        return self.label in ("LNG", "PRO", "ETH")


#: Only the substances the field trials actually used.  The flammable limits
#: are the standard ones; for the toxic gases the pair brackets the range the
#: sensors resolved, since a "level of concern" for ammonia depends on which
#: exposure standard is being applied.
SUBSTANCES: dict[str, Substance] = {
    "lng": Substance("LNG", 16.04, 111.7, ulc=0.15, llc=0.05,
                     cpk=5.6e-8, cpp=5.0),
    "methane": Substance("LNG", 16.04, 111.7, ulc=0.15, llc=0.05,
                         cpk=5.6e-8, cpp=5.0),
    "ammonia": Substance("NH3", 17.03, 239.7, ulc=0.05, llc=0.001,
                         cpk=3845.0, cpp=1.0),
    "propane": Substance("PRO", 44.10, 231.1, ulc=0.095, llc=0.021),
    "n2o4": Substance("N2O", 92.01, 294.3, ulc=0.01, llc=0.0001),
    "r12": Substance("R12", 120.91, 243.4, ulc=0.10, llc=0.005),
    "sf6": Substance("SF6", 146.06, 209.3, ulc=0.10, llc=0.001),
}


@dataclass
class TrialCase:
    """A DEGADIS case built from a trial, with its provenance."""

    trial: Trial
    case: Case
    substance: Substance
    #: Things not in ``SPECS.DAT`` that had to be supplied, and their values.
    assumptions: dict[str, str] = field(default_factory=dict)
    #: Reasons the trial cannot be modelled, if any.
    problems: list[str] = field(default_factory=list)

    @property
    def usable(self) -> bool:
        return not self.problems


def _pool_radius(trial: Trial, substance: Substance, rate: float) -> tuple[float, str]:
    """Source radius, m, and where it came from.

    A recorded pool diameter is used when present.  Otherwise the pool is
    assumed to have spread until evaporation balanced the spill rate, at the
    3e-4 kg/(m^2 s) typical of a cryogen boiling on soil -- the figure
    DEGADIS's own documentation uses for LNG on a dyke floor.
    """
    diameter = trial.value("pool diameter")
    if diameter and diameter > 0.0:
        return diameter / 2.0, f"pool diameter {diameter:.1f} m from the trial"
    flux = 3.0e-4
    radius = math.sqrt(rate / flux / PI)
    return radius, (
        f"pool radius {radius:.1f} m assumed from the spill rate at "
        f"{flux:g} kg/(m2 s) evaporation"
    )


def to_case(
    trial: Trial,
    *,
    averaging_time: float | None = None,
    surface_temperature: float | None = None,
    heat_transfer: int = 1,
    water_transfer: int = 1,
    reference_height: float | None = None,
) -> TrialCase:
    """Build a DEGADIS case from a REDIPHEM trial.

    Parameters
    ----------
    averaging_time
        Defaults to the release duration, which is what the trial's own peak
        concentrations correspond to.
    surface_temperature
        Defaults to ambient.  Over a desert playa in August the ground runs
        well above air temperature, and DEGADIS is sensitive to that through
        the ground heat flux, so this is worth setting when it is known.
    heat_transfer, water_transfer
        ``IHTFL`` and ``IWTFL``.  Both on by default: a cryogenic cloud over
        warm ground exchanges both.
    """
    assumptions: dict[str, str] = {}
    problems: list[str] = []

    key = trial.substance.lower()
    substance = SUBSTANCES.get(key)
    if substance is None:
        problems.append(f"no property data for {trial.substance!r}")
        substance = Substance("UNK", 30.0, 200.0, 0.1, 0.01)

    rate = trial.rate
    duration = trial.duration
    u0 = trial.value("site average windspeed")
    z0 = trial.wind_height
    if z0 is None and reference_height is not None:
        z0 = reference_height
        assumptions["reference height"] = (
            f"{z0:g} m, from the other trials in this series"
        )
    zr = trial.value("surface roughness")
    tamb_c = trial.value("ambient temperature")
    pamb_bar = trial.value("ambient pressure")
    relhum = trial.value("relative humidity")
    stability = trial.stability

    for name, value in (("release rate", rate), ("release duration", duration),
                        ("wind speed", u0), ("reference height", z0),
                        ("surface roughness", zr),
                        ("ambient temperature", tamb_c),
                        ("stability class", stability)):
        if value is None or (isinstance(value, float) and not math.isfinite(value)):
            problems.append(f"{name} missing")

    if problems:
        return TrialCase(trial=trial, case=None, substance=substance,
                         assumptions=assumptions, problems=problems)

    tamb = tamb_c + 273.15
    pamb = (pamb_bar or 1.013) * 1.0e5 / ATM_TO_PA
    # REDIPHEM records relative humidity as a fraction; a few trials use per cent
    relhum_pct = relhum * 100.0 if relhum is not None and relhum <= 1.0 else (relhum or 0.0)

    tsurf = surface_temperature if surface_temperature is not None else tamb
    if surface_temperature is None:
        assumptions["surface temperature"] = "ambient (not recorded)"

    avtime = averaging_time if averaging_time is not None else duration
    if averaging_time is None:
        assumptions["averaging time"] = f"{avtime:.0f} s, the release duration"

    rml = trial.value("Monin-Obukov length")
    indvel = 2 if rml is not None else 1
    if rml is None:
        rml = 0.0
        assumptions["Monin-Obukhov length"] = "from the stability class"

    # release temperature: measured for a jet, the boiling point for a pool
    exit_c = trial.value("exit temperature")
    if trial.release_type == "pool" or exit_c is None:
        gastem = substance.boiling_point
        assumptions["release temperature"] = (
            f"{gastem:.1f} K, the normal boiling point"
        )
    else:
        gastem = exit_c + 273.15

    gasrho = pamb * substance.mw / RGAS / gastem
    radius, provenance = _pool_radius(trial, substance, rate)
    if trial.release_type != "pool":
        assumptions["source"] = (
            "modelled as an evaporating pool; a pressurised jet should go "
            "through the jet/plume model instead"
        )
    assumptions["source radius"] = provenance
    assumptions["levels of concern"] = (
        f"{substance.ulc * 100:g} and {substance.llc * 100:g} mole %"
    )

    backend = LegacyBackend()
    humid, relhum_out = absolute_humidity(
        tamb, pamb, backend.water_vapour_pressure, relhum=relhum_pct
    )

    gas = GasProperties(
        name=substance.label, mw=substance.mw, temp=gastem, rho=gasrho,
        cpk=substance.cpk, cpp=substance.cpp,
        ulc=substance.ulc, llc=substance.llc, zzc=0.0,
    )
    ambient = AmbientConditions(
        tamb=tamb, pamb=pamb, humid=humid, tsurf=tsurf, isofl=0,
        ihtfl=heat_transfer, iwtfl=water_transfer, htco=0.0, wtco=0.0,
    )

    # the source table DEGINP would have written: constant rate for the
    # release, then two rows shutting it off
    times = np.array([0.0, duration, duration + 1.0, duration + 2.0])
    source = SourceTable(
        time=times,
        rate=np.array([rate, rate, 0.0, 0.0]),
        radius=np.array([radius, radius, 0.0, 0.0]),
        wc=np.ones(4),
        temp=np.full(4, gastem),
        fracv=np.ones(4),
        enthalpy=np.zeros(4),
        rho=np.zeros(4),
    )
    thermo = Thermo(gas=gas, ambient=ambient, backend=backend)
    for i in range(2):
        wa = (1.0 - source.wc[i]) / (1.0 + humid)
        h = thermo.enthalpy(source.wc[i], wa, source.temp[i])
        source.enthalpy[i] = h
        source.rho[i] = thermo.properties(
            source.wc[i], wa, h, ifl=-1, temp=source.temp[i]
        ).rho
    source.enthalpy[2:] = source.enthalpy[1]
    source.rho[2:] = source.rho[1]

    stab = stability_defaults(zr, stability, avtime)
    case = Case(
        titles=[f"{trial.series} {trial.name} -- {trial.substance}", "", "", ""],
        u0=u0, z0=z0, zr=zr, istab=stability, oodist=0.0, avtime=avtime,
        indvel=indvel, rml=rml if indvel == 2 else stab.rml,
        gas=gas, ambient=ambient, relhum=relhum_out,
        yclow=substance.llc, gmass0=0.0, source=source,
        steady_state=True, density_table=None, timestamp="",
        instantaneous=False, stability=stab,
    )
    return TrialCase(trial=trial, case=case, substance=substance,
                     assumptions=assumptions, problems=problems)


# ==========================================================================
# Pressurised flashing jets
# ==========================================================================


@dataclass(frozen=True)
class EquivalentSource:
    """A flashing jet reduced to a ground-level source.

    A pressurised release of a liquefied gas flashes at the orifice, throws
    an aerosol, and entrains air while the droplets evaporate.  None of that
    is in DEGADIS, whose jet model handles a single-phase buoyant jet and
    whose ground model starts from a pool.  The standard treatment -- and the
    one the SMEDIS evaluation used -- is to hand the model the plume state at
    the point where no liquid remains, and let it disperse from there.

    Attributes
    ----------
    distance
        Downstream distance at which the liquid has gone, m.  Becomes
        ``OODIST``, so every reported distance stays in the trial's own
        coordinate.
    velocity
        Plume velocity there, m/s.  Recorded for completeness; DEGADIS's
        ground-level source has no momentum, which is part of what this
        approximation gives up.
    molar_percent
        Contaminant mole per cent in the plume.
    temperature
        Plume temperature, K.
    half_width
        Half-width of the plume, m.
    """

    distance: float
    velocity: float
    molar_percent: float
    temperature: float
    half_width: float


#: Equivalent source terms from the SMEDIS exercise (batch 1, February 1998),
#: for the trials where they were issued.  Keys are ``(series, trial)`` in
#: REDIPHEM's naming.
SMEDIS_SOURCES: dict[tuple[str, str], EquivalentSource] = {
    ("TORTOISE", "D1"): EquivalentSource(51.0, 7.5, 13.0, 205.0, 6.40),
    ("TORTOISE", "D2"): EquivalentSource(48.3, 6.0, 13.0, 205.0, 8.40),
    ("LATHEN", "EEC36"): EquivalentSource(0.9, 22.8, 20.0, 196.0, 0.07),
    ("LATHEN", "EEC55"): EquivalentSource(5.0, 18.1, 20.0, 196.0, 0.39),
}


def computed_source(
    trial: Trial, *, fluid: str | None = None, velocity: float | None = None,
) -> EquivalentSource | None:
    """Compute an equivalent source for a flashing jet trial.

    Used when no source term was published for it.  ``velocity`` defaults to
    the site wind speed: a jet decelerating into a crossflow ends up close to
    it by the time the aerosol has evaporated, and it is the only estimate
    available from the trial specification alone.

    The distance at which the liquid disappears is not computed -- that needs
    a jet trajectory, which is the part being avoided -- so it is scaled from
    the Desert Tortoise trials, where SMEDIS placed it at about 55 nozzle
    diameters downstream.  A distance is needed only to offset the reported
    coordinates, so an error in it moves the whole profile rather than
    changing its shape.
    """
    from .flashing import equivalent_source, plume_half_width
    from ..core.fluids import resolve

    substance = SUBSTANCES.get(trial.substance.lower())
    if substance is None:
        return None
    name = fluid or resolve(substance.label, substance.mw).fluid
    if name is None:
        return None

    exit_c = trial.value("exit temperature")
    tamb_c = trial.value("ambient temperature")
    rate = trial.rate
    nozzle = trial.value("nozzle diameter")
    if exit_c is None or tamb_c is None or rate is None:
        return None
    pamb = (trial.value("ambient pressure") or 1.013) * 1.0e5
    u = velocity or trial.value("site average windspeed") or 5.0

    try:
        flash = equivalent_source(
            name,
            storage_temperature=exit_c + 273.15,
            ambient_temperature=tamb_c + 273.15,
            ambient_pressure=pamb,
            molecular_weight=substance.mw,
        )
    except (ValueError, ImportError):
        return None

    return EquivalentSource(
        distance=55.0 * (nozzle or 0.1),
        velocity=u,
        molar_percent=flash.molar_percent,
        temperature=flash.temperature,
        half_width=plume_half_width(rate, flash, u),
    )


def jet_to_case(
    trial: Trial,
    source: EquivalentSource,
    *,
    averaging_time: float | None = None,
    surface_temperature: float | None = None,
    heat_transfer: int = 1,
    water_transfer: int = 1,
    reference_height: float | None = None,
) -> TrialCase:
    """Build a DEGADIS case from a jet trial and its equivalent source.

    The plume is replaced by a circular ground-level source of the plume's
    half-width, releasing the trial's contaminant rate already diluted to the
    plume's composition and temperature -- the same substitution ``DEGBRIDG``
    makes for a jet that touches down, applied at the point the liquid
    disappears instead.

    What this gives up is the plume's momentum and its elevation: DEGADIS's
    ground source has neither. For Desert Tortoise the plume is 7.5 m/s at 51
    m, so the momentum is not negligible, and results should be read with
    that in mind.
    """
    assumptions: dict[str, str] = {}
    problems: list[str] = []

    key = trial.substance.lower()
    substance = SUBSTANCES.get(key)
    if substance is None:
        return TrialCase(trial=trial, case=None,
                         substance=Substance("UNK", 30.0, 200.0, 0.1, 0.01),
                         problems=[f"no property data for {trial.substance!r}"])

    rate = trial.rate
    duration = trial.duration
    u0 = trial.value("site average windspeed")
    z0 = trial.wind_height
    if z0 is None and reference_height is not None:
        z0 = reference_height
        assumptions["reference height"] = (
            f"{z0:g} m, from the other trials in this series"
        )
    zr = trial.value("surface roughness")
    tamb_c = trial.value("ambient temperature")
    stability = trial.stability
    for name, value in (("release rate", rate), ("release duration", duration),
                        ("wind speed", u0), ("reference height", z0),
                        ("surface roughness", zr),
                        ("ambient temperature", tamb_c),
                        ("stability class", stability)):
        if value is None:
            problems.append(f"{name} missing")
    if problems:
        return TrialCase(trial=trial, case=None, substance=substance,
                         problems=problems)

    tamb = tamb_c + 273.15
    pamb = (trial.value("ambient pressure") or 1.013) * 1.0e5 / ATM_TO_PA
    relhum = trial.value("relative humidity") or 0.0
    relhum_pct = relhum * 100.0 if relhum <= 1.0 else relhum
    tsurf = surface_temperature if surface_temperature is not None else tamb
    if surface_temperature is None:
        assumptions["surface temperature"] = "ambient (not recorded)"
    avtime = averaging_time if averaging_time is not None else duration
    rml = trial.value("Monin-Obukov length")
    indvel = 2 if rml is not None else 1

    backend = LegacyBackend()
    humid, relhum_out = absolute_humidity(
        tamb, pamb, backend.water_vapour_pressure, relhum=relhum_pct
    )

    # composition at the equivalent source: mole per cent to mass fraction
    yc = source.molar_percent / 100.0
    ya = (1.0 - yc) / (1.0 + humid * 28.96 / 18.02)
    yw = 1.0 - ya - yc
    wm = yc * substance.mw + ya * 28.96 + yw * 18.02
    wc = substance.mw / wm * yc

    gastem = substance.boiling_point
    gasrho = pamb * substance.mw / RGAS / gastem
    gas = GasProperties(
        name=substance.label, mw=substance.mw, temp=gastem, rho=gasrho,
        cpk=substance.cpk, cpp=substance.cpp,
        ulc=substance.ulc, llc=substance.llc, zzc=0.0,
    )
    ambient = AmbientConditions(
        tamb=tamb, pamb=pamb, humid=humid, tsurf=tsurf, isofl=0,
        ihtfl=heat_transfer, iwtfl=water_transfer,
    )

    times = np.array([0.0, duration, duration + 1.0, duration + 2.0])
    table = SourceTable(
        time=times,
        rate=np.array([rate, rate, 0.0, 0.0]),
        radius=np.full(4, source.half_width),
        wc=np.full(4, wc),
        temp=np.full(4, source.temperature),
        fracv=np.ones(4),
        enthalpy=np.zeros(4),
        rho=np.zeros(4),
    )
    table.radius[2:] = 0.0
    thermo = Thermo(gas=gas, ambient=ambient, backend=backend)
    for i in range(2):
        wa = (1.0 - table.wc[i]) / (1.0 + humid)
        h = thermo.enthalpy(table.wc[i], wa, table.temp[i])
        table.enthalpy[i] = h
        table.rho[i] = thermo.properties(
            table.wc[i], wa, h, ifl=-1, temp=table.temp[i]
        ).rho
    table.enthalpy[2:] = table.enthalpy[1]
    table.rho[2:] = table.rho[1]

    assumptions["source"] = (
        f"SMEDIS equivalent source: {source.molar_percent:g} mole % at "
        f"{source.temperature:g} K, half-width {source.half_width:g} m, "
        f"placed {source.distance:g} m downstream"
    )
    assumptions["momentum"] = (
        f"discarded; the plume is {source.velocity:g} m/s there and "
        "DEGADIS's ground source has no momentum"
    )
    assumptions["levels of concern"] = (
        f"{substance.ulc * 100:g} and {substance.llc * 100:g} mole %"
    )

    stab = stability_defaults(zr, stability, avtime)
    case = Case(
        titles=[f"{trial.series} {trial.name} -- {trial.substance} jet", "", "", ""],
        u0=u0, z0=z0, zr=zr, istab=stability,
        oodist=source.distance, avtime=avtime,
        indvel=indvel, rml=rml if indvel == 2 else stab.rml,
        gas=gas, ambient=ambient, relhum=relhum_out,
        yclow=substance.llc, gmass0=0.0, source=table,
        steady_state=True, density_table=None, timestamp="",
        instantaneous=False, stability=stab,
    )
    return TrialCase(trial=trial, case=case, substance=substance,
                     assumptions=assumptions, problems=problems)


def puff_to_case(
    trial: Trial,
    *,
    averaging_time: float | None = None,
    surface_temperature: float | None = None,
    heat_transfer: int = 0,
    water_transfer: int = 0,
) -> TrialCase:
    """Build a DEGADIS case from an instantaneous release.

    Thorney Island released a cylinder of Freon-air mixture by dropping its
    walls: no flux, just a volume of dense gas sitting on the ground.  DEGADIS
    handles that through ``GMASS0``, an initial mass over a source of zero
    rate, and the cloud then collapses under its own weight.

    That path exercises the van Ulden momentum balance, which the source model
    enables when the initial height-to-diameter ratio exceeds 0.1 -- a
    quasi-steady gravity current has not had time to form, so the frontal
    velocity comes from a momentum balance instead.  None of the EPA test
    cases reaches it.

    The release is a *mixture*, not pure contaminant: Thorney Island's
    cylinder held 20 to 32 mole per cent Freon in air, chosen to set the
    density.  ``initial concentration`` carries that.
    """
    assumptions: dict[str, str] = {}
    problems: list[str] = []

    substance = SUBSTANCES.get(trial.substance.lower())
    if substance is None:
        return TrialCase(trial=trial, case=None,
                         substance=Substance("UNK", 30.0, 200.0, 0.1, 0.01),
                         problems=[f"no property data for {trial.substance!r}"])

    mass = trial.rate  # for a puff the "rate" field holds the released mass
    diameter = trial.value("nozzle diameter")  # the container, for a puff
    yc0 = trial.value("initial concentration")
    u0 = trial.value("site average windspeed")
    z0 = trial.wind_height
    zr = trial.value("surface roughness")
    tamb_c = trial.value("ambient temperature")
    exit_c = trial.value("exit temperature")
    stability = trial.stability
    for name, value in (("released mass", mass), ("container diameter", diameter),
                        ("initial concentration", yc0), ("wind speed", u0),
                        ("reference height", z0), ("surface roughness", zr),
                        ("ambient temperature", tamb_c),
                        ("stability class", stability)):
        if value is None:
            problems.append(f"{name} missing")
    if problems:
        return TrialCase(trial=trial, case=None, substance=substance,
                         problems=problems)

    tamb = tamb_c + 273.15
    pamb = (trial.value("ambient pressure") or 1.013) * 1.0e5 / ATM_TO_PA
    relhum = trial.value("relative humidity") or 0.0
    relhum_pct = relhum * 100.0 if relhum <= 1.0 else relhum
    tsurf = surface_temperature if surface_temperature is not None else tamb
    gastem = (exit_c + 273.15) if exit_c is not None else tamb

    backend = LegacyBackend()
    humid, relhum_out = absolute_humidity(
        tamb, pamb, backend.water_vapour_pressure, relhum=relhum_pct
    )

    # the cylinder held a mixture; convert its mole fraction to mass fraction
    wm = yc0 * substance.mw + (1.0 - yc0) * 28.96
    wc = substance.mw / wm * yc0
    contaminant_mass = mass * wc

    gasrho = pamb * substance.mw / RGAS / gastem
    gas = GasProperties(
        name=substance.label, mw=substance.mw, temp=gastem, rho=gasrho,
        cpk=substance.cpk, cpp=substance.cpp,
        ulc=substance.ulc, llc=substance.llc, zzc=0.0,
    )
    ambient = AmbientConditions(
        tamb=tamb, pamb=pamb, humid=humid, tsurf=tsurf, isofl=0,
        ihtfl=heat_transfer, iwtfl=water_transfer,
    )

    radius = diameter / 2.0
    avtime = averaging_time if averaging_time is not None else 18.4
    # zero rate with a finite initial mass is what makes this instantaneous
    times = np.array([0.0, 1.0, 2.0, 3.0])
    table = SourceTable(
        time=times, rate=np.zeros(4),
        radius=np.array([radius, radius, 0.0, 0.0]),
        wc=np.full(4, wc), temp=np.full(4, gastem), fracv=np.ones(4),
        enthalpy=np.zeros(4), rho=np.zeros(4),
    )
    thermo = Thermo(gas=gas, ambient=ambient, backend=backend)
    for i in range(2):
        wa = (1.0 - table.wc[i]) / (1.0 + humid)
        h = thermo.enthalpy(table.wc[i], wa, table.temp[i])
        table.enthalpy[i] = h
        table.rho[i] = thermo.properties(
            table.wc[i], wa, h, ifl=-1, temp=table.temp[i]
        ).rho
    table.enthalpy[2:] = table.enthalpy[1]
    table.rho[2:] = table.rho[1]

    assumptions["source"] = (
        f"instantaneous: {contaminant_mass:.0f} kg of contaminant in "
        f"{mass:.0f} kg of mixture at {yc0 * 100:g} mole %, over a "
        f"{diameter:g} m container"
    )
    assumptions["heat transfer"] = (
        "off; the release is near ambient temperature so the ground exchanges "
        "little with it"
    )
    if surface_temperature is None:
        assumptions["surface temperature"] = "ambient (not recorded)"

    stab = stability_defaults(zr, stability, avtime)
    rml = trial.value("Monin-Obukov length")
    indvel = 2 if rml is not None else 1
    case = Case(
        titles=[f"{trial.series} {trial.name} -- {trial.substance} puff", "", "", ""],
        u0=u0, z0=z0, zr=zr, istab=stability, oodist=0.0, avtime=avtime,
        indvel=indvel, rml=rml if indvel == 2 else stab.rml,
        gas=gas, ambient=ambient, relhum=relhum_out,
        yclow=substance.llc, gmass0=contaminant_mass, source=table,
        steady_state=False, density_table=None, timestamp="",
        instantaneous=True, stability=stab,
    )
    return TrialCase(trial=trial, case=case, substance=substance,
                     assumptions=assumptions, problems=problems)
