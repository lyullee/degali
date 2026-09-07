"""The liquid hydrogen near-field comparison, computed rather than quoted.

An audit of this package found that **not one of its statistical results was
computed by the suite**. The parity claims are enforced to 1e-12 against the
original Fortran; ``MG 0.738``, its confidence interval, ``VG 1.41``, ``FAC2
0.83``, the former 0.64-to-0.94 vertical-width ratio and the 1.07-to-0.19 m
trajectory correction were all typed into documents and into an f-string. They
were produced by a working session and never pinned. A number nothing
regenerates is a number that drifts, and two of them had: the module that
prints results to a user carried a statistic superseded when the far-field arcs
were recovered.

This module closes that. It runs the model over the campaign, pairs it with the
measurements, and returns the statistics. The same call produces the as-shipped
and the corrected figures, so the effect of the five adopted corrections is a
difference between two computed numbers rather than two remembered ones.

Two paths in
------------

``from_workbooks`` reads the PRESLHY E3.5 dataset directly and is the
authoritative route. The dataset is not redistributable, so it needs
``DEGALI_E35_ROOT`` and ``DEGALI_E35_REPORT``.

``from_reduced`` reads a reduced table that ``reduce`` writes from the
workbooks: one row per sensor per trial, carrying the position, the peak and
the mean. That table is a derived product of this work rather than the dataset,
so it can travel with the repository and the statistic becomes reproducible on
a clone with no third-party data — which is the only way a reviewer can check
it.

The reduction is lossless for this comparison and lossy for every other: it
keeps what the concentration comparison uses and drops the time series. It is
not a substitute for the dataset, and ``reduce`` records which workbook and
which window each row came from.

What is compared
----------------

The arc maximum, one value per downwind distance per trial: the largest reading
across the sensors at that distance. That is what a steady model's centreline
corresponds to, and taking it is not a rejection — the reduction from every
reading to one per arc is stated in ``docs/lh2-results.md``.

**Against the model at the sensor, not at its own centreline.** This matters
more than it looks. By 6 m the modelled plume centre for a 0.5 m release sits
at 1.35 m as shipped, above the topmost sensor at 0.75 m, so the model's
centreline is a concentration the array could not have measured. Comparing
against it scores the model on a quantity that has no measurement, and the
difference is a factor of two on that arc. The prediction is therefore
evaluated at each sensor's own position and then maximised over the arc, the
same operation applied to the measurement.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..core.constants import VKC

#: Trials whose exit velocity is below this multiple of the wind speed are
#: steered by the wind rather than by their own momentum. Fixed in advance from
#: the exit-velocity ratio, not by looking at agreement.
MOMENTUM_RATIO = 10.0

# The independent-energy mechanism campaign was pre-registered on this
# population.  Alternative source-rate specifications must not admit a new
# trial, because that would mix a model change with a sample change.
INDEPENDENT_ENERGY_TRIALS = (10, 11, 12, 22, 23, 24, 25)

#: Integration step and arclength limit. ``distmx`` is the **step**, not the
#: limit -- a rebuild of this configuration passed 40.0 there, integrating the
#: whole plume in one stride, and the resulting trajectory was wrong by a
#: factor of two while still looking plausible. ``smax`` is the limit.
STEP = 0.2
REACH = 40.0

#: Columns of the reduced table, in order.
REDUCED_COLUMNS = (
    "trial", "serial", "x", "y", "z", "z_axis", "peak", "mean", "samples",
    "workbook", "window_start", "window_end",
)


@dataclass
class Pair:
    """One arc: the measured maximum and the model's, at the same positions."""

    trial: int
    x: float
    observed: float
    predicted: float
    #: Height of the sensor that carried the measured maximum, m above ground.
    observed_at: float
    #: Height of the sensor that carried the modelled maximum.
    predicted_at: float
    #: Where the model puts the plume centre at this distance.
    plume_centre: float

    @property
    def visible(self) -> bool:
        """Whether the modelled plume centre is inside the array at this arc.

        When it is not, the model's own centreline is a concentration the
        measurement could not have seen, and only the sensor-evaluated
        comparison means anything.
        """
        return self.plume_centre <= self.predicted_at + 1e-9


@dataclass
class NearField:
    """The comparison, and everything needed to audit it."""

    pairs: list[Pair] = field(default_factory=list)
    #: Trials excluded, and why. Every exclusion is recorded: a statistic a
    #: reader cannot trace back to the published dataset is not evidence.
    excluded: dict[int, str] = field(default_factory=dict)
    corrections: bool = False

    @property
    def trials(self) -> list[int]:
        return sorted({p.trial for p in self.pairs})

    def statistics(self, **kw):
        from .statistics import statistics

        return statistics(
            [p.observed for p in self.pairs],
            [p.predicted for p in self.pairs],
            **kw,
        )

    def interval(self, *, draws: int = 2000, seed: int = 0) -> tuple:
        """A 95 % interval on MG, resampling **trials** rather than points.

        Readings within one trial share a release, a wind and a source
        estimate. Resampling points instead narrows the interval by roughly
        the square root of the sensor count and overstates the evidence.
        """
        rng = np.random.default_rng(seed)
        by_trial: dict[int, list[Pair]] = {}
        for p in self.pairs:
            by_trial.setdefault(p.trial, []).append(p)
        keys = list(by_trial)
        out = []
        for _ in range(draws):
            drawn = [
                p
                for k in rng.choice(keys, size=len(keys), replace=True)
                for p in by_trial[int(k)]
            ]
            ratios = [
                math.log(p.observed / p.predicted)
                for p in drawn
                if p.observed > 0.0 and p.predicted > 0.0
            ]
            if ratios:
                out.append(math.exp(float(np.mean(ratios))))
        lo, hi = np.percentile(out, [2.5, 97.5])
        return float(lo), float(hi)

    def report(self) -> str:
        s = self.statistics()
        lo, hi = self.interval()
        state = "corrections on" if self.corrections else "as shipped"
        blind = sum(1 for p in self.pairs if not p.visible)
        return "\n".join([
            f"LH2 near field, {state}",
            f"  arcs                : {len(self.pairs)} "
            f"from {len(self.trials)} trials",
            f"  MG                  : {s.mg:.3f}  95 % CI [{lo:.3f}, {hi:.3f}]",
            f"  VG                  : {s.vg:.3f}",
            f"  FAC2                : {s.fac2:.3f}",
            f"  arcs where the model centre is above the array : {blind}",
            f"  excluded            : {len(self.excluded)} trials",
        ])


@dataclass
class CoupledFieldValidation:
    """Pre-registered PRESLHY comparison for the conserved coupled path."""

    baseline: NearField
    candidate: NearField
    baseline_vertical: list[dict]
    candidate_vertical: list[dict]
    selected_trials: list[int]
    nearfield_establishment: str
    crosswind_entrainment: str
    compatible_handoff_search: bool
    geometry_observation: str = "state"
    handoffs: dict[int, object] = field(default_factory=dict, repr=False)
    failures: dict[int, str] = field(default_factory=dict)

    @property
    def all_interfaces_accepted(self) -> bool:
        return (
            set(self.handoffs) == set(self.selected_trials)
            and not self.failures
            and all(result.accepted for result in self.handoffs.values())
        )

    @staticmethod
    def geometry(rows: list[dict]) -> dict[str, float]:
        """Summarise Gaussian width and centre-height error."""
        if not rows:
            return {
                "width_ratio": math.nan,
                "centre_mean_error": math.nan,
                "centre_mae": math.nan,
            }
        width = [
            row["modelled_sigma_z"] / row["measured_sigma_z"]
            for row in rows if row["measured_sigma_z"] > 0.0
        ]
        centre = [
            row["modelled_centre"] - row["measured_centre"] for row in rows
        ]
        return {
            "width_ratio": float(np.mean(width)),
            "centre_mean_error": float(np.mean(centre)),
            "centre_mae": float(np.mean(np.abs(centre))),
        }

    @property
    def promoted(self) -> bool:
        """Whether every frozen promotion condition is satisfied."""
        if (
            not self.all_interfaces_accepted
            or not self.candidate.pairs
            or not self.candidate_vertical
        ):
            return False
        old = self.baseline.statistics()
        new = self.candidate.statistics()
        old_geometry = self.geometry(self.baseline_vertical)
        new_geometry = self.geometry(self.candidate_vertical)
        return (
            abs(math.log(new.mg)) < abs(math.log(old.mg))
            and new.vg < old.vg
            and new.fac2 >= old.fac2
            and abs(new_geometry["width_ratio"] - 1.0)
            < abs(old_geometry["width_ratio"] - 1.0)
            and new_geometry["centre_mae"] <= old_geometry["centre_mae"]
        )

    def report(self) -> str:
        """Return the frozen baseline/candidate decision table as text."""
        lines = [
            "PRESLHY conserved near-field/crosswind validation",
            f"  Gaussian establishment: {self.nearfield_establishment}",
            f"  crosswind entrainment  : {self.crosswind_entrainment}",
            f"  compatible handoff     : {self.compatible_handoff_search}",
            f"  geometry observation   : {self.geometry_observation}",
            f"  selected trials       : {self.selected_trials}",
            f"  accepted interfaces   : "
            f"{sum(result.accepted for result in self.handoffs.values())}"
            f"/{len(self.selected_trials)}",
        ]
        if self.failures:
            lines.append("  failures:")
            lines += [
                f"    - trial {trial}: {reason}"
                for trial, reason in sorted(self.failures.items())
            ]
        for label, result in (
            ("baseline", self.baseline), ("candidate", self.candidate),
        ):
            if result.pairs:
                stats = result.statistics()
                lo, hi = result.interval()
                lines.append(
                    f"  {label:9s} concentration: n={len(result.pairs)}, "
                    f"MG={stats.mg:.3f} [{lo:.3f}, {hi:.3f}], "
                    f"VG={stats.vg:.3f}, FAC2={stats.fac2:.3f}"
                )
            else:
                lines.append(f"  {label:9s} concentration: no common arcs")
        for label, rows in (
            ("baseline", self.baseline_vertical),
            ("candidate", self.candidate_vertical),
        ):
            metric = self.geometry(rows)
            lines.append(
                f"  {label:9s} geometry     : n={len(rows)}, "
                f"sigma ratio={metric['width_ratio']:.3f}, "
                f"centre bias={metric['centre_mean_error']:.3f} m, "
                f"centre MAE={metric['centre_mae']:.3f} m"
            )
        lines.append(
            f"  pre-registered decision: "
            f"{'promote' if self.promoted else 'do not promote'}"
        )
        return "\n".join(lines)


@dataclass
class IndependentEnergyInterfaceValidation:
    """Frozen seven-trial feasibility audit for the five-flux boundary."""

    selected_trials: list[int]
    interfaces: dict[int, object] = field(default_factory=dict, repr=False)
    failures: dict[int, str] = field(default_factory=dict)

    @property
    def all_interfaces_accepted(self) -> bool:
        return (
            set(self.interfaces) == set(self.selected_trials)
            and not self.failures
            and all(result.accepted for result in self.interfaces.values())
        )

    def report(self) -> str:
        lines = [
            "PRESLHY independent-energy interface validation",
            f"  selected trials       : {self.selected_trials}",
            f"  accepted interfaces   : "
            f"{sum(result.accepted for result in self.interfaces.values())}"
            f"/{len(self.selected_trials)}",
        ]
        if self.interfaces:
            lines += [
                f"  maximum flux residual : "
                f"{max(max(result.relative_residuals.values()) for result in self.interfaces.values()):.3e}",
                f"  maximum width residual: "
                f"{max(result.halfwidth_residual for result in self.interfaces.values()):.3e}",
                f"  maximum temperature residual: "
                f"{max(result.temperature_residual for result in self.interfaces.values()):.3f} K",
            ]
        if self.failures:
            lines.append("  failures:")
            lines += [
                f"    - trial {trial}: {reason}"
                for trial, reason in sorted(self.failures.items())
            ]
        return "\n".join(lines)


# ==========================================================================
# running the model
# ==========================================================================


def hydrogen_jet(
    *, rate, diameter, wind, height, ambient_temperature,
    relative_humidity=60.0, ambient_pressure=101325.0,
    storage_pressure_barg=5.0, roughness=0.001, stability="D",
    storage_temperature=None,
    source_upstream_pressure_barg=None,
    ground_effect=False, alfa1=0.0875,
    averaging=60.0, wind_reference_height=1.5, corrections=False,
    houf_entrainment=False, source_table_consistency=True,
    source_momentum_consistency=True,
    bulk_air_phase_safe_source=False,
    condensed_air_particle_diameter=None,
    source_total_energy_consistency=False,
    source_pressure_thrust=False,
    evaporation_zone_distance=False,
    condensed_air_stationary_bound=False,
    ground_layer_entrainment=False,
    hydrogen_spin_isomer="normal",
):
    """Build the liquid hydrogen jet for one release.

    **This did not exist.** Every liquid hydrogen jet result in this project
    was produced by assembling the thermodynamics, the boundary layer, the
    flash and the coefficients by hand in a working session, and none of that
    assembly was in the package. The results were therefore not reproducible
    from the repository -- which is why none of them had a test.

    ``corrections`` switches the five adopted corrections together, so the
    as-shipped and corrected figures come from one code path with one flag
    between them. ``alfa1`` defaults to Papanicolaou and List's measured plume
    value; it remains an argument because the historical reference dumps used
    Fischer et al.'s 0.0833 and must be reproduced explicitly.

    ``bulk_air_phase_safe_source`` is an off-by-default model-form bound.  It
    advances the initial plug flow to the first all-gas N2/O2/Ar state near
    68 K.  It removes the invalid 20 K gaseous-air state, but is not part of
    ``corrections`` because the pre-registered validation worsened variance
    and centre-height error; see ``docs/prereg-bulk-air-phase-boundary.md``.

    ``condensed_air_particle_diameter`` selects the off-by-default transported
    N2/O2 source zone.  The value is a particle diameter in metres.  It starts
    from a solid-air/H2 evaporation endpoint, marches mass, momentum and
    enthalpy until retained condensate is below the declared handoff limit,
    and only then starts the single-phase plume.  It may only be used with
    ``corrections=True``.

    ``source_total_energy_consistency`` includes source kinetic energy in the
    storage-to-evaporation balance. ``source_pressure_thrust`` additionally
    uses the HyRAM+ homogeneous-equilibrium throat and Yuceil--Otugen
    pressure-thrust construction. Both are research options tied to the
    transported source; neither changes the validated default.

    ``evaporation_zone_distance`` preserves the upstream Zone-III length from
    Li et al. equation 22 instead of locating an already air-loaded
    evaporation endpoint at the physical orifice. It changes only the source
    coordinate and is off by default.

    ``condensed_air_stationary_bound`` selects Li et al.'s opposite kinematic
    limit: condensed N2/O2 has zero axial velocity and leaves the atmospheric
    stream carrying enthalpy but no axial momentum. It requires the
    total-energy source and is a diagnostic bound, not a calibrated option.

    ``storage_temperature`` separates the liquid caloric state from tanker
    driving pressure.  The historical default assumes saturation at tanker
    pressure.  PRESLHY's pressure-driven subcooled liquid can instead supply
    its independently reconstructed temperature without changing the
    mechanical-pressure input.

    ``source_upstream_pressure_barg`` is the independently measured gauge
    pressure at the same source plane as ``storage_temperature``.  Supplying
    both represents a compressed/subcooled pipe state; it is not substituted
    for the tanker pressure used elsewhere in the release description.

    ``hydrogen_spin_isomer`` selects the same normal, para or ortho CoolProp
    equation from the stored liquid through the source construction. Normal
    hydrogen remains the default; the alternatives are composition
    sensitivities unless the released spin fraction is documented.
    """
    from CoolProp.CoolProp import PropsSI

    spin_species = {
        "normal": "Hydrogen",
        "para": "ParaHydrogen",
        "ortho": "OrthoHydrogen",
    }
    try:
        hydrogen_property_species = spin_species[
            hydrogen_spin_isomer.lower()
        ]
    except (AttributeError, KeyError) as exc:
        raise ValueError(
            "hydrogen spin isomer must be 'normal', 'para' or 'ortho'"
        ) from exc
    if (
        hydrogen_property_species != "Hydrogen"
        and condensed_air_particle_diameter is not None
    ):
        raise ValueError(
            "non-normal hydrogen is not yet available with transported "
            "condensed-air particles"
        )

    if (
        source_total_energy_consistency or source_pressure_thrust
        or evaporation_zone_distance
        or condensed_air_stationary_bound
    ) and condensed_air_particle_diameter is None:
        raise ValueError(
            "source energy/thrust options require the transported "
            "condensed-air source"
        )

    from ..core.atmosphere import (
        absolute_humidity, friction_velocity, psi, stability_defaults,
    )
    from ..core.jetplume import J_X, JetCoefficients, JetPlume
    from ..core.thermo import (
        AmbientConditions, CoolPropBackend, GasProperties, Thermo,
    )
    from .flashing import equivalent_source

    backend = CoolPropBackend(
        hydrogen_property_species,
        force_contaminant_gas=corrections and source_table_consistency,
    )
    humid, _rh = absolute_humidity(
        ambient_temperature, ambient_pressure / 101325.0,
        backend.water_vapour_pressure, relhum=relative_humidity,
    )
    boil = PropsSI(
        "T", "P", ambient_pressure, "Q", 0, hydrogen_property_species
    )
    gas = GasProperties(
        name="LH2", mw=2.016, temp=boil,
        rho=PropsSI(
            "D", "T", boil, "Q", 1, hydrogen_property_species
        ),
        coolprop_name=hydrogen_property_species, ulc=0.75, llc=0.04,
    )
    ambient = AmbientConditions(
        tamb=ambient_temperature, pamb=ambient_pressure / 101325.0,
        humid=humid, tsurf=ambient_temperature, ihtfl=1, iwtfl=1,
    )
    th = Thermo(gas=gas, ambient=ambient, backend=backend,
                legacy_numerics=False)
    th.reference_enthalpies()
    th.ambient.humsrc = 0.0

    if source_upstream_pressure_barg is not None and storage_temperature is None:
        raise ValueError(
            "source upstream pressure requires an explicit source temperature"
        )
    if storage_temperature is None:
        stored = PropsSI(
            "T", "P", (storage_pressure_barg + 1.013) * 1e5,
            "Q", 0, hydrogen_property_species,
        )
    else:
        stored = float(storage_temperature)
        if not (
            PropsSI("Ttriple", hydrogen_property_species) < stored
            < PropsSI("Tcrit", hydrogen_property_species)
        ):
            raise ValueError(
                "explicit LH2 storage temperature must lie between the "
                "triple and critical temperatures"
            )
    source_upstream_pressure = (
        None
        if source_upstream_pressure_barg is None
        else ambient_pressure + float(source_upstream_pressure_barg) * 1.0e5
    )
    if source_upstream_pressure is not None and source_upstream_pressure <= ambient_pressure:
        raise ValueError("source upstream pressure must exceed ambient pressure")
    flash = equivalent_source(
        hydrogen_property_species, storage_temperature=stored,
        ambient_temperature=ambient_temperature,
        ambient_pressure=ambient_pressure, molecular_weight=2.016,
        bulk_air_phase_safe=corrections and bulk_air_phase_safe_source,
        storage_pressure=source_upstream_pressure,
    )

    transported_source = None
    source_nozzle = None
    if source_upstream_pressure is not None:
        from ..addons.notional import energy_conserving_notional_nozzle

        source_nozzle = energy_conserving_notional_nozzle(
            fluid=hydrogen_property_species,
            storage_temperature=stored,
            storage_pressure=source_upstream_pressure,
            mass_flow=rate,
            orifice_diameter=diameter,
            ambient_pressure=ambient_pressure,
        )
        if not source_nozzle.admissible:
            reasons = []
            if source_nozzle.discharge_coefficient > 1.0:
                reasons.append(
                    f"Cd={source_nozzle.discharge_coefficient:.3f} > 1"
                )
            if source_nozzle.temperature is None:
                reasons.append("no ambient-pressure total-energy state")
            raise ValueError(
                "measured-pressure source incompatible: " + ", ".join(reasons)
            )
    evaporation_distance = 0.0
    if condensed_air_particle_diameter is not None:
        if not corrections:
            raise ValueError(
                "the condensed-air source requires corrections=True"
            )
        if not source_table_consistency or not source_momentum_consistency:
            raise ValueError(
                "the condensed-air source requires consistent source-table "
                "and momentum closures"
            )
        if bulk_air_phase_safe_source:
            raise ValueError(
                "select either the transported or bulk phase-safe source, "
                "not both"
            )
        if condensed_air_particle_diameter <= 0.0:
            raise ValueError("particle diameter must be positive")
        if source_pressure_thrust and not source_total_energy_consistency:
            raise ValueError(
                "pressure thrust requires total-energy consistency"
            )
        if condensed_air_stationary_bound and not source_total_energy_consistency:
            raise ValueError(
                "the stationary-condensate bound requires total-energy "
                "consistency"
            )
        from ..addons.cryogenic_air import (
            multiphase_hydrogen_source_plane,
            transported_condensed_air_source,
        )

        specific_momentum = None
        if source_pressure_thrust:
            if source_nozzle is None:
                from ..addons.notional import energy_conserving_notional_nozzle

                source_nozzle = energy_conserving_notional_nozzle(
                    fluid=hydrogen_property_species,
                    storage_temperature=stored,
                    mass_flow=rate,
                    orifice_diameter=diameter,
                    ambient_pressure=ambient_pressure,
                )
            if not source_nozzle.admissible:
                reasons = []
                if source_nozzle.discharge_coefficient > 1.0:
                    reasons.append(
                        f"Cd={source_nozzle.discharge_coefficient:.3f} > 1"
                    )
                if source_nozzle.temperature is None:
                    reasons.append("no ambient-pressure total-energy state")
                raise ValueError(
                    "pressure-thrust source incompatible: " + ", ".join(reasons)
                )
            specific_momentum = source_nozzle.specific_momentum

        plane = multiphase_hydrogen_source_plane(
            hydrogen_flow=rate,
            orifice_diameter=diameter,
            orifice_density=flash.orifice_density,
            storage_temperature=stored,
            ambient_temperature=ambient_temperature,
            ambient_pressure=ambient_pressure,
            specific_momentum=specific_momentum,
            include_kinetic_energy=source_total_energy_consistency,
            particle_velocity_fraction=(
                0.0 if condensed_air_stationary_bound else 1.0
            ),
            storage_pressure=source_upstream_pressure,
            hydrogen_species=hydrogen_property_species,
        )
        if evaporation_zone_distance:
            evaporation_distance = plane.formation_distance
        transported_source = transported_condensed_air_source(
            hydrogen_flow=rate,
            station2_temperature=plane.endpoint.temperature,
            station2_velocity=plane.velocity,
            station2_density=plane.density,
            station2_diameter=plane.diameter,
            particle_diameter=condensed_air_particle_diameter,
            initial_nitrogen_flow=plane.nitrogen_flow,
            initial_oxygen_flow=plane.oxygen_flow,
            ambient_temperature=ambient_temperature,
            ambient_pressure=ambient_pressure,
            station2_kinetic_energy=(
                plane.kinetic_energy_flow
                if condensed_air_stationary_bound else None
            ),
            stationary_condensate=condensed_air_stationary_bound,
        )
        if not transported_source.handoff_reached:
            raise ValueError(
                "condensed-air source did not reach its single-phase "
                "handoff criterion within 2 m"
            )
        retained_air = (
            transported_source.nitrogen_gas_flow
            + transported_source.oxygen_gas_flow
            + transported_source.nitrogen_condensed_flow
            + transported_source.oxygen_condensed_flow
        )
        source_total = rate + retained_air
        source_fraction = rate / source_total
        source_air = retained_air / source_total
        source_temperature = transported_source.temperature
        th.gas.temp = source_temperature
        source_enthalpy = th.enthalpy(
            source_fraction, source_air, source_temperature
        )
        th.table = th.build_adiabatic_table(
            source_fraction, source_air, source_enthalpy
        )
    elif corrections and source_table_consistency:
        # The flash endpoint has already entrained air. It is therefore a
        # secondary source, not the pure-H2 end of the original mixing line.
        # Rebuild exactly as DEG2S does for an already mixed source; otherwise
        # the concentration lookup takes the non-monotone pure-H2 line's
        # high-concentration branch and numerically deletes the entrained air.
        total = 1.0 + flash.air_ratio * (1.0 + humid)
        source_fraction = 1.0 / total
        source_air = flash.air_ratio / total
        source_temperature = flash.temperature
        # The contaminant component is vapour at its partial pressure. Keep
        # the thermodynamic inversion's lower bracket on that gas state, not
        # on pure H2's slightly higher normal boiling point.
        th.gas.temp = source_temperature
        source_enthalpy = th.enthalpy(
            source_fraction, source_air, source_temperature
        )
        th.table = th.build_adiabatic_table(
            source_fraction, source_air, source_enthalpy
        )
    else:
        source_fraction = flash.mass_fraction if corrections else 1.0
        th.table = th.build_adiabatic_table(1.0, 0.0, th.hmrte)

    d = stability_defaults(roughness, stability, averaging)
    ustar = friction_velocity(wind, wind_reference_height, roughness, d.rml)
    coefficients = JetCoefficients(
        alfa1=alfa1 if corrections else 0.057,
        density_scaled_entrainment=corrections,
        houf_buoyant_entrainment=houf_entrainment,
        ground_layer_entrainment=ground_layer_entrainment,
    )
    jp = JetPlume(
        th, coefficients=coefficients,
        u0=wind, z0=wind_reference_height, zr=roughness, rml=d.rml,
        # rhoe is the saturated-vapour density and stays that way. Only
        # `rho_exit` in the initial conditions takes the expanded value: the
        # expansion changes the state at the starting plane, not the property
        # of the released substance.
        ustar=ustar, rhoa=th.table.rhoa, rhoe=gas.rho,
        deltay=d.deltay, betay=d.betay, deltaz=d.deltaz, betaz=d.betaz,
        gammaz=d.gammaz, yclow=1.0e-5, ground_effect=ground_effect,
    )

    # the starting plane: the orifice as shipped, the expanded plane with the
    # corrections on. That is the first source-plane correction, and the one EPA's 1991
    # evaluation records the SLAB developer objecting to.
    if corrections:
        if transported_source is not None:
            # The transport march has already closed continuity and momentum
            # at this plane.  Its sub-1%-of-peak residual condensate is folded
            # into the bulk source state as the explicit handoff approximation.
            density = transported_source.density
            fraction = source_fraction
            dia = transported_source.diameter
            jp.condensed_air_source = transported_source
        else:
            density, fraction, dia = JetPlume.expanded_source_start(
                mass_flow=rate, orifice_diameter=diameter,
                orifice_density=flash.orifice_density,
                expanded_density=(
                    th.table.rhoe if source_table_consistency else flash.density
                ),
                expanded_fraction=source_fraction,
                entrainment_momentum_conserving=source_momentum_consistency,
                source_specific_momentum=(
                    source_nozzle.specific_momentum
                    if source_nozzle is not None else None
                ),
            )
    else:
        # The orifice plane: DEGADIS as shipped starts the jet at the hole.
        #
        # This is the one input with no reference dump behind it. The session
        # that produced the published figures supplied a full parameter set
        # for the *corrected* run, which this reproduces to four significant
        # figures, and none for the baseline. Candidates -- orifice against
        # expanded density, mean against peak flow -- gave a legacy width
        # ratio between 0.49 and 0.65 against a reported 0.64, and a rise
        # between 0.20 and 0.24 m against a reported 1.07 m. Those width
        # numbers used an e-folding width mislabeled as sigma; the current
        # validation converts it to a standard deviation. The spread is
        # reachable; the rise is not, under any of them.
        density, fraction, dia = flash.orifice_density, 1.0, diameter

    ua = ustar / VKC * (
        math.log((height + roughness) / roughness) - psi(height, d.rml)
    )
    if (
        corrections and source_table_consistency
        and source_momentum_consistency
    ):
        # Retain the already-computed atmospheric source as one auditable
        # object.  The legacy JETPLU initial state stores concentration and
        # integrated widths, so reconstructing diameter and total velocity
        # from it later would lose the physical plane used above.  This is
        # also the exact input required by the conserved axisymmetric model.
        from ..addons.axisymmetric_jet import AxisymmetricJetSource

        area = math.pi * dia**2 / 4.0
        source_velocity = rate / max(fraction * density * area, 1.0e-30)
        jp.axisymmetric_source = AxisymmetricJetSource(
            diameter=dia,
            velocity=source_velocity,
            density=density,
            temperature=source_temperature,
            mass_fraction=fraction,
            theta=0.0,
            x=evaporation_distance + (
                transported_source.distance
                if transported_source is not None else 0.0
            ),
            y=height,
        )
        jp.local_source_wind = ua
    y0 = jp.initial_conditions_directed(
        erate=rate, diajet=dia, elejet=height, ua=max(ua, 0.1),
        theta0=0.0, rho_exit=density, concentration=fraction,
    )
    if transported_source is not None:
        # The source march starts at the orifice and has a real axial length;
        # the Gaussian development length begins at its downstream handoff.
        y0[J_X] += transported_source.distance + evaporation_distance
        jp.evaporation_zone_distance = evaporation_distance
    if source_nozzle is not None:
        jp.source_notional_nozzle = source_nozzle
    return jp, y0


def hydrogen_gas_jet(
    *, rate, diameter, velocity, wind, height, source_temperature,
    source_density, source_mass_fraction=1.0,
    theta=math.pi / 2.0, relative_humidity=60.0,
    ambient_temperature=288.15, ambient_pressure=101325.0,
    roughness=0.001, stability="D", averaging=60.0,
    wind_reference_height=10.0, alfa1=0.0875, sc=1.42,
    density_scaled_entrainment=True,
    momentum_entrainment_beta=0.0,
    houf_entrainment=False,
):
    """Build a cold gaseous-hydrogen jet at a measured outlet plane.

    Unlike :func:`hydrogen_jet`, this path performs no flash and no equivalent
    source expansion. ``rate``, ``diameter``, ``velocity``, ``temperature``
    and ``density`` must all describe the same atmospheric outlet plane. It is
    used for the Spadeadam ventilation mast, whose supply nozzle sits inside a
    closed room and is not the dispersion source.

    ``rate`` is hydrogen mass flow. ``source_mass_fraction`` is one at a pure
    gas outlet and below one at an already air-loaded single-phase source.
    ``velocity`` is checked against total-flow continuity so one internally
    consistent plane enters the integrator.
    """
    from ..core.atmosphere import (
        absolute_humidity, friction_velocity, psi, stability_defaults,
    )
    from ..core.jetplume import JetCoefficients, JetPlume
    from ..core.thermo import (
        AmbientConditions, CoolPropBackend, GasProperties, Thermo,
    )

    if not 0.0 < source_mass_fraction <= 1.0:
        raise ValueError("source mass fraction must lie in (0, 1]")
    area = math.pi * diameter * diameter / 4.0
    derived_velocity = (
        rate / source_mass_fraction / max(source_density * area, 1.0e-30)
    )
    if not math.isclose(derived_velocity, velocity, rel_tol=0.03):
        raise ValueError(
            "mast source rate, density, diameter and velocity are inconsistent "
            f"({derived_velocity:.3g} versus {velocity:.3g} m/s)"
        )

    backend = CoolPropBackend("Hydrogen", force_contaminant_gas=True)
    humid, _rh = absolute_humidity(
        ambient_temperature, ambient_pressure / 101325.0,
        backend.water_vapour_pressure, relhum=relative_humidity,
    )
    gas = GasProperties(
        name="H2", mw=2.016, temp=source_temperature, rho=source_density,
        coolprop_name="Hydrogen", ulc=0.75, llc=0.04,
    )
    ambient = AmbientConditions(
        tamb=ambient_temperature, pamb=ambient_pressure / 101325.0,
        humid=humid, tsurf=ambient_temperature, ihtfl=1, iwtfl=1,
    )
    th = Thermo(gas=gas, ambient=ambient, backend=backend,
                legacy_numerics=False)
    th.reference_enthalpies()
    th.ambient.humsrc = 0.0
    if source_mass_fraction == 1.0:
        th.table = th.build_adiabatic_table(1.0, 0.0, th.hmrte)
    else:
        source_air = (1.0 - source_mass_fraction) / (1.0 + humid)
        source_enthalpy = th.enthalpy(
            source_mass_fraction, source_air, source_temperature
        )
        th.table = th.build_adiabatic_table(
            source_mass_fraction, source_air, source_enthalpy
        )

    defaults = stability_defaults(roughness, stability, averaging)
    ustar = friction_velocity(
        wind, wind_reference_height, roughness, defaults.rml
    )
    jp = JetPlume(
        th,
        coefficients=JetCoefficients(
            alfa1=alfa1,
            sc=sc,
            density_scaled_entrainment=density_scaled_entrainment,
            momentum_entrainment_beta=momentum_entrainment_beta,
            houf_buoyant_entrainment=houf_entrainment,
        ),
        u0=wind, z0=wind_reference_height, zr=roughness,
        rml=defaults.rml, ustar=ustar, rhoa=th.table.rhoa,
        rhoe=source_density, deltay=defaults.deltay, betay=defaults.betay,
        deltaz=defaults.deltaz, betaz=defaults.betaz,
        gammaz=defaults.gammaz, yclow=1.0e-5, ground_effect=False,
    )
    local_wind = ustar / VKC * (
        math.log((height + roughness) / roughness)
        - psi(height, defaults.rml)
    )
    # Preserve the historical 0.1 m/s floor for atmospheric releases, while
    # allowing an explicit quiescent laboratory free jet.  The directed
    # initial-condition routine evaluates the zero-coflow limit analytically.
    initial_ambient_speed = 0.0 if wind == 0.0 else max(local_wind, 0.1)
    y0 = jp.initial_conditions_directed(
        erate=rate, diajet=diameter, elejet=height,
        ua=initial_ambient_speed, theta0=theta,
        rho_exit=source_density, concentration=source_mass_fraction,
    )
    return jp, y0


class Trajectory:
    """A finished jet run, queryable at a position.

    Wraps the ``rows`` array so that a comparison asks for the concentration
    *at a sensor* rather than on the model's own centreline. By 6 m the
    modelled centre for a 0.5 m release sits above the topmost sensor, so the
    centreline is a value the array could not have measured.
    """

    def __init__(self, table, rows):
        self.table, self.rows = table, rows

    @property
    def ok(self) -> bool:
        return len(self.rows) > 2

    def at(self, x: float):
        """Centre height, spreads and centreline concentration at ``x``."""
        if not self.ok:
            return None
        xs = self.rows[:, 0]
        if not (xs[0] <= x <= xs[-1]):
            return None
        order = np.argsort(xs)
        xs, rows = xs[order], self.rows[order]
        def take(column):
            return float(np.interp(x, xs, rows[:, column]))

        return _State(z=take(1), cc=take(2), sy=take(3), sz=take(4))

    def concentration_at(self, x: float, y: float, z: float) -> float:
        """Volume per cent at a point, with the ground image included."""
        st = self.at(x)
        if st is None or st.sz <= 0.0 or st.sy <= 0.0:
            return 0.0
        lateral = math.exp(-0.5 * (y / st.sy) ** 2)
        direct = math.exp(-0.5 * ((z - st.z) / st.sz) ** 2)
        image = math.exp(-0.5 * ((z + st.z) / st.sz) ** 2)
        cc = st.cc * lateral * (direct + image)
        if cc <= 0.0:
            return 0.0
        return 100.0 * self.table.from_concentration(cc).yc

    def temperature_at(self, x: float, y: float, z: float) -> float:
        """Temperature in K at a point, with the ground image included."""
        st = self.at(x)
        if st is None or st.sz <= 0.0 or st.sy <= 0.0:
            return float(self.table.t[0])
        lateral = math.exp(-0.5 * (y / st.sy) ** 2)
        direct = math.exp(-0.5 * ((z - st.z) / st.sz) ** 2)
        image = math.exp(-0.5 * ((z + st.z) / st.sz) ** 2)
        return float(self.table.from_concentration(
            st.cc * lateral * (direct + image)
        ).temp)


class IndependentEnergyTrajectory:
    """Query an independent-energy result at PRESLHY receptor positions."""

    def __init__(self, interface, result):
        self.model = interface.model
        self.result = result
        order = np.argsort(result.states[:, 5])
        self.x = result.states[order, 5]
        self.states = result.states[order]

    @property
    def ok(self) -> bool:
        return len(self.states) > 2

    def state_at(self, x: float) -> np.ndarray | None:
        if not self.ok or not (self.x[0] <= x <= self.x[-1]):
            return None
        return np.array([
            np.interp(x, self.x, self.states[:, column])
            for column in range(self.states.shape[1])
        ])

    def at(self, x: float):
        state = self.state_at(x)
        if state is None:
            return None
        sy, sz = self.model.section_widths(state)
        return _State(
            z=float(state[6]),
            cc=float(state[0] * state[1]),
            sy=sy,
            sz=sz,
        )

    def concentration_at(self, x: float, y: float, z: float) -> float:
        state = self.state_at(x)
        if state is None:
            return 0.0
        return 100.0 * self.model.point_mole_fraction(state, y, z)

    def temperature_at(self, x: float, y: float, z: float) -> float:
        state = self.state_at(x)
        if state is None:
            return float(self.model.thermodynamics.ambient_temperature)
        return self.model.point_temperature(state, y, z)


def _vertical_sensor_heights(trial: dict, x: float) -> np.ndarray:
    """Return the exact sensor heights admitted to a measured vertical fit."""
    return np.array(sorted({
        float(sensor["z"])
        for sensor in trial["sensors"]
        if (
            abs(float(sensor["x"]) - float(x)) < 5.0e-3
            and abs(float(sensor["y"])) < 0.01
            and float(sensor["peak"]) > 0.05
        )
    }))


def _fit_model_vertical_profile(trajectory, x: float, heights: np.ndarray) -> dict:
    """Apply the PRESLHY single-Gaussian fit to modelled sensor values.

    The plume concentration operator includes the ground image. Fitting its
    values at the actual instrument heights makes the model and measurement
    pass through the same observation operator before their fitted centre and
    width are compared.
    """
    from scipy.optimize import curve_fit

    heights = np.asarray(heights, dtype=float)
    if heights.size < 4 or not np.all(np.isfinite(heights)):
        raise ValueError("a model vertical fit requires four finite heights")
    values = np.array([
        trajectory.concentration_at(float(x), 0.0, float(height))
        for height in heights
    ])
    if not np.all(np.isfinite(values)) or np.max(values) <= 0.0:
        raise ValueError("model vertical profile has no finite positive signal")

    def gaussian(coordinate, peak, centre, sigma):
        return peak * np.exp(-0.5 * (((coordinate - centre) / sigma) ** 2))

    fitted, _covariance = curve_fit(
        gaussian,
        heights,
        values,
        p0=[
            float(np.max(values)),
            float(heights[int(np.argmax(values))]),
            0.4 / math.sqrt(2.0),
        ],
        bounds=(
            [0.0, -1.0, 0.03 / math.sqrt(2.0)],
            [110.0, 6.0, 10.0 / math.sqrt(2.0)],
        ),
        maxfev=20000,
    )
    represented = gaussian(heights, *fitted)
    residual = float(np.sum((values - represented) ** 2))
    spread = float(np.sum((values - np.mean(values)) ** 2))
    return {
        "centre": float(fitted[1]),
        "sigma_z": float(fitted[2]),
        "r_squared": float(1.0 - residual / max(spread, 1.0e-12)),
    }


@dataclass
class _State:
    z: float
    cc: float
    sy: float
    sz: float


def _plume(record: dict, *, corrections: bool):
    """The jet path for one trial. ``None`` when the source is unusable."""
    from CoolProp.CoolProp import PropsSI

    from .flashing import equivalent_source

    diameter = float(record["orifice_mm"]) / 1000.0
    rate = float(record["flow_used_gs"]) / 1000.0
    if rate <= 0.005 or diameter <= 0.0:
        return None, None

    tamb = float(record["T_mean_C"]) + 273.15
    barg = float(record["tanker_barg"])
    stored = PropsSI("T", "P", (barg + 1.013) * 1e5, "Q", 0, "Hydrogen")
    flash = equivalent_source(
        "Hydrogen", storage_temperature=stored, ambient_temperature=tamb,
        ambient_pressure=101325.0, molecular_weight=2.016,
    )
    area = math.pi * diameter * diameter / 4.0
    exit_speed = rate / (flash.orifice_density * area)

    jp, y0 = hydrogen_jet(
        rate=rate, diameter=diameter,
        wind=float(record["wind_mean_ms"]),
        height=float(record["release_height_m"]),
        ambient_temperature=tamb,
        relative_humidity=float(record["RH_mean_pct"]),
        storage_pressure_barg=barg,
        wind_reference_height=float(record.get("wind_ref_height_m") or 1.5),
        corrections=corrections,
    )
    result = jp.run(y0, distmx=STEP, smax=REACH)
    traj = Trajectory(jp.th.table, result.rows)
    return exit_speed, (traj if traj.ok else None)


def _pair_arcs(number, readings, plume, *, corrections: bool) -> list[Pair]:
    """Reduce one trial's readings to one pair per arc."""
    by_arc: dict[float, list] = {}
    for r in readings:
        by_arc.setdefault(round(r.x, 3), []).append(r)

    pairs = []
    for x, group in sorted(by_arc.items()):
        state = plume.at(x)
        if state is None:
            continue
        modelled = [(plume.concentration_at(x, r.y, r.z), r.z) for r in group]
        observed = [(r.peak, r.z) for r in group]
        o, o_at = max(observed)
        p, p_at = max(modelled)
        if o <= 0.0 or p <= 0.0:
            continue
        pairs.append(Pair(
            trial=number, x=x, observed=o, predicted=p,
            observed_at=o_at, predicted_at=max(r.z for r in group),
            plume_centre=state.z,
        ))
    return pairs


def _build(records, readings_by_trial, *, corrections: bool) -> NearField:
    out = NearField(corrections=corrections)
    for number, readings in sorted(readings_by_trial.items()):
        record = records.get(number)
        if record is None:
            out.excluded[number] = "no conditions row"
            continue
        if record.get("orientation") != "horizontal":
            out.excluded[number] = f"not horizontal ({record.get('orientation')})"
            continue
        exit_speed, plume = _plume(record, corrections=corrections)
        if plume is None:
            out.excluded[number] = "no usable source flow"
            continue
        ratio = exit_speed / float(record["wind_mean_ms"])
        if ratio < MOMENTUM_RATIO:
            out.excluded[number] = (
                f"wind-steered, exit speed {ratio:.1f}x the wind"
            )
            continue
        out.pairs += _pair_arcs(
            number, readings, plume, corrections=corrections
        )
    return out


# ==========================================================================
# the two ways in
# ==========================================================================


def _conditions(path) -> dict:
    with open(path) as fh:
        return {int(r["trial"]): r for r in csv.DictReader(fh)}


def from_workbooks(root, report, conditions, *, corrections: bool) -> NearField:
    """The authoritative route: read the dataset and compare."""
    from .preslhy import load

    records = _conditions(conditions)
    heights = {n: float(r["release_height_m"]) for n, r in records.items()}
    trials = load(Path(root), Path(report), heights=heights)
    return _build(
        records,
        {t.number: t.readings for t in trials},
        corrections=corrections,
    )


def reduce(root, report, conditions, out) -> int:
    """Write the reduced near-field table, so the statistic can be rechecked.

    One row per sensor per trial, carrying position, peak, mean and the
    provenance of each row. This is a derived product of this work, not the
    PRESLHY dataset; see ``reference/preslhy/README.md``.
    """
    from .preslhy import load

    records = _conditions(conditions)
    heights = {n: float(r["release_height_m"]) for n, r in records.items()}
    trials = load(Path(root), Path(report), heights=heights)

    rows = 0
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(REDUCED_COLUMNS)
        for t in sorted(trials, key=lambda t: t.number):
            for r in t.readings:
                writer.writerow([
                    t.number, r.serial, f"{r.x:.4g}", f"{r.y:.4g}",
                    f"{r.z:.4g}", f"{r.z_axis:.4g}", f"{r.peak:.6g}",
                    f"{r.mean:.6g}", r.samples, t.path.name,
                    t.window[0], t.window[1],
                ])
                rows += 1
    return rows


def from_reduced(
    table, *, corrections: bool, convention: str = "sensor",
    source: str = "flow_mean_gs", ground_effect: bool = False,
    liquid_temperature_mode: str = "tank_saturation",
    houf_entrainment: bool = False, source_table_consistency: bool = True,
    source_momentum_consistency: bool = True,
    bulk_air_phase_safe_source: bool = False,
    condensed_air_particle_diameter: float | None = None,
    source_total_energy_consistency: bool = False,
    source_pressure_thrust: bool = False,
    evaporation_zone_distance: bool = False,
    condensed_air_stationary_bound: bool = False,
    ground_layer_entrainment: bool = False,
) -> NearField:
    """The reproducible route: the reduced table, no third-party data needed.

    ``convention`` decides what the model is scored on, and it is not a
    detail -- see the module docstring and ``docs/lh2-results.md``.

    ``"centreline"``
        The prediction is the model's own centreline at that distance. This is
        what the reported statistics used.
    ``"sensor"``
        The prediction is evaluated at each sensor's position and maximised
        over the arc, the same operation the measurement gets.

    ``source`` selects the flow that drives the release. The reduced table
    carries a mean and a peak over the same window and they differ by up to
    fifty per cent; the reported and regression-tested statistics use the
    default window mean.
    """
    import json

    if convention not in ("centreline", "sensor"):
        raise ValueError(f"unknown convention {convention!r}")

    data = json.loads(Path(table).read_text())
    out = NearField(corrections=corrections)
    for trial in data["trials"]:
        number = trial["trial"]
        if trial["orientation"] != "horizontal":
            out.excluded[number] = f"not horizontal ({trial['orientation']})"
            continue
        rate = trial[source] / 1000.0
        if rate <= 0.005:
            out.excluded[number] = "no usable source flow"
            continue
        ratio = _exit_ratio(
            trial, rate, liquid_temperature_mode=liquid_temperature_mode
        )
        if ratio < MOMENTUM_RATIO:
            out.excluded[number] = (
                f"wind-steered, exit speed {ratio:.1f}x the wind"
            )
            continue

        try:
            jp, y0 = hydrogen_jet(
                rate=rate, diameter=trial["orifice_mm"] / 1000.0,
                wind=trial["wind_ms"], height=trial["release_height_m"],
                ambient_temperature=trial["T_C"] + 273.15,
                relative_humidity=trial["RH_pct"],
                storage_pressure_barg=trial["tanker_barg"],
                storage_temperature=_liquid_source_temperature(
                    trial, liquid_temperature_mode
                ),
                wind_reference_height=trial["wind_ref_m"],
                corrections=corrections, ground_effect=ground_effect,
                houf_entrainment=houf_entrainment,
                source_table_consistency=source_table_consistency,
                source_momentum_consistency=source_momentum_consistency,
                bulk_air_phase_safe_source=bulk_air_phase_safe_source,
                condensed_air_particle_diameter=condensed_air_particle_diameter,
                source_total_energy_consistency=source_total_energy_consistency,
                source_pressure_thrust=source_pressure_thrust,
                evaporation_zone_distance=evaporation_zone_distance,
                condensed_air_stationary_bound=condensed_air_stationary_bound,
                ground_layer_entrainment=ground_layer_entrainment,
            )
        except ValueError as error:
            if source_pressure_thrust and str(error).startswith(
                "pressure-thrust source incompatible:"
            ):
                out.excluded[number] = str(error)
                continue
            raise
        traj = Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=REACH).rows)
        if not traj.ok:
            out.excluded[number] = "the plume integration did not run"
            continue

        arcs: dict[float, list] = {}
        for sensor in trial["sensors"]:
            arcs.setdefault(round(sensor["x"], 3), []).append(sensor)

        for x, group in sorted(arcs.items()):
            state = traj.at(x)
            if state is None:
                continue
            observed = max(s["peak"] for s in group)
            if convention == "sensor":
                predicted = max(
                    traj.concentration_at(x, s["y"], s["z"]) for s in group
                )
            else:
                predicted = traj.concentration_at(x, 0.0, state.z)
            if observed <= 0.0 or predicted <= 0.0:
                continue
            out.pairs.append(Pair(
                trial=number, x=x, observed=observed, predicted=predicted,
                observed_at=max(group, key=lambda s: s["peak"])["z"],
                predicted_at=max(s["z"] for s in group),
                plume_centre=state.z,
            ))
    return out


def _liquid_source_temperature(
    trial: dict, mode: str, ambient_pressure: float = 101325.0
) -> float | None:
    """Separate mechanical tanker pressure from the LH2 caloric state."""
    if mode == "tank_saturation":
        return None
    if mode == "ambient_boiling":
        from CoolProp.CoolProp import PropsSI

        return float(PropsSI(
            "T", "P", ambient_pressure, "Q", 0, "Hydrogen"
        ))
    raise ValueError(f"unknown liquid temperature mode {mode!r}")


def _measured_pipe_source_rows(table) -> dict[int, dict]:
    """Load the frozen pressure-loss/TC3 source reconstruction."""
    if table is None:
        return {}
    import json

    data = json.loads(Path(table).read_text())
    rows = {int(row["trial"]): row for row in data.get("trials", [])}
    missing = set(INDEPENDENT_ENERGY_TRIALS) - set(rows)
    if missing:
        raise ValueError(
            "measured pipe source table is missing trials "
            f"{sorted(missing)}"
        )
    return rows


def _trial_source_inputs(
    trial: dict,
    *,
    source: str,
    liquid_temperature_mode: str,
    measured_rows: dict[int, dict],
    measured_source_mode: str = "full",
) -> tuple[float, float | None, float | None]:
    """Return rate, caloric temperature and measured nozzle gauge pressure."""
    row = measured_rows.get(int(trial["trial"]))
    if row is None:
        return (
            float(trial[source]) / 1000.0,
            _liquid_source_temperature(trial, liquid_temperature_mode),
            None,
        )
    if measured_source_mode not in {"full", "flow_only", "nozzle_only"}:
        raise ValueError(
            f"unknown measured source mode {measured_source_mode!r}"
        )
    rate = (
        float(trial[source]) / 1000.0
        if measured_source_mode == "nozzle_only"
        else float(row["pressure_loss_mass_flow_g_s"]) / 1000.0
    )
    if measured_source_mode == "flow_only":
        return rate, _liquid_source_temperature(
            trial, liquid_temperature_mode
        ), None
    if row["topology"] == "open_pipe":
        # At ambient pressure TC3 lies on a two-phase boundary, where T,P do
        # not determine quality or enthalpy. Keep the prior open-pipe state.
        return rate, _liquid_source_temperature(
            trial, liquid_temperature_mode
        ), None
    if row["topology"] != "nozzle":
        raise ValueError(f"unknown measured source topology {row['topology']!r}")
    temperature = float(row["tc3_temperature_k"]["median"])
    pressure_barg = float(row["pt2_barg"]["median"])
    if pressure_barg <= 0.05:
        raise ValueError(
            f"trial {trial['trial']} nozzle PT2 does not exceed 0.05 barg"
        )
    return rate, temperature, pressure_barg


def _momentum_dominance_distance(
    source,
    *,
    ambient_density: float,
    local_wind: float,
) -> float:
    """Coefficient-free source/crosswind momentum-balance length, m."""
    if ambient_density <= 0.0 or local_wind <= 0.0:
        raise ValueError("ambient density and local wind must be positive")
    source_momentum = (
        source.density * source.velocity**2 * source.area
    )
    if source_momentum <= 0.0:
        raise ValueError("source momentum flux must be positive")
    return math.sqrt(source_momentum / (ambient_density * local_wind**2))


def _measured_lh2_equilibrium_source(
    *,
    enabled: bool,
    setup,
    trial: dict,
    rate: float,
    source_temperature: float | None,
    source_pressure_barg: float | None,
    measured_rows: dict[int, dict],
    measured_source_mode: str,
    hydrogen_spin_isomer: str = "normal",
):
    """Replace a measured nozzle endpoint by the collective LH2 phase bound."""
    if not enabled:
        return setup.axisymmetric_source
    if not measured_rows:
        raise ValueError(
            "the measured LH2 equilibrium bound requires a source table"
        )
    if measured_source_mode != "full":
        raise ValueError(
            "the measured LH2 equilibrium bound requires full source mode"
        )
    if source_pressure_barg is None:
        # The open-pipe trials have no identifiable T,P quality state and
        # retain the pre-registered prior source without per-trial selection.
        return setup.axisymmetric_source
    if source_temperature is None:
        raise ValueError("the measured LH2 equilibrium bound requires TC3")

    from ..addons.lh2_droplets import (
        homogeneous_equilibrium_hydrogen_source,
    )

    spin_species = {
        "normal": "Hydrogen",
        "para": "ParaHydrogen",
        "ortho": "OrthoHydrogen",
    }
    try:
        hydrogen_species = spin_species[hydrogen_spin_isomer.lower()]
    except (AttributeError, KeyError) as exc:
        raise ValueError(
            "hydrogen spin isomer must be 'normal', 'para' or 'ortho'"
        ) from exc

    bound = homogeneous_equilibrium_hydrogen_source(
        mass_flow=rate,
        orifice_diameter=trial["orifice_mm"] / 1000.0,
        upstream_temperature=source_temperature,
        upstream_pressure=101325.0 + source_pressure_barg * 1.0e5,
        ambient_temperature=trial["T_C"] + 273.15,
        ambient_pressure=101325.0,
        hydrogen_species=hydrogen_species,
        theta=0.0,
        y=trial["release_height_m"],
    )
    if max(
        bound.hydrogen_mass_residual,
        bound.total_mass_residual,
        bound.momentum_residual,
        bound.energy_residual,
    ) >= 1.0e-8:
        raise ValueError("measured LH2 equilibrium source fails conservation")
    return bound.source


def _exit_ratio(
    trial: dict, rate: float,
    *, liquid_temperature_mode: str = "tank_saturation",
) -> float:
    """Exit speed over wind speed, the criterion fixed before any comparison."""
    import math

    from CoolProp.CoolProp import PropsSI

    from .flashing import equivalent_source

    stored = _liquid_source_temperature(trial, liquid_temperature_mode)
    if stored is None:
        stored = PropsSI(
            "T", "P", (trial["tanker_barg"] + 1.013) * 1e5,
            "Q", 0, "Hydrogen",
        )
    flash = equivalent_source(
        "Hydrogen", storage_temperature=stored,
        ambient_temperature=trial["T_C"] + 273.15,
        ambient_pressure=101325.0, molecular_weight=2.016,
    )
    diameter = trial["orifice_mm"] / 1000.0
    area = math.pi * diameter * diameter / 4.0
    return rate / (flash.orifice_density * area) / trial["wind_ms"]


def common_arcs(a: NearField, b: NearField) -> tuple[NearField, NearField]:
    """Restrict two comparisons to the arcs both reached.

    Source corrections can start the established plume downstream of the
    closest instruments. Comparing unlike arc sets makes part of the
    difference a change of sample rather than a change of model.
    """
    keys = {(p.trial, p.x) for p in a.pairs} & {(p.trial, p.x) for p in b.pairs}
    out = []
    for field_ in (a, b):
        trimmed = NearField(
            pairs=[p for p in field_.pairs if (p.trial, p.x) in keys],
            excluded=dict(field_.excluded), corrections=field_.corrections,
        )
        out.append(trimmed)
    return out[0], out[1]


def vertical(
    table, *, corrections: bool, source: str = "flow_mean_gs",
    liquid_temperature_mode: str = "tank_saturation",
    momentum_filter: bool = False, ground_effect: bool = False,
    houf_entrainment: bool = False, source_table_consistency: bool = True,
    source_momentum_consistency: bool = True,
    bulk_air_phase_safe_source: bool = False,
    condensed_air_particle_diameter: float | None = None,
    source_total_energy_consistency: bool = False,
    source_pressure_thrust: bool = False,
    evaporation_zone_distance: bool = False,
    condensed_air_stationary_bound: bool = False,
    ground_layer_entrainment: bool = False,
    geometry_observation: str = "state",
) -> list:
    """Measured against modelled plume centre and vertical spread, per arc.

    The measurement here is a Gaussian fitted in height at each arc, which
    gives the centre and the spread directly rather than through a
    concentration. Judging a change on concentration alone is how four of the
    five comparison errors in this project happened.

    ``momentum_filter`` decides whether the wind-steered trials are dropped.
    The default remains off for exploratory geometry, while the current
    validation claim uses the same momentum-driven population as the
    concentration comparison (23 fits).

    The reduced table declares its Gaussian-width convention. Files made
    before that metadata was added used the e-folding parameter from
    ``exp(-(z/w)**2)`` but called it ``sigma_z``. Such legacy files are
    converted on read using ``sigma = w/sqrt(2)``. New files and the bundled
    reduction store the statistical standard deviation directly, matching
    JETPLU's ``exp(-0.5*(z/sigma_z)**2)`` profile.
    """
    import json

    if geometry_observation not in {"state", "sensor_fit"}:
        raise ValueError(
            "geometry observation must be 'state' or 'sensor_fit'"
        )

    data = json.loads(Path(table).read_text())
    width_definition = data.get("gaussian_width_definition", "")
    if width_definition.startswith("standard deviation"):
        measured_width_scale = 1.0
    elif not width_definition:
        # Backward compatibility for reductions made before the convention
        # was explicit. Their stored value is an e-folding width even though
        # its key is named ``sigma_z``.
        measured_width_scale = 1.0 / math.sqrt(2.0)
    else:
        raise ValueError(
            f"unsupported Gaussian width definition: {width_definition!r}"
        )
    rows = []
    for trial in data["trials"]:
        if trial["orientation"] != "horizontal":
            continue
        rate = trial[source] / 1000.0
        if rate <= 0.005:
            continue
        if momentum_filter and _exit_ratio(
            trial, rate, liquid_temperature_mode=liquid_temperature_mode
        ) < MOMENTUM_RATIO:
            continue
        fits = [f for f in trial["vertical_fits"] if f["well_constrained"]]
        if not fits:
            continue
        try:
            jp, y0 = hydrogen_jet(
                rate=rate, diameter=trial["orifice_mm"] / 1000.0,
                wind=trial["wind_ms"], height=trial["release_height_m"],
                ambient_temperature=trial["T_C"] + 273.15,
                relative_humidity=trial["RH_pct"],
                storage_pressure_barg=trial["tanker_barg"],
                storage_temperature=_liquid_source_temperature(
                    trial, liquid_temperature_mode
                ),
                wind_reference_height=trial["wind_ref_m"],
                corrections=corrections, ground_effect=ground_effect,
                houf_entrainment=houf_entrainment,
                source_table_consistency=source_table_consistency,
                source_momentum_consistency=source_momentum_consistency,
                bulk_air_phase_safe_source=bulk_air_phase_safe_source,
                condensed_air_particle_diameter=condensed_air_particle_diameter,
                source_total_energy_consistency=source_total_energy_consistency,
                source_pressure_thrust=source_pressure_thrust,
                evaporation_zone_distance=evaporation_zone_distance,
                condensed_air_stationary_bound=condensed_air_stationary_bound,
                ground_layer_entrainment=ground_layer_entrainment,
            )
        except ValueError as error:
            if source_pressure_thrust and str(error).startswith(
                "pressure-thrust source incompatible:"
            ):
                continue
            raise
        traj = Trajectory(jp.th.table, jp.run(y0, distmx=STEP, smax=REACH).rows)
        for fit in fits:
            state = traj.at(fit["x"])
            if state is None:
                continue
            modelled_centre = state.z
            modelled_sigma_z = state.sz
            model_profile_r_squared = math.nan
            if geometry_observation == "sensor_fit":
                model_fit = _fit_model_vertical_profile(
                    traj,
                    fit["x"],
                    _vertical_sensor_heights(trial, fit["x"]),
                )
                modelled_centre = model_fit["centre"]
                modelled_sigma_z = model_fit["sigma_z"]
                model_profile_r_squared = model_fit["r_squared"]
            row = {
                "trial": trial["trial"],
                "release_height": trial["release_height_m"],
                "x": fit["x"],
                "measured_centre": fit["centre"],
                "measured_sigma_z": fit["sigma_z"] * measured_width_scale,
                "modelled_centre": modelled_centre,
                "modelled_sigma_z": modelled_sigma_z,
            }
            if geometry_observation == "sensor_fit":
                row["model_profile_r_squared"] = model_profile_r_squared
                row["internal_modelled_centre"] = state.z
                row["internal_modelled_sigma_z"] = state.sz
            rows.append(row)
    return rows


def independent_energy_interfaces_from_reduced(
    table, *, energy_transport="total", houf_width_mapping="velocity",
    equilibrium_dry_air_condensation=True,
    nearfield_establishment="entrained_mass", ground_interaction="free",
    nearfield_energy_transport="total", source="flow_mean_gs",
    liquid_temperature_mode="tank_saturation", measured_source_table=None,
    measured_source_mode="full", measured_lh2_equilibrium_bound=False,
    hydrogen_spin_isomer="normal", fit_velocity_spreading=False,
    consistent_phase_ambient=False, downstream_thermodynamic_profile="density", progress=None,
) -> IndependentEnergyInterfaceValidation:
    """Run the frozen five-flux interface audit on seven PRESLHY releases.

    This stops before downstream integration. It tests whether separating
    centre density and H2 fraction supplies the missing energy degree of
    freedom without moving the 10D boundary or changing any compatibility
    threshold.
    """
    import json

    from ..lh2 import (
        audit_lh2_independent_energy_interface,
        run_lh2_crosswind_research,
    )

    data = json.loads(Path(table).read_text())
    measured_rows = _measured_pipe_source_rows(measured_source_table)
    selected = []
    interfaces = {}
    failures = {}
    for trial in data["trials"]:
        number = trial["trial"]
        rate, source_temperature, source_pressure_barg = _trial_source_inputs(
            trial,
            source=source,
            liquid_temperature_mode=liquid_temperature_mode,
            measured_rows=measured_rows,
            measured_source_mode=measured_source_mode,
        )
        if (
            number not in INDEPENDENT_ENERGY_TRIALS
            or trial["orientation"] != "horizontal"
            or rate <= 0.005
        ):
            continue
        setup, _initial = hydrogen_jet(
            rate=rate, diameter=trial["orifice_mm"] / 1000.0,
            wind=trial["wind_ms"], height=trial["release_height_m"],
            ambient_temperature=trial["T_C"] + 273.15,
            relative_humidity=trial["RH_pct"],
            storage_pressure_barg=trial["tanker_barg"],
            storage_temperature=source_temperature,
            source_upstream_pressure_barg=source_pressure_barg,
            wind_reference_height=trial["wind_ref_m"], corrections=True,
            hydrogen_spin_isomer=hydrogen_spin_isomer,
        )
        jet_source = _measured_lh2_equilibrium_source(
            enabled=measured_lh2_equilibrium_bound,
            setup=setup,
            trial=trial,
            rate=rate,
            source_temperature=source_temperature,
            source_pressure_barg=source_pressure_barg,
            measured_rows=measured_rows,
            measured_source_mode=measured_source_mode,
            hydrogen_spin_isomer=hydrogen_spin_isomer,
        )
        source_ratio = jet_source.velocity / setup.local_source_wind
        if measured_rows:
            momentum_length = _momentum_dominance_distance(
                jet_source,
                ambient_density=setup.th.table.rhoa,
                local_wind=setup.local_source_wind,
            )
            if momentum_length < jet_source.diameter:
                failures[number] = (
                    "source/crosswind momentum length is below one source diameter"
                )
                continue
            station = min(10.0 * jet_source.diameter, momentum_length)
        else:
            if source_ratio < MOMENTUM_RATIO:
                continue
            station = 10.0 * jet_source.diameter
        selected.append(number)
        if progress is not None:
            progress(
                f"trial {number}: independent-energy interface at "
                f"{station / jet_source.diameter:.3f}D"
            )
        try:
            coupled = run_lh2_crosswind_research(
                jet_source,
                wind=trial["wind_ms"],
                height=trial["release_height_m"],
                ambient_temperature=trial["T_C"] + 273.15,
                ambient_pressure=101325.0,
                relative_humidity=trial["RH_pct"],
                roughness=0.001,
                stability="D",
                averaging=60.0,
                wind_reference_height=trial["wind_ref_m"],
                maximum_nearfield_distance=station,
                handoff_distance=station,
                minimum_mass_fraction=7.0e-4,
                radial_points=41,
                maximum_step=min(0.001, jet_source.diameter / 20.0),
                relative_tolerance=2.0e-6,
                thermodynamic_closure="phase_manifold",
                nearfield_establishment=nearfield_establishment,
                nearfield_energy_transport=nearfield_energy_transport,
                crosswind_entrainment="local_shear",
                equilibrium_dry_air_condensation=(
                    equilibrium_dry_air_condensation
                ),
                hydrogen_spin_isomer=hydrogen_spin_isomer,
                consistent_phase_ambient=consistent_phase_ambient,
            )
            interface = audit_lh2_independent_energy_interface(
                coupled.near_field,
                coupled.handoff.model,
                streamline_distance=station,
                energy_transport=energy_transport,
                houf_width_mapping=houf_width_mapping,
                ground_interaction=ground_interaction,
                fit_velocity_spreading=fit_velocity_spreading,
                downstream_thermodynamic_profile=downstream_thermodynamic_profile,
            )
        except (RuntimeError, ValueError) as error:
            failures[number] = str(error)
            continue
        interfaces[number] = interface
        if not interface.accepted:
            failures[number] = "; ".join(interface.failure_reasons)
    return IndependentEnergyInterfaceValidation(
        selected_trials=selected,
        interfaces=interfaces,
        failures=failures,
    )


def independent_energy_from_reduced(
    table, *, maximum_step=0.02, geometry_observation="state",
    energy_transport="total", houf_width_mapping="velocity",
    equilibrium_dry_air_condensation=True,
    nearfield_establishment="entrained_mass", ground_interaction="free",
    nearfield_energy_transport="total", source="flow_mean_gs",
    liquid_temperature_mode="tank_saturation", measured_source_table=None,
    measured_source_mode="full", measured_lh2_equilibrium_bound=False,
    hydrogen_spin_isomer="normal", fit_velocity_spreading=False,
    consistent_phase_ambient=False, downstream_thermodynamic_profile="density", trial_filter=None, progress=None,
) -> CoupledFieldValidation:
    """Run the frozen seven-state PRESLHY downstream comparison."""
    import json

    from ..lh2 import (
        audit_lh2_independent_energy_interface,
        run_lh2_crosswind_research,
    )

    if maximum_step <= 0.0:
        raise ValueError("independent-energy maximum step must be positive")
    if geometry_observation not in {"state", "sensor_fit"}:
        raise ValueError(
            "geometry observation must be 'state' or 'sensor_fit'"
        )
    data = json.loads(Path(table).read_text())
    measured_rows = _measured_pipe_source_rows(measured_source_table)
    if trial_filter is None:
        trial_population = INDEPENDENT_ENERGY_TRIALS
    else:
        trial_population = tuple(int(value) for value in trial_filter)
        if (
            not trial_population
            or len(set(trial_population)) != len(trial_population)
            or not set(trial_population).issubset(INDEPENDENT_ENERGY_TRIALS)
        ):
            raise ValueError(
                "trial filter must be a non-empty unique subset of the "
                "pre-registered trials"
            )
    width_definition = data.get("gaussian_width_definition", "")
    if width_definition.startswith("standard deviation"):
        measured_width_scale = 1.0
    elif not width_definition:
        measured_width_scale = 1.0 / math.sqrt(2.0)
    else:
        raise ValueError(
            f"unsupported Gaussian width definition: {width_definition!r}"
        )

    candidate = NearField(corrections=True)
    candidate_vertical = []
    selected = []
    interfaces = {}
    failures = {}
    for trial in data["trials"]:
        number = trial["trial"]
        rate, source_temperature, source_pressure_barg = _trial_source_inputs(
            trial,
            source=source,
            liquid_temperature_mode=liquid_temperature_mode,
            measured_rows=measured_rows,
            measured_source_mode=measured_source_mode,
        )
        if (
            number not in trial_population
            or trial["orientation"] != "horizontal"
            or rate <= 0.005
        ):
            continue
        setup, _initial = hydrogen_jet(
            rate=rate, diameter=trial["orifice_mm"] / 1000.0,
            wind=trial["wind_ms"], height=trial["release_height_m"],
            ambient_temperature=trial["T_C"] + 273.15,
            relative_humidity=trial["RH_pct"],
            storage_pressure_barg=trial["tanker_barg"],
            storage_temperature=source_temperature,
            source_upstream_pressure_barg=source_pressure_barg,
            wind_reference_height=trial["wind_ref_m"], corrections=True,
            hydrogen_spin_isomer=hydrogen_spin_isomer,
        )
        jet_source = _measured_lh2_equilibrium_source(
            enabled=measured_lh2_equilibrium_bound,
            setup=setup,
            trial=trial,
            rate=rate,
            source_temperature=source_temperature,
            source_pressure_barg=source_pressure_barg,
            measured_rows=measured_rows,
            measured_source_mode=measured_source_mode,
            hydrogen_spin_isomer=hydrogen_spin_isomer,
        )
        source_ratio = jet_source.velocity / setup.local_source_wind
        if measured_rows:
            momentum_length = _momentum_dominance_distance(
                jet_source,
                ambient_density=setup.th.table.rhoa,
                local_wind=setup.local_source_wind,
            )
            if momentum_length < jet_source.diameter:
                failures[number] = (
                    "source/crosswind momentum length is below one source diameter"
                )
                continue
            station = min(10.0 * jet_source.diameter, momentum_length)
        else:
            if source_ratio < MOMENTUM_RATIO:
                continue
            station = 10.0 * jet_source.diameter
        selected.append(number)
        if progress is not None:
            progress(
                f"trial {number}: build independent-energy state at "
                f"{station / jet_source.diameter:.3f}D"
            )
        try:
            coupled = run_lh2_crosswind_research(
                jet_source,
                wind=trial["wind_ms"],
                height=trial["release_height_m"],
                ambient_temperature=trial["T_C"] + 273.15,
                ambient_pressure=101325.0,
                relative_humidity=trial["RH_pct"],
                roughness=0.001,
                stability="D",
                averaging=60.0,
                wind_reference_height=trial["wind_ref_m"],
                maximum_nearfield_distance=station,
                handoff_distance=station,
                minimum_mass_fraction=7.0e-4,
                radial_points=41,
                maximum_step=min(0.001, jet_source.diameter / 20.0),
                relative_tolerance=2.0e-6,
                thermodynamic_closure="phase_manifold",
                nearfield_establishment=nearfield_establishment,
                nearfield_energy_transport=nearfield_energy_transport,
                crosswind_entrainment="local_shear",
                equilibrium_dry_air_condensation=(
                    equilibrium_dry_air_condensation
                ),
                hydrogen_spin_isomer=hydrogen_spin_isomer,
                consistent_phase_ambient=consistent_phase_ambient,
            )
            interface = audit_lh2_independent_energy_interface(
                coupled.near_field,
                coupled.handoff.model,
                streamline_distance=station,
                energy_transport=energy_transport,
                houf_width_mapping=houf_width_mapping,
                ground_interaction=ground_interaction,
                fit_velocity_spreading=fit_velocity_spreading,
                downstream_thermodynamic_profile=downstream_thermodynamic_profile,
            )
            interfaces[number] = interface
            if not interface.accepted:
                failures[number] = "; ".join(interface.failure_reasons)
                continue
            observation_x = [
                sensor["x"] for sensor in trial["sensors"]
            ] + [
                fit["x"] for fit in trial["vertical_fits"]
                if fit["well_constrained"]
            ]
            last_x = max(observation_x)
            arc_extent = max(
                0.25,
                2.0 * max(last_x - float(interface.state[5]), 0.0) + 0.25,
            )
            result = interface.run(
                maximum_distance=arc_extent,
                maximum_step=maximum_step,
                relative_tolerance=1.0e-5,
            )
            interface.downstream_result = result
            trajectory = IndependentEnergyTrajectory(interface, result)
            if not trajectory.ok:
                raise RuntimeError("independent-energy trajectory is empty")
        except (RuntimeError, ValueError) as error:
            failures[number] = str(error)
            if progress is not None:
                progress(f"trial {number}: failed ({error})")
            continue

        handoff_x = float(interface.state[5])
        arcs: dict[float, list] = {}
        for sensor in trial["sensors"]:
            arcs.setdefault(round(sensor["x"], 3), []).append(sensor)
        for x, group in sorted(arcs.items()):
            if x + 1.0e-12 < handoff_x:
                continue
            state = trajectory.at(x)
            if state is None:
                continue
            observed = max(sensor["peak"] for sensor in group)
            predicted = max(
                trajectory.concentration_at(x, sensor["y"], sensor["z"])
                for sensor in group
            )
            if observed <= 0.0 or predicted <= 0.0:
                continue
            candidate.pairs.append(Pair(
                trial=number, x=x, observed=observed, predicted=predicted,
                observed_at=max(group, key=lambda item: item["peak"])["z"],
                predicted_at=max(sensor["z"] for sensor in group),
                plume_centre=state.z,
            ))

        for fit in trial["vertical_fits"]:
            if not fit["well_constrained"] or fit["x"] < handoff_x:
                continue
            state = trajectory.at(fit["x"])
            if state is None:
                continue
            modelled_centre = state.z
            modelled_sigma_z = state.sz
            model_profile_r_squared = math.nan
            if geometry_observation == "sensor_fit":
                model_fit = _fit_model_vertical_profile(
                    trajectory,
                    fit["x"],
                    _vertical_sensor_heights(trial, fit["x"]),
                )
                modelled_centre = model_fit["centre"]
                modelled_sigma_z = model_fit["sigma_z"]
                model_profile_r_squared = model_fit["r_squared"]
            row = {
                "trial": number,
                "release_height": trial["release_height_m"],
                "x": fit["x"],
                "measured_centre": fit["centre"],
                "measured_sigma_z": fit["sigma_z"] * measured_width_scale,
                "modelled_centre": modelled_centre,
                "modelled_sigma_z": modelled_sigma_z,
            }
            if geometry_observation == "sensor_fit":
                row["model_profile_r_squared"] = model_profile_r_squared
                row["internal_modelled_centre"] = state.z
                row["internal_modelled_sigma_z"] = state.sz
            candidate_vertical.append(row)
        if progress is not None:
            progress(
                f"trial {number}: scored, balance "
                f"{result.maximum_relative_balance_residual:.3e}"
            )

    expected = list(trial_population)
    if selected != expected:
        raise ValueError(
            "the reduced table no longer selects the pre-registered trials: "
            f"{selected!r} != {expected!r}"
        )
    baseline_all = from_reduced(
        table, corrections=True, source=source,
        liquid_temperature_mode=liquid_temperature_mode,
    )
    concentration_keys = {
        (pair.trial, pair.x) for pair in candidate.pairs
    }
    baseline = NearField(
        pairs=[
            pair for pair in baseline_all.pairs
            if (pair.trial, pair.x) in concentration_keys
        ],
        excluded=dict(baseline_all.excluded),
        corrections=True,
    )
    baseline_vertical_all = vertical(
        table,
        corrections=True,
        source=source,
        liquid_temperature_mode=liquid_temperature_mode,
        momentum_filter=True,
        geometry_observation=geometry_observation,
    )
    vertical_keys = {
        (row["trial"], row["x"]) for row in candidate_vertical
    }
    baseline_vertical = [
        row for row in baseline_vertical_all
        if (row["trial"], row["x"]) in vertical_keys
    ]
    return CoupledFieldValidation(
        baseline=baseline,
        candidate=candidate,
        baseline_vertical=baseline_vertical,
        candidate_vertical=candidate_vertical,
        selected_trials=selected,
        nearfield_establishment=nearfield_establishment,
        crosswind_entrainment=(
            f"local_shear_independent_{energy_transport}_energy_"
            f"houf_{houf_width_mapping}_width_"
            f"ground_{ground_interaction}_"
            f"nearfield_{nearfield_energy_transport}_energy_"
            f"{'spreading_fitted_' if fit_velocity_spreading else ''}"
            f"{f'source_{measured_source_mode}_' if measured_rows else ''}"
            f"hydrogen_{hydrogen_spin_isomer}_"
            f"dry_air_{'equilibrium' if equilibrium_dry_air_condensation else 'metastable'}"
        ),
        compatible_handoff_search=False,
        geometry_observation=geometry_observation,
        handoffs=interfaces,
        failures=failures,
    )


def coupled_from_reduced(
    table, *, nearfield_establishment="entrained_mass",
    crosswind_entrainment="source_momentum",
    compatible_handoff_search=False, progress=None,
) -> CoupledFieldValidation:
    """Run the frozen PRESLHY test of the conserved coupled research path.

    The complete protocol was written before sensor predictions were
    calculated and is retained in
    ``docs/prereg-preslhy-coupled-crosswind-validation.md``.  ``progress`` is
    an optional one-argument callable for long-running command-line use; it
    has no effect on the calculation.
    """
    import json

    from ..core.jetplume import J_X
    from ..lh2 import run_lh2_crosswind_research

    data = json.loads(Path(table).read_text())
    width_definition = data.get("gaussian_width_definition", "")
    if width_definition.startswith("standard deviation"):
        measured_width_scale = 1.0
    elif not width_definition:
        measured_width_scale = 1.0 / math.sqrt(2.0)
    else:
        raise ValueError(
            f"unsupported Gaussian width definition: {width_definition!r}"
        )

    candidate = NearField(corrections=True)
    candidate_vertical = []
    selected = []
    handoffs = {}
    failures = {}

    for trial in data["trials"]:
        number = trial["trial"]
        rate = trial["flow_mean_gs"] / 1000.0
        if (
            trial["orientation"] != "horizontal"
            or rate <= 0.005
            or _exit_ratio(trial, rate) < MOMENTUM_RATIO
        ):
            continue

        setup, _initial = hydrogen_jet(
            rate=rate, diameter=trial["orifice_mm"] / 1000.0,
            wind=trial["wind_ms"], height=trial["release_height_m"],
            ambient_temperature=trial["T_C"] + 273.15,
            relative_humidity=trial["RH_pct"],
            storage_pressure_barg=trial["tanker_barg"],
            wind_reference_height=trial["wind_ref_m"],
            corrections=True,
        )
        source = setup.axisymmetric_source
        source_ratio = source.velocity / setup.local_source_wind
        if source_ratio < MOMENTUM_RATIO:
            continue
        selected.append(number)
        if progress is not None:
            progress(
                f"trial {number}: source ratio {source_ratio:.3f}, "
                f"handoff {10.0 * source.diameter:.5f} m"
            )

        handoff_distance = 10.0 * source.diameter
        nearfield_distance = (
            20.0 * source.diameter
            if compatible_handoff_search else handoff_distance
        )
        try:
            coupled = run_lh2_crosswind_research(
                source,
                wind=trial["wind_ms"],
                height=trial["release_height_m"],
                ambient_temperature=trial["T_C"] + 273.15,
                ambient_pressure=101325.0,
                relative_humidity=trial["RH_pct"],
                roughness=0.001,
                stability="D",
                averaging=60.0,
                wind_reference_height=trial["wind_ref_m"],
                maximum_nearfield_distance=nearfield_distance,
                handoff_distance=(
                    None if compatible_handoff_search else handoff_distance
                ),
                handoff_search_start=(
                    handoff_distance if compatible_handoff_search else None
                ),
                handoff_search_step=(
                    source.diameter if compatible_handoff_search else None
                ),
                minimum_mass_fraction=7.0e-4,
                radial_points=41,
                maximum_step=min(0.001, source.diameter / 20.0),
                relative_tolerance=2.0e-6,
                thermodynamic_closure="phase_manifold",
                nearfield_establishment=nearfield_establishment,
                crosswind_entrainment=crosswind_entrainment,
            )
        except (RuntimeError, ValueError) as error:
            failures[number] = str(error)
            if progress is not None:
                progress(f"trial {number}: failed ({error})")
            continue
        handoffs[number] = coupled
        if not coupled.accepted:
            reasons = list(coupled.applicability_failures)
            reasons += coupled.handoff.failure_reasons
            if not coupled.near_field.conservative:
                reasons.append("near-field conservation screen failed")
            failures[number] = "; ".join(reasons) or "coupled screen failed"
            if progress is not None:
                progress(f"trial {number}: interface rejected")
            continue

        result = coupled.run(distmx=0.02, tol=1.0e-5, smax=40.0)
        trajectory = Trajectory(coupled.handoff.model.th.table, result.rows)
        handoff_x = float(coupled.handoff.state[J_X])
        arcs: dict[float, list] = {}
        for sensor in trial["sensors"]:
            arcs.setdefault(round(sensor["x"], 3), []).append(sensor)
        for x, group in sorted(arcs.items()):
            if x + 1.0e-12 < handoff_x:
                continue
            state = trajectory.at(x)
            if state is None:
                continue
            observed = max(sensor["peak"] for sensor in group)
            predicted = max(
                trajectory.concentration_at(
                    x, sensor["y"], sensor["z"]
                )
                for sensor in group
            )
            if observed <= 0.0 or predicted <= 0.0:
                continue
            candidate.pairs.append(Pair(
                trial=number, x=x, observed=observed, predicted=predicted,
                observed_at=max(group, key=lambda item: item["peak"])["z"],
                predicted_at=max(sensor["z"] for sensor in group),
                plume_centre=state.z,
            ))

        for fit in trial["vertical_fits"]:
            if not fit["well_constrained"] or fit["x"] < handoff_x:
                continue
            state = trajectory.at(fit["x"])
            if state is None:
                continue
            candidate_vertical.append({
                "trial": number,
                "release_height": trial["release_height_m"],
                "x": fit["x"],
                "measured_centre": fit["centre"],
                "measured_sigma_z": (
                    fit["sigma_z"] * measured_width_scale
                ),
                "modelled_centre": state.z,
                "modelled_sigma_z": state.sz,
            })
        if progress is not None:
            progress(f"trial {number}: accepted and scored")

    expected = [10, 11, 12, 22, 23, 24, 25]
    if selected != expected:
        raise ValueError(
            "the reduced table no longer selects the pre-registered trials: "
            f"{selected!r} != {expected!r}"
        )

    baseline_all = from_reduced(table, corrections=True)
    concentration_keys = {
        (pair.trial, pair.x) for pair in candidate.pairs
    }
    baseline = NearField(
        pairs=[
            pair for pair in baseline_all.pairs
            if (pair.trial, pair.x) in concentration_keys
        ],
        excluded=dict(baseline_all.excluded),
        corrections=True,
    )
    baseline_vertical_all = vertical(
        table, corrections=True, momentum_filter=True
    )
    vertical_keys = {
        (row["trial"], row["x"]) for row in candidate_vertical
    }
    baseline_vertical = [
        row for row in baseline_vertical_all
        if (row["trial"], row["x"]) in vertical_keys
    ]
    return CoupledFieldValidation(
        baseline=baseline,
        candidate=candidate,
        baseline_vertical=baseline_vertical,
        candidate_vertical=candidate_vertical,
        selected_trials=selected,
        nearfield_establishment=nearfield_establishment,
        crosswind_entrainment=crosswind_entrainment,
        compatible_handoff_search=compatible_handoff_search,
        geometry_observation="state",
        handoffs=handoffs,
        failures=failures,
    )
