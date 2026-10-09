"""Explicit bridge from field release contracts to the LH2 flash source.

The bridge deliberately owns only the source-plane translation.  It does not
choose a pool, droplet-evaporation, or atmospheric-dispersion closure.  Those
remain separate decisions because they need different evidence and validation.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
import math
from typing import TYPE_CHECKING

from .field_contracts import (
    BoundedValue,
    FieldApplicability,
    FieldScenario,
    FieldValidationEvidence,
    PressureDrivenMassFlowBoundary,
    ReleaseSource,
    assess_field_applicability,
)
from .semi_fv_obstacle import SourceRateSchedule

if TYPE_CHECKING:
    from .lh2_droplets import FlashingHydrogenDropletSource
    from .lh2_property_table import LH2SaturationTable


_EQUILIBRIUM_FLASH_MODELS = {
    "equilibrium",
    "hem",
    "homogeneous_equilibrium",
}

# The underlying post-flash source already reports normalized mass,
# momentum, and energy residuals.  Keep the field handoff gate explicit and
# independent of any downstream transport tolerance so a source that is only
# approximately closed cannot be promoted by the scalar solver.
_FIELD_FLASH_MASS_REL_TOL = 1.0e-10
_FIELD_FLASH_MOMENTUM_REL_TOL = 1.0e-8
_FIELD_FLASH_ENERGY_REL_TOL = 1.0e-8


def _raw_flash_closure_warnings(
    flash: object,
    expected_mass_flow_kg_s: float,
) -> tuple[str, ...]:
    """Check a raw droplet flash before a schedule adapter consumes it."""
    try:
        expected = float(expected_mass_flow_kg_s)
        vapour = float(getattr(flash, "vapour_mass_flow"))
        liquid = float(getattr(flash, "liquid_mass_flow"))
        phase_residual = abs(expected - vapour - liquid) / max(expected, 1.0)
        values = {
            "phase mass-flow closure": phase_residual,
            "mass residual": float(getattr(flash, "mass_residual")),
            "momentum residual": float(getattr(flash, "momentum_residual")),
            "energy residual": float(getattr(flash, "energy_residual")),
        }
    except (AttributeError, TypeError, ValueError):
        return ("flash closure diagnostics are unavailable or non-numeric",)
    warnings = []
    if not math.isfinite(expected) or expected < 0.0:
        warnings.append("expected flash mass flow is invalid")
    if not math.isfinite(vapour) or vapour < 0.0:
        warnings.append("flash vapour mass flow is invalid or negative")
    if not math.isfinite(liquid) or liquid < 0.0:
        warnings.append("flash liquid mass flow is invalid or negative")
    if not math.isfinite(phase_residual) or phase_residual > _FIELD_FLASH_MASS_REL_TOL:
        warnings.append("flash phase mass-flow closure residual exceeds the handoff tolerance")
    if not math.isfinite(values["mass residual"]) or values["mass residual"] > _FIELD_FLASH_MASS_REL_TOL:
        warnings.append("flash mass residual exceeds the handoff tolerance")
    if not math.isfinite(values["momentum residual"]) or values["momentum residual"] > _FIELD_FLASH_MOMENTUM_REL_TOL:
        warnings.append("flash momentum residual exceeds the handoff tolerance")
    if not math.isfinite(values["energy residual"]) or values["energy residual"] > _FIELD_FLASH_ENERGY_REL_TOL:
        warnings.append("flash energy residual exceeds the handoff tolerance")
    return tuple(warnings)


@dataclass(frozen=True)
class FieldLH2FlashResult:
    """Auditable mapping from one nominal field release to a flash plane."""

    release: ReleaseSource
    flash: "FlashingHydrogenDropletSource"
    effective_area_m2: float
    effective_diameter_m: float
    model: str
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.release, ReleaseSource):
            raise TypeError("field LH2 flash result release must be a ReleaseSource")
        from .lh2_droplets import FlashingHydrogenDropletSource

        if not isinstance(self.flash, FlashingHydrogenDropletSource):
            raise TypeError(
                "field LH2 flash result flash must be a FlashingHydrogenDropletSource"
            )
        for name, value in {
            "effective_area_m2": self.effective_area_m2,
            "effective_diameter_m": self.effective_diameter_m,
        }.items():
            if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not isinstance(self.model, str) or self.model.strip().lower() not in _EQUILIBRIUM_FLASH_MODELS:
            raise ValueError("field LH2 flash result model is unsupported")
        for name, value, positive in (
            ("flash.mass_flow", self.flash.mass_flow, True),
            ("flash.vapour_mass_flow", self.flash.vapour_mass_flow, False),
            ("flash.liquid_mass_flow", self.flash.liquid_mass_flow, False),
            ("flash.mass_residual", self.flash.mass_residual, False),
            ("flash.momentum_residual", self.flash.momentum_residual, False),
            ("flash.energy_residual", self.flash.energy_residual, False),
        ):
            if isinstance(value, bool) or not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite and numeric")
            if (float(value) <= 0.0 if positive else float(value) < 0.0):
                qualifier = "positive" if positive else "non-negative"
                raise ValueError(f"{name} must be {qualifier}")
        if any(not isinstance(item, str) or not item.strip() for item in self.warnings):
            raise ValueError("field LH2 flash result warnings must be non-empty strings")

    @property
    def mass_residual_kg_s(self) -> float:
        return abs(
            self.release.mass_flow_kg_s.nominal
            - self.flash.vapour_mass_flow
            - self.flash.liquid_mass_flow
        )

    @property
    def momentum_residual(self) -> float:
        """Normalized post-flash momentum residual from the source bridge."""
        return float(self.flash.momentum_residual)

    @property
    def energy_residual(self) -> float:
        """Normalized post-flash total-specific-energy residual."""
        return float(self.flash.energy_residual)

    @property
    def phase_mass_residual(self) -> float:
        """Dimensionless mass closure residual reported by the flash source."""
        return float(self.flash.mass_residual)

    @property
    def closure_warnings(self) -> tuple[str, ...]:
        """Return explicit handoff-closure failures without inventing a fix."""
        warnings: list[str] = []
        scale = max(self.release.mass_flow_kg_s.nominal, 1.0)
        if not math.isfinite(self.mass_residual_kg_s) or self.mass_residual_kg_s > _FIELD_FLASH_MASS_REL_TOL * scale:
            warnings.append(
                "field LH2 flash mass partition residual exceeds the handoff tolerance"
            )
        if not math.isfinite(self.phase_mass_residual) or self.phase_mass_residual > _FIELD_FLASH_MASS_REL_TOL:
            warnings.append(
                "field LH2 flash phase mass-flow closure residual exceeds the handoff tolerance"
            )
        if not math.isfinite(self.momentum_residual) or self.momentum_residual > _FIELD_FLASH_MOMENTUM_REL_TOL:
            warnings.append(
                "field LH2 flash momentum residual exceeds the handoff tolerance"
            )
        if not math.isfinite(self.energy_residual) or self.energy_residual > _FIELD_FLASH_ENERGY_REL_TOL:
            warnings.append(
                "field LH2 flash energy residual exceeds the handoff tolerance"
            )
        return tuple(warnings)

    @property
    def closure_tolerances(self) -> dict[str, float]:
        """Return the explicit tolerances used by the field handoff gate."""
        return {
            "mass_partition_relative": _FIELD_FLASH_MASS_REL_TOL,
            "phase_mass_relative": _FIELD_FLASH_MASS_REL_TOL,
            "momentum_relative": _FIELD_FLASH_MOMENTUM_REL_TOL,
            "energy_relative": _FIELD_FLASH_ENERGY_REL_TOL,
        }

    @property
    def closure_diagnostics(self) -> tuple[tuple[str, float | None], ...]:
        """Return JSON-safe residual values for a blocked-source audit."""
        values = (
            ("mass_partition_residual_kg_s", self.mass_residual_kg_s),
            ("phase_mass_residual", self.phase_mass_residual),
            ("momentum_residual", self.momentum_residual),
            ("energy_residual", self.energy_residual),
        )
        return tuple(
            (name, float(value) if math.isfinite(float(value)) else None)
            for name, value in values
        )

    @property
    def conservative(self) -> bool:
        return not self.closure_warnings


@dataclass(frozen=True)
class FieldLH2SourcePreparation:
    """Failure-safe source preparation result for a field scenario.

    A blocked result deliberately contains no flash plane. Callers can persist
    the applicability reasons alongside a plant-case report without catching a
    source-physics exception and accidentally continuing with stale values.
    """

    scenario: FieldScenario
    applicability: FieldApplicability
    flash_result: FieldLH2FlashResult | None
    flash_closure_diagnostics: tuple[tuple[str, float | None], ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.scenario, FieldScenario):
            raise TypeError("field LH2 source preparation scenario must be a FieldScenario")
        if not isinstance(self.applicability, FieldApplicability):
            raise TypeError(
                "field LH2 source preparation applicability must be a FieldApplicability"
            )
        if self.flash_result is not None and not isinstance(
            self.flash_result, FieldLH2FlashResult
        ):
            raise TypeError(
                "field LH2 source preparation flash_result must be a FieldLH2FlashResult or None"
            )
        if self.applicability.status == "blocked" and self.flash_result is not None:
            raise ValueError("blocked field LH2 source preparation cannot expose a flash plane")
        diagnostics = tuple(self.flash_closure_diagnostics)
        for item in diagnostics:
            if not isinstance(item, tuple) or len(item) != 2:
                raise TypeError("flash closure diagnostics must contain name/value pairs")
            name, value = item
            if not isinstance(name, str) or not name.strip():
                raise ValueError("flash closure diagnostic names must be non-empty strings")
            if value is not None and (
                isinstance(value, bool) or not math.isfinite(float(value)) or float(value) < 0.0
            ):
                raise ValueError("flash closure diagnostic values must be finite and non-negative")
        if len({name for name, _value in diagnostics}) != len(diagnostics):
            raise ValueError("flash closure diagnostic names must be unique")
        object.__setattr__(self, "flash_closure_diagnostics", diagnostics)

    @property
    def prepared(self) -> bool:
        return self.flash_result is not None and self.applicability.status != "blocked"


@dataclass(frozen=True)
class FieldBlowdownVapourSchedule:
    """Direct post-flash vapour schedule reconstructed from a blowdown result.

    Liquid carried by a flash plane is deliberately kept out of this schedule;
    it must go through a rainout/pool path with its own atmospheric launch
    boundary. The schedule is therefore suitable only for the direct-vapour
    branch of a local field transport screen.
    """

    schedule: SourceRateSchedule
    total_discharged_mass_kg: float
    direct_vapour_mass_kg: float
    unrouted_postflash_liquid_mass_kg: float
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        for name, value in {
            "total_discharged_mass_kg": self.total_discharged_mass_kg,
            "direct_vapour_mass_kg": self.direct_vapour_mass_kg,
            "unrouted_postflash_liquid_mass_kg": self.unrouted_postflash_liquid_mass_kg,
        }.items():
            if not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
        if self.direct_vapour_mass_kg > self.total_discharged_mass_kg + 1.0e-10:
            raise ValueError("direct vapour mass cannot exceed discharged mass")
        ledger_tolerance = 1.0e-8 * max(self.total_discharged_mass_kg, 1.0)
        if abs(
            self.direct_vapour_mass_kg
            + self.unrouted_postflash_liquid_mass_kg
            - self.total_discharged_mass_kg
        ) > ledger_tolerance:
            raise ValueError(
                "direct vapour and unrouted liquid masses must close the discharged-mass ledger"
            )

    @property
    def direct_vapour_fraction(self) -> float:
        return self.direct_vapour_mass_kg / max(self.total_discharged_mass_kg, 1.0e-30)


def direct_vapour_schedule_from_cryogenic_blowdown(
    result: "CryogenicBlowdownResult",
    *,
    ambient_temperature_k: float = 295.0,
    droplet_size_coefficient: float = 15.0,
    source_id: str = "cryogenic-blowdown:direct-postflash-vapour",
) -> FieldBlowdownVapourSchedule:
    """Map every resolved blowdown interval to its direct post-flash vapour rate.

    Single-phase tank intervals use the existing pressure-thrust flash adapter.
    HEM two-phase intervals are admitted only when the original blowdown
    explicitly selected homogeneous withdrawal. Vapour-withdrawal two-phase
    states are withheld because this project has no corresponding pressure-
    thrust source adapter. No homogeneous complete-evaporation assumption is
    used here.
    """
    from .cryogenic_blowdown import (
        CryogenicBlowdownResult,
        blowdown_state_to_flashing_droplet_source,
        hem_blowdown_state_to_flashing_droplet_source,
    )

    if not isinstance(result, CryogenicBlowdownResult):
        raise TypeError("result must be a CryogenicBlowdownResult")
    if not isinstance(source_id, str) or not source_id.strip():
        raise ValueError("source_id must be non-empty")
    if not math.isfinite(float(ambient_temperature_k)) or ambient_temperature_k <= 0.0:
        raise ValueError("ambient_temperature_k must be positive and finite")
    if not math.isfinite(float(droplet_size_coefficient)) or droplet_size_coefficient <= 0.0:
        raise ValueError("droplet_size_coefficient must be positive and finite")

    times: list[float] = []
    rates: list[float] = []
    total_mass = 0.0
    used_hem = False
    states = result.states
    for state, next_state in zip(states, states[1:]):
        start, end = float(state.time_s), float(next_state.time_s)
        if end < start:
            raise ValueError("blowdown result state times are not ordered")
        if math.isclose(start, end, abs_tol=1.0e-12):
            continue
        if not times:
            times.append(start)
        elif not math.isclose(times[-1], start, abs_tol=1.0e-12):
            raise ValueError("blowdown result state intervals are discontinuous")
        flash = None
        if state.mass_flow_kg_s <= 0.0:
            vapour_rate = 0.0
        elif state.vapour_quality is None:
            flash = blowdown_state_to_flashing_droplet_source(
                result.config, state,
                ambient_temperature_k=ambient_temperature_k,
                droplet_size_coefficient=droplet_size_coefficient,
            )
            vapour_rate = flash.vapour_mass_flow
        elif result.config.two_phase_withdrawal == "homogeneous":
            flash = hem_blowdown_state_to_flashing_droplet_source(
                result.config, state,
                ambient_temperature_k=ambient_temperature_k,
                droplet_size_coefficient=droplet_size_coefficient,
            )
            vapour_rate = flash.vapour_mass_flow
            used_hem = True
        else:
            raise ValueError(
                "two-phase blowdown interval has no declared HEM pressure-thrust flash adapter"
            )
        if flash is not None:
            closure_warnings = _raw_flash_closure_warnings(
                flash, state.mass_flow_kg_s,
            )
            if closure_warnings:
                raise ValueError(
                    "blowdown flash interval is not conservatively closed: "
                    + "; ".join(closure_warnings)
                )
        rates.append(float(vapour_rate))
        times.append(end)
        total_mass += float(state.mass_flow_kg_s) * (end - start)
    if len(times) < 2 or not rates:
        raise ValueError("blowdown result contains no positive-duration source interval")
    schedule = SourceRateSchedule(
        tuple(times), tuple(rates) + (0.0,), source_id=source_id,
    )
    direct_mass = schedule.released_mass_kg
    ledger_tolerance = 1.0e-8 * max(total_mass, 1.0)
    if direct_mass > total_mass + ledger_tolerance:
        raise ValueError(
            "direct post-flash vapour mass exceeds the blowdown discharged-mass ledger"
        )
    unrouted_liquid = total_mass - direct_mass
    warnings = [
        "schedule contains direct post-flash vapour only; post-flash liquid is excluded from atmospheric scalar transport",
        "blowdown is a well-mixed source model; its declared vessel, nozzle, wall and withdrawal assumptions remain active",
    ]
    if used_hem:
        warnings.append(
            "two-phase intervals use the explicitly selected HEM pressure-thrust flash; no finite-rate evaporation was assumed",
        )
    return FieldBlowdownVapourSchedule(
        schedule=schedule,
        total_discharged_mass_kg=total_mass,
        direct_vapour_mass_kg=direct_mass,
        unrouted_postflash_liquid_mass_kg=unrouted_liquid,
        warnings=tuple(warnings),
    )


def _uncertainty_warning(release: ReleaseSource) -> tuple[str, ...]:
    names = tuple(
        name for name, value in release.uncertainty_fields().items()
        if not value.is_exact
    )
    if not names:
        return ()
    return (
        "nominal flash only; evaluate release uncertainty through "
        f"FieldScenario corners ({', '.join(names)})",
    )


def pressure_driven_lh2_mass_flow(
    release: ReleaseSource,
    *,
    ambient_pressure_pa: float = 101325.0,
    ambient_pressure_bounds_pa: tuple[float, float] | None = None,
    source_id: str = "pressure-driven-orifice",
    allow_supercritical_gas: bool = False,
) -> BoundedValue:
    """Derive a conservative mass-flow bound from an explicit leak boundary.

    The adapter uses the existing homogeneous-equilibrium isentropic throat
    closure and applies the declared discharge coefficient to the ideal
    mass-flux times geometric area.  Upstream pressure, temperature, opening
    area, discharge coefficient, and optional ambient-pressure endpoints are
    crossed as deterministic corners; no probability or fitted leak law is
    introduced.  A caller must attach the returned bound to a
    :class:`ReleaseSource` before invoking the field flash bridge.

    Supercritical gas use is opt-in and only allowed when the declared source
    contains no liquid fraction.  Mixed or liquid corners that cannot be
    represented by the throat closure fail closed.
    """
    if not isinstance(release, ReleaseSource):
        raise TypeError("release must be a ReleaseSource")
    fluid = release.fluid.strip().lower()
    if fluid not in {"hydrogen", "h2", "lh2"}:
        raise ValueError(
            "pressure-driven LH2 mass-flow adapter supports hydrogen releases only"
        )
    if release.pressure_reference != "absolute":
        raise ValueError(
            "pressure-driven LH2 mass-flow adapter requires absolute upstream pressure"
        )
    if not isinstance(source_id, str) or not source_id.strip() or source_id.strip().lower() == "unspecified":
        raise ValueError("pressure-driven mass-flow source_id must be explicitly declared")
    if not isinstance(allow_supercritical_gas, bool):
        raise TypeError("allow_supercritical_gas must be boolean")
    if allow_supercritical_gas and release.liquid_fraction.upper > 0.0:
        raise ValueError(
            "supercritical gas mass-flow derivation requires liquid_fraction upper bound zero"
        )
    if isinstance(ambient_pressure_pa, bool) or not math.isfinite(float(ambient_pressure_pa)) or ambient_pressure_pa <= 0.0:
        raise ValueError("ambient_pressure_pa must be positive and finite")
    ambient_pressure = float(ambient_pressure_pa)
    if ambient_pressure_bounds_pa is None:
        ambient_bounds = (ambient_pressure,)
    else:
        if (
            not isinstance(ambient_pressure_bounds_pa, tuple)
            or len(ambient_pressure_bounds_pa) != 2
            or any(
                isinstance(value, bool)
                or not math.isfinite(float(value))
                or float(value) <= 0.0
                for value in ambient_pressure_bounds_pa
            )
            or ambient_pressure_bounds_pa[0] > ambient_pressure_bounds_pa[1]
            or not ambient_pressure_bounds_pa[0] <= ambient_pressure <= ambient_pressure_bounds_pa[1]
        ):
            raise ValueError(
                "ambient_pressure_bounds_pa must be an ordered positive pair containing ambient_pressure_pa"
            )
        ambient_bounds = tuple(float(value) for value in ambient_pressure_bounds_pa)

    from .notional import isentropic_throat

    def evaluate(
        pressure_pa: float,
        temperature_k: float,
        area_m2: float,
        discharge_coefficient: float,
        ambient_pa: float,
    ) -> float:
        try:
            throat = isentropic_throat(
                fluid="Hydrogen",
                storage_temperature=float(temperature_k),
                storage_pressure=float(pressure_pa),
                ambient_pressure=float(ambient_pa),
                allow_supercritical_gas=allow_supercritical_gas,
            )
        except (RuntimeError, ValueError) as error:
            raise ValueError(
                "pressure-driven LH2 mass-flow corner cannot be evaluated by the declared throat closure"
            ) from error
        rate = float(discharge_coefficient) * float(area_m2) * float(throat.mass_flux)
        if not math.isfinite(rate) or rate <= 0.0:
            raise ValueError(
                "pressure-driven LH2 mass-flow corner produced a non-positive or non-finite rate"
            )
        return rate

    nominal = evaluate(
        release.upstream_pressure.nominal,
        release.upstream_temperature.nominal,
        release.opening_area_m2.nominal,
        release.discharge_coefficient.nominal,
        ambient_pressure,
    )
    corners = []
    fields = (
        release.upstream_pressure,
        release.upstream_temperature,
        release.opening_area_m2,
        release.discharge_coefficient,
    )
    for pressure_pa, temperature_k, area_m2, discharge_coefficient, ambient_pa in product(
        *(value.corners() for value in fields), ambient_bounds
    ):
        corners.append(
            evaluate(
                pressure_pa, temperature_k, area_m2,
                discharge_coefficient, ambient_pa,
            )
        )
    return BoundedValue(
        nominal,
        min(nominal, *corners),
        max(nominal, *corners),
        "kg/s",
        source_id.strip(),
    )


def release_with_pressure_driven_lh2_mass_flow(
    release: ReleaseSource,
    *,
    ambient_pressure_pa: float = 101325.0,
    ambient_pressure_bounds_pa: tuple[float, float] | None = None,
    source_id: str = "pressure-driven-orifice",
    allow_supercritical_gas: bool = False,
) -> ReleaseSource:
    """Return ``release`` with a derived pressure-driven mass-flow bound.

    A non-zero mass-flow boundary is never silently overwritten; callers must
    choose whether a measured-rate source or this pressure-driven source is
    authoritative.
    """
    if not isinstance(release, ReleaseSource):
        raise TypeError("release must be a ReleaseSource")
    if release.mass_flow_kg_s.lower != 0.0 or release.mass_flow_kg_s.upper != 0.0:
        raise ValueError(
            "pressure-driven mass-flow adapter will not overwrite a non-zero declared mass-flow boundary"
        )
    mass_flow = pressure_driven_lh2_mass_flow(
        release,
        ambient_pressure_pa=ambient_pressure_pa,
        ambient_pressure_bounds_pa=ambient_pressure_bounds_pa,
        source_id=source_id,
        allow_supercritical_gas=allow_supercritical_gas,
    )
    metadata = dict(release.metadata)
    provenance_bounds = (
        (float(ambient_pressure_bounds_pa[0]), float(ambient_pressure_bounds_pa[1]))
        if ambient_pressure_bounds_pa is not None
        else (float(ambient_pressure_pa), float(ambient_pressure_pa))
    )
    derivation = {
        "mass_flow_derivation": "pressure_driven_homogeneous_equilibrium_throat",
        "mass_flow_derivation_source_id": source_id.strip(),
        "mass_flow_derivation_ambient_pressure_pa": f"{float(ambient_pressure_pa):.17g}",
        "mass_flow_derivation_ambient_pressure_bounds_pa": ",".join(
            f"{value:.17g}" for value in provenance_bounds
        ),
        "mass_flow_derivation_allow_supercritical_gas": str(
            allow_supercritical_gas
        ).lower(),
    }
    for key, value in derivation.items():
        if key in metadata and metadata[key] != value:
            raise ValueError(
                f"release metadata {key!r} conflicts with the pressure-driven derivation"
            )
        metadata[key] = value
    return replace(
        release,
        mass_flow_kg_s=mass_flow,
        metadata=metadata,
        pressure_driven_mass_flow=PressureDrivenMassFlowBoundary(
            source_id=source_id.strip(),
            ambient_pressure_pa=BoundedValue(
                float(ambient_pressure_pa),
                provenance_bounds[0],
                provenance_bounds[1],
                "Pa",
                "pressure-driven-ambient",
            ),
            allow_supercritical_gas=allow_supercritical_gas,
        ),
    )


def lh2_flash_source_from_release(
    release: ReleaseSource,
    *,
    ambient_temperature_k: float = 295.0,
    ambient_pressure_pa: float = 101325.0,
    droplet_size_coefficient: float = 15.0,
    property_table: "LH2SaturationTable | None" = None,
) -> FieldLH2FlashResult:
    """Build the explicit equilibrium flash plane for a field LH2 release.

    ``opening_area_m2`` is the geometric opening.  The function applies the
    declared discharge coefficient exactly once by using
    ``A_eff = Cd * A_geometric``; this mirrors the existing blowdown adapter
    and prevents a common double-``Cd`` source-rate error.  A source with a
    post-flash phase fraction, gauge pressure, an unresolved phase basis or a
    non-equilibrium flash model is rejected rather than reinterpreted.
    """
    fluid = release.fluid.strip().lower()
    if fluid not in {"hydrogen", "h2", "lh2"}:
        raise ValueError("LH2 flash bridge supports hydrogen releases only")
    if release.pressure_reference != "absolute":
        raise ValueError("LH2 flash bridge requires absolute upstream pressure")
    if release.liquid_fraction_basis != "upstream":
        raise ValueError(
            "LH2 flash bridge requires an upstream liquid fraction; "
            "post-flash fraction cannot be used as the flash input"
        )
    model = release.flash_model.strip().lower()
    if model not in _EQUILIBRIUM_FLASH_MODELS:
        accepted = ", ".join(sorted(_EQUILIBRIUM_FLASH_MODELS))
        raise ValueError(
            f"unsupported flash_model={release.flash_model!r}; choose one of {accepted}"
        )
    values = {
        "ambient_temperature_k": ambient_temperature_k,
        "ambient_pressure_pa": ambient_pressure_pa,
        "droplet_size_coefficient": droplet_size_coefficient,
    }
    for name, value in values.items():
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if release.mass_flow_kg_s.nominal <= 0.0:
        raise ValueError("LH2 flash bridge requires a positive nominal mass flow")
    if release.upstream_pressure.nominal <= ambient_pressure_pa:
        raise ValueError("upstream pressure must exceed ambient pressure")

    effective_area = (
        release.opening_area_m2.nominal
        * release.discharge_coefficient.nominal
    )
    effective_diameter = math.sqrt(4.0 * effective_area / math.pi)
    liquid_fraction = release.liquid_fraction.nominal
    upstream_quality = None if liquid_fraction <= 0.0 else 1.0 - liquid_fraction

    from .lh2_droplets import flashing_hydrogen_droplet_source

    flash = flashing_hydrogen_droplet_source(
        mass_flow=release.mass_flow_kg_s.nominal,
        orifice_diameter=effective_diameter,
        upstream_temperature=release.upstream_temperature.nominal,
        upstream_pressure=release.upstream_pressure.nominal,
        ambient_temperature=ambient_temperature_k,
        ambient_pressure=ambient_pressure_pa,
        upstream_quality=upstream_quality,
        droplet_size_coefficient=droplet_size_coefficient,
        property_table=property_table,
    )
    warnings = _uncertainty_warning(release) + (
        "post-flash liquid is retained explicitly; this bridge does not assume "
        "instantaneous droplet evaporation, pool formation, or atmospheric mixing",
    )
    return FieldLH2FlashResult(
        release=release,
        flash=flash,
        effective_area_m2=effective_area,
        effective_diameter_m=effective_diameter,
        model=model,
        warnings=warnings,
    )


def prepare_field_lh2_flash(
    scenario: FieldScenario,
    *,
    ambient_temperature_k: float = 295.0,
    ambient_pressure_pa: float = 101325.0,
    droplet_size_coefficient: float = 15.0,
    lh2_validation_available: bool = False,
    validation_evidence: FieldValidationEvidence | None = None,
    property_table: "LH2SaturationTable | None" = None,
) -> FieldLH2SourcePreparation:
    """Prepare a nominal LH2 flash source with structured fail-safe output."""
    applicability = assess_field_applicability(
        scenario,
        lh2_validation_available=lh2_validation_available,
        validation_evidence=validation_evidence,
    )
    if applicability.status == "blocked":
        return FieldLH2SourcePreparation(scenario, applicability, None)
    try:
        result = lh2_flash_source_from_release(
            scenario.source,
            ambient_temperature_k=ambient_temperature_k,
            ambient_pressure_pa=ambient_pressure_pa,
            droplet_size_coefficient=droplet_size_coefficient,
            property_table=property_table,
        )
    except ValueError as error:
        blocked = FieldApplicability(
            "blocked",
            reasons=applicability.reasons + (str(error),),
            warnings=applicability.warnings,
            uncertainty_complete=False,
        )
        return FieldLH2SourcePreparation(scenario, blocked, None)
    if not result.conservative:
        blocked = FieldApplicability(
            "blocked",
            reasons=applicability.reasons + result.closure_warnings,
            warnings=applicability.warnings + result.warnings,
            uncertainty_complete=False,
        )
        # Do not expose an unclosed flash plane to downstream callers.  The
        # residual reasons remain in the applicability record for audit, while
        # the preparation object follows the same no-plane contract as every
        # other blocked source path.
        return FieldLH2SourcePreparation(
            scenario, blocked, None, result.closure_diagnostics,
        )
    return FieldLH2SourcePreparation(scenario, applicability, result)


def build_lh2_saturation_table_for_temperature_bounds(
    minimum_temperature_k: float,
    maximum_temperature_k: float,
    *,
    ambient_pressure_pa: float = 101325.0,
    ambient_pressure_bounds_pa: tuple[float, float] | None = None,
    nodes: int = 161,
    margin_k: float = 0.25,
) -> "LH2SaturationTable":
    """Build a bounded table covering declared LH2 source-temperature bounds.

    This table is intended for a repeated uncertainty sweep.  Its construction
    still uses CoolProp once, but subsequent saturated LH2 calls in
    :func:`lh2_flash_source_from_release` use interpolation. The caller must
    supply every temperature that could be queried: an out-of-domain query is
    still an error rather than a fallback to CoolProp. When ambient pressure
    bounds are declared, both endpoint saturation temperatures are included in
    the table span.
    """
    from CoolProp.CoolProp import PropsSI

    if not all(math.isfinite(float(value)) for value in (
        minimum_temperature_k, maximum_temperature_k,
    )):
        raise ValueError("table source temperatures must be finite")
    if minimum_temperature_k > maximum_temperature_k:
        raise ValueError("minimum_temperature_k cannot exceed maximum_temperature_k")
    if not math.isfinite(float(ambient_pressure_pa)) or ambient_pressure_pa <= 0.0:
        raise ValueError("ambient_pressure_pa must be positive and finite")
    if ambient_pressure_bounds_pa is None:
        pressure_values = (float(ambient_pressure_pa),)
    else:
        if (
            not isinstance(ambient_pressure_bounds_pa, tuple)
            or len(ambient_pressure_bounds_pa) != 2
            or any(
                not math.isfinite(float(value)) or float(value) <= 0.0
                for value in ambient_pressure_bounds_pa
            )
            or ambient_pressure_bounds_pa[0] > ambient_pressure_bounds_pa[1]
        ):
            raise ValueError(
                "ambient_pressure_bounds_pa must be an ordered positive finite pair"
            )
        pressure_values = tuple(float(value) for value in ambient_pressure_bounds_pa)
    if not math.isfinite(float(margin_k)) or margin_k <= 0.0:
        raise ValueError("margin_k must be positive and finite")
    ambient_saturation_temperatures = tuple(
        float(PropsSI("T", "P", pressure, "Q", 0, "Hydrogen"))
        for pressure in pressure_values
    )
    triple = float(PropsSI("Ttriple", "Hydrogen"))
    critical = float(PropsSI("Tcrit", "Hydrogen"))
    lower = min(minimum_temperature_k, *ambient_saturation_temperatures) - margin_k
    upper = max(maximum_temperature_k, *ambient_saturation_temperatures) + margin_k
    lower = max(lower, triple + 1.0e-4)
    upper = min(upper, critical - 1.0e-4)
    if lower >= upper:
        raise ValueError("release temperatures cannot define a valid LH2 saturation table")
    from .lh2_property_table import LH2SaturationTable

    return LH2SaturationTable.build(lower, upper, nodes=nodes)


def build_lh2_saturation_table_for_release(
    release: ReleaseSource,
    *,
    ambient_pressure_pa: float = 101325.0,
    ambient_pressure_bounds_pa: tuple[float, float] | None = None,
    nodes: int = 161,
    margin_k: float = 0.25,
) -> "LH2SaturationTable":
    """Build a bounded saturation table covering every release-temperature corner.

    Only the normal-hydrogen field bridge is currently table-enabled;
    para/ortho cases retain their direct, explicit CoolProp path.
    """
    if release.fluid.strip().lower() not in {"hydrogen", "h2", "lh2"}:
        raise ValueError("LH2 saturation table supports hydrogen releases only")
    return build_lh2_saturation_table_for_temperature_bounds(
        release.upstream_temperature.lower,
        release.upstream_temperature.upper,
        ambient_pressure_pa=ambient_pressure_pa,
        ambient_pressure_bounds_pa=ambient_pressure_bounds_pa,
        nodes=nodes,
        margin_k=margin_k,
    )


__all__ = [
    "FieldLH2FlashResult", "FieldLH2SourcePreparation", "FieldBlowdownVapourSchedule",
    "pressure_driven_lh2_mass_flow", "release_with_pressure_driven_lh2_mass_flow",
    "lh2_flash_source_from_release", "prepare_field_lh2_flash",
    "build_lh2_saturation_table_for_temperature_bounds",
    "build_lh2_saturation_table_for_release",
    "direct_vapour_schedule_from_cryogenic_blowdown",
]
