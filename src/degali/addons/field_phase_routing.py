"""Field-contract bridge to the existing conservative LH2 rainout/pool ledger."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from itertools import product
import math
from types import MappingProxyType
from typing import TYPE_CHECKING, Mapping

from .dynamic_pool import ConstantHeatFluxSurface
from .field_contracts import (
    BoundedValue,
    FieldApplicability,
    FieldCoordinateReference,
    FieldScenario,
    FieldValidationEvidence,
    assess_field_applicability,
)
from .droplet_rainout import DropletClass
from .field_meteorology import FieldStabilityAlternatives
from .field_lh2 import (
    FieldLH2SourcePreparation,
    build_lh2_saturation_table_for_release,
    prepare_field_lh2_flash,
)

if TYPE_CHECKING:
    from .lh2_property_table import LH2SaturationTable


_STABILITY_MAP = {
    "very_unstable": "A",
    "unstable": "B",
    "neutral": "D",
    "stable": "E",
    "very_stable": "F",
}


@dataclass(frozen=True)
class FieldPhaseRoutingConfig:
    """Declared boundaries required for liquid rainout/pool routing.

    Pool footprint and the d-squared flight coefficient are not inferred from
    the release; they must be supplied as visible engineering assumptions.
    ``post_release_duration_s`` controls how long deposited liquid is retained
    in the existing dynamic-pool ledger after the release ends.
    """

    post_release_duration_s: float
    puff_duration_s: float
    pool_area_m2: float
    pool_time_step_s: float
    evaporation_coefficient_m2_s: float
    relative_humidity_pct: float = 0.0
    roughness_m: float = 0.001
    averaging_time_s: float = 60.0
    wind_reference_height_m: float = 10.0
    pool_model: str = "dynamic"
    gas_model_options: Mapping[str, object] = field(default_factory=dict)
    phase_model_options: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        positive = {
            "post_release_duration_s": self.post_release_duration_s,
            "puff_duration_s": self.puff_duration_s,
            "pool_area_m2": self.pool_area_m2,
            "pool_time_step_s": self.pool_time_step_s,
            "evaporation_coefficient_m2_s": self.evaporation_coefficient_m2_s,
            "roughness_m": self.roughness_m,
            "averaging_time_s": self.averaging_time_s,
            "wind_reference_height_m": self.wind_reference_height_m,
        }
        for name, value in positive.items():
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not math.isfinite(float(self.relative_humidity_pct)) or not 0.0 <= self.relative_humidity_pct <= 100.0:
            raise ValueError("relative_humidity_pct must lie in [0, 100]")
        if self.pool_model not in {"dynamic", "fixed"}:
            raise ValueError("pool_model must be dynamic or fixed")
        for name, options in (
            ("gas_model_options", self.gas_model_options),
            ("phase_model_options", self.phase_model_options),
        ):
            if not isinstance(options, Mapping):
                raise TypeError(f"{name} must be a mapping")
            object.__setattr__(self, name, MappingProxyType(dict(options)))


@dataclass(frozen=True)
class FieldPhaseRoutingUncertainty:
    """Explicit bounded phase/pool inputs propagated as deterministic corners.

    The bounds are optional because legacy callers may still supply a nominal
    configuration.  When present, each interval must be backed by its own
    ``BoundedValue.source``; this object never infers a range from a nominal
    value.  Droplet population alternatives are accepted only as complete,
    evidence-backed simplex corners; no mass-fraction renormalization or
    probability distribution is inferred. Categorical model choices remain
    outside this envelope.
    """

    post_release_duration_s: BoundedValue | None = None
    puff_duration_s: BoundedValue | None = None
    pool_area_m2: BoundedValue | None = None
    pool_time_step_s: BoundedValue | None = None
    evaporation_coefficient_m2_s: BoundedValue | None = None
    droplet_population_corners: tuple[tuple[DropletClass, ...], ...] | None = None
    droplet_population_uncertainty_evidence_id: str | None = None

    def __post_init__(self) -> None:
        names = (
            "post_release_duration_s", "puff_duration_s", "pool_area_m2",
            "pool_time_step_s", "evaporation_coefficient_m2_s",
        )
        for name in names:
            value = getattr(self, name)
            if value is not None and not isinstance(value, BoundedValue):
                raise TypeError(f"{name} must be BoundedValue or None")
            if value is not None and (
                not isinstance(value.source, str)
                or not value.source.strip()
                or value.source.lower() == "unspecified"
            ):
                raise ValueError(f"{name} requires an explicit bounded-value source")
            if value is not None and (value.lower <= 0.0 or value.upper <= 0.0):
                raise ValueError(f"{name} bounds must be positive")
        corners = self.droplet_population_corners
        if corners is not None:
            if not corners:
                raise ValueError("droplet_population_corners must not be empty")
            if len({tuple(population) for population in corners}) != len(corners):
                raise ValueError("droplet_population_corners must be unique")
            for index, population in enumerate(corners):
                if not population or not all(
                    isinstance(item, DropletClass) for item in population
                ):
                    raise TypeError(
                        "droplet_population_corners must contain non-empty "
                        "DropletClass tuples"
                    )
                if any(
                    not math.isfinite(float(item.diameter_m))
                    or item.diameter_m <= 0.0
                    or not math.isfinite(float(item.mass_fraction))
                    or item.mass_fraction <= 0.0
                    for item in population
                ):
                    raise ValueError(
                        f"droplet population corner {index} must contain positive finite values"
                    )
                if not math.isclose(
                    sum(item.mass_fraction for item in population),
                    1.0,
                    rel_tol=1.0e-10,
                    abs_tol=1.0e-12,
                ):
                    raise ValueError(
                        f"droplet population corner {index} mass fractions must sum to one"
                    )
            evidence = self.droplet_population_uncertainty_evidence_id
            if (
                not isinstance(evidence, str)
                or not evidence.strip()
                or evidence.lower() == "unspecified"
            ):
                raise ValueError(
                    "droplet population uncertainty requires an explicit evidence_id"
                )
        elif self.droplet_population_uncertainty_evidence_id is not None:
            raise ValueError(
                "droplet_population_uncertainty_evidence_id requires population corners"
            )

    def items(self) -> tuple[tuple[str, BoundedValue], ...]:
        """Return only the declared bounded fields in stable schema order."""
        names = (
            "post_release_duration_s", "puff_duration_s", "pool_area_m2",
            "pool_time_step_s", "evaporation_coefficient_m2_s",
        )
        return tuple(
            (name, value)
            for name in names
            if (value := getattr(self, name)) is not None
        )

    def as_record(self) -> dict[str, object]:
        record = {name: value.as_dict() for name, value in self.items()}
        if self.droplet_population_corners is not None:
            record["droplet_population"] = {
                "evidence_id": self.droplet_population_uncertainty_evidence_id,
                "corners": [
                    {
                        "classes": [
                            {
                                "diameter_m": item.diameter_m,
                                "mass_fraction": item.mass_fraction,
                            }
                            for item in population
                        ]
                    }
                    for population in self.droplet_population_corners
                ],
            }
        return record

    def max_value(self, config: "FieldPhaseRoutingConfig", name: str) -> float:
        """Return the largest declared value for a config field."""
        declared = getattr(self, name)
        return float(getattr(config, name) if declared is None else declared.upper)

    def config_corners(
        self, config: "FieldPhaseRoutingConfig",
    ) -> tuple[tuple[FieldPhaseRoutingConfig, tuple[tuple[str, float], ...]], ...]:
        """Build phase-config corners and verify each bound's nominal anchor."""
        choices: list[tuple[str, tuple[float, ...]]] = []
        for name, bounded in self.items():
            nominal = float(getattr(config, name))
            if not math.isclose(nominal, bounded.nominal, rel_tol=1.0e-12, abs_tol=1.0e-15):
                raise ValueError(
                    f"phase uncertainty {name}.nominal must match phase config nominal"
                )
            if bounded.lower <= 0.0 or bounded.upper <= 0.0:
                raise ValueError(f"phase uncertainty {name} bounds must be positive")
            choices.append((name, (bounded.nominal,) if bounded.is_exact else bounded.corners()))
        configured_population = tuple(config.phase_model_options.get("droplet_classes", ()))
        if self.droplet_population_corners is None:
            population_choices = ((configured_population, ()),)
        else:
            if not configured_population:
                raise ValueError(
                    "droplet population uncertainty requires configured droplet_classes"
                )
            if not any(
                tuple(population) == configured_population
                for population in self.droplet_population_corners
            ):
                raise ValueError(
                    "droplet population uncertainty corners must include the nominal population"
                )
            population_choices = tuple(
                (
                    tuple(population),
                    tuple(
                        item
                        for class_index, droplet in enumerate(population)
                        for item in (
                            (f"phase_droplet_class_{class_index}_diameter_m", float(droplet.diameter_m)),
                            (f"phase_droplet_class_{class_index}_mass_fraction", float(droplet.mass_fraction)),
                        )
                    ),
                )
                for population in self.droplet_population_corners
            )
        if not choices:
            choices = []
        names = tuple(name for name, _ in choices)
        return tuple(
            (
                replace(
                    config,
                    **dict(zip(names, values)),
                ) if population == configured_population else replace(
                    config,
                    **{
                        **dict(zip(names, values)),
                        "phase_model_options": {
                            **dict(config.phase_model_options),
                            "droplet_classes": population,
                        },
                    },
                ),
                tuple((f"phase_{name}", float(value)) for name, value in zip(names, values))
                + population_selection,
            )
            for values in product(*(values for _, values in choices))
            for population, population_selection in population_choices
        )


@dataclass(frozen=True)
class FieldPhaseRoutingResult:
    """Field source preparation plus the existing multiphase routing result."""

    scenario: FieldScenario
    config: FieldPhaseRoutingConfig
    applicability: FieldApplicability
    source_preparation: FieldLH2SourcePreparation | None
    heat_flux_w_m2: float | None
    coupled_result: object | None

    def __post_init__(self) -> None:
        if not isinstance(self.scenario, FieldScenario):
            raise TypeError("phase-routing scenario must be a FieldScenario")
        if not isinstance(self.config, FieldPhaseRoutingConfig):
            raise TypeError("phase-routing config must be a FieldPhaseRoutingConfig")
        if not isinstance(self.applicability, FieldApplicability):
            raise TypeError("phase-routing applicability must be a FieldApplicability")
        if self.source_preparation is not None and not isinstance(
            self.source_preparation, FieldLH2SourcePreparation
        ):
            raise TypeError(
                "phase-routing source_preparation must be a "
                "FieldLH2SourcePreparation or None"
            )
        if self.heat_flux_w_m2 is not None:
            if isinstance(self.heat_flux_w_m2, bool):
                raise TypeError("phase-routing heat_flux_w_m2 cannot be boolean")
            if not math.isfinite(float(self.heat_flux_w_m2)) or self.heat_flux_w_m2 < 0.0:
                raise ValueError(
                    "phase-routing heat_flux_w_m2 must be finite and non-negative"
                )

    @property
    def completed(self) -> bool:
        return self.coupled_result is not None and self.applicability.status != "blocked"


@dataclass(frozen=True)
class FieldPhaseRoutingEnvelopeCase:
    """One deterministic phase-routing input corner and its separate ledger."""

    values: tuple[tuple[str, float | str], ...]
    result: FieldPhaseRoutingResult

    def __post_init__(self) -> None:
        if not isinstance(self.values, tuple):
            raise TypeError("phase-routing corner values must be a tuple")
        seen: set[str] = set()
        for item in self.values:
            if not isinstance(item, tuple) or len(item) != 2:
                raise TypeError("phase-routing corner values must be key/value pairs")
            key, value = item
            if not isinstance(key, str) or not key.strip():
                raise ValueError("phase-routing corner keys must be non-empty strings")
            if key in seen:
                raise ValueError("phase-routing corner keys must be unique")
            seen.add(key)
            if isinstance(value, bool):
                raise TypeError("phase-routing corner values cannot be boolean")
            if isinstance(value, (int, float)):
                if not math.isfinite(float(value)):
                    raise ValueError("phase-routing corner numeric values must be finite")
            elif not isinstance(value, str) or not value.strip():
                raise TypeError(
                    "phase-routing corner values must be finite numbers or non-empty strings"
                )
        if not isinstance(self.result, FieldPhaseRoutingResult):
            raise TypeError("phase-routing envelope result must be a FieldPhaseRoutingResult")

    @property
    def corner(self) -> dict[str, float | str]:
        return dict(self.values)


@dataclass(frozen=True)
class FieldPhaseRoutingEnvelope:
    """Untruncated phase-routing sensitivity envelope over physical inputs."""

    scenario: FieldScenario
    config: FieldPhaseRoutingConfig
    cases: tuple[FieldPhaseRoutingEnvelopeCase, ...]
    property_table_used: bool
    warnings: tuple[str, ...] = ()
    phase_uncertainty: FieldPhaseRoutingUncertainty | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.scenario, FieldScenario):
            raise TypeError("phase-routing envelope scenario must be a FieldScenario")
        if not isinstance(self.config, FieldPhaseRoutingConfig):
            raise TypeError("phase-routing envelope config must be a FieldPhaseRoutingConfig")
        if not isinstance(self.cases, tuple) or not self.cases:
            raise ValueError("phase-routing envelope requires at least one case")
        if any(not isinstance(case, FieldPhaseRoutingEnvelopeCase) for case in self.cases):
            raise TypeError(
                "phase-routing envelope cases must contain only "
                "FieldPhaseRoutingEnvelopeCase values"
            )
        if len({case.values for case in self.cases}) != len(self.cases):
            raise ValueError("phase-routing envelope corner selections must be unique")
        if not isinstance(self.property_table_used, bool):
            raise TypeError("phase-routing envelope property_table_used must be boolean")
        if not isinstance(self.warnings, tuple) or any(
            not isinstance(item, str) or not item.strip() for item in self.warnings
        ):
            raise TypeError("phase-routing envelope warnings must be non-empty strings")
        if self.phase_uncertainty is not None and not isinstance(
            self.phase_uncertainty, FieldPhaseRoutingUncertainty
        ):
            raise TypeError("phase_uncertainty must be FieldPhaseRoutingUncertainty or None")

    @property
    def completed_case_count(self) -> int:
        return sum(case.result.completed for case in self.cases)

    @property
    def blocked_case_count(self) -> int:
        return sum(case.result.applicability.status == "blocked" for case in self.cases)


def _merge(*items: FieldApplicability, warnings: tuple[str, ...] = ()) -> FieldApplicability:
    reasons = tuple(reason for item in items for reason in item.reasons)
    merged_warnings = tuple(warning for item in items for warning in item.warnings) + warnings
    if reasons:
        return FieldApplicability("blocked", reasons, merged_warnings, False)
    return FieldApplicability(
        "conditional" if any(item.status == "conditional" for item in items) or merged_warnings else "accepted",
        warnings=merged_warnings,
        uncertainty_complete=all(item.uncertainty_complete for item in items),
    )


def run_field_phase_routing(
    scenario: FieldScenario,
    config: FieldPhaseRoutingConfig,
    *,
    ambient_temperature_k: float = 295.0,
    ambient_pressure_pa: float = 101325.0,
    lh2_validation_available: bool = False,
    validation_evidence: FieldValidationEvidence | None = None,
    property_table: "LH2SaturationTable | None" = None,
    coordinate_reference: FieldCoordinateReference | None = None,
) -> FieldPhaseRoutingResult:
    """Route an explicit field LH2 flash through droplet/rainout/pool ledgers.

    A declared convective surface boundary is converted to a constant heat
    flux only for the existing pool model: ``q'' = h(T_surface - T_sat,H2)``.
    It is not promoted to a universal soil, concrete, water or ice model.
    Atmospheric source terms from in-flight and pool evaporation remain
    explicitly unresolved unless a separate pool-vapour launch model is given.
    """
    if coordinate_reference is not None and not isinstance(
        coordinate_reference, FieldCoordinateReference
    ):
        raise TypeError("coordinate_reference must be a FieldCoordinateReference or None")
    gate = assess_field_applicability(
        scenario,
        lh2_validation_available=lh2_validation_available,
        validation_evidence=validation_evidence,
    )
    if scenario.model_family != "degali":
        gate = _merge(gate, FieldApplicability(
            "blocked", reasons=("field phase routing is a DEGALI path, not a SLABx solver",)
        ))
    source = scenario.source
    if scenario.temporal_mode != "transient" or source.duration_s is None:
        gate = _merge(gate, FieldApplicability(
            "blocked", reasons=("rainout/pool routing requires a finite release duration",)
        ))
    if source.location_m[2] <= 0.0:
        gate = _merge(gate, FieldApplicability(
            "blocked", reasons=("coupled phase routing requires a positive release height",)
        ))
    surface = scenario.surface
    if surface.substrate == "unspecified" or surface.evidence_id == "unspecified":
        gate = _merge(gate, FieldApplicability(
            "blocked", reasons=(
                "pool heat boundary requires declared substrate and evidence_id",
            )
        ))
    if gate.status == "blocked":
        return FieldPhaseRoutingResult(scenario, config, gate, None, None, None)

    preparation = prepare_field_lh2_flash(
        scenario,
        ambient_temperature_k=ambient_temperature_k,
        ambient_pressure_pa=ambient_pressure_pa,
        lh2_validation_available=lh2_validation_available,
        validation_evidence=validation_evidence,
        property_table=property_table,
    )
    applicability = _merge(gate, preparation.applicability)
    if applicability.status == "blocked" or preparation.flash_result is None:
        return FieldPhaseRoutingResult(
            scenario, config, applicability, preparation, None, None
        )

    from CoolProp.CoolProp import PropsSI

    saturation_temperature = float(PropsSI(
        "T", "P", ambient_pressure_pa, "Q", 0, "Hydrogen"
    ))
    heat_flux = surface.heat_transfer_w_m2_k.nominal * max(
        0.0, surface.surface_temperature_k.nominal - saturation_temperature
    )
    if heat_flux <= 0.0:
        blocked = _merge(applicability, FieldApplicability(
            "blocked", reasons=(
                "declared surface temperature does not provide positive heat flux to an LH2 pool",
            )
        ))
        return FieldPhaseRoutingResult(
            scenario, config, blocked, preparation, heat_flux, None
        )
    substrate = ConstantHeatFluxSurface(
        heat_flux_w_m2=heat_flux,
        calibration_label=surface.evidence_id,
    )
    direction = source.direction_unit
    horizontal = math.hypot(direction[0], direction[1])
    azimuth = math.atan2(direction[1], direction[0])
    elevation = math.atan2(direction[2], horizontal)
    if preparation.flash_result.flash.vapour_mass_flow > 0.0 and not math.isclose(
        elevation, 0.0, abs_tol=1.0e-12
    ):
        blocked = _merge(applicability, FieldApplicability(
            "blocked", reasons=(
                "direct-vapour coupled branch requires a horizontal release direction",
            )
        ))
        return FieldPhaseRoutingResult(
            scenario, config, blocked, preparation, heat_flux, None
        )
    from ..lh2 import run_lh2_coupled_transient_research

    wind_to_earth_rad = scenario.weather.wind_to_math_radians()
    wind_to_site_rad = (
        wind_to_earth_rad if coordinate_reference is None
        else coordinate_reference.earth_to_local_math_radians(wind_to_earth_rad)
    )
    try:
        coupled = run_lh2_coupled_transient_research(
            preparation.flash_result.flash,
            release_duration_s=source.duration_s,
            post_release_duration_s=config.post_release_duration_s,
            puff_duration_s=config.puff_duration_s,
            release_position_m=source.location_m,
            release_azimuth_rad=azimuth,
            release_elevation_rad=elevation,
            wind_speed_m_s=scenario.weather.speed_m_s.nominal,
            wind_to_angle_rad=wind_to_site_rad,
            evaporation_coefficient_m2_s=config.evaporation_coefficient_m2_s,
            pool_area_m2=config.pool_area_m2,
            pool_time_step_s=config.pool_time_step_s,
            substrate=substrate,
            ambient_temperature_k=ambient_temperature_k,
            ambient_pressure_pa=ambient_pressure_pa,
            relative_humidity=config.relative_humidity_pct,
            roughness=config.roughness_m,
            stability=_STABILITY_MAP[scenario.weather.stability],
            averaging=config.averaging_time_s,
            wind_reference_height=config.wind_reference_height_m,
            pool_model=config.pool_model,
            gas_model_options=dict(config.gas_model_options),
            phase_model_options=dict(config.phase_model_options),
        )
    except (ValueError, RuntimeError) as error:
        blocked = _merge(applicability, FieldApplicability(
            "blocked", reasons=(str(error),)
        ))
        return FieldPhaseRoutingResult(
            scenario, config, blocked, preparation, heat_flux, None
        )

    warnings = [
        "surface boundary is a declared constant heat flux derived from h(Tsurface-Tsat), not a resolved substrate-conduction model",
    ]
    if not coupled.atmospherically_complete:
        warnings.append(
            "pool and in-flight vapour masses are conserved but retain unresolved atmospheric launch boundaries",
        )
    if not coupled.conservative:
        blocked = _merge(applicability, FieldApplicability(
            "blocked", reasons=("coupled rainout/pool mass conservation screen failed",)
        ), warnings=tuple(warnings))
        return FieldPhaseRoutingResult(
            scenario, config, blocked, preparation, heat_flux, coupled)
    final = _merge(applicability, FieldApplicability(
        "conditional" if not coupled.accepted else "accepted",
        warnings=tuple(warnings),
        uncertainty_complete=True,
    ))
    return FieldPhaseRoutingResult(
        scenario, config, final, preparation, heat_flux, coupled
    )


def _phase_routing_corner_cases(
    scenario: FieldScenario,
    *,
    max_cases: int,
) -> tuple[dict[str, float], ...]:
    """Vary source/weather/surface bounds, not irrelevant sensor settings."""
    if max_cases < 1:
        raise ValueError("max_cases must be positive")
    all_fields = scenario.uncertainty_fields()
    phase_names = tuple(
        list(scenario.source.uncertainty_fields())
        + list(scenario.weather.uncertainty_fields())
        + list(scenario.surface.uncertainty_fields())
    )
    choices = {
        name: ((all_fields[name].nominal,) if all_fields[name].is_exact else all_fields[name].corners())
        for name in phase_names
    }
    count = math.prod(len(values) for values in choices.values())
    if count > max_cases:
        raise ValueError(f"{count} phase-routing uncertainty corners exceed max_cases={max_cases}")
    nominal = {name: field.nominal for name, field in all_fields.items()}
    return tuple(
        {**nominal, **dict(zip(phase_names, values))}
        for values in product(*(choices[name] for name in phase_names))
    )


def run_field_phase_routing_envelope(
    scenario: FieldScenario,
    config: FieldPhaseRoutingConfig,
    *,
    phase_uncertainty: FieldPhaseRoutingUncertainty | None = None,
    ambient_temperature_k: float = 295.0,
    ambient_pressure_pa: float = 101325.0,
    lh2_validation_available: bool = False,
    validation_evidence: FieldValidationEvidence | None = None,
    property_table: "LH2SaturationTable | None" = None,
    coordinate_reference: FieldCoordinateReference | None = None,
    stability_alternatives: FieldStabilityAlternatives | None = None,
    max_cases: int = 128,
    table_nodes: int = 161,
) -> FieldPhaseRoutingEnvelope:
    """Run all declared physical corners as separate ledgers.

    This includes the optional evidence-backed scalar phase/pool bounds in
    ``phase_uncertainty``.  Each selected phase configuration receives its own
    ledger; no corner histories are merged.

    Sensor gain, bias and response do not affect rainout or pool physics, so
    they are fixed at nominal values here rather than multiplying identical
    mass-ledger cases. Each output remains a separate deterministic scenario,
    not a probabilistic interval or a merged pool history.
    """
    if stability_alternatives is not None and not isinstance(
        stability_alternatives, FieldStabilityAlternatives
    ):
        raise TypeError("stability_alternatives must be a FieldStabilityAlternatives")
    if coordinate_reference is not None and not isinstance(
        coordinate_reference, FieldCoordinateReference
    ):
        raise TypeError("coordinate_reference must be a FieldCoordinateReference or None")
    if phase_uncertainty is not None and not isinstance(
        phase_uncertainty, FieldPhaseRoutingUncertainty
    ):
        raise TypeError("phase_uncertainty must be FieldPhaseRoutingUncertainty or None")
    corners = _phase_routing_corner_cases(scenario, max_cases=max_cases)
    config_corners = (
        ((config, ()),) if phase_uncertainty is None
        else phase_uncertainty.config_corners(config)
    )
    stability_classes = (
        (scenario.weather.stability,)
        if stability_alternatives is None
        else stability_alternatives.classes_for(scenario.weather.stability)
    )
    total_cases = len(corners) * len(stability_classes) * len(config_corners)
    if total_cases > max_cases:
        raise ValueError(
            f"{total_cases} phase-routing uncertainty corners including stability alternatives "
            f"exceed max_cases={max_cases}"
        )
    table = property_table
    warnings: list[str] = []
    if table is None:
        try:
            table = build_lh2_saturation_table_for_release(
                scenario.source,
                ambient_pressure_pa=ambient_pressure_pa,
                nodes=table_nodes,
            )
        except (ValueError, RuntimeError) as error:
            warnings.append(
                "LH2 saturation table was not constructed; direct property path used: "
                f"{error}"
            )
    cases = []
    for values in corners:
        for stability in stability_classes:
            corner_scenario = scenario.at_corner(values)
            if stability_alternatives is not None:
                corner_scenario = replace(
                    corner_scenario,
                    weather=replace(corner_scenario.weather, stability=stability),
                )
            base_selection: tuple[tuple[str, float | str], ...] = tuple(
                (name, float(value)) for name, value in values.items()
            )
            if stability_alternatives is not None:
                base_selection += (("weather_stability", stability),)
            for corner_config, phase_selection in config_corners:
                cases.append(FieldPhaseRoutingEnvelopeCase(
                    base_selection + phase_selection,
                    run_field_phase_routing(
                        corner_scenario,
                        corner_config,
                        ambient_temperature_k=ambient_temperature_k,
                        ambient_pressure_pa=ambient_pressure_pa,
                        lh2_validation_available=lh2_validation_available,
                        validation_evidence=validation_evidence,
                        property_table=table,
                        coordinate_reference=coordinate_reference,
                    ),
                ))
    return FieldPhaseRoutingEnvelope(
        scenario, config, tuple(cases), property_table_used=table is not None,
        warnings=tuple(warnings), phase_uncertainty=phase_uncertainty,
    )


__all__ = [
    "FieldPhaseRoutingConfig", "FieldPhaseRoutingResult",
    "FieldPhaseRoutingUncertainty",
    "FieldPhaseRoutingEnvelopeCase", "FieldPhaseRoutingEnvelope",
    "run_field_phase_routing", "run_field_phase_routing_envelope",
]
