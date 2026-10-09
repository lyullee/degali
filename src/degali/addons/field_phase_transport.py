"""Conservative field transport handoff for phase-routed LH2 pool vapour.

The phase-routing model owns the rainout and dynamic-pool mass ledger.  This
module is the deliberately narrow bridge that can add only its timed,
conservative pool evaporation to the local semi-FV transport alongside the
separate direct-flash vapour branch.  It never converts a pool ledger into a
nozzle source or invents an atmospheric closure for unresolved in-flight
droplet evaporation.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
import math

from .field_contracts import BoundedValue
from .field_phase_routing import (
    FieldPhaseRoutingConfig,
    FieldPhaseRoutingEnvelope,
    FieldPhaseRoutingEnvelopeCase,
    FieldPhaseRoutingResult,
    FieldPhaseRoutingUncertainty,
    run_field_phase_routing_envelope,
)
from .field_pool_launch import (
    FieldPoolVapourSchedule,
    PoolVapourLaunchBoundary,
    distributed_source_from_pool_vapour_schedule,
    pool_vapour_schedule_from_phase_routing,
)
from .field_workflow import (
    FieldSemiFVRequest,
    FieldSemiFVScreeningResult,
    _obstacle_geometry_corner_cases,
    run_field_semi_fv_screening,
)


FIELD_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA = (
    "degali.field-phase-routing-transport-envelope.v1"
)


@dataclass(frozen=True)
class FieldPhaseRoutingTransportEnvelopeCase:
    """One phase/surface corner and its optional conservative field screen."""

    selection: tuple[tuple[str, float | str], ...]
    phase_routing: FieldPhaseRoutingResult
    pool_schedule: FieldPoolVapourSchedule | None
    screening: FieldSemiFVScreeningResult | None
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.selection, tuple):
            raise TypeError("phase-transport selection must be a tuple")
        seen: set[str] = set()
        for item in self.selection:
            if not isinstance(item, tuple) or len(item) != 2:
                raise TypeError("phase-transport selection entries must be key/value pairs")
            key, value = item
            if not isinstance(key, str) or not key.strip():
                raise ValueError("phase-transport selection keys must be non-empty strings")
            if key in seen:
                raise ValueError("phase-transport selection keys must be unique")
            seen.add(key)
            if isinstance(value, bool):
                raise TypeError("phase-transport selection values cannot be boolean")
            if isinstance(value, (int, float)):
                if not math.isfinite(float(value)):
                    raise ValueError("phase-transport selection numeric values must be finite")
            elif not isinstance(value, str) or not value.strip():
                raise TypeError(
                    "phase-transport selection values must be finite numbers or non-empty strings"
                )
        if not isinstance(self.phase_routing, FieldPhaseRoutingResult):
            raise TypeError("phase-transport phase_routing must be a FieldPhaseRoutingResult")
        if self.pool_schedule is not None and not isinstance(
            self.pool_schedule, FieldPoolVapourSchedule
        ):
            raise TypeError(
                "phase-transport pool_schedule must be a FieldPoolVapourSchedule or None"
            )
        if self.screening is not None and not isinstance(
            self.screening, FieldSemiFVScreeningResult
        ):
            raise TypeError(
                "phase-transport screening must be a FieldSemiFVScreeningResult or None"
            )
        if not isinstance(self.warnings, tuple) or any(
            not isinstance(item, str) or not item.strip() for item in self.warnings
        ):
            raise TypeError("phase-transport warnings must be non-empty strings")

    @property
    def completed(self) -> bool:
        return self.screening is not None and self.screening.completed


@dataclass(frozen=True)
class FieldPhaseRoutingTransportEnvelope:
    """Separate, unmerged field screens for every phase-routing corner.

    The object intentionally has no operational decision.  It propagates the
    liquid/pool surface boundary to transport, but unresolved launch physics
    and any later sensor/refinement envelope must remain visible to the caller.
    """

    request: FieldSemiFVRequest
    phase_envelope: FieldPhaseRoutingEnvelope
    pool_launch: PoolVapourLaunchBoundary
    pool_vertical_sigma_m: float
    cases: tuple[FieldPhaseRoutingTransportEnvelopeCase, ...]
    pool_vertical_sigma_uncertainty: BoundedValue | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.request, FieldSemiFVRequest):
            raise TypeError("phase-transport request must be a FieldSemiFVRequest")
        if not isinstance(self.phase_envelope, FieldPhaseRoutingEnvelope):
            raise TypeError("phase-transport phase_envelope must be a FieldPhaseRoutingEnvelope")
        if not isinstance(self.pool_launch, PoolVapourLaunchBoundary):
            raise TypeError("phase-transport pool_launch must be a PoolVapourLaunchBoundary")
        if not math.isfinite(float(self.pool_vertical_sigma_m)) or self.pool_vertical_sigma_m <= 0.0:
            raise ValueError("pool_vertical_sigma_m must be positive and finite")
        if not isinstance(self.cases, tuple) or not self.cases:
            raise ValueError("phase-transport envelope requires at least one case")
        if any(not isinstance(case, FieldPhaseRoutingTransportEnvelopeCase) for case in self.cases):
            raise TypeError(
                "phase-transport cases must contain only "
                "FieldPhaseRoutingTransportEnvelopeCase values"
            )
        sigma_values, include_sigma_selection = _pool_sigma_values(
            self.pool_vertical_sigma_m, self.pool_vertical_sigma_uncertainty,
        )
        obstacle_cases = _obstacle_geometry_corner_cases(
            self.request, max_cases=1_000_000,
        )
        expected = {
            phase_case.values
            + (("pool_vertical_sigma_m", sigma),) * int(include_sigma_selection)
            + obstacle_selection
            for phase_case in self.phase_envelope.cases
            for sigma in sigma_values
            for _obstacles, obstacle_selection in obstacle_cases
        }
        actual = {case.selection for case in self.cases}
        if len(self.cases) != len(expected):
            raise ValueError("phase-transport envelope must retain every phase/launch corner")
        if expected != actual:
            raise ValueError("phase-transport selections do not cover every phase/launch corner")

    @property
    def completed_case_count(self) -> int:
        return sum(case.completed for case in self.cases)

    @property
    def withheld_case_count(self) -> int:
        return sum(not case.completed for case in self.cases)


def _declared_pool_vapour_mass(result: FieldPhaseRoutingResult) -> float | None:
    """Read an available pool ledger mass only to distinguish zero from missing."""
    coupled = result.coupled_result
    phase = None if coupled is None else getattr(coupled, "phase_routing", None)
    coupling = None if phase is None else getattr(phase, "pool_coupling", None)
    value = None if coupling is None else getattr(coupling, "pool_vapour_mass_kg", None)
    if value is None:
        return None
    try:
        mass = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("phase-routing pool_vapour_mass_kg is not numeric") from error
    if not math.isfinite(mass) or mass < 0.0:
        raise ValueError("phase-routing pool_vapour_mass_kg is invalid")
    return mass


def _pool_sigma_values(
    nominal: float,
    uncertainty: BoundedValue | None,
) -> tuple[tuple[float, ...], bool]:
    """Validate and enumerate an optional bounded pool-source width."""
    if not math.isfinite(float(nominal)) or nominal <= 0.0:
        raise ValueError("pool_vertical_sigma_m must be positive and finite")
    if uncertainty is None:
        return (float(nominal),), False
    if not isinstance(uncertainty, BoundedValue):
        raise TypeError("pool_vertical_sigma_uncertainty must be BoundedValue or None")
    if (
        not isinstance(uncertainty.source, str)
        or not uncertainty.source.strip()
        or uncertainty.source.lower() == "unspecified"
    ):
        raise ValueError("pool_vertical_sigma_uncertainty requires an explicit source")
    if uncertainty.lower <= 0.0 or uncertainty.upper <= 0.0:
        raise ValueError("pool_vertical_sigma_uncertainty bounds must be positive")
    if not math.isclose(
        float(nominal), uncertainty.nominal, rel_tol=1.0e-12, abs_tol=1.0e-15,
    ):
        raise ValueError(
            "pool_vertical_sigma_uncertainty.nominal must match pool_vertical_sigma_m"
        )
    values = (uncertainty.nominal,) if uncertainty.is_exact else uncertainty.corners()
    return tuple(float(value) for value in values), True


def _phase_request(
    request: FieldSemiFVRequest,
    phase_result: FieldPhaseRoutingResult,
    pool_schedule: FieldPoolVapourSchedule | None,
    launch: PoolVapourLaunchBoundary,
    *,
    pool_vertical_sigma_m: float,
    post_release_duration_s: float,
) -> FieldSemiFVRequest:
    """Build direct-flash plus pool-vapour transport without schedule replacement.

    The phase-routing continuation is authoritative for this derived request:
    it converts the primary source to a finite zero-after-shutoff schedule and
    gives the distributed pool ledger the same post-release time window.
    """
    scenario = phase_result.scenario
    if pool_schedule is None:
        return replace(
            request, scenario=scenario,
            post_release_duration_s=post_release_duration_s,
        )
    pool_source = distributed_source_from_pool_vapour_schedule(
        pool_schedule, launch, vertical_sigma_m=pool_vertical_sigma_m,
    )
    return replace(
        request,
        scenario=scenario,
        post_release_duration_s=post_release_duration_s,
        distributed_vapour_sources=(*request.distributed_vapour_sources, pool_source),
    )


def _ambient_corner_cases(
    request: FieldSemiFVRequest,
    *,
    max_cases: int,
) -> tuple[tuple[dict[str, float], tuple[tuple[str, float], ...]], ...]:
    """Enumerate declared ambient boundary corners for phase-to-field transport.

    The phase ledger consumes temperature and pressure directly, while the
    local sensor conversion also consumes air density.  Keep all three as
    explicit deterministic selections so a nominal phase ledger cannot be
    mistaken for a resolved ambient boundary.
    """
    fields = request.ambient_uncertainty_fields()
    choices = {
        name: ((bound.nominal,) if bound.is_exact else bound.corners())
        for name, bound in fields.items()
    }
    names = tuple(choices)
    count = math.prod(len(values) for values in choices.values()) if choices else 1
    if count > max_cases:
        raise ValueError(
            f"{count} ambient boundary uncertainty corners exceed max_cases={max_cases}"
        )
    return tuple(
        (
            {name: float(value) for name, value in zip(names, selected)},
            tuple((name, float(value)) for name, value in zip(names, selected)),
        )
        for selected in (product(*(choices[name] for name in names)) if names else ((),))
    )


def run_field_phase_routing_transport_envelope(
    request: FieldSemiFVRequest,
    phase_config: FieldPhaseRoutingConfig,
    pool_launch: PoolVapourLaunchBoundary,
    *,
    pool_vertical_sigma_m: float,
    pool_vertical_sigma_uncertainty: BoundedValue | None = None,
    phase_uncertainty: FieldPhaseRoutingUncertainty | None = None,
    max_cases: int = 128,
    table_nodes: int = 161,
) -> FieldPhaseRoutingTransportEnvelope:
    """Propagate phase-routing source/weather/surface corners into semi-FV transport.

    A request cannot already own an atmospheric schedule or supplemental source:
    this adapter needs one unambiguous primary flash plus one pool-ledger branch.
    The local transport duration must cover release plus the declared pool
    continuation; otherwise a late pool source is withheld instead of clipped.
    """
    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(phase_config, FieldPhaseRoutingConfig):
        raise TypeError("phase_config must be a FieldPhaseRoutingConfig")
    if phase_uncertainty is not None and not isinstance(
        phase_uncertainty, FieldPhaseRoutingUncertainty
    ):
        raise TypeError("phase_uncertainty must be FieldPhaseRoutingUncertainty or None")
    if not isinstance(pool_launch, PoolVapourLaunchBoundary):
        raise TypeError("pool_launch must be a PoolVapourLaunchBoundary")
    sigma_values, include_sigma_selection = _pool_sigma_values(
        pool_vertical_sigma_m, pool_vertical_sigma_uncertainty,
    )
    if request.direct_vapour_schedule is not None:
        raise ValueError("phase-routing transport requires no pre-existing direct-vapour schedule")
    if request.distributed_vapour_sources:
        raise ValueError("phase-routing transport requires no pre-existing distributed vapour source")
    release_duration = request.scenario.source.duration_s
    if release_duration is None:
        raise ValueError("phase-routing transport requires a finite source duration")
    post_release_duration_s = (
        phase_config.post_release_duration_s
        if phase_uncertainty is None
        else phase_uncertainty.max_value(phase_config, "post_release_duration_s")
    )
    required_duration = release_duration + post_release_duration_s
    if request.transport.duration_s + 1.0e-12 < required_duration:
        raise ValueError(
            "transport.duration_s must cover source.duration_s plus phase-routing "
            "post_release_duration_s"
        )

    obstacle_cases = _obstacle_geometry_corner_cases(
        request, max_cases=max_cases,
    )
    ambient_cases = _ambient_corner_cases(request, max_cases=max_cases)
    phase_case_limit = max_cases // (
        len(ambient_cases) * len(sigma_values) * len(obstacle_cases)
    )
    if phase_case_limit < 1:
        raise ValueError(
            "ambient, pool-launch and obstacle uncertainty corners leave no room under max_cases"
        )
    phase_envelopes: list[
        tuple[dict[str, float], tuple[tuple[str, float], ...], FieldPhaseRoutingEnvelope]
    ] = []
    for ambient_values, ambient_selection in ambient_cases:
        phase_envelopes.append((
            ambient_values,
            ambient_selection,
            run_field_phase_routing_envelope(
                request.scenario,
                phase_config,
                ambient_temperature_k=ambient_values.get(
                    "ambient_temperature_k", request.ambient_temperature_k,
                ),
                ambient_pressure_pa=ambient_values.get(
                    "ambient_pressure_pa", request.ambient_pressure_pa,
                ),
                lh2_validation_available=request.lh2_validation_available,
                validation_evidence=request.validation_evidence,
                property_table=request.property_table,
                coordinate_reference=request.coordinate_reference,
                stability_alternatives=request.stability_alternatives,
                phase_uncertainty=phase_uncertainty,
                max_cases=phase_case_limit,
                table_nodes=table_nodes,
            ),
        ))
    phase_cases = tuple(
        FieldPhaseRoutingEnvelopeCase(
            phase_case.values + ambient_selection,
            phase_case.result,
        )
        for _ambient_values, ambient_selection, envelope in phase_envelopes
        for phase_case in envelope.cases
    )
    phase_warnings = tuple(
        warning
        for _values, _selection, envelope in phase_envelopes
        for warning in envelope.warnings
    )
    if any(selection for _values, selection, _envelope in phase_envelopes):
        phase_warnings += (
            "ambient temperature, pressure and air-density boundary corners were "
            "propagated into phase-to-field transport",
        )
    phase_envelope = FieldPhaseRoutingEnvelope(
        request.scenario,
        phase_config,
        phase_cases,
        property_table_used=all(
            envelope.property_table_used
            for _values, _selection, envelope in phase_envelopes
        ),
        warnings=phase_warnings,
        phase_uncertainty=phase_uncertainty,
    )
    total_cases = len(phase_envelope.cases) * len(sigma_values) * len(obstacle_cases)
    if total_cases > max_cases:
        raise ValueError(
            f"{total_cases} phase-routing transport corners exceed max_cases={max_cases}"
        )
    cases = []
    for phase_case in phase_envelope.cases:
        for sigma in sigma_values:
            for obstacle_geometry, obstacle_selection in obstacle_cases:
                phase_result = phase_case.result
                warnings = list(phase_result.applicability.warnings)
                selection = phase_case.values + (
                    (("pool_vertical_sigma_m", sigma),) if include_sigma_selection else ()
                ) + obstacle_selection
                selected_values = dict(phase_case.values)
                if not phase_result.completed:
                    cases.append(FieldPhaseRoutingTransportEnvelopeCase(
                        selection, phase_result, None, None,
                        tuple(warnings) + ("phase-routing ledger did not complete; field transport withheld",),
                    ))
                    continue
                try:
                    reported_pool_mass = _declared_pool_vapour_mass(phase_result)
                    if reported_pool_mass is not None and math.isclose(
                        reported_pool_mass, 0.0, abs_tol=1.0e-12,
                    ):
                        pool_schedule = None
                        warnings.append("phase-routing ledger reports zero pool vapour; only direct-flash transport is present")
                    else:
                        pool_schedule = pool_vapour_schedule_from_phase_routing(phase_result)
                    if (
                        pool_schedule is not None
                        and pool_schedule.schedule.duration_s > request.transport.duration_s + 1.0e-12
                    ):
                        raise ValueError("pool-vapour schedule exceeds declared transport duration")
                    field_request = _phase_request(
                        request, phase_result, pool_schedule, pool_launch,
                        pool_vertical_sigma_m=float(sigma),
                        post_release_duration_s=phase_result.config.post_release_duration_s,
                    )
                    field_request = replace(
                        field_request,
                        ambient_temperature_k=selected_values.get(
                            "ambient_temperature_k", request.ambient_temperature_k,
                        ),
                        ambient_pressure_pa=selected_values.get(
                            "ambient_pressure_pa", request.ambient_pressure_pa,
                        ),
                        ambient_air_density_kg_m3=selected_values.get(
                            "ambient_air_density_kg_m3", request.ambient_air_density_kg_m3,
                        ),
                        ambient_temperature_uncertainty_k=None,
                        ambient_pressure_uncertainty_pa=None,
                        ambient_air_density_uncertainty_kg_m3=None,
                        obstacle=(
                            obstacle_geometry[0]
                            if request.obstacle is not None else None
                        ),
                        obstacles=(
                            () if request.obstacle is not None else obstacle_geometry
                        ),
                        obstacle_geometry_uncertainty=(),
                    )
                    warnings.append(
                        "pool vapour is an explicit conservative internal scalar source; unresolved in-flight droplet vapour is not injected"
                    )
                    cases.append(FieldPhaseRoutingTransportEnvelopeCase(
                        selection,
                        phase_result,
                        pool_schedule,
                        run_field_semi_fv_screening(field_request),
                        tuple(warnings),
                    ))
                except (TypeError, ValueError, RuntimeError) as error:
                    cases.append(FieldPhaseRoutingTransportEnvelopeCase(
                        selection,
                        phase_result,
                        None,
                        None,
                        tuple(warnings) + (
                            "pool-vapour field transport withheld: " + str(error),
                        ),
                    ))
    return FieldPhaseRoutingTransportEnvelope(
        request, phase_envelope, pool_launch, float(pool_vertical_sigma_m), tuple(cases),
        pool_vertical_sigma_uncertainty=pool_vertical_sigma_uncertainty,
    )


def field_phase_routing_transport_envelope_report(
    result: FieldPhaseRoutingTransportEnvelope,
) -> dict[str, object]:
    """Serialise every phase-to-field handoff and each explicit withholding."""
    if not isinstance(result, FieldPhaseRoutingTransportEnvelope):
        raise TypeError("result must be a FieldPhaseRoutingTransportEnvelope")
    from .field_report import field_screening_report

    return {
        "schema": FIELD_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA,
        "pool_launch": {
            "closure_id": result.pool_launch.closure_id,
            "source_height_m": result.pool_launch.source_height_m,
            "evidence_id": result.pool_launch.evidence_id,
            "vertical_sigma_m": result.pool_vertical_sigma_m,
            "vertical_sigma_uncertainty": (
                None if result.pool_vertical_sigma_uncertainty is None
                else result.pool_vertical_sigma_uncertainty.as_dict()
            ),
        },
        "completed_case_count": result.completed_case_count,
        "withheld_case_count": result.withheld_case_count,
        "property_table_used": result.phase_envelope.property_table_used,
        "phase_envelope_warnings": list(result.phase_envelope.warnings),
        "phase_uncertainty": (
            None if result.phase_envelope.phase_uncertainty is None
            else result.phase_envelope.phase_uncertainty.as_record()
        ),
        "scope": (
            "direct-flash plus conservative dynamic-pool evaporation only; "
            "unresolved in-flight vapour, pool thermal/momentum launch and 3-D "
            "footprint are not supplied by this reduced-order transport handoff"
        ),
        "cases": [
            {
                "selection": dict(case.selection),
                "phase_routing": {
                    "status": case.phase_routing.applicability.status,
                    "reasons": list(case.phase_routing.applicability.reasons),
                    "warnings": list(case.phase_routing.applicability.warnings),
                    "heat_flux_w_m2": case.phase_routing.heat_flux_w_m2,
                },
                "pool_vapour": None if case.pool_schedule is None else {
                    "position_m": list(case.pool_schedule.position_m),
                    "duration_s": case.pool_schedule.schedule.duration_s,
                    "evaporated_mass_kg": case.pool_schedule.evaporated_mass_kg,
                    "maximum_wet_area_m2": case.pool_schedule.maximum_wet_area_m2,
                    "warnings": list(case.pool_schedule.warnings),
                },
                "field_screening": (
                    None if case.screening is None else field_screening_report(case.screening)
                ),
                "warnings": list(case.warnings),
            }
            for case in result.cases
        ],
    }


__all__ = [
    "FIELD_PHASE_ROUTING_TRANSPORT_ENVELOPE_SCHEMA",
    "FieldPhaseRoutingTransportEnvelopeCase", "FieldPhaseRoutingTransportEnvelope",
    "run_field_phase_routing_transport_envelope",
    "field_phase_routing_transport_envelope_report",
]
