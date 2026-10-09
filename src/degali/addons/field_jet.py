"""Conservative axisymmetric-jet handoff to the local field scalar screen.

The semi-FV operator has a fixed vertical wind plane; it cannot resolve the
three-dimensional jet.  This module makes the reduction explicit: only a
conservation-screened near-field station may provide the location, vertical
Gaussian scale and direct-hydrogen rate for that operator.  The result keeps
the original source contract intact and is re-positioned with every declared
wind-direction corner.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import TYPE_CHECKING

import numpy as np

from .field_contracts import FieldScenario, ReleaseSource, WeatherState

if TYPE_CHECKING:
    from ..lh2 import LH2NearFieldResearchResult
    from .field_workflow import FieldSemiFVRequest


@dataclass(frozen=True)
class FieldJetScalarHandoff:
    """Auditable source boundary derived from one axisymmetric jet station.

    ``downwind_offset_m`` and ``elevation_offset_m`` are relative to the
    original field release location.  The first is reprojected along the
    scenario's wind-*to* direction, rather than frozen at the nominal wind.
    This keeps deterministic weather-direction corners physically aligned.
    """

    source_id: str
    evidence_id: str
    streamline_distance_m: float
    downwind_offset_m: float
    elevation_offset_m: float
    vertical_sigma_m: float
    hydrogen_mass_flow_kg_s: float
    species_flux_relative_residual: float
    maximum_flash_rate_relative_residual: float = 1.0e-5
    maximum_wind_misalignment_deg: float = 10.0
    model_id: str = "axisymmetric_jet_vertical_marginal_scalar_handoff"
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name, value in {
            "source_id": self.source_id,
            "evidence_id": self.evidence_id,
            "model_id": self.model_id,
        }.items():
            if (
                not isinstance(value, str)
                or not value.strip()
                or value.strip().lower() == "unspecified"
            ):
                raise ValueError(f"jet scalar handoff requires a declared {name}")
        positive = {
            "streamline_distance_m": self.streamline_distance_m,
            "vertical_sigma_m": self.vertical_sigma_m,
            "hydrogen_mass_flow_kg_s": self.hydrogen_mass_flow_kg_s,
            "maximum_wind_misalignment_deg": self.maximum_wind_misalignment_deg,
            "maximum_flash_rate_relative_residual": self.maximum_flash_rate_relative_residual,
        }
        for name, value in positive.items():
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if not math.isfinite(float(self.downwind_offset_m)) or self.downwind_offset_m < 0.0:
            raise ValueError("downwind_offset_m must be finite and non-negative")
        if not math.isfinite(float(self.elevation_offset_m)):
            raise ValueError("elevation_offset_m must be finite")
        if (
            not math.isfinite(float(self.species_flux_relative_residual))
            or self.species_flux_relative_residual < 0.0
        ):
            raise ValueError("species_flux_relative_residual must be finite and non-negative")
        if any(not isinstance(item, str) or not item.strip() for item in self.warnings):
            raise ValueError("jet scalar handoff warnings must contain non-empty strings")


def _finite_increasing(values: object, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or len(array) < 2 or not np.all(np.isfinite(array)):
        raise ValueError(f"near-field {name} must be a finite one-dimensional station array")
    if np.any(np.diff(array) <= 0.0):
        raise ValueError("near-field streamline stations must be strictly increasing")
    return array


def field_jet_scalar_handoff_from_near_field(
    result: "LH2NearFieldResearchResult",
    *,
    streamline_distance_m: float,
    source_id: str,
    evidence_id: str,
    maximum_species_flux_relative_residual: float = 2.0e-4,
    maximum_wind_misalignment_deg: float = 10.0,
) -> FieldJetScalarHandoff:
    """Project one verified near-field station into a vertical scalar source.

    For the axisymmetric species Gaussian
    ``rho*Y ~ exp(-r**2/(lambda_Y*B)**2)``, the vertical marginal has
    ``sigma_z = lambda_Y*B/sqrt(2)``.  The total direct-H2 rate is retained;
    this projection does *not* represent lateral dilution, jet momentum or a
    three-dimensional obstacle interaction.
    """
    from ..lh2 import LH2NearFieldResearchResult

    if not isinstance(result, LH2NearFieldResearchResult):
        raise TypeError("result must be an LH2NearFieldResearchResult")
    if not result.conservative:
        raise ValueError("jet scalar handoff requires a conservation-screened near-field result")
    if not math.isfinite(float(streamline_distance_m)) or streamline_distance_m < 0.0:
        raise ValueError("streamline_distance_m must be finite and non-negative")
    if (
        not math.isfinite(float(maximum_species_flux_relative_residual))
        or maximum_species_flux_relative_residual <= 0.0
    ):
        raise ValueError("maximum_species_flux_relative_residual must be positive and finite")

    solution = result.solution
    station = _finite_increasing(solution.S, "streamline distance")
    if not station[0] <= streamline_distance_m <= station[-1]:
        raise ValueError(
            "jet handoff streamline distance lies outside the near-field solution"
        )
    x = np.asarray(solution.x, dtype=float)
    y = np.asarray(solution.y, dtype=float)
    width = np.asarray(solution.width, dtype=float)
    species_flux = np.asarray(solution.species_flux, dtype=float)
    if any(
        array.shape != station.shape or not np.all(np.isfinite(array))
        for array in (x, y, width, species_flux)
    ):
        raise ValueError("near-field handoff arrays must align with streamline stations")
    if np.any(width <= 0.0) or np.any(species_flux <= 0.0):
        raise ValueError("near-field handoff requires positive width and species flux")
    direct_hydrogen_rate = float(result.source.fuel_mass_flow)
    if not math.isfinite(direct_hydrogen_rate) or direct_hydrogen_rate <= 0.0:
        raise ValueError("near-field source must contain a positive hydrogen mass flow")
    station_flux = float(np.interp(streamline_distance_m, station, species_flux))
    residual = abs(station_flux - direct_hydrogen_rate) / direct_hydrogen_rate
    if residual > maximum_species_flux_relative_residual:
        raise ValueError(
            "near-field H2 species flux does not close at the requested jet handoff station"
        )
    spreading_ratio = float(result.model.spreading_ratio)
    if not math.isfinite(spreading_ratio) or spreading_ratio <= 0.0:
        raise ValueError("near-field jet spreading ratio must be positive and finite")
    handoff_width = float(np.interp(streamline_distance_m, station, width))
    warnings = (
        "axisymmetric jet species profile is projected to a vertical scalar marginal; "
        "the local semi-FV operator does not resolve lateral dilution or jet momentum",
        "jet handoff remains conditional on field-relevant near-field validation",
    ) + tuple(str(item) for item in result.warnings)
    return FieldJetScalarHandoff(
        source_id=source_id,
        evidence_id=evidence_id,
        streamline_distance_m=float(streamline_distance_m),
        downwind_offset_m=float(np.interp(streamline_distance_m, station, x) - result.source.x),
        elevation_offset_m=float(np.interp(streamline_distance_m, station, y) - result.source.y),
        vertical_sigma_m=spreading_ratio * handoff_width / math.sqrt(2.0),
        hydrogen_mass_flow_kg_s=direct_hydrogen_rate,
        species_flux_relative_residual=residual,
        maximum_wind_misalignment_deg=maximum_wind_misalignment_deg,
        warnings=warnings,
    )


def jet_scalar_handoff_location(
    source: ReleaseSource,
    weather: WeatherState,
    handoff: FieldJetScalarHandoff,
) -> tuple[float, float, float]:
    """Locate a jet handoff in the current wind plane or reject misalignment."""
    if not isinstance(handoff, FieldJetScalarHandoff):
        raise TypeError("handoff must be a FieldJetScalarHandoff")
    release_direction = source.direction_unit
    horizontal_norm = math.hypot(release_direction[0], release_direction[1])
    if horizontal_norm <= 1.0e-12:
        raise ValueError("jet scalar handoff requires a release direction with horizontal component")
    release_x = release_direction[0] / horizontal_norm
    release_y = release_direction[1] / horizontal_norm
    wind_angle = weather.wind_to_math_radians()
    wind_x, wind_y = math.cos(wind_angle), math.sin(wind_angle)
    dot = min(1.0, max(-1.0, release_x * wind_x + release_y * wind_y))
    misalignment = math.degrees(math.acos(dot))
    if misalignment > handoff.maximum_wind_misalignment_deg + 1.0e-12:
        raise ValueError(
            "jet direction differs from the local wind plane by "
            f"{misalignment:.6g} deg, above the declared "
            f"{handoff.maximum_wind_misalignment_deg:.6g} deg limit"
        )
    location = (
        float(source.location_m[0]) + handoff.downwind_offset_m * wind_x,
        float(source.location_m[1]) + handoff.downwind_offset_m * wind_y,
        float(source.location_m[2]) + handoff.elevation_offset_m,
    )
    if location[2] < 0.0:
        raise ValueError("jet handoff lies below ground level")
    return location


def request_with_jet_scalar_handoff(
    request: "FieldSemiFVRequest",
    handoff: FieldJetScalarHandoff,
) -> "FieldSemiFVRequest":
    """Attach an audited jet handoff without modifying the source contract."""
    from .field_workflow import FieldSemiFVRequest

    if not isinstance(request, FieldSemiFVRequest):
        raise TypeError("request must be a FieldSemiFVRequest")
    if not isinstance(handoff, FieldJetScalarHandoff):
        raise TypeError("handoff must be a FieldJetScalarHandoff")
    if request.jet_scalar_handoff is not None:
        raise ValueError("request already declares a jet scalar handoff")
    if request.direct_vapour_schedule is not None:
        raise ValueError(
            "a scalar jet handoff has no time-resolved jet boundary; "
            "do not combine it with a direct-vapour schedule"
        )
    # Validate nominal alignment now. Direction corners are re-evaluated by
    # the workflow so they can fail closed independently.
    jet_scalar_handoff_location(request.scenario.source, request.scenario.weather, handoff)
    return replace(request, jet_scalar_handoff=handoff)


__all__ = [
    "FieldJetScalarHandoff", "field_jet_scalar_handoff_from_near_field",
    "jet_scalar_handoff_location", "request_with_jet_scalar_handoff",
]
