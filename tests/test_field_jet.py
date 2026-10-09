import math
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from degali.addons.axisymmetric_jet import AxisymmetricJetSource
from degali.addons.field_contracts import (
    BoundedValue,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    SurfaceBoundary,
    WeatherState,
)
from degali.addons.field_jet import (
    FieldJetScalarHandoff,
    field_jet_scalar_handoff_from_near_field,
    jet_scalar_handoff_location,
    request_with_jet_scalar_handoff,
)
from degali.addons.field_lh2 import prepare_field_lh2_flash
from degali.addons.field_report import field_screening_report
from degali.addons.field_workflow import (
    FieldSemiFVRequest,
    run_field_semi_fv_envelope,
    run_field_semi_fv_screening,
)
from degali.addons.semi_fv_obstacle import SemiFVConfig
from degali.lh2 import LH2NearFieldResearchResult


def _near_field_result():
    rate = 0.2
    source = AxisymmetricJetSource(
        diameter=math.sqrt(4.0 * rate / math.pi), velocity=1.0, density=1.0,
        temperature=290.0, theta=0.0, x=2.0, y=0.5,
    )
    solution = SimpleNamespace(
        S=np.array((0.0, 1.0)),
        x=np.array((2.0, 3.0)),
        y=np.array((0.5, 0.7)),
        width=np.array((0.2, 0.4)),
        species_flux=np.array((rate, rate)),
    )
    return LH2NearFieldResearchResult(
        source=source,
        solution=solution,
        model=SimpleNamespace(spreading_ratio=1.16),
        maximum_boundary_residual=0.0,
        maximum_species_drift=0.0,
        maximum_energy_drift=0.0,
        warnings=["near-field test qualification"],
    )


def _scenario():
    return FieldScenario(
        source=ReleaseSource(
            fluid="lh2", location_m=(0.0, 0.0, 0.5), direction_m=(1.0, 0.0, 0.0),
            upstream_pressure=BoundedValue(0.4e6),
            upstream_temperature=BoundedValue(26.084),
            mass_flow_kg_s=BoundedValue(0.265),
            opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8),
            discharge_coefficient=BoundedValue(0.8),
            liquid_fraction=BoundedValue(0.922),
            flash_model="homogeneous_equilibrium", duration_s=1.0,
        ),
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        sensor=SensorModel((1.6, 0.0, 0.5)),
        temporal_mode="transient",
    )


def _transport():
    return SemiFVConfig(
        length_m=4.0, height_m=2.0, nx=20, nz=10,
        duration_s=1.0, time_step_s=0.005, source_sigma_m=0.2,
    )


def test_near_field_jet_station_becomes_conservative_vertical_scalar_boundary():
    handoff = field_jet_scalar_handoff_from_near_field(
        _near_field_result(),
        streamline_distance_m=0.5,
        source_id="jet-case-01",
        evidence_id="nearfield-audit-01",
    )

    assert handoff.downwind_offset_m == pytest.approx(0.5)
    assert handoff.elevation_offset_m == pytest.approx(0.1)
    assert handoff.vertical_sigma_m == pytest.approx(1.16 * 0.3 / math.sqrt(2.0))
    assert handoff.hydrogen_mass_flow_kg_s == pytest.approx(0.2)
    assert handoff.species_flux_relative_residual == pytest.approx(0.0)
    assert any("lateral dilution" in warning for warning in handoff.warnings)


def test_jet_handoff_is_relocated_with_wind_and_rejects_crosswind_jet():
    handoff = FieldJetScalarHandoff(
        source_id="jet-case-01", evidence_id="nearfield-audit-01",
        streamline_distance_m=0.5, downwind_offset_m=0.5, elevation_offset_m=0.1,
        vertical_sigma_m=0.2, hydrogen_mass_flow_kg_s=0.2,
        species_flux_relative_residual=0.0, maximum_wind_misalignment_deg=10.0,
    )
    source = ReleaseSource(location_m=(1.0, 2.0, 0.5), direction_m=(1.0, 0.0, 0.0))
    aligned = WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0))

    assert jet_scalar_handoff_location(source, aligned, handoff) == pytest.approx((1.5, 2.0, 0.6))
    crosswind = WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(180.0))
    with pytest.raises(ValueError, match="differs from the local wind plane"):
        jet_scalar_handoff_location(source, crosswind, handoff)


def test_field_screen_uses_jet_boundary_only_when_fresh_flash_rate_closes():
    scenario = _scenario()
    preparation = prepare_field_lh2_flash(scenario)
    assert preparation.flash_result is not None
    handoff = FieldJetScalarHandoff(
        source_id="jet-case-01", evidence_id="nearfield-audit-01",
        streamline_distance_m=0.5, downwind_offset_m=0.5, elevation_offset_m=0.1,
        vertical_sigma_m=0.15,
        hydrogen_mass_flow_kg_s=preparation.flash_result.flash.vapour_mass_flow,
        species_flux_relative_residual=0.0,
    )
    request = request_with_jet_scalar_handoff(
        FieldSemiFVRequest(scenario, transport=_transport()), handoff,
    )
    result = run_field_semi_fv_screening(request)

    assert result.completed
    assert result.wind_frame is not None
    assert result.wind_frame.origin_x_m == pytest.approx(0.5)
    assert result.transport is not None
    assert result.transport.diagnostics.mass_injected_kg == pytest.approx(
        handoff.hydrogen_mass_flow_kg_s
    )
    report = field_screening_report(result)
    assert report["jet_scalar_handoff"]["resolved_location_m"] == pytest.approx((0.5, 0.0, 0.6))
    assert report["jet_scalar_handoff"]["vertical_sigma_m"] == pytest.approx(0.15)

    stale = FieldJetScalarHandoff(
        source_id="jet-case-01", evidence_id="nearfield-audit-01",
        streamline_distance_m=0.5, downwind_offset_m=0.5, elevation_offset_m=0.1,
        vertical_sigma_m=0.15,
        hydrogen_mass_flow_kg_s=0.1,
        species_flux_relative_residual=0.0,
    )
    blocked = run_field_semi_fv_screening(
        request_with_jet_scalar_handoff(
            FieldSemiFVRequest(scenario, transport=_transport()), stale,
        )
    )
    assert not blocked.completed
    assert any("does not match" in reason for reason in blocked.applicability.reasons)


def test_jet_direction_corners_are_rechecked_during_field_envelope():
    base = _scenario()
    scenario = replace(
        base,
        source=replace(
            base.source,
            direction_uncertainty_m=(
                BoundedValue(1.0, 1.0, 1.0, "1", "orientation-review-01"),
                BoundedValue(0.0, 0.0, 1.0, "1", "orientation-review-01"),
                BoundedValue(0.0, 0.0, 0.0, "1", "orientation-review-01"),
            ),
        ),
        surface=SurfaceBoundary(heat_transfer_w_m2_k=BoundedValue(10.0)),
    )
    preparation = prepare_field_lh2_flash(base)
    assert preparation.flash_result is not None
    handoff = FieldJetScalarHandoff(
        source_id="jet-case-01", evidence_id="nearfield-audit-01",
        streamline_distance_m=0.5, downwind_offset_m=0.5, elevation_offset_m=0.1,
        vertical_sigma_m=0.15,
        hydrogen_mass_flow_kg_s=preparation.flash_result.flash.vapour_mass_flow,
        species_flux_relative_residual=0.0,
    )
    request = request_with_jet_scalar_handoff(
        FieldSemiFVRequest(scenario, transport=_transport()), handoff,
    )

    envelope = run_field_semi_fv_envelope(request, max_cases=2)

    assert len(envelope.cases) == 2
    assert envelope.completed_case_count == 1
    blocked = next(case for case in envelope.cases if case.result.applicability.status == "blocked")
    assert blocked.corner["source_direction_y"] == pytest.approx(1.0)
    assert any("differs from the local wind plane" in reason for reason in blocked.result.applicability.reasons)
