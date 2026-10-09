import math
from types import SimpleNamespace

import degali.lh2 as lh2
import degali.addons.field_phase_routing as field_phase_routing
import pytest
from dataclasses import replace

from degali.addons.dynamic_pool import ConstantHeatFluxSurface
from degali.addons.droplet_rainout import DropletClass
from degali.addons.field_contracts import (
    BoundedValue,
    FieldCoordinateReference,
    FieldScenario,
    ReleaseSource,
    SensorModel,
    SurfaceBoundary,
    WeatherState,
)
from degali.addons.field_phase_routing import (
    FieldPhaseRoutingConfig,
    FieldPhaseRoutingResult,
    FieldPhaseRoutingUncertainty,
    run_field_phase_routing,
    run_field_phase_routing_envelope,
)
from degali.addons.field_contracts import FieldApplicability
from degali.addons.field_meteorology import FieldStabilityAlternatives
from degali.addons.field_pool_launch import pool_vapour_schedule_from_phase_routing


def _scenario(*, surface=None):
    return FieldScenario(
        source=ReleaseSource(
            fluid="lh2",
            location_m=(0.0, 0.0, 0.5),
            upstream_pressure=BoundedValue(0.4e6),
            upstream_temperature=BoundedValue(26.084),
            mass_flow_kg_s=BoundedValue(0.265),
            opening_area_m2=BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8),
            discharge_coefficient=BoundedValue(0.8),
            liquid_fraction=BoundedValue(0.922),
            flash_model="homogeneous_equilibrium",
            duration_s=1.0,
        ),
        weather=WeatherState(speed_m_s=BoundedValue(2.0), direction_deg=BoundedValue(270.0)),
        surface=surface or SurfaceBoundary(
            heat_transfer_w_m2_k=BoundedValue(10.0),
            surface_temperature_k=BoundedValue(293.15),
            substrate="concrete",
            evidence_id="facility-surface-log-01",
        ),
        temporal_mode="transient",
    )


def _config():
    return FieldPhaseRoutingConfig(
        post_release_duration_s=1.0,
        puff_duration_s=1.0,
        pool_area_m2=0.5,
        pool_time_step_s=0.1,
        evaporation_coefficient_m2_s=1.0e-7,
    )


def test_phase_routing_model_options_are_immutable_snapshots():
    gas_options = {"radial_points": 41}
    phase_options = {"maximum_droplet_time_s": 0.1}
    config = replace(
        _config(),
        gas_model_options=gas_options,
        phase_model_options=phase_options,
    )
    gas_options["radial_points"] = 5
    phase_options["maximum_droplet_time_s"] = 9.0

    assert config.gas_model_options["radial_points"] == 41
    assert config.phase_model_options["maximum_droplet_time_s"] == pytest.approx(0.1)
    with pytest.raises(TypeError):
        config.gas_model_options["new"] = 1


def test_field_phase_routing_uses_declared_surface_heat_boundary(monkeypatch):
    calls = {}

    def fake_coupled(source, **kwargs):
        calls["source"] = source
        calls.update(kwargs)
        return SimpleNamespace(
            atmospherically_complete=False,
            conservative=True,
            accepted=True,
        )

    monkeypatch.setattr(lh2, "run_lh2_coupled_transient_research", fake_coupled)
    result = run_field_phase_routing(
        _scenario(), _config(),
        coordinate_reference=FieldCoordinateReference(
            "plant-grid-rev-A", "transfer-skid-origin", 90.0, "site-grade", "layout-drawing-A",
        ),
    )
    assert result.completed
    assert result.heat_flux_w_m2 is not None and result.heat_flux_w_m2 > 0.0
    assert isinstance(calls["substrate"], ConstantHeatFluxSurface)
    assert calls["substrate"].calibration_label == "facility-surface-log-01"
    assert calls["release_duration_s"] == 1.0
    assert calls["wind_to_angle_rad"] == pytest.approx(1.5 * math.pi)
    assert any("unresolved atmospheric launch" in item for item in result.applicability.warnings)


def test_phase_routing_direction_corners_fail_closed_for_nonhorizontal_vapour(monkeypatch):
    def fake_coupled(source, **kwargs):
        return SimpleNamespace(
            atmospherically_complete=False,
            conservative=True,
            accepted=True,
        )

    monkeypatch.setattr(lh2, "run_lh2_coupled_transient_research", fake_coupled)
    base = _scenario()
    scenario = replace(
        base,
        source=replace(
            base.source,
            direction_uncertainty_m=(
                BoundedValue(1.0, 1.0, 1.0, "1", "orientation-review-01"),
                BoundedValue(0.0, 0.0, 0.0, "1", "orientation-review-01"),
                BoundedValue(0.0, 0.0, 0.2, "1", "orientation-review-01"),
            ),
        ),
    )

    envelope = run_field_phase_routing_envelope(scenario, _config(), max_cases=2)

    assert len(envelope.cases) == 2
    assert envelope.completed_case_count == 1
    blocked = next(case for case in envelope.cases if case.result.applicability.status == "blocked")
    assert blocked.corner["source_direction_z"] == pytest.approx(0.2)
    assert any("horizontal release direction" in reason for reason in blocked.result.applicability.reasons)


def test_field_phase_routing_blocks_undefined_surface_evidence():
    result = run_field_phase_routing(
        _scenario(surface=SurfaceBoundary()), _config()
    )
    assert not result.completed
    assert result.applicability.status == "blocked"
    assert "substrate and evidence_id" in result.applicability.reasons[0]


def test_phase_routing_envelope_varies_surface_ledger_inputs_but_not_sensor_calibration(monkeypatch):
    calls = []
    table = object()

    def fake_table(*args, **kwargs):
        return table

    def fake_routing(scenario, config, **kwargs):
        calls.append((scenario, kwargs["property_table"]))
        return FieldPhaseRoutingResult(
            scenario, config, FieldApplicability("accepted", uncertainty_complete=True),
            None, None, object(),
        )

    monkeypatch.setattr(field_phase_routing, "build_lh2_saturation_table_for_release", fake_table)
    monkeypatch.setattr(field_phase_routing, "run_field_phase_routing", fake_routing)
    base = _scenario(surface=SurfaceBoundary(
        heat_transfer_w_m2_k=BoundedValue(10.0, 5.0, 20.0),
        surface_temperature_k=BoundedValue(293.15),
        substrate="concrete", evidence_id="facility-surface-log-01",
    ))
    scenario = replace(
        base,
        sensor=SensorModel(
            (1.0, 0.0, 0.5), response_time_s=BoundedValue(1.0, 0.5, 2.0),
        ),
    )
    envelope = run_field_phase_routing_envelope(scenario, _config(), max_cases=2)

    assert len(envelope.cases) == 2
    assert envelope.property_table_used
    assert envelope.completed_case_count == 2
    assert {case.corner["surface_heat_transfer_w_m2_k"] for case in envelope.cases} == {5.0, 20.0}
    assert {case.corner["sensor_response_time_s"] for case in envelope.cases} == {1.0}
    assert all(received_table is table for _scenario, received_table in calls)
    with pytest.raises(ValueError, match="phase-routing uncertainty corners"):
        run_field_phase_routing_envelope(scenario, _config(), max_cases=1)


def test_phase_routing_envelope_includes_declared_stability_alternatives(monkeypatch):
    calls = []

    def fake_routing(scenario, config, **_kwargs):
        calls.append(scenario.weather.stability)
        return FieldPhaseRoutingResult(
            scenario, config, FieldApplicability("accepted", uncertainty_complete=True),
            None, None, object(),
        )

    monkeypatch.setattr(field_phase_routing, "run_field_phase_routing", fake_routing)
    envelope = run_field_phase_routing_envelope(
        _scenario(), _config(),
        stability_alternatives=FieldStabilityAlternatives(
            ("stable",), "met-stability-review-A",
        ),
        max_cases=2,
    )

    assert len(envelope.cases) == 2
    assert calls == ["neutral", "stable"]
    assert {case.corner["weather_stability"] for case in envelope.cases} == {
        "neutral", "stable",
    }


def test_phase_routing_envelope_propagates_only_declared_phase_scalar_bounds(monkeypatch):
    calls = []

    def fake_routing(scenario, config, **_kwargs):
        calls.append(config)
        return FieldPhaseRoutingResult(
            scenario, config, FieldApplicability("accepted", uncertainty_complete=True),
            None, None, object(),
        )

    monkeypatch.setattr(field_phase_routing, "run_field_phase_routing", fake_routing)
    uncertainty = FieldPhaseRoutingUncertainty(
        pool_area_m2=BoundedValue(0.5, 0.4, 0.6, unit="m2", source="pool-footprint-A"),
        evaporation_coefficient_m2_s=BoundedValue(
            1.0e-7, 5.0e-8, 2.0e-7, unit="m2/s", source="evaporation-A",
        ),
    )
    envelope = run_field_phase_routing_envelope(
        _scenario(), _config(), phase_uncertainty=uncertainty, max_cases=4,
    )

    assert len(envelope.cases) == 4
    assert {
        (case.corner["phase_pool_area_m2"], case.corner["phase_evaporation_coefficient_m2_s"])
        for case in envelope.cases
    } == {(0.4, 5.0e-8), (0.4, 2.0e-7), (0.6, 5.0e-8), (0.6, 2.0e-7)}
    assert {config.pool_area_m2 for config in calls} == {0.4, 0.6}
    assert {config.evaporation_coefficient_m2_s for config in calls} == {5.0e-8, 2.0e-7}
    assert envelope.phase_uncertainty is uncertainty


def test_phase_routing_uncertainty_requires_nominal_anchor_and_explicit_source():
    with pytest.raises(ValueError, match="nominal must match"):
        FieldPhaseRoutingUncertainty(
            pool_area_m2=BoundedValue(0.6, 0.5, 0.7, unit="m2", source="pool-footprint-A"),
        ).config_corners(_config())
    with pytest.raises(ValueError, match="explicit bounded-value source"):
        FieldPhaseRoutingUncertainty(
            pool_area_m2=BoundedValue(0.5, 0.4, 0.6, unit="m2"),
        )


def test_phase_routing_population_corners_preserve_the_declared_simplex(monkeypatch):
    calls = []

    def fake_routing(scenario, config, **_kwargs):
        calls.append(tuple(config.phase_model_options["droplet_classes"]))
        return FieldPhaseRoutingResult(
            scenario, config, FieldApplicability("accepted", uncertainty_complete=True),
            None, None, object(),
        )

    monkeypatch.setattr(field_phase_routing, "run_field_phase_routing", fake_routing)
    nominal = (
        DropletClass(1.0e-4, 0.25), DropletClass(1.0e-3, 0.75),
    )
    config = replace(
        _config(), phase_model_options={"droplet_classes": nominal},
    )
    uncertainty = FieldPhaseRoutingUncertainty(
        droplet_population_corners=(
            nominal,
            (DropletClass(8.0e-5, 0.5), DropletClass(1.2e-3, 0.5)),
        ),
        droplet_population_uncertainty_evidence_id="spray-population-envelope-A",
    )
    envelope = run_field_phase_routing_envelope(
        _scenario(), config, phase_uncertainty=uncertainty, max_cases=2,
    )

    assert len(envelope.cases) == 2
    assert len(calls) == 2
    assert {case.corner["phase_droplet_class_0_mass_fraction"] for case in envelope.cases} == {
        0.25, 0.5,
    }
    assert all(
        sum(item.mass_fraction for item in population) == pytest.approx(1.0)
        for population in calls
    )
    assert envelope.phase_uncertainty is uncertainty
    assert uncertainty.as_record()["droplet_population"]["evidence_id"] == (
        "spray-population-envelope-A"
    )


def test_phase_routing_population_corners_require_the_nominal_simplex():
    nominal = (
        DropletClass(1.0e-4, 0.25), DropletClass(1.0e-3, 0.75),
    )
    uncertainty = FieldPhaseRoutingUncertainty(
        droplet_population_corners=(
            (DropletClass(8.0e-5, 0.5), DropletClass(1.2e-3, 0.5)),
        ),
        droplet_population_uncertainty_evidence_id="spray-population-envelope-A",
    )
    config = replace(_config(), phase_model_options={"droplet_classes": nominal})
    with pytest.raises(ValueError, match="include the nominal population"):
        uncertainty.config_corners(config)
    with pytest.raises(ValueError, match="must be unique"):
        FieldPhaseRoutingUncertainty(
            droplet_population_corners=(nominal, nominal),
            droplet_population_uncertainty_evidence_id="spray-population-envelope-A",
        )


@pytest.mark.slow
def test_field_phase_routing_real_short_flash_preserves_mass_ledger():
    base = _scenario()
    scenario = FieldScenario(
        source=replace(base.source, duration_s=0.002),
        weather=WeatherState(speed_m_s=BoundedValue(2.5), direction_deg=BoundedValue(270.0)),
        surface=base.surface,
        temporal_mode="transient",
    )
    config = FieldPhaseRoutingConfig(
        post_release_duration_s=0.1,
        puff_duration_s=0.05,
        pool_area_m2=0.5,
        pool_time_step_s=0.05,
        evaporation_coefficient_m2_s=1.0e-7,
        gas_model_options={
            "maximum_nearfield_distance": 0.2,
            "radial_points": 41,
            "nearfield_maximum_step": 0.001,
            "nearfield_relative_tolerance": 2.0e-6,
            "crosswind_maximum_distance": 2.0,
            "crosswind_maximum_step": 0.02,
            "puff_time_step": 0.01,
        },
        phase_model_options={"maximum_droplet_time_s": 0.1},
    )
    result = run_field_phase_routing(scenario, config)
    assert result.completed
    assert result.coupled_result is not None
    assert result.coupled_result.conservative
    assert result.coupled_result.hydrogen_mass_residual_kg == pytest.approx(0.0, abs=1.0e-12)
    assert not result.coupled_result.atmospherically_complete
    # This very short release has no resolved timed pool ledger. The pool
    # handoff must refuse it rather than manufacture a pool-vapour schedule.
    with pytest.raises(ValueError, match="time-resolved dynamic-pool ledger"):
        pool_vapour_schedule_from_phase_routing(result)
