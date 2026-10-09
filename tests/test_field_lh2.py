import math
from dataclasses import replace
from types import SimpleNamespace

import pytest

from degali.addons.field_contracts import BoundedValue, FieldScenario, ReleaseSource
from degali.addons.cryogenic_blowdown import CryogenicBlowdownConfig, run_cryogenic_blowdown
import degali.addons.cryogenic_blowdown as cryogenic_blowdown
from degali.addons.field_lh2 import (
    build_lh2_saturation_table_for_release,
    direct_vapour_schedule_from_cryogenic_blowdown,
    lh2_flash_source_from_release,
    pressure_driven_lh2_mass_flow,
    prepare_field_lh2_flash,
    release_with_pressure_driven_lh2_mass_flow,
)
import degali.addons.field_lh2 as field_lh2


def _lh2_release(**changes):
    values = {
        "fluid": "lh2",
        "upstream_pressure": BoundedValue(0.4e6, unit="Pa"),
        "upstream_temperature": BoundedValue(26.084, unit="K"),
        "mass_flow_kg_s": BoundedValue(0.265, unit="kg/s"),
        # The geometric area is chosen so Cd*A gives a 12 mm effective throat.
        "opening_area_m2": BoundedValue(math.pi * 0.012**2 / 4.0 / 0.8, unit="m2"),
        "discharge_coefficient": BoundedValue(0.8),
        "liquid_fraction": BoundedValue(0.922),
        "flash_model": "homogeneous_equilibrium",
    }
    values.update(changes)
    return ReleaseSource(**values)


def test_field_lh2_bridge_applies_cd_once_and_closes_phase_mass():
    result = lh2_flash_source_from_release(
        _lh2_release(), ambient_temperature_k=293.15
    )
    assert result.effective_diameter_m == pytest.approx(0.012)
    assert result.flash.mass_flow == pytest.approx(0.265)
    assert result.mass_residual_kg_s == pytest.approx(0.0, abs=1.0e-12)
    assert result.conservative
    assert result.flash.liquid_mass_flow > 0.0
    assert result.phase_mass_residual == pytest.approx(0.0, abs=1.0e-12)
    assert result.momentum_residual < 1.0e-8
    assert result.energy_residual < 1.0e-8
    assert result.closure_warnings == ()
    assert result.closure_tolerances == {
        "mass_partition_relative": 1.0e-10,
        "phase_mass_relative": 1.0e-10,
        "momentum_relative": 1.0e-8,
        "energy_relative": 1.0e-8,
    }


def test_pressure_driven_lh2_mass_flow_crosses_release_and_ambient_bounds():
    release = _lh2_release(
        mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"),
        upstream_pressure=BoundedValue(0.4e6, 0.38e6, 0.42e6, "Pa", "pressure-review"),
    )
    bound = pressure_driven_lh2_mass_flow(
        release,
        ambient_pressure_pa=101325.0,
        ambient_pressure_bounds_pa=(100000.0, 102000.0),
        source_id="orifice-review-A",
    )

    assert bound.unit == "kg/s"
    assert bound.source == "orifice-review-A"
    assert 0.0 < bound.lower < bound.nominal < bound.upper
    attached = release_with_pressure_driven_lh2_mass_flow(
        release,
        ambient_pressure_pa=101325.0,
        ambient_pressure_bounds_pa=(100000.0, 102000.0),
        source_id="orifice-review-A",
    )
    assert attached.mass_flow_kg_s == bound
    assert attached.pressure_driven_mass_flow is not None
    assert attached.pressure_driven_mass_flow.source_id == "orifice-review-A"
    assert attached.pressure_driven_mass_flow.ambient_pressure_pa.lower == pytest.approx(100000.0)
    assert attached.pressure_driven_mass_flow.ambient_pressure_pa.upper == pytest.approx(102000.0)
    # The derived rate is not an independent corner dimension: the scenario
    # recomputes it from the typed ambient/source boundary.
    assert "mass_flow_kg_s" not in FieldScenario(source=attached).uncertainty_fields()
    assert "pressure_driven_ambient_pressure_pa" in FieldScenario(source=attached).uncertainty_fields()
    assert attached.metadata["mass_flow_derivation"] == (
        "pressure_driven_homogeneous_equilibrium_throat"
    )
    assert attached.metadata["mass_flow_derivation_source_id"] == "orifice-review-A"


def test_pressure_driven_lh2_mass_flow_does_not_overwrite_a_measured_rate():
    with pytest.raises(ValueError, match="will not overwrite"):
        release_with_pressure_driven_lh2_mass_flow(_lh2_release())


@pytest.mark.parametrize(
    "changes, message",
    [
        ({"fluid": "methane", "mass_flow_kg_s": BoundedValue(0.0, unit="kg/s")}, "hydrogen releases only"),
        ({"pressure_reference": "gauge", "mass_flow_kg_s": BoundedValue(0.0, unit="kg/s")}, "absolute upstream pressure"),
    ],
)
def test_pressure_driven_lh2_mass_flow_rejects_ambiguous_source_boundaries(changes, message):
    with pytest.raises(ValueError, match=message):
        pressure_driven_lh2_mass_flow(_lh2_release(**changes))


def test_pressure_driven_lh2_mass_flow_rejects_an_ambient_corner_without_driving_pressure():
    release = _lh2_release(mass_flow_kg_s=BoundedValue(0.0, unit="kg/s"))
    with pytest.raises(ValueError, match="corner cannot be evaluated"):
        pressure_driven_lh2_mass_flow(
            release,
            ambient_pressure_pa=101325.0,
            ambient_pressure_bounds_pa=(101325.0, 500000.0),
        )


def test_field_lh2_preparation_blocks_a_nonconservative_flash_handoff(monkeypatch):
    release = _lh2_release()
    nominal = lh2_flash_source_from_release(release)
    invalid = replace(
        nominal,
        flash=replace(nominal.flash, energy_residual=1.0e-3),
    )
    monkeypatch.setattr(
        field_lh2,
        "lh2_flash_source_from_release",
        lambda *args, **kwargs: invalid,
    )

    prepared = prepare_field_lh2_flash(FieldScenario(source=release))

    assert prepared.flash_result is None
    assert dict(prepared.flash_closure_diagnostics)["energy_residual"] == pytest.approx(1.0e-3)
    assert prepared.applicability.status == "blocked"
    assert any("energy residual" in reason for reason in prepared.applicability.reasons)


def test_field_lh2_flash_result_rejects_negative_phase_flow_or_residual():
    nominal = lh2_flash_source_from_release(_lh2_release())
    with pytest.raises(ValueError, match="liquid_mass_flow"):
        replace(
            nominal,
            flash=replace(nominal.flash, liquid_mass_flow=-1.0e-6),
        )
    with pytest.raises(ValueError, match="energy_residual"):
        replace(
            nominal,
            flash=replace(nominal.flash, energy_residual=-1.0e-6),
        )


def test_blocked_field_lh2_preparation_cannot_expose_a_flash_plane():
    nominal = lh2_flash_source_from_release(_lh2_release())
    with pytest.raises(ValueError, match="cannot expose a flash plane"):
        field_lh2.FieldLH2SourcePreparation(
            FieldScenario(source=_lh2_release()),
            field_lh2.FieldApplicability(
                "blocked", reasons=("test block",), uncertainty_complete=False,
            ),
            nominal,
        )


@pytest.mark.parametrize(
    "changes, message",
    [
        ({"pressure_reference": "gauge"}, "absolute"),
        ({"liquid_fraction_basis": "post_flash"}, "upstream liquid fraction"),
        ({"flash_model": "declared"}, "unsupported flash_model"),
    ],
)
def test_field_lh2_bridge_refuses_ambiguous_source_interpretations(changes, message):
    with pytest.raises(ValueError, match=message):
        lh2_flash_source_from_release(_lh2_release(**changes))


def test_field_lh2_preparation_returns_blocked_result_without_flash_plane():
    prepared = prepare_field_lh2_flash(
        FieldScenario(source=_lh2_release(flash_model="declared")),
    )
    assert not prepared.prepared
    assert prepared.flash_result is None
    assert prepared.applicability.status == "blocked"
    assert "unsupported flash_model" in prepared.applicability.reasons[0]


def test_field_lh2_preparation_keeps_conditional_status_with_explicit_flash():
    prepared = prepare_field_lh2_flash(FieldScenario(source=_lh2_release()))
    assert prepared.prepared
    assert prepared.flash_result is not None
    assert prepared.applicability.status == "conditional"


def test_lh2_saturation_table_matches_reference_flash_within_declared_tolerance():
    release = _lh2_release()
    table = build_lh2_saturation_table_for_release(release, nodes=161)
    direct = lh2_flash_source_from_release(release, ambient_temperature_k=293.15)
    tabular = lh2_flash_source_from_release(
        release, ambient_temperature_k=293.15, property_table=table
    )
    for name in (
        "postflash_quality", "postflash_velocity", "postflash_density",
        "droplet_diameter", "jet_weber",
    ):
        assert getattr(tabular.flash, name) == pytest.approx(
            getattr(direct.flash, name), rel=1.0e-4
        )
    with pytest.raises(ValueError, match="outside the LH2 saturation-table domain"):
        table.props("H", "T", table.maximum_temperature_k + 0.01, "Q", 0.0, "Hydrogen")


def test_lh2_saturation_table_covers_declared_ambient_pressure_bounds():
    coolprop = pytest.importorskip("CoolProp.CoolProp")
    release = _lh2_release()
    lower_pressure, upper_pressure = 80_000.0, 120_000.0
    table = build_lh2_saturation_table_for_release(
        release,
        ambient_pressure_pa=101325.0,
        ambient_pressure_bounds_pa=(lower_pressure, upper_pressure),
        nodes=161,
    )

    for pressure in (lower_pressure, upper_pressure):
        saturation_temperature = float(
            coolprop.PropsSI("T", "P", pressure, "Q", 0, "Hydrogen")
        )
        assert table.minimum_temperature_k <= saturation_temperature <= table.maximum_temperature_k


def test_blowdown_adapter_exports_only_direct_postflash_vapour_to_a_schedule():
    blowdown = run_cryogenic_blowdown(CryogenicBlowdownConfig(
        vessel_volume_m3=2.815e-3,
        nozzle_diameter_m=4.0e-3,
        initial_temperature_k=80.0,
        initial_pressure_pa=20.0e6,
        discharge_coefficient=0.7,
        duration_s=0.04,
        time_step_s=0.02,
    ))
    adapted = direct_vapour_schedule_from_cryogenic_blowdown(blowdown)

    assert adapted.schedule.time_s == pytest.approx((0.0, 0.02, 0.04))
    assert adapted.schedule.rate_kg_s[-1] == 0.0
    assert adapted.direct_vapour_mass_kg == pytest.approx(
        adapted.schedule.released_mass_kg
    )
    assert adapted.unrouted_postflash_liquid_mass_kg > 0.0
    assert adapted.direct_vapour_mass_kg + adapted.unrouted_postflash_liquid_mass_kg == pytest.approx(
        adapted.total_discharged_mass_kg, rel=1.0e-8
    )
    assert "liquid is excluded" in adapted.warnings[0]


def test_blowdown_adapter_blocks_an_unclosed_flash_interval(monkeypatch):
    blowdown = run_cryogenic_blowdown(CryogenicBlowdownConfig(
        vessel_volume_m3=2.815e-3,
        nozzle_diameter_m=4.0e-3,
        initial_temperature_k=80.0,
        initial_pressure_pa=20.0e6,
        discharge_coefficient=0.7,
        duration_s=0.04,
        time_step_s=0.02,
    ))
    state = blowdown.states[0]
    nominal = cryogenic_blowdown.blowdown_state_to_flashing_droplet_source(
        blowdown.config, state, ambient_temperature_k=295.0,
        droplet_size_coefficient=15.0,
    )
    invalid = replace(nominal, energy_residual=1.0e-3)
    monkeypatch.setattr(
        cryogenic_blowdown,
        "blowdown_state_to_flashing_droplet_source",
        lambda *args, **kwargs: invalid,
    )

    with pytest.raises(ValueError, match="not conservatively closed"):
        direct_vapour_schedule_from_cryogenic_blowdown(blowdown)


def test_blowdown_adapter_rejects_negative_flash_phase_flow(monkeypatch):
    blowdown = run_cryogenic_blowdown(CryogenicBlowdownConfig(
        vessel_volume_m3=2.815e-3,
        nozzle_diameter_m=4.0e-3,
        initial_temperature_k=80.0,
        initial_pressure_pa=20.0e6,
        discharge_coefficient=0.7,
        duration_s=0.04,
        time_step_s=0.02,
    ))
    invalid = SimpleNamespace(
        vapour_mass_flow=0.0,
        liquid_mass_flow=-1.0e-6,
        mass_residual=0.0,
        momentum_residual=0.0,
        energy_residual=0.0,
    )
    monkeypatch.setattr(
        cryogenic_blowdown,
        "blowdown_state_to_flashing_droplet_source",
        lambda *args, **kwargs: invalid,
    )

    with pytest.raises(ValueError, match="negative"):
        direct_vapour_schedule_from_cryogenic_blowdown(blowdown)


def test_blowdown_adapter_does_not_clip_an_overreleased_vapour_schedule(monkeypatch):
    blowdown = run_cryogenic_blowdown(CryogenicBlowdownConfig(
        vessel_volume_m3=2.815e-3,
        nozzle_diameter_m=4.0e-3,
        initial_temperature_k=80.0,
        initial_pressure_pa=20.0e6,
        discharge_coefficient=0.7,
        duration_s=0.04,
        time_step_s=0.02,
    ))
    state = blowdown.states[0]
    invalid = SimpleNamespace(
        vapour_mass_flow=2.0 * state.mass_flow_kg_s,
        liquid_mass_flow=0.0,
        mass_residual=0.0,
        momentum_residual=0.0,
        energy_residual=0.0,
    )
    monkeypatch.setattr(
        cryogenic_blowdown,
        "blowdown_state_to_flashing_droplet_source",
        lambda *args, **kwargs: invalid,
    )
    # Isolate the aggregate ledger guard from the per-interval closure guard.
    monkeypatch.setattr(field_lh2, "_raw_flash_closure_warnings", lambda *args: ())

    with pytest.raises(ValueError, match="discharged-mass ledger"):
        direct_vapour_schedule_from_cryogenic_blowdown(blowdown)


def test_blowdown_adapter_handles_only_the_explicit_hem_two_phase_branch():
    blowdown = run_cryogenic_blowdown(CryogenicBlowdownConfig(
        vessel_volume_m3=2.815e-3,
        nozzle_diameter_m=4.0e-3,
        initial_temperature_k=80.0,
        initial_pressure_pa=20.0e6,
        discharge_coefficient=0.7,
        duration_s=1.2,
        time_step_s=0.02,
        two_phase_withdrawal="homogeneous",
    ))
    adapted = direct_vapour_schedule_from_cryogenic_blowdown(blowdown)

    assert any(state.vapour_quality is not None for state in blowdown.states)
    assert adapted.direct_vapour_mass_kg > 0.0
    assert "HEM pressure-thrust flash" in adapted.warnings[-1]
