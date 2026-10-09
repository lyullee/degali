import json

import pytest

from degali.addons.source_state_ledger import SourceState, SourceStateLedger


def _ledger(liquid_fraction=0.0, final_rate=0.0):
    state = lambda time, rate, liquid=liquid_fraction: SourceState(
        time_s=time,
        h2_rate_kg_s=rate,
        h2_mass_fraction=0.9,
        temperature_k=25.0,
        density_kg_m3=2.0,
        area_m2=0.02,
        velocity_m_s=10.0,
        liquid_fraction=liquid,
        height_m=0.5,
        bearing_to_deg=350.0,
    )
    return SourceStateLedger(
        substance="hydrogen",
        stage="post_flash_atmospheric",
        states=(state(0.0, 1.0), state(1.0, 0.5), state(2.0, final_rate)),
        duration_s=2.0,
        rate_operator="piecewise_constant",
        observation_operator="common-clock-275s",
        metadata={"event_id": "ffi-test4"},
    )


def test_ledger_closes_h2_and_carrier_mass_and_converts_to_source_history():
    ledger = _ledger()

    assert ledger.has_zero_endpoint
    assert ledger.released_h2_mass_kg == pytest.approx(1.5)
    assert ledger.released_carrier_mass_kg == pytest.approx(1.5 / 0.9)
    history = ledger.to_source_history()
    assert history.released_mass_kg() == pytest.approx(1.5)
    assert history.source_direction_to_deg == (350.0, 350.0, 350.0)
    assert ledger.as_record()["schema"] == "degali.source-state-ledger.v1"
    assert ledger.as_record()["transport_ready"] is True


def test_liquid_bearing_ledger_cannot_bypass_phase_routing():
    with pytest.raises(ValueError, match="flash/rainout/pool"):
        _ledger(liquid_fraction=0.2).to_source_history()
    with pytest.raises(ValueError, match="flash/rainout/pool"):
        _ledger(liquid_fraction=0.2).to_source_rate_schedule()


def test_to_source_rate_schedule_preserves_hydrogen_mass_and_source_id():
    schedule = _ledger().to_source_rate_schedule(source_id="impact-01")
    assert schedule.source_id == "impact-01"
    assert schedule.time_s == (0.0, 1.0, 2.0)
    assert schedule.rate_kg_s == (1.0, 0.5, 0.0)
    assert schedule.released_mass_kg == pytest.approx(1.5)
    assert schedule.has_zero_endpoint


def test_ledger_requires_explicit_zero_endpoint_and_multiple_clock_states():
    with pytest.raises(ValueError, match="zero H2 rate"):
        _ledger(final_rate=0.1).to_source_history()
    state = SourceState(
        time_s=0.0, h2_rate_kg_s=0.1, h2_mass_fraction=0.9,
        temperature_k=25.0, density_kg_m3=2.0, area_m2=0.02,
    )
    with pytest.raises(ValueError, match="at least two"):
        SourceStateLedger("hydrogen", "static", (state,), 1.0)


def test_ledger_json_roundtrip(tmp_path):
    path = tmp_path / "ledger.json"
    original = _ledger()
    original.write_json(path)
    loaded = SourceStateLedger.read_json(path)
    assert loaded.as_record() == json.loads(path.read_text(encoding="utf-8"))
    assert loaded.metadata["event_id"] == "ffi-test4"
