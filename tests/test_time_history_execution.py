import json

import pytest

from degali.addons.field_evidence_audit import audit_field_evidence
from degali.addons.field_evidence_manifest import FieldEvidenceManifest
from degali.addons.time_history_execution import run_manifest_bound_time_history
from degali.addons.transient_receptor import (
    FixedReceptor,
    SourceHistory,
    WindHistory,
)


def _manifest(root):
    root.mkdir(parents=True, exist_ok=True)
    paths = {
        "source_boundary": root / "source.csv",
        "weather": root / "weather.csv",
        "obstacle_geometry": root / "obstacle.json",
        "receptor_observations": root / "receptors.csv",
        "common_clock": root / "clock.csv",
    }
    paths["source_boundary"].write_text(
        "source_boundary_id,mass_flow_kg_s\nsource-a,0.1\n", encoding="utf-8"
    )
    paths["weather"].write_text(
        "weather_id,wind_speed_m_s\nweather-a,1.0\n", encoding="utf-8"
    )
    paths["obstacle_geometry"].write_text(
        json.dumps({"obstacle_geometry_id": "obstacle-a"}), encoding="utf-8"
    )
    paths["receptor_observations"].write_text(
        "sensor,x_downwind_m,y_crosswind_m,height_m,observed_mole_fraction,common_clock_id\n"
        "s1,1,0,0,0.1,clock-a\n", encoding="utf-8"
    )
    paths["common_clock"].write_text("common_clock_id\nclock-a\n", encoding="utf-8")
    audit = audit_field_evidence(root)
    assert audit.status == "candidate_complete"
    return FieldEvidenceManifest.from_audit(
        audit,
        manifest_id="manifest-a",
        event_id="event-a",
        selected_paths={key: str(value) for key, value in paths.items()},
        dataset_id="dataset-a",
        observed_row_count=1,
        source_boundary_id="source-a",
        weather_id="weather-a",
        obstacle_geometry_id="obstacle-a",
        receptor_geometry_id="receptors-a",
        temporal_operator_id="instantaneous-test",
        common_clock_id="clock-a",
    )


def _histories():
    source = SourceHistory([0.0, 1.0, 2.0], [1.0, 0.0, 0.0])
    wind = WindHistory([0.0, 1.0, 2.0], [1.0, 1.0, 1.0], [270.0, 270.0, 270.0])
    receptors = [FixedReceptor("s1", 1.0, 0.0, 0.0)]
    return source, wind, receptors


def _run(manifest, **overrides):
    source, wind, receptors = _histories()
    values = dict(
        source_event_id="event-a",
        weather_event_id="event-a",
        source_common_clock_id="clock-a",
        weather_common_clock_id="clock-a",
        receptor_geometry_id="receptors-a",
        temporal_operator_id="instantaneous-test",
        source_history=source,
        wind_history=wind,
        receptors=receptors,
        packet_sampler=lambda packet, age, x, y, z: 0.1 if abs(x) < 1.0e-12 else 0.0,
        observation_time_s=[0.0, 0.5, 1.5, 2.0],
    )
    values.update(overrides)
    return run_manifest_bound_time_history(manifest, **values)


def test_manifest_bound_history_runs_but_retains_conditional_metadata_gate(tmp_path):
    result = _run(_manifest(tmp_path / "evidence"))

    assert result.gate.status == "conditional"
    assert result.gate.executable
    assert result.packet_count == 1
    assert result.released_mass_kg == pytest.approx(1.0)
    assert len(result.traces) == 1
    assert "manifest_metadata_conditional" in result.gate.gate_codes
    assert "time_history_valid" in result.gate.gate_codes
    assert result.as_record()["schema"] == "degali.source-weather-receptor-time-history-execution.v1"


def test_manifest_bound_history_withholds_id_mismatch_without_traces(tmp_path):
    result = _run(_manifest(tmp_path / "evidence"), weather_event_id="other-event")

    assert result.gate.status == "withheld"
    assert not result.gate.executable
    assert result.traces == ()
    assert "event_id_mismatch" in result.gate.gate_codes
    assert "execution_withheld" in result.gate.gate_codes


def test_manifest_bound_history_withholds_file_digest_drift(tmp_path):
    root = tmp_path / "evidence"
    manifest = _manifest(root)
    (root / "weather.csv").write_text(
        "weather_id,wind_speed_m_s\nweather-a,2.0\n", encoding="utf-8"
    )

    result = _run(manifest)

    assert result.gate.status == "withheld"
    assert result.traces == ()
    assert result.packet_count == 0
    assert result.gate.gate_codes == ("manifest_files_withheld", "execution_withheld")


def test_manifest_bound_history_withholds_uncovered_common_clock(tmp_path):
    source, _wind, receptors = _histories()
    short_wind = WindHistory([0.0, 0.25], [1.0, 1.0], [270.0, 270.0])
    manifest = _manifest(tmp_path / "evidence")
    result = run_manifest_bound_time_history(
        manifest,
        source_event_id="event-a", weather_event_id="event-a",
        source_common_clock_id="clock-a", weather_common_clock_id="clock-a",
        receptor_geometry_id="receptors-a", temporal_operator_id="instantaneous-test",
        source_history=source, wind_history=short_wind, receptors=receptors,
        packet_sampler=lambda packet, age, x, y, z: 0.0,
    )

    assert result.gate.status == "withheld"
    assert "time_history_withheld" in result.gate.gate_codes
    assert result.traces == ()
