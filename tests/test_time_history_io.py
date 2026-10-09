import hashlib

import pytest

from degali.addons.time_history_io import (
    SourceHistoryCsvMap,
    WindHistoryCsvMap,
    read_source_history_csv,
    read_wind_history_csv,
)


def test_source_history_csv_imports_si_channels_and_pins_digest(tmp_path):
    path = tmp_path / "source.csv"
    path.write_text(
        "time_s,mass_rate_kg_s,pressure_pa,temperature_k,source_direction_to_deg\n"
        "0,1.0,100000,25.0,350\n"
        "1,0.5,99000,25.5,10\n"
        "2,0.0,98000,26.0,20\n",
        encoding="utf-8",
    )

    imported = read_source_history_csv(
        path,
        SourceHistoryCsvMap(
            time_s_column="time_s",
            mass_rate_kg_s_column="mass_rate_kg_s",
            pressure_pa_column="pressure_pa",
            temperature_k_column="temperature_k",
            source_direction_to_deg_column="source_direction_to_deg",
            event_id="ffi-test4",
            evidence_id="source-history-001",
            rate_operator="linear",
        ),
    )

    assert imported.row_count == 3
    assert imported.sha256 == hashlib.sha256(path.read_bytes()).hexdigest()
    assert imported.event_id == "ffi-test4"
    assert imported.history.released_mass_kg() == pytest.approx(1.0)
    packets = imported.history.to_packets(subdivisions_per_interval=2)
    assert packets[0].source_direction_to_deg == pytest.approx(355.0)
    assert imported.as_record()["released_mass_kg"] == pytest.approx(1.0)


def test_wind_history_csv_imports_and_normalizes_direction(tmp_path):
    path = tmp_path / "wind.csv"
    path.write_text(
        "time_s,speed_m_s,direction_from_deg\n"
        "0,2.0,270\n"
        "1,2.5,-90\n"
        "2,3.0,450\n",
        encoding="utf-8",
    )

    imported = read_wind_history_csv(
        path,
        WindHistoryCsvMap(
            time_s_column="time_s",
            speed_m_s_column="speed_m_s",
            direction_from_deg_column="direction_from_deg",
            event_id="ffi-test4",
            evidence_id="wind-history-001",
        ),
    )

    _time, _speed, direction = imported.history.arrays()
    assert direction.tolist() == pytest.approx([270.0, 270.0, 90.0])
    assert imported.as_record()["speed_range_m_s"] == pytest.approx([2.0, 3.0])


def test_source_history_csv_rejects_static_source_state_without_time_history(tmp_path):
    path = tmp_path / "source_state.csv"
    path.write_text(
        "reported_rate_kg_s,P04_pressure_pa,P04_temperature_K\n"
        "0.1,200000,25\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="at least two data rows"):
        read_source_history_csv(
            path,
            SourceHistoryCsvMap(
                time_s_column="time_s",
                mass_rate_kg_s_column="reported_rate_kg_s",
                event_id="ffi-test4",
                evidence_id="source-state-only",
            ),
        )


@pytest.mark.parametrize(
    ("csv_text", "message"),
    [
        (
            "time_s,mass_rate_kg_s\n0,1\n1,0,unexpected\n",
            "more values",
        ),
        (
            "time_s,time_s,mass_rate_kg_s\n0,0,1\n1,1,0\n",
            "unique",
        ),
        (
            "time_s\n0\n1\n",
            "missing columns",
        ),
    ],
)
def test_history_csv_importer_rejects_ambiguous_or_incomplete_rows(
    tmp_path, csv_text, message
):
    path = tmp_path / "invalid.csv"
    path.write_text(csv_text, encoding="utf-8")
    mapping = SourceHistoryCsvMap(
        time_s_column="time_s",
        mass_rate_kg_s_column="mass_rate_kg_s",
        event_id="event",
        evidence_id="evidence",
    )
    with pytest.raises(ValueError, match=message):
        read_source_history_csv(path, mapping)
