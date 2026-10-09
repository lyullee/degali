import pytest

from degali.validation.open_channel_h2 import (
    read_open_channel_h2_csv,
    threshold_window,
)


def _write_run(path):
    path.write_text(
        "flow time [s],mass flow meter 1 [g/s],h2 sensor time,"
        "sensor 1 h2 concentration [%],sensor 2 h2 concentration [%]\n"
        "0,0.1,0,0.1,0.2\n"
        "0.01,0.2,0.5,1.2,0.3\n"
        "0.02,0.3,1.0,2.1,0.4\n"
        "0.03,0.2,1.5,0.8,0.3\n",
        encoding="utf-8",
    )


def test_open_channel_h2_local_intake_keeps_independent_clocks_and_scope(tmp_path):
    path = tmp_path / "23_FFI_P101_T00014.csv"
    _write_run(path)
    run = read_open_channel_h2_csv(path)

    assert run.flow_time_s.tolist() == pytest.approx([0, .01, .02, .03])
    assert run.sensor_time_s.tolist() == pytest.approx([0, .5, 1, 1.5])
    assert run.sensor_names == ("sensor 1", "sensor 2")
    assert run.flow_column == "mass flow meter 1 [g/s]"
    assert not run.quantitative_lh2_pool_validation_allowed
    assert run.validation_scope == "hydrogen_source_and_sensor_timing_in_channel_only"

    window = threshold_window(run, "sensor 1", threshold_percent=1.0)
    assert window.arrival_s == pytest.approx(.5)
    assert window.peak_time_s == pytest.approx(1.0)
    assert window.departure_s == pytest.approx(1.0)
    assert window.peak_percent == pytest.approx(2.1)
    assert window.duration_s == pytest.approx(.5)
    assert threshold_window(run, "sensor 2", threshold_percent=1.0) is None


def test_open_channel_h2_accepts_the_public_meter_2_variant(tmp_path):
    path = tmp_path / "23_FFI_P101_T00026.csv"
    path.write_text(
        "flow time [s],mass flow meter 2 [g/s],h2 sensor time,"
        "sensor 1 h2 concentration [%]\n"
        "0,0.1,0,0.1\n"
        "0.01,0.2,0.5,1.2\n"
        "0.02,0.3,1.0,2.1\n",
        encoding="utf-8",
    )

    run = read_open_channel_h2_csv(path)

    assert run.flow_column == "mass flow meter 2 [g/s]"
    assert run.mass_flow_g_s.tolist() == pytest.approx([.1, .2, .3])


def test_open_channel_h2_rejects_missing_schema_and_bad_threshold(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("flow time [s],h2 sensor time\n0,0\n1,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing required columns"):
        read_open_channel_h2_csv(bad)

    path = tmp_path / "valid.csv"
    _write_run(path)
    with pytest.raises(ValueError, match="threshold_percent"):
        threshold_window(read_open_channel_h2_csv(path), "sensor 1", threshold_percent=0.0)
