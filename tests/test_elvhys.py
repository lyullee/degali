import math

import pytest

from degali.validation.elvhys import read_elvhys_stream, read_elvhys_test


def _write_stream(directory, test, stream, body):
    path = directory / f"ELE402HSE{test:03d}{stream}20250101.csv"
    path.write_text(body, encoding="utf-8")
    return path


def test_elvhys_local_stream_preserves_clock_channels_and_scope(tmp_path):
    _write_stream(tmp_path, 10, "CONC", "Time,H2Bottom1,H2Side1\n0,0.1,\n0.05,1.2,2.3\n")
    _write_stream(tmp_path, 10, "TEMP", "Time,Bottom1Temp\n0,280\n0.05,270\n")
    _write_stream(tmp_path, 10, "PRES", "Time,PT2Nozzle\n0,0\n0.05,2\n")
    _write_stream(tmp_path, 10, "FLMT", "Time,FanFlowMeter\n0,500\n0.05,500\n")
    test = read_elvhys_test(tmp_path, 10)

    assert test.stream("CONC").time_s.tolist() == [0.0, 0.05]
    assert math.isnan(test.stream("CONC").channels["H2Side1"][0])
    assert test.stream("FLMT").channel_names == ("FanFlowMeter",)
    assert not test.hydrogen_mass_flow_measured
    assert not test.quantitative_outdoor_lh2_pool_validation_allowed
    assert test.validation_scope == "confined_cryogenic_h2_observation_only"


def test_elvhys_local_intake_requires_declared_streams_and_valid_time(tmp_path):
    path = _write_stream(tmp_path, 10, "CONC", "Time,H2Bottom1\n1,1\n0,2\n")
    with pytest.raises(ValueError, match="nondecreasing"):
        read_elvhys_stream(path, stream="CONC")

    _write_stream(tmp_path, 11, "CONC", "Time,H2Bottom1\n0,1\n0.05,2\n")
    with pytest.raises(ValueError, match="required TEMP"):
        read_elvhys_test(tmp_path, 11)

    with pytest.raises(ValueError, match="positive integer"):
        read_elvhys_test(tmp_path, 0)


def test_elvhys_file_match_does_not_confuse_padded_test_numbers(tmp_path):
    # Real ELVHYS exports contain both HSE0024 (test 24) and HSE0043
    # (test 43).  Test 4 must not consume the latter through a prefix match.
    for stream in ("CONC", "TEMP", "PRES", "FLMT"):
        _write_stream(tmp_path, 4, stream, "Time,Value\n0,1\n1,2\n")
    path = tmp_path / "ELE402HSE0043TEMP20250101.csv"
    path.write_text("Time,Value\n0,3\n1,4\n", encoding="utf-8")
    test = read_elvhys_test(tmp_path, 4)
    assert test.stream("TEMP").source_path.name == "ELE402HSE004TEMP20250101.csv"


def test_elvhys_accepts_extra_zero_padded_public_filename(tmp_path):
    for stream in ("CONC", "TEMP", "PRES", "FLMT"):
        path = tmp_path / f"ELE402HSE0024{stream}20250101.csv"
        path.write_text("Time,Value\n0,1\n1,2\n", encoding="utf-8")
    test = read_elvhys_test(tmp_path, 24)
    assert test.stream("CONC").source_path.name.startswith("ELE402HSE0024")


def test_elvhys_trims_reserved_trailing_blank_rows(tmp_path):
    path = _write_stream(
        tmp_path, 12, "CONC", "Time,H2Bottom1\n0,1\n0.05,2\n,\n,\n"
    )
    stream = read_elvhys_stream(path, stream="CONC")
    assert stream.time_s.tolist() == [0.0, 0.05]
    assert stream.channels["H2Bottom1"].tolist() == [1.0, 2.0]


def test_elvhys_optional_empty_flow_placeholder_is_omitted(tmp_path):
    for stream in ("CONC", "TEMP", "PRES"):
        _write_stream(tmp_path, 13, stream, "Time,Value\n0,1\n1,2\n")
    (tmp_path / "ELE402HSE013FLMT20250101.csv").write_text(
        "Time\n0\n1\n", encoding="utf-8"
    )
    test = read_elvhys_test(tmp_path, 13, required_streams=("CONC", "TEMP", "PRES"))
    assert "FLMT" not in test.streams
