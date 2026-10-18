import csv
import hashlib
import json

from tools import audit_elvhys_tcs_component as audit_module


def _csv(path, header, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def test_elvhys_component_audit_pins_sensor_traces_and_withholds_wake_promotion(
    tmp_path, monkeypatch,
):
    declarations = "\n".join((
        "Hydrogen sensors are calibrated at 5 points between 0 - 100 vol%H2",
        "PT1ReleaseRig: Pressure in LH2 release rig (barg)",
        "PT2Nozzle: Pressure at LH2 nozzle (barg)",
        "FanFlowMeter: Flow rate of ventilation fan (L/min)",
        "Hydrogen Sensors | Xensor | +- 1 vol% H2",
    ))
    declarations += "\nSensor\tX (mm)\tY (mm)\tZ (mm)\tSensors\n"
    declarations += "Bottom1_Front\t100\t500\t200\tT,C\n"
    declarations += "Bottom2\t200\t500\t200\tT\n"
    declarations += "Side1_Top\t500\t500\t900\tT,C\n\n"
    (tmp_path / "00_ELVHYSWP4_2_README.txt").write_text(declarations)
    (tmp_path / "ELE402HSEMETA.csv").write_text("metadata\n")
    (tmp_path / "ELVHYS_WP4.2_TCS_Sensor_Details.pdf").write_bytes(b"fixture pdf")
    for test_id, interval in ((29, 0.05), (31, 0.002)):
        prefix = f"test{test_id:03d}"
        _csv(
            tmp_path / f"{prefix}_CONC.csv",
            ["Time", "H2Bottom1", "H2Bottom2", "H2Side1"],
            [[0, -0.5, 1.0, 4.0], [interval, 2.5, 2.0, 6.0]],
        )
        _csv(
            tmp_path / f"{prefix}_FLMT.csv",
            ["Time", "FanFlowMeter"],
            [[0, 10], [interval, 10]],
        )
        _csv(
            tmp_path / f"{prefix}_MISC.csv",
            ["Time", "AmbientTemperature", "Humidity", "Pressure"],
            [[0, 20, 50, 1000], [interval, 20, 50, 1000]],
        )
        _csv(
            tmp_path / f"{prefix}_PRES.csv",
            ["Time", "PT1ReleaseRig", "PT2Nozzle"],
            [[0, 2.0, 1.9], [interval, 1.8, 1.7]],
        )

    monkeypatch.setattr(
        audit_module,
        "EXPECTED_SHA256",
        {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in tmp_path.iterdir()
        },
    )

    record = audit_module.audit(tmp_path)

    assert record["status"] == "conditional_cryogenic_h2_enclosure_vent_restriction_component"
    assert record["promotion_allowed"] is False
    assert record["tests"]["29"]["concentration"]["h2_sensor_count"] == 3
    assert record["tests"]["31"]["concentration"]["sample_interval_s"]["minimum"] == 0.002
    assert record["tests"]["31"]["concentration"]["h2_sensor_range_vol_percent"]["H2Bottom1"] == [-0.5, 2.5]
    assert record["tests"]["31"]["release_mass_flow_measured"] is False
    assert record["tests"]["31"]["common_clock_verified"] is True
    assert record["tests"]["31"]["sensor_registration"]["observed_channel_without_concentration_registry"] == ["H2Bottom2"]
    assert record["tests"]["31"]["sensor_registration"]["all_observed_channels_have_coordinates"] is True
    assert record["tests"]["31"]["field_evidence_manifest"]["evidence_readiness"]["status"] == "conditional"
    assert record["geometry_registry"]["quantitative_geometry_allowed"] is False
    geometry = json.loads(
        (tmp_path / "ELVHYS_WP4.2_geometry_registry_v1.json").read_text()
    )
    assert geometry["pressure_peaking_test_vent_configurations"]["31"][
        "reported_total_open_area_m2"
    ] == 0.0104
    assert any("not reconciled" in item for item in geometry["unresolved_items"])
    obstacle_artifact = next(
        item for item in record["tests"]["31"]["field_evidence_manifest"]["channel_artifacts"]
        if item["channel"] == "obstacle_geometry"
    )
    assert obstacle_artifact["sha256"] == record["geometry_registry"]["sha256"]
    assert any("not an open-field obstacle wake" in limitation for limitation in record["limitations"])
