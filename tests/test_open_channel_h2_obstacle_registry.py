import json
import hashlib
import math
from pathlib import Path

import pytest

from degali.validation.open_channel_h2_obstacle import (
    conditional_ffi_o1_crossbeam_assembly,
    conditional_ffi_o1_thin_top_crossbeam_assembly,
    load_registered_ffi_obstacle_geometry,
    load_open_channel_h2_obstacle_benchmark,
)


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = (
    ROOT / "reference" / "hydrogen" / "open_channel_usn"
    / "obstacle_benchmark_registry_v1.json"
)


def test_ffi_obstacle_registry_is_complete_but_fails_closed_on_geometry():
    record = json.loads(REGISTRY.read_text(encoding="utf-8"))

    sensors = record["sensors"]
    assert record["schema"] == "degali.ffi-open-channel-obstacle-benchmark.v1"
    assert len(sensors) == record["time_operator"]["sensor_count"] == 29
    assert len({sensor["id"] for sensor in sensors}) == 29
    assert {sensor["id"] for sensor in sensors} == {
        f"sensor {number}" for number in range(1, 31) if number != 16
    }
    for sensor in sensors:
        assert 0.0 <= sensor["x_m"] <= 0.9
        assert 0.0 <= sensor["y_m"] <= 5.8
        assert 0.0 <= sensor["z_m"] <= 0.8

    obstacle = record["obstacle"]
    assert obstacle["reported_blockage_ratio"] == 0.368
    configuration = obstacle["article_configuration_evidence"]
    assert "two top cross beams" in configuration["first_configuration"]
    assert configuration["reported_configuration_blockage_ratios"] == [0.37, 0.45, 0.53]
    assert "not established as a solid rectangular prism" in configuration["interpretation_limit"]
    projection = obstacle["conditional_rectangular_projection_check"]
    assert math.isclose(projection["projected_obstacle_area_m2"], 0.8 * 0.75)
    assert math.isclose(projection["channel_cross_section_area_m2"], 0.9 * 0.8)
    assert math.isclose(
        projection["implied_blockage_ratio"],
        projection["projected_obstacle_area_m2"]
        / projection["channel_cross_section_area_m2"],
    )
    assert obstacle["quantitative_geometry_allowed"] is False
    assert record["promotion"]["obstacle_h2_validation_allowed"] is False
    assert "obstacle_assembly_dimensions_unresolved" in record["promotion"]["gate_codes"]
    assert "source_flow_history_ends_before_sensor_record" in record["promotion"]["gate_codes"]
    assert "source-flow record ends before the sensor record" in record["time_operator"]["observed_coverage_issue"]
    sensor_uncertainty = record["sensor_instrumentation"]
    assert sensor_uncertainty["model_reported_by_article"] == "Xensor XEN-5320"
    assert sensor_uncertainty["article_reported_accuracy_text"] == "±1%"
    assert "1% full scale" in sensor_uncertainty["manufacturer_datasheet_inaccuracy"]
    assert sensor_uncertainty["event_specific_calibration_certificate_available"] is False
    assert "do not assert event-level accuracy" in sensor_uncertainty["scoring_policy"]
    relative_uncertainty = sensor_uncertainty[
        "derived_relative_uncertainty_if_that_interpretation_applies"
    ]
    assert relative_uncertainty == {
        "at_2_vol_percent": "50% of reading",
        "at_4_vol_percent": "25% of reading",
    }


def test_ffi_registry_pairs_each_obstacle_condition_to_the_article_control():
    record = json.loads(REGISTRY.read_text(encoding="utf-8"))
    cases = {case["article_case"]: case for case in record["case_crosswalk"]}

    assert cases["A"]["dataset_files"] == [
        "23_FFI_P101_T00012.csv", "23_FFI_P101_T00013.csv"
    ]
    assert cases["O1"]["dataset_files"] == [
        "23_FFI_P101_T00027.csv", "23_FFI_P101_T00028.csv"
    ]
    assert cases["E"]["dataset_files"] == [
        "23_FFI_P101_T00003.csv", "23_FFI_P101_T00004.csv"
    ]
    assert cases["O2"]["dataset_files"] == [
        "23_FFI_P101_T00029.csv", "23_FFI_P101_T00030.csv"
    ]
    assert cases["A"]["obstacle"] is False and cases["O1"]["obstacle"] is True
    assert cases["E"]["obstacle"] is False and cases["O2"]["obstacle"] is True


def _materialize_synthetic_ffi_bundle(tmp_path):
    record = json.loads(REGISTRY.read_text(encoding="utf-8"))
    data_root = tmp_path / "extracted"
    data_root.mkdir()
    sensor_ids = [sensor["id"] for sensor in record["sensors"]]
    headers = [
        "flow time [s]", "mass flow meter 1 [g/s]", "h2 sensor time",
        *(f"{sensor} h2 concentration [%]" for sensor in sensor_ids),
    ]
    rows = [
        ["0", "0.1", "0", *("0.1" for _ in sensor_ids)],
        ["1", "0.2", "1", *("0.2" for _ in sensor_ids)],
    ]
    payload = ",".join(headers) + "\n" + "\n".join(
        ",".join(row) for row in rows
    ) + "\n"
    for case in record["case_crosswalk"]:
        hashes = []
        for filename in case["dataset_files"]:
            path = data_root / filename
            path.write_text(payload, encoding="utf-8")
            hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
        case["dataset_csv_sha256"] = hashes
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(record), encoding="utf-8")
    return registry_path, data_root


def test_ffi_obstacle_bundle_loader_verifies_hashes_and_keeps_score_withheld(tmp_path):
    registry_path, data_root = _materialize_synthetic_ffi_bundle(tmp_path)

    benchmark = load_open_channel_h2_obstacle_benchmark(registry_path, data_root)

    assert len(benchmark.sensor_receptors) == 29
    assert [case.article_case for case in benchmark.cases] == ["A", "O1", "E", "O2"]
    assert all(len(case.runs) == 2 for case in benchmark.cases)
    assert benchmark.geometry_quantitative_allowed is False
    assert benchmark.obstacle_h2_validation_allowed is False
    assert benchmark.degali_score_allowed is False
    assert benchmark.cases[0].runs[0].sensor_percent["sensor 1"].tolist() == pytest.approx(
        [0.1, 0.2]
    )
    run = benchmark.cases[0].runs[0]
    assert run.flow_time_s is not run.sensor_time_s
    assert run.flow_time_s.tolist() == run.sensor_time_s.tolist()


def test_ffi_obstacle_bundle_loader_rejects_changed_raw_csv_hash(tmp_path):
    registry_path, data_root = _materialize_synthetic_ffi_bundle(tmp_path)
    changed = data_root / "23_FFI_P101_T00027.csv"
    changed.write_text(changed.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        load_open_channel_h2_obstacle_benchmark(registry_path, data_root)


def test_ffi_obstacle_score_gate_uses_manifest_field_and_requires_zero_gate_codes(tmp_path):
    registry_path, data_root = _materialize_synthetic_ffi_bundle(tmp_path)
    record = json.loads(registry_path.read_text(encoding="utf-8"))
    record["obstacle"]["quantitative_geometry_allowed"] = True
    record["promotion"]["obstacle_h2_validation_allowed"] = True
    record["promotion"]["degali_score_available"] = True
    registry_path.write_text(json.dumps(record), encoding="utf-8")

    gated = load_open_channel_h2_obstacle_benchmark(registry_path, data_root)
    assert gated.obstacle_h2_validation_allowed is False
    assert gated.degali_score_allowed is False

    record["promotion"]["gate_codes"] = []
    registry_path.write_text(json.dumps(record), encoding="utf-8")
    benchmark = load_open_channel_h2_obstacle_benchmark(registry_path, data_root)

    assert benchmark.obstacle_h2_validation_allowed is True
    assert benchmark.degali_score_allowed is True
    assert benchmark.gate_codes == ()


def test_conditional_ffi_crossbeam_surrogate_preserves_reported_blockage():
    beams = conditional_ffi_o1_crossbeam_assembly()

    assert len(beams) == 2
    blocked_area = sum(
        (beam.x_max_m - beam.x_min_m) * (beam.z_max_m - beam.z_min_m)
        for beam in beams
    )
    assert blocked_area / (0.9 * 0.8) == pytest.approx(0.368, rel=1.0e-12)
    assert all(0.0 < beam.x_min_m < beam.x_max_m < 0.9 for beam in beams)
    assert all(2.1 <= beam.y_min_m < beam.y_max_m <= 3.1 for beam in beams)
    assert all(beam.z_max_m < 0.8 for beam in beams)


def test_conditional_ffi_thin_top_crossbeam_hypothesis_preserves_blockage():
    beams = conditional_ffi_o1_thin_top_crossbeam_assembly()

    blocked_area = sum(
        (beam.x_max_m - beam.x_min_m) * (beam.z_max_m - beam.z_min_m)
        for beam in beams
    )
    assert blocked_area / (0.9 * 0.8) == pytest.approx(0.368, rel=1.0e-12)
    assert all(beam.z_max_m == pytest.approx(0.79) for beam in beams)
    assert all(beam.z_max_m - beam.z_min_m == pytest.approx(0.15) for beam in beams)


def test_registered_ffi_geometry_requires_exact_registration_and_preserves_bounds(tmp_path):
    path = tmp_path / "registered-o1.json"
    path.write_text(json.dumps({
        "schema": "degali.ffi-obstacle-geometry.v1",
        "registration_status": "machine_registered_exact",
        "article_case": "O1",
        "coordinate_frame": "ffi-open-channel-v1",
        "obstacles": [{
            "shape": "axis_aligned_cuboid",
            "label": "beam-A",
            "x_min_m": 0.1, "x_max_m": 0.8,
            "y_min_m": 2.1, "y_max_m": 2.5,
            "z_min_m": 0.55, "z_max_m": 0.75,
        }],
    }), encoding="utf-8")

    beams = load_registered_ffi_obstacle_geometry(path)
    assert len(beams) == 1
    assert beams[0].label == "beam-A"
    assert beams[0].bounding_box_m == pytest.approx(
        (0.1, 0.8, 2.1, 2.5, 0.55, 0.75)
    )

    record = json.loads(path.read_text(encoding="utf-8"))
    record["registration_status"] = "conditional_surrogate"
    path.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="machine_registered_exact"):
        load_registered_ffi_obstacle_geometry(path)
