from __future__ import annotations

import json
import hashlib
from pathlib import Path

from tools.audit_exjet_benchmark import audit


def _benchmark(tmp_path):
    source = json.loads(
        (Path(__file__).resolve().parents[1] / "reference" / "obstacle_wake" / "exjet_obstructed_free_jet_benchmark.json").read_text()
    )
    source["observations"]["paper_plot_digitisation_sha256"] = None
    path = tmp_path / "benchmark.json"
    path.write_text(json.dumps(source), encoding="utf-8")
    return path


def test_exjet_audit_fails_closed_without_machine_data(tmp_path):
    record = audit(_benchmark(tmp_path), None, None, None, None)
    assert record["status"] == "withheld_spatial_profile_inputs"
    assert record["quantitative_promotion_allowed"] is False
    assert "machine_registered_cylinder_array_geometry" in record["missing_requirements"]
    assert "common_release_time_definition" in record["missing_requirements"]


def test_obstructed_case_does_not_inherit_unobstructed_12mm_nozzle():
    benchmark = json.loads(
        (Path(__file__).resolve().parents[1] / "reference" / "obstacle_wake" / "exjet_obstructed_free_jet_benchmark.json").read_text()
    )
    assert benchmark["release"]["direct_orifice_diameter_mm"] is None
    assert "do_not_borrow_12_mm" in benchmark["release"]["direct_orifice_diameter_status"]
    pseudo_sources = benchmark["release"]["published_model_specific_obstructed_pseudo_sources"]
    assert len(pseudo_sources) == 4
    assert benchmark["release"]["pseudo_sources_are_model_inputs_not_measured_release_truth"] is True
    assert benchmark["observations"]["temporal_observation_history_available"] is False


def test_published_plot_digitisation_is_counted_but_not_promoted():
    root = Path(__file__).resolve().parents[1]
    benchmark = root / "reference" / "obstacle_wake" / "exjet_obstructed_free_jet_benchmark.json"
    observations = root / "reference" / "obstacle_wake" / "exjet_obstructed_h2_plot_digitisation_v1.csv"
    record = audit(benchmark, observations, None, None, None)
    assert record["observation_count"] == 14
    assert record["observation_quantities"] == [
        "centreline_hydrogen_molar_fraction",
        "radial_hydrogen_molar_fraction",
    ]
    assert record["observation_stations_m"] == [1.4, 2.5, 4.4]
    assert record["status"] == "withheld_spatial_profile_inputs"
    assert record["quantitative_promotion_allowed"] is False
    assert record["observations_have_time_resolved_history"] is False
    assert record["observation_integrity_verified"] is True
    assert "common_release_time_definition" in record["missing_transient_requirements"]


def test_exjet_audit_verifies_pinned_source_pdf_hash():
    root = Path(__file__).resolve().parents[1]
    benchmark = root / "reference" / "obstacle_wake" / "exjet_obstructed_free_jet_benchmark.json"
    source_pdf = root / "reference" / "obstacle_wake" / "exjet_obstructed_free_jet_id252.pdf"
    observations = root / "reference" / "obstacle_wake" / "exjet_obstructed_h2_plot_digitisation_v1.csv"
    expected = json.loads(benchmark.read_text(encoding="utf-8"))["source"]["local_pdf_sha256"]
    actual = hashlib.sha256(source_pdf.read_bytes()).hexdigest()
    assert actual == expected
    record = audit(benchmark, observations, None, None, None, source_pdf)
    assert record["source_pdf_integrity_verified"] is True
    assert record["source_pdf_sha256"] == expected


def test_exjet_audit_accepts_only_complete_machine_package(tmp_path):
    benchmark = _benchmark(tmp_path)
    observations = tmp_path / "observations.csv"
    observations.write_text(
        "station_x_m,radial_coordinate_m,h2_molar_fraction,quantity_kind,coordinate_uncertainty_m,concentration_uncertainty_mole_fraction,observation_time_s\n"
        "1.4,0.0,0.25,centreline_hydrogen_molar_fraction,0.03,0.01,1.0\n"
        "1.4,0.1,0.10,radial_hydrogen_molar_fraction,0.03,0.015,1.0\n"
        "2.5,0.0,0.15,centreline_hydrogen_molar_fraction,0.03,0.01,1.0\n"
        "4.4,0.0,0.05,radial_hydrogen_molar_fraction,0.03,0.015,1.0\n",
        encoding="utf-8",
    )
    geometry = tmp_path / "geometry.json"
    geometry.write_text(
        '{"shape":"cylinder_array","registration_status":"machine_registered_exact"}',
        encoding="utf-8",
    )
    clock = tmp_path / "clock.json"
    clock.write_text(
        '{"release_time_s":0.0,"time_basis_id":"synthetic-test-event-clock"}',
        encoding="utf-8",
    )
    source_term = tmp_path / "source-term.json"
    source_term.write_text('{"provenance":"reported"}', encoding="utf-8")
    record = audit(benchmark, observations, geometry, clock, source_term)
    assert record["status"] == "accepted_for_spatial_and_transient_quantitative_scoring"
    assert record["quantitative_promotion_allowed"] is True
    assert record["observation_count"] == 4


def test_exact_steady_profile_package_can_be_spatial_ready_without_clock(tmp_path):
    benchmark = _benchmark(tmp_path)
    observations = tmp_path / "steady-observations.csv"
    observations.write_text(
        "station_x_m,radial_coordinate_m,h2_molar_fraction,quantity_kind,coordinate_uncertainty_m,concentration_uncertainty_mole_fraction\n"
        "1.4,0.0,0.25,centreline_hydrogen_molar_fraction,0.03,0.01\n"
        "1.4,0.1,0.10,radial_hydrogen_molar_fraction,0.03,0.015\n"
        "2.5,0.0,0.15,centreline_hydrogen_molar_fraction,0.03,0.01\n"
        "4.4,0.0,0.05,radial_hydrogen_molar_fraction,0.03,0.015\n",
        encoding="utf-8",
    )
    geometry = tmp_path / "geometry.json"
    geometry.write_text(
        '{"shape":"cylinder_array","registration_status":"machine_registered_exact"}',
        encoding="utf-8",
    )
    source_term = tmp_path / "source-term.json"
    source_term.write_text('{"provenance":"reported"}', encoding="utf-8")
    record = audit(benchmark, observations, geometry, None, source_term)
    assert record["status"] == "spatial_profile_ready_transient_withheld"
    assert record["spatial_profile_comparison_allowed"] is True
    assert record["transient_validation_allowed"] is False
    assert record["quantitative_promotion_allowed"] is False


def test_schematic_geometry_does_not_satisfy_exact_registration_gate(tmp_path):
    root = Path(__file__).resolve().parents[1]
    benchmark = root / "reference" / "obstacle_wake" / "exjet_obstructed_free_jet_benchmark.json"
    observations = root / "reference" / "obstacle_wake" / "exjet_obstructed_h2_plot_digitisation_v1.csv"
    geometry = root / "reference" / "obstacle_wake" / "exjet_obstructed_schematic_geometry.json"
    record = audit(benchmark, observations, geometry, None, None)
    assert "machine_registered_cylinder_array_geometry" in record["missing_requirements"]
    assert record["spatial_profile_comparison_allowed"] is False
    assert record["quantitative_promotion_allowed"] is False
