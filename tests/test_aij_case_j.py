from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from degali.validation.aij_case_j import (
    AIJ_CASE_J_DATASET_SCHEMA,
    read_aij_case_j,
)
from tools.audit_aij_case_j_dataset import build
from degali.validation.aij_case_j_comparison import (
    ETHYLENE_MOLAR_MASS_KG_MOL,
    compare_aij_case_j,
)
from degali.validation.aij_case_j_setup import build_aij_case_j_resolved_setup


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reference" / "obstacle_wake" / "aij_case_j_zenodo_15447888"


def test_official_case_j_package_is_hash_pinned_and_units_are_explicit():
    dataset = read_aij_case_j(DATA)

    assert dataset.schema == AIJ_CASE_J_DATASET_SCHEMA
    assert dataset.record_doi == "10.5281/zenodo.15447888"
    assert (len(dataset.approach_points), len(dataset.surface_points), len(dataset.observations)) == (20, 121, 444)
    assert dataset.approach_sha256 == "1ecfda6e3e99f889e92b5c2899bdbf5f2a13a46e37b77f7f0da436dc0f10c958"
    assert dataset.surface_temperature_sha256 == "8874909df799fac347b3e1c7335cc7ebc025a8c596a037071e2661270d8bdcba"
    assert dataset.observations_sha256 == "39ffbbbbca133b370bcde3bbd52a258fddebd599b2c531d69efd6699f3d4f729"
    peak = max(dataset.observations, key=lambda item: item.concentration_ppm)
    assert peak.concentration_ppm == pytest.approx(5370.794)
    assert peak.concentration_vol_percent == pytest.approx(0.5370794)
    assert dataset.tracer == "ethylene"
    assert any("not hydrogen" in warning for warning in dataset.warnings)
    assert any("not a transient" in warning for warning in dataset.warnings)
    velocity_z, velocity_u = dataset.inlet_velocity_profile
    temperature_z, temperature_k = dataset.inlet_temperature_profile
    assert len(velocity_z) == len(velocity_u) == 21
    assert velocity_z[0] == velocity_u[0] == 0.0
    assert len(temperature_z) == len(temperature_k) == 20
    assert temperature_z == tuple(sorted(temperature_z))
    assert temperature_k[0] == pytest.approx(294.5)
    assert dataset.mean_building_surface_temperature_k == pytest.approx(297.6223140495868)
    assert dataset.building_surface_temperature_range_k == pytest.approx((295.15, 302.05))


def test_case_j_loader_rejects_non_officially_checksummed_files(tmp_path):
    for name in ("readme_caseJ.md", "AF_caseJ.csv", "BD-ST_caseJ.csv", "RS_caseJ.csv"):
        (tmp_path / name).write_bytes((DATA / name).read_bytes())
    with (tmp_path / "AF_caseJ.csv").open("ab") as handle:
        handle.write(b"\n")

    with pytest.raises(ValueError, match="MD5 mismatch"):
        read_aij_case_j(tmp_path)


def test_case_j_dataset_audit_does_not_claim_model_or_hydrogen_validation():
    record = build()

    assert record["status"] == "scored_development_component_comparison_validation_withheld"
    assert record["row_counts"] == {
        "approach_flow": 20,
        "building_surface_temperature": 121,
        "wake_observations": 444,
    }
    assert record["wake_observation_ranges"]["ethylene_vol_percent"][1] == pytest.approx(0.5370794)
    assert record["model_comparison"]["status"] == "scored_development_component_comparison_validation_withheld"
    assert record["model_comparison"]["point_count"] == 440
    assert record["model_comparison"]["experimental_statistic_match"] is False
    assert record["model_comparison"]["ethylene_vol_percent_mae_vol_percentage_points"] == pytest.approx(0.004859552711164247)
    station = {item["x_experiment_m"]: item for item in record["model_comparison"]["concentration_by_downstream_station"]}
    assert station[0.04]["measured_max_vol_percent"] == pytest.approx(0.5370794)
    assert station[0.04]["predicted_max_vol_percent"] == pytest.approx(0.1735, abs=1e-3)
    assert station[0.02]["predicted_max_vol_percent"] > station[0.04]["predicted_max_vol_percent"]
    assert record["model_comparison"]["pointwise_scoring_implementation_available"] is True
    assert record["model_comparison"]["hydrogen_validation"] is False


def test_case_j_comparator_interpolates_velocity_temperature_and_tracer_as_vol_percent():
    dataset = read_aij_case_j(DATA)
    original = dataset.observations[0]
    x_axis = np.asarray((original.x_m - 0.01, original.x_m + 0.01))
    y_axis = np.asarray((-0.01, 0.01))
    z_axis = np.asarray((original.z_m - 0.01, original.z_m + 0.01))
    zz, yy, xx = np.meshgrid(z_axis, y_axis, x_axis, indexing="ij")
    u = 1.0 + xx + 2.0 * yy + 3.0 * zz
    v = -0.5 * xx + yy - zz
    w = 0.25 * xx - 0.5 * yy + 2.0 * zz
    temperature = 280.0 + 2.0 * xx + yy + 3.0 * zz
    predicted = (
        float(1.0 + original.x_m + 3.0 * original.z_m),
        float(-0.5 * original.x_m - original.z_m),
        float(0.25 * original.x_m + 2.0 * original.z_m),
        float(280.0 + 2.0 * original.x_m + 3.0 * original.z_m),
    )
    total_molar_density = 101325.0 / (8.31446261815324 * predicted[3])
    ethylene_mass = (
        original.concentration_vol_percent / 100.0
        * total_molar_density * ETHYLENE_MOLAR_MASS_KG_MOL
    )
    synthetic_observation = replace(
        original,
        u_m_s=predicted[0], v_m_s=predicted[1], w_m_s=predicted[2],
        temperature_k=predicted[3],
    )
    one_point_dataset = replace(dataset, observations=(synthetic_observation,))
    result = SimpleNamespace(
        z_m=z_axis, y_m=y_axis, x_m=x_axis,
        velocity_u_m_s=u, velocity_v_m_s=v, velocity_w_m_s=w,
        temperature_k=temperature, concentration_kg_m3=np.full_like(u, ethylene_mass),
        solid_mask=np.zeros_like(u, dtype=bool),
        ambient_temperature_k=285.32, ambient_pressure_pa=101325.0,
    )

    comparison = compare_aij_case_j(
        one_point_dataset, result,
        coordinate_translation_model_minus_experiment_m=(0.0, 0.0, 0.0),
    )

    assert comparison.status == "scored_component_comparison_not_validation"
    assert comparison.metrics["point_count"] == 1
    assert comparison.metrics["u_m_s_rmse"] == pytest.approx(0.0, abs=1.0e-12)
    assert comparison.metrics["temperature_k_rmse"] == pytest.approx(0.0, abs=1.0e-12)
    assert comparison.metrics["ethylene_vol_percent_rmse"] == pytest.approx(0.0, abs=1.0e-12)


def test_case_j_resolved_setup_uses_physical_fv_properties_and_source_mapping():
    dataset = read_aij_case_j(DATA)
    setup = build_aij_case_j_resolved_setup(
        dataset, grid_spacing_m=0.02, duration_s=0.02, time_step_s=0.001,
    )
    config = setup.config

    assert config.active_momentum_solver == "finite_volume_projection"
    assert config.carrier_flow_initialization == "projected_inlet"
    assert config.kinematic_viscosity_m2_s == pytest.approx(1.46e-5, rel=0.08)
    assert setup.reynolds_from_properties == pytest.approx(16_000, rel=0.12)
    assert setup.reference_reynolds_number == 16_000
    assert config.cells == (44, 20, 20)
    assert config.obstacles[0].x_min_m == pytest.approx(0.32)
    assert config.obstacles[0].x_max_m == pytest.approx(0.40)
    assert config.source_position_m == pytest.approx((0.44, 0.20, 0.0025))
    assert setup.source_hole_resolved_cells_across == pytest.approx(0.25)
    assert config.source_total_mass_to_species_mass_ratio > 20.0
    assert config.source_schedule.released_mass_kg == pytest.approx(
        setup.source_ethylene_mass_rate_kg_s * config.duration_s
    )
    assert config.transient_average_start_s == pytest.approx(0.01)
    assert config.inlet_temperature_profile_z_m is not None
    assert config.obstacle_wall_temperature_k == pytest.approx(
        dataset.mean_building_surface_temperature_k
    )
    assert config.ground_temperature_k == pytest.approx(dataset.floor_temperature_k)
