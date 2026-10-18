from dataclasses import replace

import pytest

from degali.addons.resolved_obstacle_flow_2d import (
    ResolvedObstacleFlow2DConfig,
    solve_resolved_obstacle_flow_2d,
)
from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.semi_fv_obstacle import (
    RectangularObstacle2D,
    SourceRateSchedule,
)
from degali.validation.aij_case_h import (
    AIJ_CASE_H_DATASET_SCHEMA,
    AijCaseHDataset,
    AijCaseHObservation,
    aij_case_h_flow_comparison_record,
    compare_aij_case_h_centerplane_flow,
    compare_aij_case_h_centerplane_flow_3d,
    read_aij_case_h,
)


def _write_aij_fixture(root):
    approach_header = (
        "x (m), y (m), z (m),U (m/s),V (m/s),W (m/s),u_rms (m/s),"
        "v_rms (m/s),w_rms (m/s),uw (m2/s2),uv (m2/s2),k (m2/s2)\n"
    )
    approach_rows = [
        f"-0.25,0,{0.01 * (index + 1):.2f},{0.2 * (index + 1):.3f},"
        "0,0,0.1,0.1,0.1,-0.01,0,0.015\n"
        for index in range(17)
    ]
    (root / "AF_caseH.csv").write_text(
        approach_header + "".join(approach_rows), encoding="utf-8"
    )

    result_header = (
        "No.,x (m),y (m),z (m),U (m/s),V (m/s),W (m/s),u_rms (m/s),"
        "v_rms (m/s),w_rsm (m/s),k (m2/s2),<c> (ppm),c_rms (ppm)\n"
    )
    result_rows = []
    for index in range(876):
        missing = "-999.9" if index == 0 else "0.1"
        result_rows.append(
            f"{index + 1},{0.001 * index:.3f},0,0.25,{missing},0,0,"
            "0.1,0.1,0.1,0.015,2,0.2\n"
        )
    (root / "RS_caseH.csv").write_text(
        result_header + "".join(result_rows), encoding="utf-8"
    )
    (root / "readme_caseH.md").write_text("fixture\n", encoding="utf-8")


def _flow_result():
    config = ResolvedObstacleFlow2DConfig(
        length_m=4.0,
        height_m=2.0,
        nx=40,
        nz=20,
        wind_speed_m_s=2.0,
        kinematic_viscosity_m2_s=0.3,
        flow_steps=100,
        minimum_flow_steps=100,
        lattice_inlet_speed=0.03,
        duration_s=0.01,
        scalar_time_step_s=0.001,
        scalar_diffusivity_m2_s=0.03,
        source_schedule=SourceRateSchedule((0.0, 0.01), (0.0, 0.0)),
        source_x_m=0.5,
        source_height_m=0.5,
        source_sigma_m=0.15,
        obstacle=RectangularObstacle2D(2.0, 0.5, 1.0, label="block"),
    )
    return solve_resolved_obstacle_flow_2d(config)


def _observation(point_id, x_m, z_m, u_m_s, w_m_s):
    return AijCaseHObservation(
        point_id=point_id,
        x_m=x_m,
        y_m=0.0,
        z_m=z_m,
        u_m_s=u_m_s,
        v_m_s=0.0,
        w_m_s=w_m_s,
        u_rms_m_s=0.1,
        v_rms_m_s=0.1,
        w_rms_m_s=0.1,
        turbulent_kinetic_energy_m2_s2=0.015,
        concentration_ppm=2.0,
        concentration_rms_ppm=0.2,
    )


def test_aij_case_h_loader_hashes_counts_and_maps_missing_values(tmp_path):
    _write_aij_fixture(tmp_path)
    dataset = read_aij_case_h(tmp_path)

    assert dataset.schema == AIJ_CASE_H_DATASET_SCHEMA
    assert len(dataset.approach_points) == 17
    assert len(dataset.observations) == 876
    assert dataset.observations[0].u_m_s is None
    assert len(dataset.approach_flow_sha256) == 64
    assert len(dataset.result_sha256) == 64
    assert len(dataset.readme_sha256) == 64


def test_aij_comparison_is_velocity_only_conditional_component_evidence(tmp_path):
    _write_aij_fixture(tmp_path)
    base = read_aij_case_h(tmp_path)
    dataset = replace(
        base,
        observations=(
            _observation(1, 1.0, 1.25, 1.8, 0.0),
            _observation(2, 3.0, 1.25, 2.1, 0.05),
            _observation(3, 2.0, 0.5, -0.2, 0.0),  # inside the obstacle
        ),
    )

    comparison = compare_aij_case_h_centerplane_flow(
        dataset, _flow_result(), model_x_origin_in_experiment_m=0.0
    )
    record = aij_case_h_flow_comparison_record(comparison)

    assert comparison.status == "component_evidence_conditional"
    assert len(comparison.points) == 2
    assert comparison.skipped_point_ids == (3,)
    assert record["diagnostics"]["concentration_validation_withheld"] is True
    assert record["diagnostics"]["quantitative_hydrogen_validation"] is False
    assert "u_rmse_normalized" in record["metrics"]


def test_aij_case_h_loader_rejects_an_incomplete_result_package(tmp_path):
    _write_aij_fixture(tmp_path)
    lines = (tmp_path / "RS_caseH.csv").read_text(encoding="utf-8").splitlines()
    (tmp_path / "RS_caseH.csv").write_text(
        "\n".join(lines[:-1]) + "\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match="expected 876"):
        read_aij_case_h(tmp_path)


def test_aij_comparison_accepts_the_finite_width_3d_flow_contract(tmp_path):
    _write_aij_fixture(tmp_path)
    base = read_aij_case_h(tmp_path)
    dataset = replace(
        base,
        observations=(
            _observation(1, 1.0, 1.25, 1.8, 0.0),
            _observation(2, 3.0, 1.25, 2.1, 0.05),
        ),
    )
    result = solve_fast_obstacle_transport_3d(FastObstacleTransport3DConfig(
        domain_m=(4.0, 2.4, 2.0),
        cells=(40, 24, 20),
        duration_s=0.001,
        time_step_s=0.001,
        wind_speed_m_s=2.0,
        diffusivity_m2_s=0.03,
        source_schedule=SourceRateSchedule((0.0, 0.001), (0.0, 0.0)),
        source_position_m=(0.5, 1.2, 0.5),
        source_sigma_m=0.15,
    ))

    comparison = compare_aij_case_h_centerplane_flow_3d(
        dataset,
        result,
        model_x_origin_in_experiment_m=0.0,
        model_y_origin_in_experiment_m=-1.2,
    )

    assert comparison.status == "component_evidence_conditional"
    assert len(comparison.points) == 2
    assert comparison.diagnostics["model_dimension"] == "finite-width 3-D"
    assert comparison.diagnostics["concentration_validation_withheld"] is True
