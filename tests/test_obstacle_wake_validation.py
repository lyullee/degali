import math

import pytest

from degali.addons.obstacle_wake_validation import (
    DenseGasFenceObservation,
    NeutralObstacleWakePoint,
    pair_smedis_fence_trials,
    read_aij_case_h,
)
from degali.validation.smedis import Sensor, SmedisTrial


def _write_case_h(path, rows):
    header = (
        "No.,x (m),y (m),z (m),U (m/s),V (m/s),W (m/s),"
        "u_rms (m/s),v_rms (m/s),w_rsm (m/s),k (m2/s2),<c> (ppm),c_rms (ppm)"
    )
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")


def test_local_aij_case_h_intake_preserves_measurements_and_scope(tmp_path):
    path = tmp_path / "RS_caseH.csv"
    _write_case_h(path, [
        "1,0,0,0.1,1,0,0,0.1,0.1,0.1,0.015,20,3",
        "2,0.1,0,0.1,-999.9,0,0,0.2,0.2,0.2,0.06,-999.9,4",
    ])
    benchmark = read_aij_case_h(path)

    assert benchmark.source_path == path
    assert benchmark.building_height_m == pytest.approx(0.20)
    assert benchmark.points[0].concentration_measured
    assert benchmark.points[0].velocity_measured
    assert not benchmark.points[1].concentration_measured
    assert not benchmark.points[1].velocity_measured
    assert math.isnan(benchmark.points[1].concentration_ppm)
    assert benchmark.validation_scope == "neutral_gas_single_cuboid_only"
    assert not benchmark.quantitative_wake_prediction_allowed


def test_local_aij_case_h_intake_preserves_signed_background_subtracted_mean(tmp_path):
    path = tmp_path / "RS_caseH.csv"
    _write_case_h(path, [
        "1,0,0,0.1,1,0,0,0.1,0.1,0.1,0.015,-0.211,3",
    ])

    benchmark = read_aij_case_h(path)

    assert benchmark.points[0].concentration_ppm == pytest.approx(-0.211)
    assert benchmark.points[0].concentration_measured


def test_local_aij_case_h_intake_rejects_schema_or_duplicate_identifier(tmp_path):
    incomplete = tmp_path / "incomplete.csv"
    incomplete.write_text("No.,x (m)\n1,0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing required columns"):
        read_aij_case_h(incomplete)

    duplicate = tmp_path / "duplicate.csv"
    row = "1,0,0,0.1,1,0,0,0.1,0.1,0.1,0.015,20,3"
    _write_case_h(duplicate, [row, row])
    with pytest.raises(ValueError, match="identifiers must be unique"):
        read_aij_case_h(duplicate)


def test_obstacle_observation_boundaries_reject_corrupt_physical_values():
    with pytest.raises(ValueError, match="kinetic energy cannot be negative"):
        NeutralObstacleWakePoint(
            identifier=1,
            x_m=0.0,
            y_m=0.0,
            z_m=0.1,
            velocity_m_s=(1.0, 0.0, 0.0),
            velocity_rms_m_s=(0.1, 0.1, 0.1),
            turbulent_kinetic_energy_m2_s2=-1.0e-6,
            concentration_ppm=20.0,
            concentration_rms_ppm=3.0,
        )

    with pytest.raises(ValueError, match="concentration_rms_ppm cannot be negative"):
        NeutralObstacleWakePoint(
            identifier=1,
            x_m=0.0,
            y_m=0.0,
            z_m=0.1,
            velocity_m_s=(1.0, 0.0, 0.0),
            velocity_rms_m_s=(0.1, 0.1, 0.1),
            turbulent_kinetic_energy_m2_s2=0.015,
            concentration_ppm=1.0,
            concentration_rms_ppm=-3.0,
        )

    with pytest.raises(ValueError, match="fence_std_percent cannot be negative"):
        DenseGasFenceObservation(
            position_m=(1.0, 0.0, 0.05),
            control_percent=1.0,
            fence_percent=0.5,
            control_std_percent=0.1,
            fence_std_percent=-0.01,
        )


def test_obstacle_observation_boundaries_keep_missing_nan_values_as_missing():
    point = NeutralObstacleWakePoint(
        identifier=1,
        x_m=0.0,
        y_m=0.0,
        z_m=0.1,
        velocity_m_s=(math.nan, 0.0, 0.0),
        velocity_rms_m_s=(0.1, math.nan, 0.1),
        turbulent_kinetic_energy_m2_s2=math.nan,
        concentration_ppm=math.nan,
        concentration_rms_ppm=math.nan,
    )
    assert not point.velocity_measured
    assert not point.concentration_measured

    observation = DenseGasFenceObservation(
        position_m=(1.0, 0.0, 0.05),
        control_percent=1.0,
        fence_percent=0.5,
        control_std_percent=math.nan,
        fence_std_percent=math.nan,
    )
    assert observation.fence_to_control_ratio == pytest.approx(0.5)


def _fence_trial(path, *, fences, rate=0.11, sensors=None):
    return SmedisTrial(
        path=path,
        substance="Propane",
        conditions={
            "number of fences": fences,
            "release rate": rate,
            "release duration": 50.0,
            "release point x": 0.0,
            "release point y": 0.0,
            "release point z": 0.05,
            "site average windspeed at zref": 3.4,
            "ideal wind direction": 0.0,
            "surface roughness": 0.1,
        },
        sensors=sensors or [
            Sensor(1.0, 0.0, 0.05, 1.0, 0.1),
            Sensor(1.0, 0.0, 0.05, 2.0, 0.2),
            Sensor(2.0, 0.0, 0.05, 0.0, 0.1),
        ],
    )


def test_paired_smedis_fence_intake_preserves_nonuniform_observations(tmp_path):
    control = _fence_trial(tmp_path / "control.xls", fences=0.0)
    fence = _fence_trial(tmp_path / "fence.xls", fences=1.0, sensors=[
        Sensor(1.0, 0.0, 0.05, 0.5, 0.05),
        Sensor(1.0, 0.0, 0.05, 4.0, 0.4),
        Sensor(2.0, 0.0, 0.05, 0.2, 0.02),
    ])
    benchmark = pair_smedis_fence_trials(control, fence)

    assert len(benchmark.observations) == 3
    assert benchmark.observations[0].fence_to_control_ratio == pytest.approx(0.5)
    assert benchmark.observations[1].fence_to_control_ratio == pytest.approx(2.0)
    assert math.isnan(benchmark.observations[2].fence_to_control_ratio)
    assert benchmark.validation_scope == "paired_dense_gas_fence_observation_only"
    assert not benchmark.quantitative_lh2_prediction_allowed


def test_paired_smedis_fence_intake_rejects_mismatched_conditions_or_layout(tmp_path):
    control = _fence_trial(tmp_path / "control.xls", fences=0.0)
    changed_rate = _fence_trial(tmp_path / "fence.xls", fences=1.0, rate=0.12)
    with pytest.raises(ValueError, match="release rate"):
        pair_smedis_fence_trials(control, changed_rate)

    missing_sensor = _fence_trial(
        tmp_path / "fence.xls", fences=1.0,
        sensors=[Sensor(1.0, 0.0, 0.05, 1.0, 0.1)],
    )
    with pytest.raises(ValueError, match="sensor positions"):
        pair_smedis_fence_trials(control, missing_sensor)
