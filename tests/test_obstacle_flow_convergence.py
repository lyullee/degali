import json

import pytest

from degali.addons.fast_obstacle_transport import FastObstacleReceptor
from degali.addons.resolved_obstacle_flow_2d import ResolvedObstacleFlow2DConfig
from degali.addons.semi_fv_obstacle import RectangularObstacle2D, SourceRateSchedule
from degali.validation.obstacle_flow_convergence import (
    RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA,
    resolved_obstacle_grid_convergence_record,
    run_resolved_obstacle_grid_convergence,
)


def _configs():
    schedule = SourceRateSchedule((0.0, 0.5, 3.0), (0.01, 0.0, 0.0))
    obstacle = RectangularObstacle2D(3.0, 0.8, 1.5, label="block")
    configs = []
    for nx, nz in ((40, 20), (60, 30), (80, 40)):
        dx = 8.0 / nx
        flow_steps = int(round(3.0 / (0.03 * dx / 2.0)))
        configs.append(ResolvedObstacleFlow2DConfig(
            length_m=8.0,
            height_m=4.0,
            nx=nx,
            nz=nz,
            wind_speed_m_s=2.0,
            kinematic_viscosity_m2_s=0.3,
            flow_steps=flow_steps,
            lattice_inlet_speed=0.03,
            duration_s=3.0,
            scalar_time_step_s=0.02 * dx,
            scalar_diffusivity_m2_s=0.03,
            source_schedule=schedule,
            source_x_m=0.6,
            source_height_m=0.5,
            source_sigma_m=0.2,
            obstacle=obstacle,
            minimum_flow_steps=flow_steps,
        ))
    return tuple(configs)


def _receptors():
    return (
        FastObstacleReceptor("windward", 2.4, 0.8),
        FastObstacleReceptor("roof", 3.0, 1.7),
        FastObstacleReceptor("lee-low", 4.0, 0.8),
        FastObstacleReceptor("lee-shear", 4.0, 1.7),
    )


@pytest.fixture(scope="module")
def convergence_result():
    return run_resolved_obstacle_grid_convergence(
        _configs(), receptors=_receptors(), relative_tolerance=0.10,
    )


def test_grid_convergence_runs_three_matched_physical_grids(convergence_result):
    result = convergence_result

    assert result.schema == RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA
    assert [case.nx for case in result.cases] == [40, 60, 80]
    assert [case.cell_size_m for case in result.cases] == pytest.approx(
        [0.2, 8.0 / 60.0, 0.1]
    )
    assert len({case.flow_physical_duration_s for case in result.cases}) == 1
    assert all(case.maximum_absolute_divergence_s_1 < 1.0e-9 for case in result.cases)
    assert all(case.maximum_mass_residual_per_span_kg_m < 1.0e-12 for case in result.cases)
    assert all(len(case.receptor_metrics) == 4 for case in result.cases)


def test_grid_convergence_record_is_json_native_and_fail_closed(convergence_result):
    record = resolved_obstacle_grid_convergence_record(convergence_result)

    assert record["schema"] == RESOLVED_OBSTACLE_GRID_CONVERGENCE_SCHEMA
    assert record["diagnostics"]["quantitative_prediction_allowed"] is False
    assert "result" not in record["cases"][0]
    json.dumps(record, allow_nan=False)


def test_grid_convergence_rejects_unmatched_or_underpopulated_studies():
    configs = _configs()
    with pytest.raises(ValueError, match="at least three"):
        run_resolved_obstacle_grid_convergence(
            configs[:2], receptors=_receptors(),
        )
    mismatched = list(configs)
    mismatched[-1] = ResolvedObstacleFlow2DConfig(**{
        **mismatched[-1].__dict__, "wind_speed_m_s": 2.1,
    })
    with pytest.raises(ValueError, match="one physical scenario"):
        run_resolved_obstacle_grid_convergence(
            mismatched, receptors=_receptors(),
        )
