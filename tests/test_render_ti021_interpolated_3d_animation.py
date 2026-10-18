import numpy as np

from tools.render_ti021_interpolated_3d_animation import (
    _apply_plan,
    _frame_times,
    _interpolate_time,
    _inverse_distance_plan,
    reconstruct_fields,
)


def test_time_interpolation_and_early_detail_warp():
    source_time = np.array([0.0, 5.0, 10.0])
    traces = np.array([[0.0, 10.0, 20.0], [4.0, 2.0, 0.0]])
    frame_time = _frame_times(
        source_time, frames=5, end_time_s=10.0, time_warp="early_detail"
    )
    assert np.allclose(frame_time, [0.0, 0.625, 2.5, 5.625, 10.0])
    values = _interpolate_time(source_time, traces, frame_time)
    assert values.shape == (5, 2)
    assert np.allclose(values[[0, -1], 0], [0.0, 20.0])


def test_idw_exact_sensor_and_nonnegative_reconstruction():
    sensors = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    targets = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    plan = _inverse_distance_plan(
        sensors,
        targets,
        neighbours=2,
        power=2.0,
        vertical_scale=10.0,
        support_radius_m=100.0,
    )
    values = _apply_plan(plan, np.array([10.0, 0.0]))
    assert np.isclose(values[0], 10.0)
    assert 0.0 < values[1] < 10.0

    fields = reconstruct_fields(
        np.array([[10.0, 0.0], [0.0, 8.0]]),
        plan,
        shape=(1, 1, 2),
        smooth_sigma_cells=0.0,
    )
    assert fields.shape == (2, 1, 1, 2)
    assert np.all(fields >= 0.0)
