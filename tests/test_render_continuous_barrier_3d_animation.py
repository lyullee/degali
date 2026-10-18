from argparse import Namespace

import numpy as np

from tools.render_continuous_barrier_3d_animation import (
    _barrier_mask,
    _cloud_isosurface_triangles,
    _config,
    _display_field,
    _temporal_upsample,
)


def test_temporal_upsample_adds_smooth_display_only_frames():
    snapshots = [
        np.zeros((1, 1, 2), dtype=np.float32),
        np.full((1, 1, 2), 2.0, dtype=np.float32),
    ]
    fields, times = _temporal_upsample(snapshots, np.array([0.0, 10.0]), 2)
    assert len(fields) == 3
    assert np.allclose(times, [0.0, 5.0, 10.0])
    assert np.allclose(fields[1], 1.0)
    assert fields[0] is not snapshots[1]


def test_display_interpolation_reapplies_impermeable_barrier_mask():
    args = Namespace(
        duration_s=10.0,
        release_duration_s=10.0,
        release_rate_kg_s=0.1,
        wind_speed_m_s=2.2,
    )
    config, barrier = _config(args)
    raw = np.ones((config.cells[2], config.cells[1], config.cells[0])) * 1.0e-3
    shown = _display_field(
        raw,
        config,
        upscale=2.0,
        smooth_sigma_cells=0.65,
        barrier=barrier,
    )
    nz, ny, nx = shown.shape
    x = np.linspace(0.5 * config.domain_m[0] / nx, config.domain_m[0] - 0.5 * config.domain_m[0] / nx, nx)
    y = np.linspace(0.5 * config.domain_m[1] / ny, config.domain_m[1] - 0.5 * config.domain_m[1] / ny, ny)
    z = np.linspace(0.5 * config.domain_m[2] / nz, config.domain_m[2] - 0.5 * config.domain_m[2] / nz, nz)
    iz = np.flatnonzero(z <= barrier["height_m"])
    iy = np.flatnonzero((y >= barrier["y_min_m"]) & (y <= barrier["y_max_m"]))
    ix = np.flatnonzero((x >= barrier["x_m"]) & (x <= barrier["x_m"] + barrier["thickness_m"]))
    assert np.all(shown[np.ix_(iz, iy, ix)] == 0.0)


def test_split_barrier_display_mask_keeps_central_gap_open():
    args = Namespace(
        duration_s=10.0,
        release_duration_s=10.0,
        release_rate_kg_s=0.1,
        wind_speed_m_s=2.2,
        barrier_layout="split-gap",
    )
    config, barrier = _config(args)
    raw = np.ones((config.cells[2], config.cells[1], config.cells[0])) * 1.0e-3
    shown = _display_field(
        raw,
        config,
        upscale=2.0,
        smooth_sigma_cells=0.65,
        barrier=barrier,
    )
    nz, ny, nx = shown.shape
    x = np.linspace(0.5 * config.domain_m[0] / nx, config.domain_m[0] - 0.5 * config.domain_m[0] / nx, nx)
    y = np.linspace(0.5 * config.domain_m[1] / ny, config.domain_m[1] - 0.5 * config.domain_m[1] / ny, ny)
    z = np.linspace(0.5 * config.domain_m[2] / nz, config.domain_m[2] - 0.5 * config.domain_m[2] / nz, nz)
    iz = np.flatnonzero(z <= barrier["height_m"])
    ix = np.flatnonzero((x >= barrier["x_m"]) & (x <= barrier["x_m"] + barrier["thickness_m"]))
    segments = barrier["segments"]
    lower = np.flatnonzero(
        (y >= segments[0]["y_min_m"]) & (y <= segments[0]["y_max_m"])
    )
    gap = np.flatnonzero(
        (y > barrier["gap_y_min_m"]) & (y < barrier["gap_y_max_m"])
    )
    upper = np.flatnonzero(
        (y >= segments[1]["y_min_m"]) & (y <= segments[1]["y_max_m"])
    )
    assert np.all(shown[np.ix_(iz, lower, ix)] == 0.0)
    assert np.all(shown[np.ix_(iz, upper, ix)] == 0.0)
    assert np.all(shown[np.ix_(iz, gap, ix)] > 0.0)


def test_display_interpolation_does_not_create_false_transport_through_wall():
    args = Namespace(
        duration_s=10.0,
        release_duration_s=10.0,
        release_rate_kg_s=0.1,
        wind_speed_m_s=2.2,
        barrier_layout="blocking-wall",
    )
    config, barrier = _config(args)
    nz, ny, nx = config.cells[2], config.cells[1], config.cells[0]
    raw = np.zeros((nz, ny, nx), dtype=float)
    x_raw = (np.arange(nx) + 0.5) * config.domain_m[0] / nx
    raw[:, :, x_raw < float(barrier["x_m"])] = 1.0e-3
    solid = _barrier_mask(config.domain_m, config.cells, **barrier)
    raw[solid] = 0.0
    shown = _display_field(
        raw,
        config,
        upscale=2.0,
        smooth_sigma_cells=0.65,
        barrier=barrier,
    )
    dz, dy, dx = shown.shape
    x = (np.arange(dx) + 0.5) * config.domain_m[0] / dx
    z = (np.arange(dz) + 0.5) * config.domain_m[2] / dz
    below_wall = z <= float(barrier["height_m"])
    downstream = x > float(barrier["x_m"]) + float(barrier["thickness_m"])
    y = (np.arange(dy) + 0.5) * config.domain_m[1] / dy
    within_wall = (
        (y >= float(barrier["y_min_m"]))
        & (y <= float(barrier["y_max_m"]))
    )
    assert np.max(shown[np.ix_(below_wall, within_wall, downstream)]) == 0.0


def test_wall_shadow_display_blur_is_capped_by_linear_raw_field():
    args = Namespace(
        duration_s=10.0,
        release_duration_s=10.0,
        release_rate_kg_s=0.1,
        wind_speed_m_s=2.2,
        barrier_layout="blocking-wall",
    )
    config, barrier = _config(args)
    nz, ny, nx = config.cells[2], config.cells[1], config.cells[0]
    raw = np.zeros((nz, ny, nx), dtype=float)
    x_raw = (np.arange(nx) + 0.5) * config.domain_m[0] / nx
    raw[:, :, x_raw < float(barrier["x_m"])] = 1.0e-3
    first_downstream = np.flatnonzero(
        x_raw >= float(barrier["x_m"]) + float(barrier["thickness_m"])
    )[0]
    raw[:, :, first_downstream] = 1.0e-7
    solid = _barrier_mask(config.domain_m, config.cells, **barrier)
    raw[solid] = 0.0
    linear = _display_field(
        raw, config, upscale=2.0, smooth_sigma_cells=0.0, barrier=barrier
    )
    smoothed = _display_field(
        raw, config, upscale=2.0, smooth_sigma_cells=0.65, barrier=barrier
    )
    dz, dy, dx = smoothed.shape
    x = (np.arange(dx) + 0.5) * config.domain_m[0] / dx
    y = (np.arange(dy) + 0.5) * config.domain_m[1] / dy
    z = (np.arange(dz) + 0.5) * config.domain_m[2] / dz
    below_wall = z <= float(barrier["height_m"])
    first_cell_shadow = (
        (x > float(barrier["x_m"]) + float(barrier["thickness_m"]))
        & (
            x <= float(barrier["x_m"])
            + float(barrier["thickness_m"])
            + config.domain_m[0] / config.cells[0]
        )
    )
    within_wall = (
        (y >= float(barrier["y_min_m"]))
        & (y <= float(barrier["y_max_m"]))
    )
    selection = np.ix_(below_wall, within_wall, first_cell_shadow)
    assert np.all(smoothed[selection] <= linear[selection])


def test_building_cluster_display_mask_keeps_every_building_empty():
    args = Namespace(
        duration_s=10.0,
        release_duration_s=10.0,
        release_rate_kg_s=0.1,
        wind_speed_m_s=2.2,
        barrier_layout="building-cluster",
    )
    config, geometry = _config(args)
    raw = np.ones((config.cells[2], config.cells[1], config.cells[0])) * 1.0e-3
    shown = _display_field(
        raw,
        config,
        upscale=2.0,
        smooth_sigma_cells=0.65,
        barrier=geometry,
    )
    solid = _barrier_mask(
        config.domain_m,
        (shown.shape[2], shown.shape[1], shown.shape[0]),
        **geometry,
    )
    assert np.count_nonzero(solid) > 0
    assert np.all(shown[solid] == 0.0)


def test_cloud_isosurface_is_closed_and_uses_physical_coordinates():
    field = np.zeros((6, 8, 10), dtype=float)
    field[2:5, 2:6, 3:8] = 1.0
    triangles = _cloud_isosurface_triangles(
        field,
        (50.0, 40.0, 30.0),
        level_vol_percent=0.5,
        vertical_exaggeration=1.0,
    )
    assert triangles.ndim == 3
    assert triangles.shape[1:] == (3, 3)
    assert triangles.shape[0] > 0
    assert np.min(triangles[..., 0]) >= 0.0
    assert np.max(triangles[..., 0]) <= 50.0
    assert np.min(triangles[..., 1]) >= 0.0
    assert np.max(triangles[..., 1]) <= 40.0
    assert np.min(triangles[..., 2]) >= 0.0
    assert np.max(triangles[..., 2]) <= 30.0
