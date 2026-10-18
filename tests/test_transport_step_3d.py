import numpy as np
import pytest

from degali.addons.fast_obstacle_transport_3d import (
    _mass_consistent_face_projection_3d,
    _transport_step,
)


def _uniform_x_faces(shape, speed):
    nz, ny, nx = shape
    return (
        np.full((nz, ny, nx + 1), speed),
        np.zeros((nz, ny + 1, nx)),
        np.zeros((nz + 1, ny, nx)),
    )


def _advance(initial, *, scheme, speed, dx, dt, steps):
    field = initial.copy()
    solid = np.zeros(field.shape, dtype=bool)
    faces = _uniform_x_faces(field.shape, speed)
    outflow = 0.0
    for _ in range(steps):
        field, step_outflow = _transport_step(
            field,
            u_face=faces[0], v_face=faces[1], w_face=faces[2],
            solid=solid, diffusivity=0.0,
            dt=dt, spacing=np.asarray((dx, dx, dx)),
            advection_scheme=scheme,
        )
        outflow += step_outflow
    return field, outflow


def test_muscl_van_leer_reduces_smooth_translation_error_and_conserves_mass():
    shape = (4, 4, 160)
    dx = 0.05
    speed = 0.5
    dt = 0.4 * dx / speed
    steps = 25
    x = (np.arange(shape[2]) + 0.5) * dx
    initial_line = np.exp(-0.5 * ((x - 2.0) / 0.20) ** 2)
    initial = np.broadcast_to(initial_line, shape).copy()
    exact_line = np.exp(-0.5 * (
        (x - (2.0 + speed * dt * steps)) / 0.20
    ) ** 2)
    exact = np.broadcast_to(exact_line, shape)

    upwind, upwind_outflow = _advance(
        initial, scheme="upwind", speed=speed, dx=dx, dt=dt, steps=steps,
    )
    muscl, muscl_outflow = _advance(
        initial, scheme="muscl_van_leer", speed=speed, dx=dx, dt=dt,
        steps=steps,
    )
    upwind_error = float(np.sqrt(np.mean((upwind - exact) ** 2)))
    muscl_error = float(np.sqrt(np.mean((muscl - exact) ** 2)))
    cell_volume = dx**3
    initial_mass = float(np.sum(initial) * cell_volume)

    assert muscl_error < 0.60 * upwind_error
    assert np.min(muscl) >= -1.0e-14
    assert np.max(muscl) <= 1.0 + 1.0e-12
    assert float(np.sum(upwind) * cell_volume + upwind_outflow) == pytest.approx(
        initial_mass, abs=1.0e-12,
    )
    assert float(np.sum(muscl) * cell_volume + muscl_outflow) == pytest.approx(
        initial_mass, abs=1.0e-12,
    )


def test_muscl_van_leer_keeps_translated_step_bounded():
    shape = (3, 3, 100)
    dx = 0.05
    initial = np.zeros(shape)
    initial[:, :, 20:35] = 1.0
    result, _ = _advance(
        initial, scheme="muscl_van_leer", speed=0.5, dx=dx,
        dt=0.02, steps=30,
    )

    assert np.min(result) >= -1.0e-14
    assert np.max(result) <= 1.0 + 1.0e-12


def test_muscl_ssprk2_improves_multidimensional_translation_and_conserves_mass():
    shape = (32, 32, 32)
    dx = 0.05
    velocity = (0.20, 0.15, 0.10)
    dt = 0.4 * dx / sum(velocity)
    steps = 10
    x = (np.arange(shape[2]) + 0.5) * dx
    y = (np.arange(shape[1]) + 0.5) * dx
    z = (np.arange(shape[0]) + 0.5) * dx
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    center = (0.45, 0.55, 0.65)
    sigma = 0.15
    initial = np.exp(-0.5 * (
        (xx - center[0])**2 + (yy - center[1])**2 + (zz - center[2])**2
    ) / sigma**2)
    exact = np.exp(-0.5 * (
        (xx - center[0] - velocity[0] * dt * steps)**2
        + (yy - center[1] - velocity[1] * dt * steps)**2
        + (zz - center[2] - velocity[2] * dt * steps)**2
    ) / sigma**2)
    faces = (
        np.full((shape[0], shape[1], shape[2] + 1), velocity[0]),
        np.full((shape[0], shape[1] + 1, shape[2]), velocity[1]),
        np.full((shape[0] + 1, shape[1], shape[2]), velocity[2]),
    )
    solid = np.zeros(shape, dtype=bool)

    def advance(scheme):
        field = initial.copy()
        outflow = 0.0
        for _ in range(steps):
            field, step_outflow = _transport_step(
                field, u_face=faces[0], v_face=faces[1], w_face=faces[2],
                solid=solid, diffusivity=0.0, dt=dt,
                spacing=np.full(3, dx), advection_scheme=scheme,
            )
            outflow += step_outflow
        return field, outflow

    hancock, _ = advance("muscl_van_leer")
    ssprk2, outflow = advance("muscl_ssprk2")
    hancock_error = float(np.sqrt(np.mean((hancock - exact)**2)))
    ssprk2_error = float(np.sqrt(np.mean((ssprk2 - exact)**2)))
    cell_volume = dx**3

    assert ssprk2_error < hancock_error
    assert np.min(ssprk2) >= -1.0e-14
    assert np.max(ssprk2) <= 1.0 + 1.0e-12
    assert float(np.sum(ssprk2) * cell_volume + outflow) == pytest.approx(
        float(np.sum(initial) * cell_volume), abs=1.0e-12,
    )


def test_transport_rejects_unknown_advection_scheme():
    shape = (3, 3, 5)
    with pytest.raises(ValueError, match="advection_scheme"):
        _transport_step(
            np.zeros(shape),
            u_face=_uniform_x_faces(shape, 0.0)[0],
            v_face=_uniform_x_faces(shape, 0.0)[1],
            w_face=_uniform_x_faces(shape, 0.0)[2],
            solid=np.zeros(shape, dtype=bool),
            diffusivity=0.0, dt=0.1, spacing=np.ones(3),
            advection_scheme="unknown",
        )


def test_lateral_pressure_open_face_advects_scalar_out_and_closes_ledger():
    shape = (3, 4, 5)
    state = np.zeros(shape)
    state[:, -1, :] = 2.0
    solid = np.zeros(shape, dtype=bool)
    u_face, v_face, w_face = _uniform_x_faces(shape, 0.0)
    v_face[:, -1, :] = 0.25
    dt = 0.1

    result, outflow = _transport_step(
        state,
        u_face=u_face, v_face=v_face, w_face=w_face,
        solid=solid, diffusivity=0.0, dt=dt, spacing=np.ones(3),
        open_boundary_faces=("y_max",),
    )

    initial_mass = float(np.sum(state))
    final_mass = float(np.sum(result))
    assert outflow > 0.0
    assert final_mass + outflow == pytest.approx(initial_mass, abs=1.0e-13)
    assert np.all(result[:, -1, :] < state[:, -1, :])


def test_closed_lateral_face_does_not_leak_scalar():
    shape = (3, 4, 5)
    state = np.zeros(shape)
    state[:, -1, :] = 2.0
    solid = np.zeros(shape, dtype=bool)
    u_face, v_face, w_face = _uniform_x_faces(shape, 0.0)
    v_face[:, -1, :] = 0.25

    result, outflow = _transport_step(
        state,
        u_face=u_face, v_face=v_face, w_face=w_face,
        solid=solid, diffusivity=0.0, dt=0.1, spacing=np.ones(3),
        signed_boundary_flux=True,
    )

    assert outflow == 0.0
    assert np.array_equal(result, state)


def test_embedded_blocked_face_suppresses_advective_and_diffusive_scalar_flux():
    state = np.zeros((3, 3, 3), dtype=float)
    state[:, :, 0] = 1.0
    solid = np.zeros_like(state, dtype=bool)
    u_face, v_face, w_face = _uniform_x_faces(state.shape, 0.2)
    blocked_x = np.zeros_like(u_face, dtype=bool)
    blocked_y = np.zeros_like(v_face, dtype=bool)
    blocked_z = np.zeros_like(w_face, dtype=bool)
    blocked_x[:, :, 1] = True

    isolated, outflow = _transport_step(
        state.copy(),
        u_face=u_face, v_face=v_face, w_face=w_face,
        solid=solid, diffusivity=0.1, dt=0.1, spacing=np.ones(3),
        blocked_faces=(blocked_x, blocked_y, blocked_z),
    )
    permeable, _ = _transport_step(
        state.copy(),
        u_face=u_face, v_face=v_face, w_face=w_face,
        solid=solid, diffusivity=0.1, dt=0.1, spacing=np.ones(3),
    )

    assert np.array_equal(isolated, state)
    assert outflow == 0.0
    assert not np.array_equal(permeable, state)
    assert float(isolated.sum()) == pytest.approx(float(state.sum()))


def test_embedded_blocked_face_masks_require_staggered_shapes():
    state = np.zeros((2, 3, 4), dtype=float)
    faces = _uniform_x_faces(state.shape, 0.0)
    with pytest.raises(ValueError, match="staggered face arrays"):
        _transport_step(
            state,
            u_face=faces[0], v_face=faces[1], w_face=faces[2],
            solid=np.zeros_like(state, dtype=bool), diffusivity=0.0,
            dt=0.1, spacing=np.ones(3),
            blocked_faces=(np.zeros_like(state, dtype=bool), faces[1], faces[2]),
        )


def test_pressure_projection_keeps_embedded_face_closed_and_mass_consistent():
    shape = (1, 3, 4)
    solid = np.zeros(shape, dtype=bool)
    u_face, v_face, w_face = (
        np.ones((1, 3, 5)), np.zeros((1, 4, 4)), np.zeros((2, 3, 4)),
    )
    blocked_x = np.zeros_like(u_face, dtype=bool)
    blocked_y = np.zeros_like(v_face, dtype=bool)
    blocked_z = np.zeros_like(w_face, dtype=bool)
    blocked_x[0, 1, 2] = True

    projected_u, projected_v, projected_w, before, after, _ = (
        _mass_consistent_face_projection_3d(
            u_face, v_face, w_face, solid, np.ones(3),
            blocked_faces=(blocked_x, blocked_y, blocked_z),
        )
    )
    divergence = (
        projected_u[:, :, 1:] - projected_u[:, :, :-1]
        + projected_v[:, 1:, :] - projected_v[:, :-1, :]
        + projected_w[1:, :, :] - projected_w[:-1, :, :]
    )

    assert before > 0.0
    assert after < before
    assert np.max(np.abs(divergence[:, :, :-1])) < 1.0e-8
    assert projected_u[0, 1, 2] == 0.0


def test_vectorized_pressure_projection_closes_irregular_obstacle_faces():
    rng = np.random.default_rng(20261013)
    shape = (4, 5, 6)
    solid = np.zeros(shape, dtype=bool)
    solid[1:3, 1:4, 2] = True
    solid[2, 2, 3:5] = True
    u_face = rng.normal(size=(4, 5, 7))
    v_face = rng.normal(size=(4, 6, 6))
    w_face = rng.normal(size=(5, 5, 6))
    blocked_x = np.zeros_like(u_face, dtype=bool)
    blocked_y = np.zeros_like(v_face, dtype=bool)
    blocked_z = np.zeros_like(w_face, dtype=bool)
    blocked_x[0, 2, 3] = True
    blocked_y[1, 3, 4] = True
    blocked_z[2, 2, 1] = True

    projected = _mass_consistent_face_projection_3d(
        u_face, v_face, w_face, solid, np.asarray((0.7, 1.1, 0.9)),
        blocked_faces=(blocked_x, blocked_y, blocked_z),
    )
    projected_u, projected_v, projected_w, before, after, iterations = projected
    divergence = (
        (projected_u[:, :, 1:] - projected_u[:, :, :-1]) / 0.7
        + (projected_v[:, 1:, :] - projected_v[:, :-1, :]) / 1.1
        + (projected_w[1:, :, :] - projected_w[:-1, :, :]) / 0.9
    )

    assert before > 0.0
    assert after < before
    assert after < 1.0e-8
    assert iterations > 0
    assert projected_u[blocked_x] == pytest.approx(0.0)
    assert projected_v[blocked_y] == pytest.approx(0.0)
    assert projected_w[blocked_z] == pytest.approx(0.0)
    unknown = ~solid
    unknown[:, :, -1] = False
    assert np.max(np.abs(divergence[unknown])) < 1.0e-8


def test_spatial_diffusivity_uses_conservative_face_average():
    shape = (1, 1, 4)
    state = np.asarray([[[1.0, 0.0, 0.0, 0.0]]])
    diffusivity = np.asarray([[[0.0, 2.0, 4.0, 4.0]]])
    solid = np.zeros(shape, dtype=bool)
    u_face, v_face, w_face = _uniform_x_faces(shape, 0.0)

    result, outflow = _transport_step(
        state,
        u_face=u_face, v_face=v_face, w_face=w_face,
        solid=solid, diffusivity=diffusivity, dt=0.1,
        spacing=np.ones(3),
    )

    assert result == pytest.approx(np.asarray([[[0.9, 0.1, 0.0, 0.0]]]))
    assert np.sum(result) == pytest.approx(np.sum(state), abs=1.0e-14)


def test_prescribed_inlet_state_enters_conservatively_and_rejects_bad_shape():
    shape = (1, 1, 2)
    state = np.zeros(shape)
    solid = np.zeros(shape, dtype=bool)
    u_face, v_face, w_face = _uniform_x_faces(shape, 1.0)
    result, signed_outflow = _transport_step(
        state,
        u_face=u_face, v_face=v_face, w_face=w_face,
        solid=solid, diffusivity=0.0, dt=0.1, spacing=np.ones(3),
        signed_boundary_flux=True,
        open_boundary_faces=("x_min",),
        prescribed_x_min_inflow_state=np.asarray([[10.0]]),
    )

    assert np.sum(result) == pytest.approx(1.0)
    assert signed_outflow == pytest.approx(-1.0)
    assert np.sum(result) + signed_outflow == pytest.approx(0.0)
    with pytest.raises(ValueError, match=r"shape \(nz, ny\)"):
        _transport_step(
            state,
            u_face=u_face, v_face=v_face, w_face=w_face,
            solid=solid, diffusivity=0.0, dt=0.1, spacing=np.ones(3),
            prescribed_x_min_inflow_state=np.zeros((2, 1)),
        )


def test_spatial_diffusivity_rejects_wrong_shape_and_negative_values():
    shape = (2, 2, 3)
    state = np.zeros(shape)
    solid = np.zeros(shape, dtype=bool)
    u_face, v_face, w_face = _uniform_x_faces(shape, 0.0)
    arguments = dict(
        concentration=state,
        u_face=u_face, v_face=v_face, w_face=w_face,
        solid=solid, dt=0.1, spacing=np.ones(3),
    )

    with pytest.raises(ValueError, match="same shape"):
        _transport_step(diffusivity=np.zeros((2, 2, 2)), **arguments)
    invalid = np.zeros(shape)
    invalid[0, 0, 0] = -1.0
    with pytest.raises(ValueError, match="finite and non-negative"):
        _transport_step(diffusivity=invalid, **arguments)
