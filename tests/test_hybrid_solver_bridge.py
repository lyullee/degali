from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from degali.addons.hybrid_flux_packets import CartesianInterfacePlane, boundary_cell_mask
from degali.addons.hybrid_solver_bridge import SynchronizedSolverBridge


def _faces():
    return (
        np.zeros((2, 2, 3)),
        np.zeros((2, 3, 2)),
        np.zeros((3, 2, 2)),
    )


def test_synchronized_bridge_exchanges_face_packet_on_next_shared_step():
    bridge = SynchronizedSolverBridge(time_step_s=0.1, timeout_s=2.0)
    shape = (2, 2, 2)
    bridge.register_side(
        name="outer",
        plane=CartesianInterfacePlane(0, -1, 1.0, 1.0),
        target_mask=boundary_cell_mask(shape, axis=2, side="min"),
        cell_volume_m3=1.0,
    )
    bridge.register_side(
        name="patch",
        plane=CartesianInterfacePlane(0, 0, -1.0, 1.0),
        target_mask=boundary_cell_mask(shape, axis=2, side="max"),
        cell_volume_m3=1.0,
        momentum_application="characteristic",
    )
    outer = bridge.callbacks_for("outer")
    patch = bridge.callbacks_for("patch")
    u, v, w = [np.zeros(shape) for _ in range(3)]
    ux, vx, wx = _faces()
    ux[:, :, -1] = 1.0
    up, vp, wp = [np.zeros(shape) for _ in range(3)]
    uxp, vxp, wxp = _faces()
    uxp[:, :, 0] = -1.0

    outer["interface_momentum_callback"](u, v, w, ux, vx, wx, 0.0, 0.1)
    patch["interface_momentum_callback"](up, vp, wp, uxp, vxp, wxp, 0.0, 0.1)
    species_outer = _faces()
    species_outer[0][:, :, -1] = 0.2
    species_patch = _faces()
    species_patch[0][:, :, 0] = -0.1
    outer["face_flux_callback"](*species_outer)
    patch["face_flux_callback"](*species_patch)

    def outer_state():
        return outer["interface_state_callback"](
            np.zeros(shape), 0.1, 0.1,
        )

    def patch_state():
        return patch["interface_state_callback"](
            np.zeros(shape), 0.1, 0.1,
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(lambda callback: callback(), (outer_state, patch_state)))

    assert bridge.exchanged_steps == 1
    assert np.sum(results[0]) == pytest.approx(0.04)
    assert np.sum(results[1]) == pytest.approx(0.08)
    assert outer["interface_state_mass_callback"](0.1, 0.1) == pytest.approx(0.04)
    assert patch["interface_state_mass_callback"](0.1, 0.1) == pytest.approx(0.08)

    # The characteristic mode reconstructs u_in=p_in/m_in at the patch face;
    # it does not add the same impulse again on every subsequent step.
    updated = patch["interface_momentum_callback"](
        np.zeros(shape), np.zeros(shape), np.zeros(shape),
        np.zeros_like(uxp), np.zeros_like(vxp), np.zeros_like(wxp),
        0.1, 0.1,
    )
    np.testing.assert_allclose(updated[0][..., -1], 1.0)
    np.testing.assert_allclose(updated[3][..., 0], 1.0)

