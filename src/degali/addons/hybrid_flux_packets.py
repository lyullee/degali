"""Turn face-resolved solver fluxes into conservative hybrid packets."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from .hybrid_conservative_exchange import ConservativeInterfacePacket


@dataclass(frozen=True)
class CartesianInterfacePlane:
    """A solver face plane and its outward orientation."""

    axis: int
    face_index: int
    outward_sign: float
    face_area_m2: float

    def __post_init__(self) -> None:
        if int(self.axis) not in (0, 1, 2):
            raise ValueError("axis must be 0, 1, or 2")
        if float(self.outward_sign) not in (-1.0, 1.0):
            raise ValueError("outward_sign must be -1 or +1")
        area = float(self.face_area_m2)
        if not math.isfinite(area) or area <= 0.0:
            raise ValueError("face_area_m2 must be positive and finite")


def _face_values(field: np.ndarray, plane: CartesianInterfacePlane) -> np.ndarray:
    array = np.asarray(field, dtype=float)
    if array.ndim != 3 or not np.all(np.isfinite(array)):
        raise ValueError("face flux arrays must be finite 3-D arrays")
    # Flux arrays are stored as (z, y, x-face), (z, y-face, x),
    # (z-face, y, x), while ``plane.axis`` is the physical normal x/y/z.
    array_axis = (2, 1, 0)[int(plane.axis)]
    return np.take(array, int(plane.face_index), axis=array_axis)


def outgoing_packet_from_face_fluxes(
    *,
    total_mass_fluxes: Sequence[np.ndarray],
    species_mass_fluxes: Sequence[np.ndarray],
    plane: CartesianInterfacePlane,
    time_step_s: float,
    sensible_enthalpy_fluxes: Sequence[np.ndarray] | None = None,
    momentum_fluxes: Sequence[Sequence[np.ndarray]] | None = None,
) -> ConservativeInterfacePacket:
    """Integrate positive outward flux on one plane into one packet.

    ``total_mass_fluxes`` and ``species_mass_fluxes`` are ordered ``(x, y,
    z)`` and use the same staggered face shape.  Optional enthalpy and
    momentum fluxes are integrated only over faces carrying positive outward
    total mass, preventing a reverse-flow face from creating an unrelated
    energy or momentum packet.
    """

    dt = float(time_step_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("time_step_s must be positive and finite")
    if len(total_mass_fluxes) != 3 or len(species_mass_fluxes) != 3:
        raise ValueError("flux sequences must contain x, y, and z arrays")
    total = plane.outward_sign * _face_values(
        total_mass_fluxes[plane.axis], plane,
    )
    species = plane.outward_sign * _face_values(
        species_mass_fluxes[plane.axis], plane,
    )
    outgoing = total > 0.0
    total_mass = float(np.sum(np.maximum(total, 0.0)) * plane.face_area_m2 * dt)
    species_mass = float(
        np.sum(np.maximum(species, 0.0) * outgoing)
        * plane.face_area_m2
        * dt
    )
    if species_mass > total_mass + 1.0e-10 * max(total_mass, 1.0):
        raise ValueError("outgoing species mass exceeds outgoing total mass")
    enthalpy = 0.0
    if sensible_enthalpy_fluxes is not None:
        if len(sensible_enthalpy_fluxes) != 3:
            raise ValueError("enthalpy flux sequence must contain x, y, and z arrays")
        enthalpy_values = plane.outward_sign * _face_values(
            sensible_enthalpy_fluxes[plane.axis], plane,
        )
        enthalpy = float(
            np.sum(enthalpy_values * outgoing) * plane.face_area_m2 * dt
        )
    momentum = np.zeros(3, dtype=float)
    if momentum_fluxes is not None:
        if len(momentum_fluxes) != 3 or any(len(item) != 3 for item in momentum_fluxes):
            raise ValueError("momentum_fluxes must be three vector components of x/y/z arrays")
        for component in range(3):
            values = plane.outward_sign * _face_values(
                momentum_fluxes[component][plane.axis], plane,
            )
            momentum[component] = float(
                np.sum(values * outgoing) * plane.face_area_m2 * dt
            )
    return ConservativeInterfacePacket(
        total_mass_kg=total_mass,
        species_mass_kg=max(0.0, species_mass),
        sensible_enthalpy_j=enthalpy,
        momentum_kg_m_s=tuple(float(value) for value in momentum),
    )


def boundary_cell_mask(
    shape: tuple[int, int, int], *, axis: int, side: str, layers: int = 1,
) -> np.ndarray:
    """Return a cell strip receiving an incoming interface packet."""

    if len(shape) != 3 or any(int(value) < 1 for value in shape):
        raise ValueError("shape must contain three positive cell counts")
    if int(axis) not in (0, 1, 2) or side not in {"min", "max"}:
        raise ValueError("axis must be 0/1/2 and side must be min/max")
    count = int(layers)
    if count < 1 or count > int(shape[axis]):
        raise ValueError("layers must fit within the selected axis")
    mask = np.zeros(shape, dtype=bool)
    selector = [slice(None)] * 3
    selector[axis] = slice(0, count) if side == "min" else slice(-count, None)
    mask[tuple(selector)] = True
    return mask


def apply_incoming_packet_to_fields(
    concentration_kg_m3: np.ndarray,
    *,
    packet: ConservativeInterfacePacket,
    target_mask: np.ndarray,
    cell_volume_m3: float,
    sensible_enthalpy_j_m3: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Distribute incoming species/enthalpy over a boundary strip."""

    concentration = np.asarray(concentration_kg_m3, dtype=float).copy()
    mask = np.asarray(target_mask, dtype=bool)
    volume = float(cell_volume_m3)
    if concentration.ndim != 3 or mask.shape != concentration.shape:
        raise ValueError("concentration and target_mask must be matching 3-D arrays")
    if not np.all(np.isfinite(concentration)) or np.any(concentration < 0.0):
        raise ValueError("concentration must be finite and non-negative")
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("cell_volume_m3 must be positive and finite")
    count = int(np.count_nonzero(mask))
    if count == 0 and (packet.species_mass_kg > 0.0 or packet.sensible_enthalpy_j != 0.0):
        raise ValueError("positive incoming packet requires target cells")
    if count:
        concentration[mask] += packet.species_mass_kg / (count * volume)
    if sensible_enthalpy_j_m3 is None:
        return concentration, None
    enthalpy = np.asarray(sensible_enthalpy_j_m3, dtype=float).copy()
    if enthalpy.shape != concentration.shape or not np.all(np.isfinite(enthalpy)):
        raise ValueError("sensible enthalpy field must be finite and solver-shaped")
    if count:
        enthalpy[mask] += packet.sensible_enthalpy_j / (count * volume)
    return concentration, enthalpy


def apply_incoming_momentum_to_velocity(
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    w_m_s: np.ndarray,
    *,
    packet: ConservativeInterfacePacket,
    target_mask: np.ndarray,
    density_kg_m3: np.ndarray,
    cell_volume_m3: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply an integrated momentum packet to a density-weighted cell strip."""

    velocity = [np.asarray(value, dtype=float).copy() for value in (u_m_s, v_m_s, w_m_s)]
    mask = np.asarray(target_mask, dtype=bool)
    density = np.asarray(density_kg_m3, dtype=float)
    if any(value.shape != mask.shape for value in velocity) or density.shape != mask.shape:
        raise ValueError("velocity, density, and target_mask must have matching shapes")
    if not np.all(np.isfinite(density)) or np.any(density <= 0.0):
        raise ValueError("density must be finite and positive")
    volume = float(cell_volume_m3)
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("cell_volume_m3 must be positive and finite")
    target_mass = float(np.sum(density[mask]) * volume)
    if target_mass <= 0.0 and any(abs(value) > 0.0 for value in packet.momentum_kg_m_s):
        raise ValueError("positive incoming momentum requires target mass")
    if target_mass > 0.0:
        delta = np.asarray(packet.momentum_kg_m_s, dtype=float) / target_mass
        for index in range(3):
            velocity[index][mask] += delta[index]
    return tuple(velocity)  # type: ignore[return-value]


def apply_incoming_momentum_characteristic_to_velocity(
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    w_m_s: np.ndarray,
    u_face_m_s: np.ndarray,
    v_face_m_s: np.ndarray,
    w_face_m_s: np.ndarray,
    *,
    packet: ConservativeInterfacePacket,
    target_mask: np.ndarray,
    density_kg_m3: np.ndarray,
    cell_volume_m3: float,
    plane: CartesianInterfacePlane,
    relaxation: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Map an integrated packet to a characteristic velocity boundary.

    The additive helper above is a conservative momentum impulse, but applying
    that impulse repeatedly to the same coarse boundary strip can accumulate
    an unphysical velocity.  This mode instead reconstructs the incoming
    carrier velocity from packet momentum divided by packet total mass and
    relaxes the boundary cells toward that state.  The normal staggered face
    is updated on the declared Cartesian plane; tangential faces are left to
    the local solver reconstruction.
    """
    velocity = [
        np.asarray(value, dtype=float).copy()
        for value in (u_m_s, v_m_s, w_m_s)
    ]
    faces = [
        np.asarray(value, dtype=float).copy()
        for value in (u_face_m_s, v_face_m_s, w_face_m_s)
    ]
    mask = np.asarray(target_mask, dtype=bool)
    density = np.asarray(density_kg_m3, dtype=float)
    if any(value.shape != mask.shape for value in velocity) or density.shape != mask.shape:
        raise ValueError("velocity, density, and target_mask must have matching shapes")
    if any(value.ndim != 3 or not np.all(np.isfinite(value)) for value in (*velocity, *faces)):
        raise ValueError("velocity and face fields must be finite 3-D arrays")
    if not np.all(np.isfinite(density)) or np.any(density <= 0.0):
        raise ValueError("density must be finite and positive")
    alpha = float(relaxation)
    if not math.isfinite(alpha) or not 0.0 <= alpha <= 1.0:
        raise ValueError("relaxation must lie in [0, 1]")
    volume = float(cell_volume_m3)
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("cell_volume_m3 must be positive and finite")
    if np.count_nonzero(mask) == 0:
        return (*velocity, *faces)  # type: ignore[return-value]
    incoming_mass = float(packet.total_mass_kg)
    momentum = np.asarray(packet.momentum_kg_m_s, dtype=float)
    if incoming_mass <= 1.0e-18:
        # No carrier mass crossed the face, so there is no characteristic
        # velocity to impose.  In particular, do not divide a numerical
        # momentum round-off by an empty packet.
        return (*velocity, *faces)  # type: ignore[return-value]
    target_velocity = momentum / incoming_mass
    for component in range(3):
        velocity[component][mask] = (
            (1.0 - alpha) * velocity[component][mask]
            + alpha * target_velocity[component]
        )
    normal = int(plane.axis)
    array_axis = (2, 1, 0)[normal]
    face_index = int(plane.face_index)
    selector = [slice(None)] * 3
    selector[array_axis] = face_index
    faces[normal][tuple(selector)] = (
        (1.0 - alpha) * faces[normal][tuple(selector)]
        + alpha * target_velocity[normal]
    )
    return (*velocity, *faces)  # type: ignore[return-value]


__all__ = [
    "CartesianInterfacePlane",
    "outgoing_packet_from_face_fluxes",
    "boundary_cell_mask",
    "apply_incoming_packet_to_fields",
    "apply_incoming_momentum_to_velocity",
    "apply_incoming_momentum_characteristic_to_velocity",
]
