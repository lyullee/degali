"""Shared-clock callback bridge for a fast outer and resolved local solver.

The bridge is intentionally explicit about its limitations: it couples only
declared Cartesian interface planes and requires both solvers to advance on
the same scalar clock.  Any missing flux channel or timeout is an error, not a
silent one-way fallback.
"""
from __future__ import annotations

from dataclasses import dataclass
import threading
from typing import Callable

import numpy as np

from .hybrid_conservative_exchange import ConservativeInterfacePacket
from .hybrid_flux_packets import (
    CartesianInterfacePlane,
    apply_incoming_momentum_to_velocity,
    apply_incoming_momentum_characteristic_to_velocity,
    apply_incoming_packet_to_fields,
    outgoing_packet_from_face_fluxes,
)


@dataclass(frozen=True)
class SolverBridgeSide:
    name: str
    plane: CartesianInterfacePlane
    target_mask: np.ndarray
    cell_volume_m3: float
    density_provider: Callable[[tuple[int, int, int]], np.ndarray]
    nominal_density_kg_m3: float
    expects_enthalpy: bool
    momentum_application: str


@dataclass
class _StepRecord:
    total_mass_fluxes: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
    species_mass_fluxes: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
    sensible_enthalpy_fluxes: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
    momentum_fluxes: tuple[tuple[np.ndarray, ...], ...] | None = None
    enthalpy_seen: bool = False


class SynchronizedSolverBridge:
    """Exchange one-step interface packets between two concurrent solvers."""

    def __init__(self, *, time_step_s: float, timeout_s: float = 30.0):
        self.time_step_s = float(time_step_s)
        self.timeout_s = float(timeout_s)
        if self.time_step_s <= 0.0 or self.timeout_s <= 0.0:
            raise ValueError("time_step_s and timeout_s must be positive")
        self._condition = threading.Condition()
        self._sides: dict[str, SolverBridgeSide] = {}
        self._records: dict[int, dict[str, _StepRecord]] = {}
        self._incoming: dict[int, dict[str, ConservativeInterfacePacket]] = {}
        self._exchanged: set[int] = set()
        self._errors: list[BaseException] = []

    def register_side(
        self,
        *,
        name: str,
        plane: CartesianInterfacePlane,
        target_mask: np.ndarray,
        cell_volume_m3: float,
        density_provider: Callable[[tuple[int, int, int]], np.ndarray] | None = None,
        nominal_density_kg_m3: float = 1.2,
        expects_enthalpy: bool = False,
        momentum_application: str = "additive",
    ) -> None:
        if name in self._sides or len(self._sides) >= 2:
            raise ValueError("bridge requires exactly two uniquely named sides")
        mask = np.asarray(target_mask, dtype=bool)
        if mask.ndim != 3 or not np.any(mask):
            raise ValueError("target_mask must be a non-empty 3-D mask")
        volume = float(cell_volume_m3)
        density = float(nominal_density_kg_m3)
        if volume <= 0.0 or density <= 0.0:
            raise ValueError("cell volume and nominal density must be positive")
        if momentum_application not in {"additive", "characteristic"}:
            raise ValueError(
                "momentum_application must be 'additive' or 'characteristic'"
            )
        provider = density_provider or (
            lambda shape: np.full(shape, density, dtype=float)
        )
        self._sides[name] = SolverBridgeSide(
            name=name,
            plane=plane,
            target_mask=mask,
            cell_volume_m3=volume,
            density_provider=provider,
            nominal_density_kg_m3=density,
            expects_enthalpy=bool(expects_enthalpy),
            momentum_application=momentum_application,
        )

    def _step_index(self, start: float) -> int:
        index = int(round(float(start) / self.time_step_s))
        if abs(index * self.time_step_s - float(start)) > 1.0e-8:
            raise ValueError("solver callback time is not on the shared bridge clock")
        return index

    def _record(self, step: int, side: str) -> _StepRecord:
        if side not in self._sides:
            raise KeyError(f"unregistered bridge side: {side}")
        records = self._records.setdefault(step, {})
        return records.setdefault(side, _StepRecord())

    def _wait_for_exchange(self, step: int) -> None:
        with self._condition:
            end = __import__("time").monotonic() + self.timeout_s
            while step not in self._exchanged:
                if self._errors:
                    raise RuntimeError("synchronized solver bridge failed") from self._errors[0]
                if len(self._records.get(step, {})) == len(self._sides) and all(
                    record.species_mass_fluxes is not None
                    and record.total_mass_fluxes is not None
                    and (
                        not self._sides[name].expects_enthalpy
                        or record.enthalpy_seen
                    )
                    for name, record in self._records[step].items()
                ):
                    names = tuple(self._sides)
                    packets = {}
                    for name in names:
                        record = self._records[step][name]
                        side = self._sides[name]
                        packets[name] = outgoing_packet_from_face_fluxes(
                            total_mass_fluxes=record.total_mass_fluxes,
                            species_mass_fluxes=record.species_mass_fluxes,
                            plane=side.plane,
                            time_step_s=self.time_step_s,
                            sensible_enthalpy_fluxes=record.sensible_enthalpy_fluxes,
                            momentum_fluxes=record.momentum_fluxes,
                        )
                    self._incoming[step] = {
                        names[0]: packets[names[1]],
                        names[1]: packets[names[0]],
                    }
                    self._exchanged.add(step)
                    self._condition.notify_all()
                    return
                remaining = end - __import__("time").monotonic()
                if remaining <= 0.0:
                    record_state = {
                        name: {
                            "total": record.total_mass_fluxes is not None,
                            "species": record.species_mass_fluxes is not None,
                            "momentum": record.momentum_fluxes is not None,
                        }
                        for name, record in self._records.get(step, {}).items()
                    }
                    raise TimeoutError(
                        f"solver bridge timed out waiting for step {step}; "
                        f"records={record_state}"
                    )
                self._condition.wait(timeout=remaining)

    def _incoming_for(self, side: str, step: int) -> ConservativeInterfacePacket:
        if step == 0:
            return ConservativeInterfacePacket.zero()
        self._wait_for_exchange(step - 1)
        return self._incoming[step - 1][side]

    def _momentum_fluxes_from_velocity(
        self,
        side: SolverBridgeSide,
        u: np.ndarray,
        v: np.ndarray,
        w: np.ndarray,
        u_face: np.ndarray,
        v_face: np.ndarray,
        w_face: np.ndarray,
    ) -> tuple[tuple[np.ndarray, ...], ...]:
        density = np.asarray(side.density_provider(u.shape), dtype=float)
        if density.shape != u.shape or np.any(~np.isfinite(density)) or np.any(density <= 0.0):
            raise ValueError("density_provider returned an invalid cell density field")
        rho_x = np.pad(
            0.5 * (density[:, :, :-1] + density[:, :, 1:]),
            ((0, 0), (0, 0), (1, 1)), mode="edge",
        )
        rho_y = np.pad(
            0.5 * (density[:, :-1, :] + density[:, 1:, :]),
            ((0, 0), (1, 1), (0, 0)), mode="edge",
        )
        rho_z = np.pad(
            0.5 * (density[:-1, :, :] + density[1:, :, :]),
            ((1, 1), (0, 0), (0, 0)), mode="edge",
        )
        return (
            (rho_x * u_face * u_face,
             np.zeros_like(v_face),
             np.zeros_like(w_face)),
            (np.zeros_like(u_face),
             rho_y * v_face * v_face,
             np.zeros_like(w_face)),
            (np.zeros_like(u_face),
             np.zeros_like(v_face),
             rho_z * w_face * w_face),
        )

    def callbacks_for(self, side_name: str) -> dict[str, Callable]:
        if side_name not in self._sides:
            raise KeyError(f"unregistered bridge side: {side_name}")
        side = self._sides[side_name]
        state: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

        def conserved_state_callback(concentration, sensible_enthalpy, start, dt):
            step = self._step_index(start)
            incoming = self._incoming_for(side_name, step)
            updated, enthalpy = apply_incoming_packet_to_fields(
                concentration,
                packet=incoming,
                target_mask=side.target_mask,
                cell_volume_m3=side.cell_volume_m3,
                sensible_enthalpy_j_m3=sensible_enthalpy,
            )
            state[step] = (updated, enthalpy, incoming)
            return (
                updated, enthalpy
            )

        def conserved_mass_callback(start, dt):
            step = self._step_index(start)
            incoming = self._incoming_for(side_name, step)
            return incoming.species_mass_kg, incoming.sensible_enthalpy_j

        def scalar_state_callback(concentration, start, dt):
            step = self._step_index(start)
            incoming = self._incoming_for(side_name, step)
            updated, _ = apply_incoming_packet_to_fields(
                concentration,
                packet=incoming,
                target_mask=side.target_mask,
                cell_volume_m3=side.cell_volume_m3,
            )
            state[step] = (updated, np.empty(0), incoming)
            return updated

        def scalar_state_mass_callback(start, dt):
            step = self._step_index(start)
            return self._incoming_for(side_name, step).species_mass_kg

        def momentum_callback(u, v, w, u_face, v_face, w_face, start, dt):
            step = self._step_index(start)
            incoming = self._incoming_for(side_name, step)
            density = np.asarray(side.density_provider(u.shape), dtype=float)
            if side.momentum_application == "characteristic":
                u, v, w, u_face, v_face, w_face = (
                    apply_incoming_momentum_characteristic_to_velocity(
                        u, v, w, u_face, v_face, w_face,
                        packet=incoming,
                        target_mask=side.target_mask,
                        density_kg_m3=density,
                        cell_volume_m3=side.cell_volume_m3,
                        plane=side.plane,
                    )
                )
            else:
                u, v, w = apply_incoming_momentum_to_velocity(
                    u, v, w,
                    packet=incoming,
                    target_mask=side.target_mask,
                    density_kg_m3=density,
                    cell_volume_m3=side.cell_volume_m3,
                )
            # Use the post-injection velocities to form the total and momentum
            # flux schedule for this step.  The solver's scalar callback supplies
            # the species and (for thermal runs) enthalpy flux later.
            normal = side.plane.axis
            velocity_faces = (u_face, v_face, w_face)
            rho = side.nominal_density_kg_m3
            total = tuple(rho * np.asarray(value, dtype=float) for value in velocity_faces)
            record = self._record(step, side_name)
            record.total_mass_fluxes = total
            record.momentum_fluxes = self._momentum_fluxes_from_velocity(
                side, u, v, w, u_face, v_face, w_face,
            )
            return u, v, w, u_face, v_face, w_face

        def face_flux_callback(flux_x, flux_y, flux_z):
            step = max(state.keys(), default=0)
            record = self._record(step, side_name)
            record.species_mass_fluxes = (
                np.asarray(flux_x, dtype=float).copy(),
                np.asarray(flux_y, dtype=float).copy(),
                np.asarray(flux_z, dtype=float).copy(),
            )

        def enthalpy_flux_callback(flux_x, flux_y, flux_z):
            step = max(state.keys(), default=0)
            record = self._record(step, side_name)
            record.sensible_enthalpy_fluxes = (
                np.asarray(flux_x, dtype=float).copy(),
                np.asarray(flux_y, dtype=float).copy(),
                np.asarray(flux_z, dtype=float).copy(),
            )
            record.enthalpy_seen = True

        return {
            "interface_state_callback": scalar_state_callback,
            "interface_state_mass_callback": scalar_state_mass_callback,
            "interface_conserved_state_callback": conserved_state_callback,
            "interface_conserved_mass_callback": conserved_mass_callback,
            "interface_momentum_callback": momentum_callback,
            "face_flux_callback": face_flux_callback,
            "enthalpy_flux_callback": enthalpy_flux_callback,
        }

    @property
    def exchanged_steps(self) -> int:
        return len(self._exchanged)

    def finalize_step(self, step: int) -> None:
        """Wait for and ledger the final step that has no next receiver step."""
        self._wait_for_exchange(int(step))


__all__ = ["SolverBridgeSide", "SynchronizedSolverBridge"]
