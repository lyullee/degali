"""Fast, opt-in source histories for cryo-compressed hydrogen blowdown.

This is deliberately upstream of the dispersion solver.  A stationary
integral jet cannot turn a tank pressure record into a physically consistent
transient source on its own.  The model below is a well-mixed, fixed-volume
control volume with the real-gas properties used elsewhere in DEGALI.  It
implements the mass and internal-energy update used by Cirrone et al. (2023)
for their PRESLHY cryo-compressed hydrogen blowdown model:

``d(m u)/dt = Qdot - h * mdot`` and ``dm/dt = -mdot``.

The optional ``tank_wall_ua_w_k`` is a *lumped prescribed-wall* heat path;
it is not a replacement for that paper's resolved tank-wall and discharge
pipe conduction model.  In particular, a caller must provide pipe geometry
and thermal data before claiming a pipe-heating prediction.  The default
has zero heat transfer and is useful as a transparent adiabatic lower layer
for a source-bound sensitivity study.

An independently prescribed discharge-pipe wall temperature and ``UA`` can
also be supplied.  That option applies the pipe energy balance before the
critical nozzle calculation by a fixed-point solve.  It is intentionally not
enabled from a gas thermocouple record: a welded flow thermocouple is not a
pipe-wall boundary condition.

The optional two-phase continuation is deliberately narrower still.  It is a
well-mixed equilibrium tank with **vapour withdrawal only**.  It needs the
tank diameter and outlet elevation, and stops when its calculated liquid
level reaches that outlet.  It therefore cannot silently reinterpret a
liquid-withdrawal release as a gas release.

No core DEGADIS calculation selects this model.  It is only an explicit
research boundary for cryo-compressed, single-component hydrogen releases.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING

import numpy as np

from .notional import isentropic_throat

if TYPE_CHECKING:
    from ..lh2 import LH2ExpandedSource
    from .lh2_droplets import (
        FlashingHydrogenDropletSource,
        HomogeneousEquilibriumHydrogenSource,
    )


@dataclass(frozen=True)
class TransientTankWall:
    """One-dimensional planar tank-wall thermal boundary.

    ``inner_heat_transfer_w_m2_k`` and ``outer_heat_transfer_w_m2_k`` are
    deliberately explicit closures.  The wall solves conduction and stores
    heat; it does not silently select a natural/forced-convection correlation
    for a geometry it cannot observe.
    """

    thickness_m: float
    area_m2: float
    density_kg_m3: float
    specific_heat_j_kg_k: float
    conductivity_w_m_k: float
    initial_temperature_k: float
    external_temperature_k: float
    inner_heat_transfer_w_m2_k: float
    outer_heat_transfer_w_m2_k: float
    nodes: int = 20


@dataclass(frozen=True)
class TransientPipeWall:
    """Thin-wall one-dimensional pipe thermal store with explicit convection."""

    length_m: float
    inner_diameter_m: float
    thickness_m: float
    density_kg_m3: float
    specific_heat_j_kg_k: float
    conductivity_w_m_k: float
    initial_temperature_k: float
    external_temperature_k: float
    inner_heat_transfer_w_m2_k: float
    outer_heat_transfer_w_m2_k: float
    nodes: int = 5


@dataclass(frozen=True)
class CryogenicBlowdownConfig:
    """Explicit inputs for a well-mixed cryo-compressed H2 vessel."""

    vessel_volume_m3: float
    nozzle_diameter_m: float
    initial_temperature_k: float
    initial_pressure_pa: float
    discharge_coefficient: float
    duration_s: float
    time_step_s: float = 0.01
    ambient_pressure_pa: float = 101325.0
    tank_wall_temperature_k: float | None = None
    tank_wall_ua_w_k: float = 0.0
    transient_tank_wall: TransientTankWall | None = None
    pipe_wall_temperature_k: float | None = None
    pipe_wall_ua_w_k: float = 0.0
    pipe_inner_diameter_m: float | None = None
    transient_pipe_wall: TransientPipeWall | None = None
    stop_pressure_ratio: float = 1.003
    two_phase_withdrawal: str = "stop"
    tank_internal_diameter_m: float | None = None
    outlet_height_from_bottom_m: float | None = None


@dataclass(frozen=True)
class CryogenicBlowdownState:
    """One tank state immediately before its indicated source flux."""

    time_s: float
    mass_kg: float
    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    specific_internal_energy_j_kg: float
    specific_enthalpy_j_kg: float
    mass_flow_kg_s: float
    heat_transfer_w: float
    pipe_heat_transfer_w: float
    nozzle_inlet_temperature_k: float
    tank_wall_inner_temperature_k: float | None
    vapour_quality: float | None
    liquid_height_m: float | None


@dataclass(frozen=True)
class _HomogeneousEquilibriumThroat:
    """Private HEM critical state used consistently by tank and source paths."""

    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    specific_enthalpy_j_kg: float
    vapour_quality: float
    mass_flux_kg_m2_s: float


@dataclass(frozen=True)
class CryogenicBlowdownResult:
    """A conservative source history and the reason integration stopped."""

    config: CryogenicBlowdownConfig
    states: tuple[CryogenicBlowdownState, ...]
    termination: str
    maximum_discrete_energy_residual_j: float
    final_tank_wall_temperature_profile_k: tuple[float, ...] | None
    final_pipe_wall_temperature_profile_k: tuple[float, ...] | None

    @property
    def released_mass_kg(self) -> float:
        return self.states[0].mass_kg - self.states[-1].mass_kg


def _validate(config: CryogenicBlowdownConfig) -> None:
    positive = {
        "vessel_volume_m3": config.vessel_volume_m3,
        "nozzle_diameter_m": config.nozzle_diameter_m,
        "initial_temperature_k": config.initial_temperature_k,
        "initial_pressure_pa": config.initial_pressure_pa,
        "duration_s": config.duration_s,
        "time_step_s": config.time_step_s,
        "ambient_pressure_pa": config.ambient_pressure_pa,
        "stop_pressure_ratio": config.stop_pressure_ratio,
    }
    for name, value in positive.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite")
    if not np.isfinite(config.discharge_coefficient) or not (
        0.0 < config.discharge_coefficient <= 1.0
    ):
        raise ValueError("discharge_coefficient must lie in (0, 1]")
    if config.initial_pressure_pa <= config.ambient_pressure_pa:
        raise ValueError("initial_pressure_pa must exceed ambient_pressure_pa")
    if not np.isfinite(config.tank_wall_ua_w_k) or config.tank_wall_ua_w_k < 0.0:
        raise ValueError("tank_wall_ua_w_k must be finite and non-negative")
    if config.tank_wall_ua_w_k > 0.0 and (
        config.tank_wall_temperature_k is None
        or not np.isfinite(config.tank_wall_temperature_k)
        or config.tank_wall_temperature_k <= 0.0
    ):
        raise ValueError("positive tank_wall_ua_w_k requires tank_wall_temperature_k")
    if config.transient_tank_wall is not None and config.tank_wall_ua_w_k > 0.0:
        raise ValueError("choose either transient_tank_wall or tank_wall_ua_w_k")
    if config.transient_tank_wall is not None:
        _validate_transient_wall(config.transient_tank_wall)
    if not np.isfinite(config.pipe_wall_ua_w_k) or config.pipe_wall_ua_w_k < 0.0:
        raise ValueError("pipe_wall_ua_w_k must be finite and non-negative")
    if config.pipe_wall_ua_w_k > 0.0 and (
        config.pipe_wall_temperature_k is None
        or not np.isfinite(config.pipe_wall_temperature_k)
        or config.pipe_wall_temperature_k <= 0.0
        or config.pipe_inner_diameter_m is None
        or not np.isfinite(config.pipe_inner_diameter_m)
        or config.pipe_inner_diameter_m <= 0.0
    ):
        raise ValueError(
            "positive pipe_wall_ua_w_k requires pipe_wall_temperature_k and pipe_inner_diameter_m"
        )
    if config.transient_pipe_wall is not None and config.pipe_wall_ua_w_k > 0.0:
        raise ValueError("choose either transient_pipe_wall or pipe_wall_ua_w_k")
    if config.transient_pipe_wall is not None:
        _validate_transient_pipe_wall(config.transient_pipe_wall)
    if config.two_phase_withdrawal not in {"stop", "vapour", "homogeneous"}:
        raise ValueError(
            "two_phase_withdrawal must be 'stop', 'vapour', or 'homogeneous'"
        )
    if config.two_phase_withdrawal == "vapour":
        for name, value in {
            "tank_internal_diameter_m": config.tank_internal_diameter_m,
            "outlet_height_from_bottom_m": config.outlet_height_from_bottom_m,
        }.items():
            if value is None or not np.isfinite(value) or value < 0.0:
                raise ValueError(f"vapour withdrawal requires {name}")
        assert config.tank_internal_diameter_m is not None
        assert config.outlet_height_from_bottom_m is not None
        if config.outlet_height_from_bottom_m <= 0.0:
            raise ValueError("vapour withdrawal requires a positive outlet height")
    if config.two_phase_withdrawal in {"vapour", "homogeneous"} and (
        config.pipe_wall_ua_w_k > 0.0 or config.transient_pipe_wall is not None
    ):
        raise ValueError(
            "two-phase withdrawal does not yet support a pipe thermal boundary"
        )


def _validate_transient_wall(wall: TransientTankWall) -> None:
    positive = {
        "thickness_m": wall.thickness_m,
        "area_m2": wall.area_m2,
        "density_kg_m3": wall.density_kg_m3,
        "specific_heat_j_kg_k": wall.specific_heat_j_kg_k,
        "conductivity_w_m_k": wall.conductivity_w_m_k,
        "initial_temperature_k": wall.initial_temperature_k,
        "external_temperature_k": wall.external_temperature_k,
        "inner_heat_transfer_w_m2_k": wall.inner_heat_transfer_w_m2_k,
    }
    for name, value in positive.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"transient tank wall {name} must be positive and finite")
    if not np.isfinite(wall.outer_heat_transfer_w_m2_k) or wall.outer_heat_transfer_w_m2_k < 0.0:
        raise ValueError("transient tank wall outer_heat_transfer_w_m2_k must be finite and non-negative")
    if not isinstance(wall.nodes, int) or wall.nodes < 3:
        raise ValueError("transient tank wall nodes must be an integer of at least three")


def _validate_transient_pipe_wall(wall: TransientPipeWall) -> None:
    _validate_transient_wall(
        TransientTankWall(
            thickness_m=wall.thickness_m,
            area_m2=math.pi * (wall.inner_diameter_m + wall.thickness_m) * wall.length_m,
            density_kg_m3=wall.density_kg_m3,
            specific_heat_j_kg_k=wall.specific_heat_j_kg_k,
            conductivity_w_m_k=wall.conductivity_w_m_k,
            initial_temperature_k=wall.initial_temperature_k,
            external_temperature_k=wall.external_temperature_k,
            inner_heat_transfer_w_m2_k=wall.inner_heat_transfer_w_m2_k,
            outer_heat_transfer_w_m2_k=wall.outer_heat_transfer_w_m2_k,
            nodes=wall.nodes,
        )
    )
    if not np.isfinite(wall.length_m) or wall.length_m <= 0.0:
        raise ValueError("transient pipe wall length_m must be positive and finite")
    if not np.isfinite(wall.inner_diameter_m) or wall.inner_diameter_m <= 0.0:
        raise ValueError("transient pipe wall inner_diameter_m must be positive and finite")


def _make_pipe_wall_state(wall: TransientPipeWall) -> _TankWallState:
    return _TankWallState(
        TransientTankWall(
            thickness_m=wall.thickness_m,
            area_m2=math.pi * (wall.inner_diameter_m + wall.thickness_m) * wall.length_m,
            density_kg_m3=wall.density_kg_m3,
            specific_heat_j_kg_k=wall.specific_heat_j_kg_k,
            conductivity_w_m_k=wall.conductivity_w_m_k,
            initial_temperature_k=wall.initial_temperature_k,
            external_temperature_k=wall.external_temperature_k,
            inner_heat_transfer_w_m2_k=wall.inner_heat_transfer_w_m2_k,
            outer_heat_transfer_w_m2_k=wall.outer_heat_transfer_w_m2_k,
            nodes=wall.nodes,
        )
    )


class _TankWallState:
    """Finite-volume wall state; positive flux is into the hydrogen tank."""

    def __init__(self, wall: TransientTankWall) -> None:
        self.wall = wall
        self.temperature = np.full(wall.nodes, wall.initial_temperature_k, dtype=float)
        self.dx = wall.thickness_m / wall.nodes
        self.cell_heat_capacity = (
            wall.density_kg_m3 * wall.specific_heat_j_kg_k * wall.area_m2 * self.dx
        )

    def heat_to_tank_w(self, tank_temperature_k: float) -> float:
        return self.wall.inner_heat_transfer_w_m2_k * self.wall.area_m2 * (
            self.temperature[-1] - tank_temperature_k
        )

    def advance(self, duration_s: float, tank_temperature_k: float) -> None:
        """Advance explicit finite volumes under a conservative CFL substep."""
        alpha = self.wall.conductivity_w_m_k / (
            self.wall.density_kg_m3 * self.wall.specific_heat_j_kg_k
        )
        stable_step = 0.45 * self.dx**2 / alpha
        steps = max(1, int(math.ceil(duration_s / stable_step)))
        dt = duration_s / steps
        conductance = self.wall.conductivity_w_m_k * self.wall.area_m2 / self.dx
        for _ in range(steps):
            flux = conductance * (self.temperature[:-1] - self.temperature[1:])
            outer = self.wall.outer_heat_transfer_w_m2_k * self.wall.area_m2 * (
                self.wall.external_temperature_k - self.temperature[0]
            )
            inner = self.heat_to_tank_w(tank_temperature_k)
            rate = np.empty_like(self.temperature)
            rate[0] = outer - flux[0]
            rate[1:-1] = flux[:-1] - flux[1:]
            rate[-1] = flux[-1] - inner
            self.temperature += dt * rate / self.cell_heat_capacity
        if not np.all(np.isfinite(self.temperature)) or np.any(self.temperature <= 0.0):
            raise RuntimeError("transient tank wall produced a non-physical temperature")


def _saturated_vapour_throat_mass_flux(
    *, pressure_pa: float, ambient_pressure_pa: float
) -> float:
    """Return the critical mass flux from a saturated-vapour stagnation state."""
    from CoolProp.CoolProp import PropsSI
    from scipy.optimize import minimize_scalar

    h0 = float(PropsSI("H", "P", pressure_pa, "Q", 1, "Hydrogen"))
    s0 = float(PropsSI("S", "P", pressure_pa, "Q", 1, "Hydrogen"))

    def negative_mass_flux(log_pressure: float) -> float:
        pressure = math.exp(log_pressure)
        enthalpy = float(PropsSI("H", "P", pressure, "S", s0, "Hydrogen"))
        density = float(PropsSI("D", "P", pressure, "S", s0, "Hydrogen"))
        velocity = math.sqrt(max(2.0 * (h0 - enthalpy), 0.0))
        return -(density * velocity)

    low, high = math.log(ambient_pressure_pa), math.log(pressure_pa)
    optimum = minimize_scalar(
        negative_mass_flux, bounds=(low, high), method="bounded", options={"xatol": 1.0e-12}
    )
    return max(-float(optimum.fun), 0.0)


def _homogeneous_equilibrium_throat(
    *, density_kg_m3: float, specific_internal_energy_j_kg: float,
    storage_pressure_pa: float, ambient_pressure_pa: float,
) -> _HomogeneousEquilibriumThroat:
    """Return HEM critical mass flux and throat enthalpy for a tank mixture.

    The liquid and vapour are assumed to remain in local equilibrium and to
    share one velocity through the throat.  This is an explicit homogeneous-
    equilibrium model, not a separated-flow or flashing-kinetics closure.
    """
    from CoolProp.CoolProp import PropsSI
    from scipy.optimize import minimize_scalar

    h0 = float(PropsSI(
        "H", "D", density_kg_m3, "U", specific_internal_energy_j_kg, "Hydrogen"
    ))
    entropy = float(PropsSI(
        "S", "D", density_kg_m3, "U", specific_internal_energy_j_kg, "Hydrogen"
    ))

    def state(pressure_pa: float) -> tuple[float, float, float, float, float]:
        enthalpy = float(PropsSI("H", "P", pressure_pa, "S", entropy, "Hydrogen"))
        density = float(PropsSI("D", "P", pressure_pa, "S", entropy, "Hydrogen"))
        temperature = float(PropsSI("T", "P", pressure_pa, "S", entropy, "Hydrogen"))
        quality = float(PropsSI("Q", "P", pressure_pa, "S", entropy, "Hydrogen"))
        velocity = math.sqrt(max(2.0 * (h0 - enthalpy), 0.0))
        return density * velocity, enthalpy, density, temperature, quality

    lower, upper = math.log(ambient_pressure_pa), math.log(storage_pressure_pa)
    optimum = minimize_scalar(
        lambda log_pressure: -state(math.exp(log_pressure))[0],
        bounds=(lower, upper), method="bounded", options={"xatol": 1.0e-12},
    )
    candidates = (lower, float(optimum.x), upper)
    selected = max(candidates, key=lambda value: state(math.exp(value))[0])
    pressure = math.exp(selected)
    mass_flux, enthalpy, density, temperature, quality = state(pressure)
    if not 0.0 <= quality <= 1.0:
        raise RuntimeError("HEM critical state is not two-phase")
    return _HomogeneousEquilibriumThroat(
        pressure_pa=pressure,
        temperature_k=temperature,
        density_kg_m3=density,
        specific_enthalpy_j_kg=enthalpy,
        vapour_quality=quality,
        mass_flux_kg_m2_s=mass_flux,
    )


def _liquid_height(
    *, mass_kg: float, vapour_quality: float, pressure_pa: float, diameter_m: float
) -> float:
    """Equilibrium liquid level in a cylindrical, well-mixed tank."""
    from CoolProp.CoolProp import PropsSI

    liquid_density = float(PropsSI("D", "P", pressure_pa, "Q", 0, "Hydrogen"))
    area = math.pi * diameter_m**2 / 4.0
    return mass_kg * (1.0 - vapour_quality) / max(liquid_density * area, 1.0e-30)


def _pipe_heated_mass_flow(
    *,
    storage_temperature_k: float,
    storage_pressure_pa: float,
    storage_enthalpy_j_kg: float,
    discharge_coefficient: float,
    nozzle_area_m2: float,
    ambient_pressure_pa: float,
    pipe_wall_temperature_k: float | None,
    pipe_wall_ua_w_k: float,
    pipe_inner_diameter_m: float | None,
    prescribed_pipe_heat_transfer_w: float | None = None,
) -> tuple[float, float, float]:
    """Fixed-point version of the published pipe energy balance (Eqs. 5--7)."""
    from CoolProp.CoolProp import PropsSI
    from scipy.optimize import brentq

    throat = isentropic_throat(
        fluid="Hydrogen",
        storage_temperature=storage_temperature_k,
        storage_pressure=storage_pressure_pa,
        ambient_pressure=ambient_pressure_pa,
        allow_supercritical_gas=True,
    )
    mass_flow = discharge_coefficient * throat.mass_flux * nozzle_area_m2
    if pipe_wall_ua_w_k == 0.0 and prescribed_pipe_heat_transfer_w is None:
        return mass_flow, storage_temperature_k, 0.0
    assert pipe_inner_diameter_m is not None
    pipe_area = math.pi * pipe_inner_diameter_m**2 / 4.0
    if prescribed_pipe_heat_transfer_w is None:
        assert pipe_wall_temperature_k is not None
        heat_transfer = pipe_wall_ua_w_k * (pipe_wall_temperature_k - storage_temperature_k)
    else:
        heat_transfer = float(prescribed_pipe_heat_transfer_w)

    def mapped_mass_flow(current_flow: float) -> tuple[float, float]:
        specific_heat = heat_transfer / max(current_flow, 1.0e-30)
        h2 = storage_enthalpy_j_kg + specific_heat
        density2 = float(PropsSI("D", "H", h2, "P", storage_pressure_pa, "Hydrogen"))
        velocity2 = current_flow / max(pipe_area * density2, 1.0e-30)
        h2 -= 0.5 * velocity2**2
        inlet_temperature = float(
            PropsSI("T", "H", h2, "P", storage_pressure_pa, "Hydrogen")
        )
        throat = isentropic_throat(
            fluid="Hydrogen",
            storage_temperature=inlet_temperature,
            storage_pressure=storage_pressure_pa,
            ambient_pressure=ambient_pressure_pa,
            allow_supercritical_gas=True,
        )
        return discharge_coefficient * throat.mass_flux * nozzle_area_m2, inlet_temperature

    def residual(current_flow: float) -> float:
        return mapped_mass_flow(current_flow)[0] - current_flow

    lower, upper = mass_flow * 0.1, mass_flow * 10.0
    try:
        solved_flow = float(brentq(residual, lower, upper, xtol=1.0e-10, rtol=1.0e-8))
    except ValueError as error:
        raise RuntimeError(
            "pipe heat-flow balance has no positive bracket; revise the prescribed pipe boundary"
        ) from error
    _mapped, inlet_temperature = mapped_mass_flow(solved_flow)
    return solved_flow, inlet_temperature, heat_transfer


def blowdown_state_to_lh2_source(
    config: CryogenicBlowdownConfig,
    state: CryogenicBlowdownState,
    *,
    theta: float = math.pi / 2.0,
    x: float = 0.0,
    y: float = 0.0,
) -> "LH2ExpandedSource":
    """Convert one single-phase blowdown state to a conserved gas source plane.

    The adapter uses the same pressure-preserving discharge-line convention as
    the blowdown calculation and expands the reconstructed critical throat to
    ambient pressure through the existing mass, momentum, and total-energy
    conservation mapping.  ``Cd`` is represented as an effective throat area
    (``A_eff = Cd * A_orifice``), so the supplied state mass flow remains
    consistent with its throat density and velocity.

    It is deliberately a snapshot adapter, not a transient atmospheric solver.
    A two-phase tank, a zero-flow terminal state, or an ambient two-phase
    expansion has no single gas source plane under this model and is rejected.
    """
    if state.mass_flow_kg_s <= 0.0:
        raise ValueError("a positive-flow blowdown state is required for a source plane")
    if state.vapour_quality is not None:
        raise ValueError(
            "a two-phase blowdown state has no single-phase gas source plane; "
            "supply an explicit multiphase source model"
        )
    throat = isentropic_throat(
        fluid="Hydrogen",
        storage_temperature=state.nozzle_inlet_temperature_k,
        storage_pressure=state.pressure_pa,
        ambient_pressure=config.ambient_pressure_pa,
        allow_supercritical_gas=True,
    )
    effective_diameter = config.nozzle_diameter_m * math.sqrt(
        config.discharge_coefficient
    )
    from ..lh2 import lh2_source_from_measured_throat

    expanded = lh2_source_from_measured_throat(
        throat_diameter=effective_diameter,
        throat_pressure=throat.pressure,
        throat_temperature=throat.temperature,
        throat_density=throat.density,
        throat_velocity=throat.velocity,
        throat_enthalpy=throat.enthalpy,
        mass_flow=state.mass_flow_kg_s,
        ambient_pressure=config.ambient_pressure_pa,
        theta=theta,
        x=x,
        y=y,
    )
    from CoolProp.CoolProp import PropsSI

    atmospheric_quality = float(PropsSI(
        "Q", "P", config.ambient_pressure_pa, "H", expanded.expansion.enthalpy,
        "Hydrogen",
    ))
    if 0.0 <= atmospheric_quality < 1.0 - 1.0e-10:
        raise ValueError(
            "the ambient-pressure expansion is two-phase; use "
            "blowdown_state_to_homogeneous_evaporation_source() or an "
            "explicit finite-rate multiphase source"
        )
    return expanded


def blowdown_state_to_flashing_droplet_source(
    config: CryogenicBlowdownConfig,
    state: CryogenicBlowdownState,
    *,
    ambient_temperature_k: float = 295.0,
    droplet_size_coefficient: float = 15.0,
) -> "FlashingHydrogenDropletSource":
    """Map one single-phase tank state to its pressure-thrust flash plane.

    The returned source retains its ambient-pressure liquid and vapour mass
    flows, post-flash velocity, and the explicitly selected atomisation
    diagnostic.  It is not yet a single-phase atmospheric jet.  Use the
    homogeneous-evaporation adapter only when its fast complete-evaporation
    bound is appropriate, or use the phase/slip research path separately.

    The declared ``Cd`` is expressed once as an effective flow area.  The tank
    itself must still be single phase; a two-phase withdrawal needs an
    independently specified source law.
    """
    if state.mass_flow_kg_s <= 0.0:
        raise ValueError("a positive-flow blowdown state is required for a source plane")
    if state.vapour_quality is not None:
        raise ValueError(
            "a two-phase tank state requires an explicit withdrawal source law"
        )
    effective_diameter = config.nozzle_diameter_m * math.sqrt(
        config.discharge_coefficient
    )
    from .lh2_droplets import flashing_hydrogen_droplet_source

    return flashing_hydrogen_droplet_source(
        mass_flow=state.mass_flow_kg_s,
        orifice_diameter=effective_diameter,
        upstream_temperature=state.nozzle_inlet_temperature_k,
        upstream_pressure=state.pressure_pa,
        ambient_temperature=ambient_temperature_k,
        ambient_pressure=config.ambient_pressure_pa,
        droplet_size_coefficient=droplet_size_coefficient,
    )


def blowdown_state_to_homogeneous_evaporation_source(
    config: CryogenicBlowdownConfig,
    state: CryogenicBlowdownState,
    *,
    ambient_temperature_k: float = 295.0,
    theta: float = 0.0,
    x: float = 0.0,
    y: float = 0.0,
    droplet_size_coefficient: float = 15.0,
) -> "HomogeneousEquilibriumHydrogenSource":
    """Map a tank state to the explicit fast complete-evaporation bound.

    This first calls :func:`blowdown_state_to_flashing_droplet_source` and
    then applies the common-velocity, heat-limited complete-evaporation
    source.  It is a fast homogeneous-equilibrium bound, not a prediction of
    finite droplet evaporation or gas--liquid slip.
    """
    postflash = blowdown_state_to_flashing_droplet_source(
        config,
        state,
        ambient_temperature_k=ambient_temperature_k,
        droplet_size_coefficient=droplet_size_coefficient,
    )
    from .lh2_droplets import homogeneous_equilibrium_hydrogen_source_from_postflash

    return homogeneous_equilibrium_hydrogen_source_from_postflash(
        postflash,
        ambient_temperature=ambient_temperature_k,
        ambient_pressure=config.ambient_pressure_pa,
        theta=theta,
        x=x,
        y=y,
    )


def hem_blowdown_state_to_flashing_droplet_source(
    config: CryogenicBlowdownConfig,
    state: CryogenicBlowdownState,
    *,
    ambient_temperature_k: float = 295.0,
    droplet_size_coefficient: float = 15.0,
) -> "FlashingHydrogenDropletSource":
    """Map an explicit HEM tank state to its downstream flash plane.

    This adapter is available only for the opt-in homogeneous-equilibrium
    withdrawal path.  It reconstructs its equilibrium two-phase critical
    state, preserves the declared effective flow area, and then applies the
    pressure-thrust flash.  It does not reinterpret a vapour-withdrawal or
    separated-flow state as HEM.
    """
    if config.two_phase_withdrawal != "homogeneous":
        raise ValueError(
            "HEM source adapter requires two_phase_withdrawal='homogeneous'"
        )
    if state.mass_flow_kg_s <= 0.0 or state.vapour_quality is None:
        raise ValueError("a positive-flow two-phase HEM state is required")
    throat = _homogeneous_equilibrium_throat(
        density_kg_m3=state.density_kg_m3,
        specific_internal_energy_j_kg=state.specific_internal_energy_j_kg,
        storage_pressure_pa=state.pressure_pa,
        ambient_pressure_pa=config.ambient_pressure_pa,
    )
    effective_diameter = config.nozzle_diameter_m * math.sqrt(
        config.discharge_coefficient
    )
    from .lh2_droplets import flashing_hydrogen_droplet_source

    return flashing_hydrogen_droplet_source(
        mass_flow=state.mass_flow_kg_s,
        orifice_diameter=effective_diameter,
        upstream_temperature=throat.temperature_k,
        upstream_pressure=throat.pressure_pa,
        upstream_quality=throat.vapour_quality,
        ambient_temperature=ambient_temperature_k,
        ambient_pressure=config.ambient_pressure_pa,
        droplet_size_coefficient=droplet_size_coefficient,
    )


def hem_blowdown_state_to_homogeneous_evaporation_source(
    config: CryogenicBlowdownConfig,
    state: CryogenicBlowdownState,
    *,
    ambient_temperature_k: float = 295.0,
    theta: float = 0.0,
    x: float = 0.0,
    y: float = 0.0,
    droplet_size_coefficient: float = 15.0,
) -> "HomogeneousEquilibriumHydrogenSource":
    """Map an HEM tank state to the separate fast evaporation bound.

    This is the HEM counterpart of the single-phase homogeneous-evaporation
    adapter.  It retains the same limitations: common phase velocity and
    complete evaporation after the flash are an explicit fast bound, not a
    finite-rate liquid or slip prediction.
    """
    postflash = hem_blowdown_state_to_flashing_droplet_source(
        config,
        state,
        ambient_temperature_k=ambient_temperature_k,
        droplet_size_coefficient=droplet_size_coefficient,
    )
    from .lh2_droplets import homogeneous_equilibrium_hydrogen_source_from_postflash

    return homogeneous_equilibrium_hydrogen_source_from_postflash(
        postflash,
        ambient_temperature=ambient_temperature_k,
        ambient_pressure=config.ambient_pressure_pa,
        theta=theta,
        x=x,
        y=y,
    )


def run_cryogenic_blowdown(
    config: CryogenicBlowdownConfig,
) -> CryogenicBlowdownResult:
    """Advance a real-gas, well-mixed H2 tank by mass and energy balance.

    The source flux is the published one-dimensional critical-flow bound
    multiplied by the caller-declared discharge coefficient.  It does not
    infer an effective area or select ``Cd``.  ``Cd`` therefore belongs to an
    independently declared source experiment or sensitivity envelope.  The
    default stops at an equilibrium two-phase tank state.  ``vapour`` is an
    explicit, geometry-limited continuation for a top-of-liquid outlet only.
    """
    _validate(config)
    try:
        from CoolProp.CoolProp import PropsSI
    except ImportError as error:  # pragma: no cover - environment dependent
        raise ImportError("cryogenic blowdown requires CoolProp") from error

    area = math.pi * config.nozzle_diameter_m**2 / 4.0
    density = float(
        PropsSI(
            "D", "T", config.initial_temperature_k, "P", config.initial_pressure_pa,
            "Hydrogen",
        )
    )
    mass = density * config.vessel_volume_m3
    temperature = float(config.initial_temperature_k)
    pressure = float(config.initial_pressure_pa)
    internal_energy = float(PropsSI("U", "T", temperature, "P", pressure, "Hydrogen"))
    wall_state = (
        _TankWallState(config.transient_tank_wall)
        if config.transient_tank_wall is not None
        else None
    )
    pipe_wall_state = (
        _make_pipe_wall_state(config.transient_pipe_wall)
        if config.transient_pipe_wall is not None
        else None
    )
    states: list[CryogenicBlowdownState] = []
    maximum_residual = 0.0
    time = 0.0
    termination = "duration"

    while time < config.duration_s - 1.0e-12:
        if pressure <= config.stop_pressure_ratio * config.ambient_pressure_pa:
            termination = "near_ambient_pressure"
            break
        quality_value = float(PropsSI("Q", "D", density, "U", internal_energy, "Hydrogen"))
        two_phase = 0.0 <= quality_value <= 1.0
        liquid_height: float | None = None
        if two_phase:
            if config.two_phase_withdrawal == "stop":
                termination = "two_phase_tank_boundary"
                break
            if config.two_phase_withdrawal == "vapour":
                assert config.tank_internal_diameter_m is not None
                assert config.outlet_height_from_bottom_m is not None
                liquid_height = _liquid_height(
                    mass_kg=mass,
                    vapour_quality=quality_value,
                    pressure_pa=pressure,
                    diameter_m=config.tank_internal_diameter_m,
                )
                if liquid_height >= config.outlet_height_from_bottom_m:
                    termination = "liquid_reaches_outlet"
                    break
                enthalpy = float(PropsSI("H", "P", pressure, "Q", 1, "Hydrogen"))
                mass_flow = (
                    config.discharge_coefficient * area
                    * _saturated_vapour_throat_mass_flux(
                        pressure_pa=pressure,
                        ambient_pressure_pa=config.ambient_pressure_pa,
                    )
                )
            else:
                hem_throat = _homogeneous_equilibrium_throat(
                    density_kg_m3=density,
                    specific_internal_energy_j_kg=internal_energy,
                    storage_pressure_pa=pressure,
                    ambient_pressure_pa=config.ambient_pressure_pa,
                )
                mass_flow = (
                    config.discharge_coefficient * area
                    * hem_throat.mass_flux_kg_m2_s
                )
                # The tank control-volume balance loses the stagnation
                # enthalpy of the mixture at the tank boundary.  The HEM
                # throat enthalpy determines the critical flux, but has
                # already exchanged internal energy for kinetic energy inside
                # the nozzle and must not be used as the tank outflow term.
                enthalpy = float(PropsSI(
                    "H", "D", density, "U", internal_energy, "Hydrogen"
                ))
            nozzle_inlet_temperature = temperature
            pipe_heat_transfer = 0.0
        else:
            enthalpy = float(PropsSI("H", "D", density, "U", internal_energy, "Hydrogen"))
            pipe_heat_guess = (
                pipe_wall_state.heat_to_tank_w(temperature)
                if pipe_wall_state is not None else None
            )
            mass_flow, nozzle_inlet_temperature, pipe_heat_transfer = _pipe_heated_mass_flow(
                storage_temperature_k=temperature,
                storage_pressure_pa=pressure,
                storage_enthalpy_j_kg=enthalpy,
                discharge_coefficient=config.discharge_coefficient,
                nozzle_area_m2=area,
                ambient_pressure_pa=config.ambient_pressure_pa,
                pipe_wall_temperature_k=config.pipe_wall_temperature_k,
                pipe_wall_ua_w_k=config.pipe_wall_ua_w_k,
                pipe_inner_diameter_m=(
                    config.transient_pipe_wall.inner_diameter_m
                    if config.transient_pipe_wall is not None
                    else config.pipe_inner_diameter_m
                ),
                prescribed_pipe_heat_transfer_w=pipe_heat_guess,
            )
        heat_transfer = 0.0
        tank_wall_inner_temperature = None
        if wall_state is not None:
            heat_transfer = wall_state.heat_to_tank_w(temperature)
            tank_wall_inner_temperature = float(wall_state.temperature[-1])
        elif config.tank_wall_ua_w_k > 0.0:
            assert config.tank_wall_temperature_k is not None
            heat_transfer = config.tank_wall_ua_w_k * (
                config.tank_wall_temperature_k - temperature
            )
        states.append(
            CryogenicBlowdownState(
                time_s=time,
                mass_kg=mass,
                pressure_pa=pressure,
                temperature_k=temperature,
                density_kg_m3=density,
                specific_internal_energy_j_kg=internal_energy,
                specific_enthalpy_j_kg=enthalpy,
                mass_flow_kg_s=mass_flow,
                heat_transfer_w=heat_transfer,
                pipe_heat_transfer_w=pipe_heat_transfer,
                nozzle_inlet_temperature_k=nozzle_inlet_temperature,
                tank_wall_inner_temperature_k=tank_wall_inner_temperature,
                vapour_quality=quality_value if two_phase else None,
                liquid_height_m=liquid_height,
            )
        )
        step = min(config.time_step_s, config.duration_s - time)
        removed_mass = mass_flow * step
        if removed_mass >= mass:
            termination = "source_exhausted"
            break
        next_mass = mass - removed_mass
        next_internal_energy = (
            mass * internal_energy + step * (heat_transfer - enthalpy * mass_flow)
        ) / next_mass
        discrete_residual = (
            next_mass * next_internal_energy
            - (mass * internal_energy + step * (heat_transfer - enthalpy * mass_flow))
        )
        maximum_residual = max(maximum_residual, abs(discrete_residual))
        next_density = next_mass / config.vessel_volume_m3
        try:
            next_temperature = float(
                PropsSI("T", "D", next_density, "U", next_internal_energy, "Hydrogen")
            )
            next_pressure = float(
                PropsSI("P", "D", next_density, "U", next_internal_energy, "Hydrogen")
            )
        except ValueError as error:
            raise RuntimeError(
                "blowdown state leaves the single-fluid CoolProp domain; "
                "a two-phase tank model is required"
            ) from error
        if not all(np.isfinite(value) and value > 0.0 for value in (
            next_temperature, next_pressure, next_density, next_internal_energy,
        )):
            raise RuntimeError("non-positive or non-finite blowdown state")
        time += step
        if wall_state is not None:
            wall_state.advance(step, temperature)
        if pipe_wall_state is not None:
            pipe_wall_state.advance(step, 0.5 * (temperature + nozzle_inlet_temperature))
        mass, density = next_mass, next_density
        temperature, pressure, internal_energy = (
            next_temperature, next_pressure, next_internal_energy
        )

    enthalpy = float(PropsSI("H", "D", density, "U", internal_energy, "Hydrogen"))
    states.append(
        CryogenicBlowdownState(
            time_s=time,
            mass_kg=mass,
            pressure_pa=pressure,
            temperature_k=temperature,
            density_kg_m3=density,
            specific_internal_energy_j_kg=internal_energy,
            specific_enthalpy_j_kg=enthalpy,
            mass_flow_kg_s=0.0,
            heat_transfer_w=0.0,
            pipe_heat_transfer_w=0.0,
            nozzle_inlet_temperature_k=temperature,
            tank_wall_inner_temperature_k=(
                float(wall_state.temperature[-1]) if wall_state is not None else None
            ),
            vapour_quality=(
                float(PropsSI("Q", "D", density, "U", internal_energy, "Hydrogen"))
                if 0.0 <= float(PropsSI("Q", "D", density, "U", internal_energy, "Hydrogen")) <= 1.0
                else None
            ),
            liquid_height_m=None,
        )
    )
    return CryogenicBlowdownResult(
        config=config,
        states=tuple(states),
        termination=termination,
        maximum_discrete_energy_residual_j=maximum_residual,
        final_tank_wall_temperature_profile_k=(
            tuple(float(value) for value in wall_state.temperature)
            if wall_state is not None else None
        ),
        final_pipe_wall_temperature_profile_k=(
            tuple(float(value) for value in pipe_wall_state.temperature)
            if pipe_wall_state is not None else None
        ),
    )
