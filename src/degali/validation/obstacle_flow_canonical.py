"""Canonical numerical checks for the resolved obstacle-flow tier.

These checks verify numerical building blocks independently of hydrogen data.
They do not validate facility dispersion.  The Poiseuille case exercises the
same D2Q9 equilibrium, BGK collision and link bounce-back used by the Tier-2
obstacle solver against an analytic pressure-gradient-driven channel profile.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import time

import numpy as np

from degali.addons.resolved_obstacle_flow_2d import (
    _CS2,
    _CX,
    _CZ,
    _OPPOSITE,
    _WEIGHTS,
    _equilibrium,
)
from degali.addons.resolved_obstacle_flow_3d import (
    _C3,
    _CX3,
    _CY3,
    _CZ3,
    _OPPOSITE3,
    _CS2_3D,
    _equilibrium_d3q19,
    _guo_force_d3q19,
    _moments_d3q19,
)


POISEUILLE_CHANNEL_VERIFICATION_SCHEMA = (
    "degali.poiseuille-channel-verification.v1"
)
D3Q19_POISEUILLE_CHANNEL_VERIFICATION_SCHEMA = (
    "degali.d3q19-poiseuille-channel-verification.v1"
)


@dataclass(frozen=True)
class PoiseuilleChannelVerificationConfig:
    """Lattice-unit inputs and predeclared analytic error gates."""

    nx: int = 48
    nz: int = 24
    relaxation_time: float = 0.8
    body_acceleration_lattice: float = 2.0e-6
    steps: int = 20_000
    convergence_check_interval: int = 500
    relative_l2_tolerance: float = 0.03
    relative_linf_tolerance: float = 0.05
    streamwise_uniformity_tolerance: float = 1.0e-6

    def __post_init__(self) -> None:
        for name, minimum in (
            ("nx", 8), ("nz", 8), ("steps", 1),
            ("convergence_check_interval", 1),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer of at least {minimum}")
        for name in (
            "relaxation_time", "body_acceleration_lattice",
            "relative_l2_tolerance", "relative_linf_tolerance",
            "streamwise_uniformity_tolerance",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if self.relaxation_time <= 0.5:
            raise ValueError("relaxation_time must exceed the BGK stability limit 0.5")


@dataclass(frozen=True)
class PoiseuilleChannelVerificationResult:
    schema: str
    wall_normal_coordinate_lattice: np.ndarray
    numerical_velocity_lattice: np.ndarray
    analytic_velocity_lattice: np.ndarray
    relative_l2_error: float
    relative_linf_error: float
    relative_mass_flux_error: float
    maximum_streamwise_nonuniformity: float
    final_convergence_residual: float
    maximum_relative_density_variation: float
    runtime_s: float
    passed: bool
    diagnostics: dict[str, object]


@dataclass(frozen=True)
class D3Q19PoiseuilleVerificationConfig:
    """Periodic x/y channel with two no-slip z walls in lattice units."""

    nx: int = 16
    ny: int = 8
    nz: int = 18
    relaxation_time: float = 0.8
    body_acceleration_lattice: float = 2.0e-6
    steps: int = 15_000
    convergence_check_interval: int = 500
    relative_l2_tolerance: float = 0.03
    relative_linf_tolerance: float = 0.05
    horizontal_uniformity_tolerance: float = 1.0e-6

    def __post_init__(self) -> None:
        for name, minimum in (
            ("nx", 8), ("ny", 4), ("nz", 8), ("steps", 1),
            ("convergence_check_interval", 1),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer of at least {minimum}")
        for name in (
            "relaxation_time", "body_acceleration_lattice",
            "relative_l2_tolerance", "relative_linf_tolerance",
            "horizontal_uniformity_tolerance",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if self.relaxation_time <= 0.5:
            raise ValueError("relaxation_time must exceed 0.5")


@dataclass(frozen=True)
class D3Q19PoiseuilleVerificationResult:
    schema: str
    wall_normal_coordinate_lattice: np.ndarray
    numerical_velocity_lattice: np.ndarray
    analytic_velocity_lattice: np.ndarray
    relative_l2_error: float
    relative_linf_error: float
    relative_mass_flux_error: float
    maximum_horizontal_nonuniformity: float
    final_convergence_residual: float
    maximum_relative_density_variation: float
    runtime_s: float
    passed: bool
    diagnostics: dict[str, object]


def _moments_with_force(
    distribution: np.ndarray,
    solid: np.ndarray,
    acceleration: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rho = np.maximum(np.sum(distribution, axis=0), 1.0e-14)
    force_x = rho * acceleration
    ux = (
        np.tensordot(_CX, distribution, axes=(0, 0)) + 0.5 * force_x
    ) / rho
    uz = np.tensordot(_CZ, distribution, axes=(0, 0)) / rho
    ux[solid] = 0.0
    uz[solid] = 0.0
    return rho, ux, uz


def _guo_force(
    rho: np.ndarray,
    ux: np.ndarray,
    uz: np.ndarray,
    acceleration: float,
    relaxation_time: float,
) -> np.ndarray:
    force_x = rho * acceleration
    result = np.empty((9, *rho.shape), dtype=float)
    prefactor = 1.0 - 0.5 / relaxation_time
    for direction, (cx, cz, weight) in enumerate(zip(_CX, _CZ, _WEIGHTS)):
        cu = cx * ux + cz * uz
        # F=(rho*a, 0); Guo et al. forcing in lattice units.
        projection = (3.0 * (cx - ux) + 9.0 * cu * cx) * force_x
        result[direction] = weight * prefactor * projection
    return result


def verify_poiseuille_channel(
    config: PoiseuilleChannelVerificationConfig = (
        PoiseuilleChannelVerificationConfig()
    ),
) -> PoiseuilleChannelVerificationResult:
    """Compare periodic body-force channel flow with the analytic parabola."""
    if not isinstance(config, PoiseuilleChannelVerificationConfig):
        raise TypeError("config must be PoiseuilleChannelVerificationConfig")
    started = time.perf_counter()
    solid = np.zeros((config.nz, config.nx), dtype=bool)
    solid[0, :] = True
    solid[-1, :] = True
    rho = np.ones((config.nz, config.nx), dtype=float)
    ux = np.zeros_like(rho)
    uz = np.zeros_like(rho)
    distribution = _equilibrium(rho, ux, uz)
    previous_profile = np.zeros(config.nz - 2, dtype=float)
    residual = math.inf

    for step in range(1, config.steps + 1):
        rho, ux, uz = _moments_with_force(
            distribution, solid, config.body_acceleration_lattice,
        )
        equilibrium = _equilibrium(rho, ux, uz)
        post_collision = (
            distribution
            - (distribution - equilibrium) / config.relaxation_time
            + _guo_force(
                rho, ux, uz, config.body_acceleration_lattice,
                config.relaxation_time,
            )
        )
        streamed = np.empty_like(distribution)
        for direction in range(9):
            streamed[direction] = np.roll(
                post_collision[direction],
                shift=(_CZ[direction], _CX[direction]),
                axis=(0, 1),
            )
        for direction in range(1, 9):
            destination_solid = np.roll(
                solid,
                shift=(-_CZ[direction], -_CX[direction]),
                axis=(0, 1),
            )
            bounce = (~solid) & destination_solid
            streamed[_OPPOSITE[direction], bounce] = post_collision[direction, bounce]
        distribution = streamed
        if step % config.convergence_check_interval == 0 or step == config.steps:
            _, current_ux, _ = _moments_with_force(
                distribution, solid, config.body_acceleration_lattice,
            )
            profile = np.mean(current_ux[1:-1, :], axis=1)
            scale = max(float(np.max(np.abs(profile))), 1.0e-15)
            residual = float(np.max(np.abs(profile - previous_profile)) / scale)
            previous_profile = profile

    rho, ux, uz = _moments_with_force(
        distribution, solid, config.body_acceleration_lattice,
    )
    numerical = np.mean(ux[1:-1, :], axis=1)
    # Link bounce-back puts the hydrodynamic walls halfway between the first
    # fluid row and the solid row.
    channel_height = float(config.nz - 2)
    coordinate = np.arange(config.nz - 2, dtype=float) + 0.5
    viscosity = (config.relaxation_time - 0.5) * _CS2
    analytic = (
        config.body_acceleration_lattice
        * coordinate * (channel_height - coordinate)
        / (2.0 * viscosity)
    )
    difference = numerical - analytic
    l2 = float(np.linalg.norm(difference) / np.linalg.norm(analytic))
    linf = float(np.max(np.abs(difference)) / np.max(np.abs(analytic)))
    mass_flux_error = float(
        abs(np.sum(numerical) - np.sum(analytic)) / abs(np.sum(analytic))
    )
    nonuniformity = float(
        np.max(np.ptp(ux[1:-1, :], axis=1))
        / max(float(np.max(np.abs(numerical))), 1.0e-15)
    )
    density_variation = float(
        (np.max(rho[~solid]) - np.min(rho[~solid])) / np.mean(rho[~solid])
    )
    passed = bool(
        l2 <= config.relative_l2_tolerance
        and linf <= config.relative_linf_tolerance
        and nonuniformity <= config.streamwise_uniformity_tolerance
    )
    runtime = time.perf_counter() - started
    return PoiseuilleChannelVerificationResult(
        schema=POISEUILLE_CHANNEL_VERIFICATION_SCHEMA,
        wall_normal_coordinate_lattice=coordinate,
        numerical_velocity_lattice=numerical,
        analytic_velocity_lattice=analytic,
        relative_l2_error=l2,
        relative_linf_error=linf,
        relative_mass_flux_error=mass_flux_error,
        maximum_streamwise_nonuniformity=nonuniformity,
        final_convergence_residual=residual,
        maximum_relative_density_variation=density_variation,
        runtime_s=runtime,
        passed=passed,
        diagnostics={
            "verification_only": True,
            "facility_validation": False,
            "flow_model": "periodic D2Q9 BGK with Guo body forcing",
            "wall_boundary": "link bounce-back",
            "kinematic_viscosity_lattice": viscosity,
            "channel_height_lattice": channel_height,
            "steps": config.steps,
            "acceptance": {
                "relative_l2_tolerance": config.relative_l2_tolerance,
                "relative_linf_tolerance": config.relative_linf_tolerance,
                "streamwise_uniformity_tolerance": (
                    config.streamwise_uniformity_tolerance
                ),
            },
        },
    )


def verify_d3q19_poiseuille_channel(
    config: D3Q19PoiseuilleVerificationConfig = (
        D3Q19PoiseuilleVerificationConfig()
    ),
) -> D3Q19PoiseuilleVerificationResult:
    """Verify the D3Q19 equilibrium, forcing and link bounce-back together."""
    if not isinstance(config, D3Q19PoiseuilleVerificationConfig):
        raise TypeError("config must be D3Q19PoiseuilleVerificationConfig")
    started = time.perf_counter()
    solid = np.zeros((config.nz, config.ny, config.nx), dtype=bool)
    solid[0, :, :] = True
    solid[-1, :, :] = True
    rho = np.ones_like(solid, dtype=float)
    ux = np.zeros_like(rho)
    uy = np.zeros_like(rho)
    uz = np.zeros_like(rho)
    distribution = _equilibrium_d3q19(rho, ux, uy, uz)
    acceleration = (config.body_acceleration_lattice, 0.0, 0.0)
    previous_profile = np.zeros(config.nz - 2, dtype=float)
    residual = math.inf
    tau = np.full_like(rho, config.relaxation_time)

    for step in range(1, config.steps + 1):
        rho, ux, uy, uz = _moments_d3q19(distribution, solid, acceleration)
        equilibrium = _equilibrium_d3q19(rho, ux, uy, uz)
        post_collision = (
            distribution
            - (distribution - equilibrium) / config.relaxation_time
            + _guo_force_d3q19(rho, ux, uy, uz, acceleration, tau)
        )
        streamed = np.empty_like(distribution)
        for direction, (cx, cy, cz) in enumerate(_C3):
            streamed[direction] = np.roll(
                post_collision[direction],
                shift=(cz, cy, cx),
                axis=(0, 1, 2),
            )
        for direction in range(1, 19):
            cx, cy, cz = _C3[direction]
            destination_solid = np.roll(
                solid,
                shift=(-cz, -cy, -cx),
                axis=(0, 1, 2),
            )
            bounce = (~solid) & destination_solid
            streamed[_OPPOSITE3[direction], bounce] = post_collision[direction, bounce]
        distribution = streamed
        if step % config.convergence_check_interval == 0 or step == config.steps:
            _, current_u, _current_v, _current_w = _moments_d3q19(
                distribution, solid, acceleration,
            )
            profile = np.mean(current_u[1:-1, :, :], axis=(1, 2))
            scale = max(float(np.max(np.abs(profile))), 1.0e-15)
            residual = float(np.max(np.abs(profile - previous_profile)) / scale)
            previous_profile = profile

    rho, ux, uy, uz = _moments_d3q19(distribution, solid, acceleration)
    numerical = np.mean(ux[1:-1, :, :], axis=(1, 2))
    channel_height = float(config.nz - 2)
    coordinate = np.arange(config.nz - 2, dtype=float) + 0.5
    viscosity = (config.relaxation_time - 0.5) * _CS2_3D
    analytic = (
        config.body_acceleration_lattice
        * coordinate * (channel_height - coordinate)
        / (2.0 * viscosity)
    )
    difference = numerical - analytic
    l2 = float(np.linalg.norm(difference) / np.linalg.norm(analytic))
    linf = float(np.max(np.abs(difference)) / np.max(np.abs(analytic)))
    mass_flux_error = float(
        abs(np.sum(numerical) - np.sum(analytic)) / abs(np.sum(analytic))
    )
    horizontal_mean = np.mean(ux[1:-1, :, :], axis=(1, 2))
    nonuniformity = float(
        np.max(np.abs(ux[1:-1, :, :] - horizontal_mean[:, None, None]))
        / max(float(np.max(np.abs(horizontal_mean))), 1.0e-15)
    )
    fluid = ~solid
    density_variation = float(
        (np.max(rho[fluid]) - np.min(rho[fluid])) / np.mean(rho[fluid])
    )
    passed = bool(
        l2 <= config.relative_l2_tolerance
        and linf <= config.relative_linf_tolerance
        and nonuniformity <= config.horizontal_uniformity_tolerance
    )
    return D3Q19PoiseuilleVerificationResult(
        schema=D3Q19_POISEUILLE_CHANNEL_VERIFICATION_SCHEMA,
        wall_normal_coordinate_lattice=coordinate,
        numerical_velocity_lattice=numerical,
        analytic_velocity_lattice=analytic,
        relative_l2_error=l2,
        relative_linf_error=linf,
        relative_mass_flux_error=mass_flux_error,
        maximum_horizontal_nonuniformity=nonuniformity,
        final_convergence_residual=residual,
        maximum_relative_density_variation=density_variation,
        runtime_s=time.perf_counter() - started,
        passed=passed,
        diagnostics={
            "verification_only": True,
            "facility_validation": False,
            "flow_model": "periodic D3Q19 BGK with Guo body forcing",
            "wall_boundary": "halfway link bounce-back",
            "periodic_axes": ["x", "y"],
            "kinematic_viscosity_lattice": viscosity,
            "channel_height_lattice": channel_height,
            "steps": config.steps,
            "acceptance": {
                "relative_l2_tolerance": config.relative_l2_tolerance,
                "relative_linf_tolerance": config.relative_linf_tolerance,
                "horizontal_uniformity_tolerance": config.horizontal_uniformity_tolerance,
            },
        },
    )


__all__ = [
    "POISEUILLE_CHANNEL_VERIFICATION_SCHEMA",
    "PoiseuilleChannelVerificationConfig",
    "PoiseuilleChannelVerificationResult",
    "verify_poiseuille_channel",
    "D3Q19_POISEUILLE_CHANNEL_VERIFICATION_SCHEMA",
    "D3Q19PoiseuilleVerificationConfig",
    "D3Q19PoiseuilleVerificationResult",
    "verify_d3q19_poiseuille_channel",
]
