"""Conditional 3-D obstacle flow and local hydrogen transport.

The flow tier supports both D3Q19 BGK and a conservative variable-density
finite-volume projection path. D3Q19 uses interpolated bounce-back for
zero-thickness circular arcs; the finite-volume path uses blocked normal
faces and a smooth-wall immersed shear model based on cell-integrated surface
geometry. These wall and subgrid treatments are research closures, not yet
qualified facility models. Ground and outer-boundary conditions are declared
separately from obstacle-wall conditions.

Hydrogen is advanced by a local conservative 3-D finite-volume operator.  The
default neutral mode preserves the original one-way isothermal comparison.
An opt-in active-scalar mode also conserves signed sensible enthalpy, recovers
temperature and ideal-gas mixture density, and feeds vertical buoyancy back to
the D3Q19 momentum equation through the Guo forcing operator.  A fixed-point
corrector enforces the discrete variable-density continuity equation on face
mass fluxes and reconciles the pressure-corrected velocity with D3Q19 momentum.
The carrier thermodynamics remains a bounded active-scalar approximation, so
strong-density-ratio cells are diagnosed and quantitative facility use remains
prohibited.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
import os
from pathlib import Path
import time
from typing import Callable, Sequence

import numpy as np

from .fast_obstacle_transport import HYDROGEN_MOLAR_MASS_KG_MOL, R_GAS_J_MOL_K
from .fast_obstacle_transport_3d import (
    FastObstacleReceptor3D,
    _mass_consistent_face_projection_3d,
    _receptor_stencil,
    _source_weights,
    _source_weights_xyz,
    _source_weights_cosine_tapered_box,
    _source_weights_variance_matched_box,
    _transport_step,
)
from .hydrogen_air_thermodynamics import HydrogenAirActiveScalarClosure
from .binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
from .binary_mixture_corrector_3d import (
    solve_binary_mixture_initial_cloud_step_3d,
)
from .binary_mixture_handoff_plane_3d import (
    BinaryMixtureHandoffPlane3DFrame,
    BinaryMixtureHandoffPlane3DRequest,
    extract_binary_mixture_handoff_plane_3d,
)
from .finite_volume_initial_cloud_3d import (
    RegularPolygonPrismInitialCloud,
    initialise_regular_polygon_prism_cloud,
)
from .low_mach_thermal_corrector_3d import (
    _face_density,
    solve_low_mach_thermal_corrector_3d,
)
from .variable_density_momentum_3d import (
    advance_variable_density_momentum_3d,
    smooth_wall_log_law_force_density_3d,
)
from .thermal_wall_boundary_3d import isothermal_wall_heat_source_3d
from .semi_fv_obstacle import SourceRateSchedule
from .site_geometry import (
    AxisAlignedCuboid,
    CircularArcWall,
    CrosswindCylinder,
    SolidObstacle3D,
    circular_arc_surface_cell_metrics,
    circular_arc_surface_face_masks,
    circular_arc_surface_link_fractions,
    solid_obstacle_bounds_m,
    solid_obstacle_contains_point,
    solid_obstacle_voxel_mask,
)
from .subgrid_viscosity_3d import (
    subgrid_eddy_viscosity_physical,
    subgrid_relaxation_time_lattice,
)
from .atmospheric_surface_layer_mixing import NeutralSurfaceLayerMixing
from .dense_gas_stratification_3d import (
    diagnose_gradient_richardson_3d,
    stratification_damped_diffusivity_tensor_3d,
)


RESOLVED_OBSTACLE_FLOW_3D_SCHEMA = "degali.resolved-obstacle-flow-3d.v25"
RESOLVED_OBSTACLE_FLOW_3D_RESTART_SCHEMA = (
    "degali.resolved-obstacle-flow-3d-restart.v1"
)

_TURBULENT_DIFFUSIVITY_MODELS = {
    "legacy_total_kinematic_viscosity_ratio",
    "molecular_plus_subgrid_eddy_viscosity",
}
_CS2_3D = 1.0 / 3.0
_C3 = np.asarray((
    (0, 0, 0),
    (1, 0, 0), (-1, 0, 0),
    (0, 1, 0), (0, -1, 0),
    (0, 0, 1), (0, 0, -1),
    (1, 1, 0), (-1, 1, 0), (1, -1, 0), (-1, -1, 0),
    (1, 0, 1), (-1, 0, 1), (1, 0, -1), (-1, 0, -1),
    (0, 1, 1), (0, -1, 1), (0, 1, -1), (0, -1, -1),
), dtype=int)
_CX3 = _C3[:, 0]
_CY3 = _C3[:, 1]
_CZ3 = _C3[:, 2]
_W3 = np.asarray((
    1.0 / 3.0,
    *(1.0 / 18.0 for _ in range(6)),
    *(1.0 / 36.0 for _ in range(12)),
), dtype=float)
_OPPOSITE3 = np.asarray((
    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15,
), dtype=int)
_DIRECTION_INDEX = {tuple(direction): index for index, direction in enumerate(_C3)}


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _source_kernel_fwhm_cells(
    sigma_xyz_m: tuple[float, float, float], spacing_m: float,
) -> tuple[float, float, float]:
    """Return Gaussian FWHM in cells, a source-shape resolution diagnostic."""
    sigma = tuple(_finite(value, f"source sigma[{index}]") for index, value in enumerate(sigma_xyz_m))
    spacing = _finite(spacing_m, "grid spacing")
    if len(sigma) != 3 or any(value <= 0.0 for value in sigma) or spacing <= 0.0:
        raise ValueError("source sigmas and grid spacing must be positive")
    return tuple(2.0 * math.sqrt(2.0 * math.log(2.0)) * value / spacing for value in sigma)


def _effective_transport_diffusivity(
    *,
    molecular_diffusivity_m2_s: float,
    effective_dimensionless_number: float | None,
    base_kinematic_viscosity_m2_s: float,
    subgrid_eddy_viscosity_m2_s: np.ndarray,
    model: str,
) -> float | np.ndarray:
    """Return the declared molecular or explicitly selected turbulent closure.

    The legacy option is retained solely so old v14 results can be reproduced.
    It treats the configured base momentum viscosity as a uniform total
    viscosity.  The SGS option instead applies the gradient-diffusion closure
    to the local subgrid eddy viscosity and keeps molecular diffusion separate.
    """
    if effective_dimensionless_number is None:
        return float(molecular_diffusivity_m2_s)
    if model == "legacy_total_kinematic_viscosity_ratio":
        return float(
            base_kinematic_viscosity_m2_s / effective_dimensionless_number
        )
    if model == "molecular_plus_subgrid_eddy_viscosity":
        return (
            float(molecular_diffusivity_m2_s)
            + np.asarray(subgrid_eddy_viscosity_m2_s, dtype=float)
            / effective_dimensionless_number
        )
    raise ValueError(f"unsupported turbulent diffusivity model {model!r}")


def _apply_vertical_stratification_sensitivity(
    effective_scalar_diffusivity_m2_s: float | np.ndarray,
    effective_thermal_diffusivity_m2_s: float | np.ndarray,
    *,
    molecular_scalar_diffusivity_m2_s: float,
    molecular_thermal_diffusivity_m2_s: float,
    density_kg_m3: np.ndarray,
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    solid_mask: np.ndarray,
    grid_spacing_m: float,
    config: "ResolvedObstacleFlow3DConfig",
) -> tuple[
    float | np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray],
    float | np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray],
    float,
]:
    """Apply the opt-in, lagged-Ri_g anisotropic mixing sensitivity."""
    if config.stratification_damping_alpha is None:
        return (
            effective_scalar_diffusivity_m2_s,
            effective_thermal_diffusivity_m2_s,
            0.0,
        )
    assert config.stratification_shear_regularization_s_inv is not None
    diagnostic = diagnose_gradient_richardson_3d(
        density_kg_m3, u_m_s, v_m_s, solid_mask,
        grid_spacing_m=grid_spacing_m,
        shear_regularization_s_inv=(
            config.stratification_shear_regularization_s_inv
        ),
        gravitational_acceleration_m_s2=config.gravity_m_s2,
    )
    fluid = ~np.asarray(solid_mask, dtype=bool)

    def build_tensor(effective_value, molecular_value):
        molecular = float(molecular_value)
        effective = np.asarray(effective_value, dtype=float)
        if effective.ndim == 0:
            effective = np.full_like(density_kg_m3, float(effective), dtype=float)
        turbulent = np.maximum(effective - molecular, 0.0)
        return stratification_damped_diffusivity_tensor_3d(
            molecular,
            turbulent,
            diagnostic.gradient_richardson_number,
            fluid,
            alpha=config.stratification_damping_alpha,
            exponent=config.stratification_damping_exponent,
        )

    return (
        build_tensor(
            effective_scalar_diffusivity_m2_s,
            molecular_scalar_diffusivity_m2_s,
        ),
        build_tensor(
            effective_thermal_diffusivity_m2_s,
            molecular_thermal_diffusivity_m2_s,
        ),
        float(np.max(diagnostic.gradient_richardson_number)),
    )


def _diffusivity_statistics(
    diffusivity_m2_s: float | np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray],
    fluid: np.ndarray,
) -> tuple[float, float, float]:
    if isinstance(diffusivity_m2_s, tuple):
        components = [
            np.asarray(component, dtype=float)[fluid]
            for component in diffusivity_m2_s
        ]
        return (
            float(min(np.min(values) for values in components)),
            float(max(np.max(values) for values in components)),
            float(np.mean([np.mean(values) for values in components])),
        )
    if np.ndim(diffusivity_m2_s) == 0:
        value = float(diffusivity_m2_s)
        return value, value, value
    values = np.asarray(diffusivity_m2_s, dtype=float)[fluid]
    return float(np.min(values)), float(np.max(values)), float(np.mean(values))


@dataclass(frozen=True)
class ResolvedInletVelocityHistory3D:
    """Measured or declared vector wind on the solver event clock.

    Values are linearly interpolated and must cover the complete transient.
    The streamwise component stays strictly positive because the current
    x-min boundary is an inlet, not a bidirectional open boundary.
    """

    time_s: tuple[float, ...]
    u_m_s: tuple[float, ...]
    v_m_s: tuple[float, ...]
    w_m_s: tuple[float, ...]
    history_id: str
    interpolation: str = "linear"

    def __post_init__(self) -> None:
        if not isinstance(self.history_id, str) or not self.history_id.strip():
            raise ValueError("inlet velocity history_id must be non-empty")
        if self.interpolation not in {"linear", "piecewise_constant"}:
            raise ValueError(
                "inlet velocity interpolation must be linear or piecewise_constant"
            )
        values = {
            "time_s": self.time_s,
            "u_m_s": self.u_m_s,
            "v_m_s": self.v_m_s,
            "w_m_s": self.w_m_s,
        }
        if any(not isinstance(item, tuple) for item in values.values()):
            raise TypeError("inlet velocity history channels must be tuples")
        lengths = {len(item) for item in values.values()}
        if len(lengths) != 1 or next(iter(lengths)) < 2:
            raise ValueError(
                "inlet velocity history requires equal channels of at least two values"
            )
        converted = {
            name: tuple(_finite(value, name) for value in channel)
            for name, channel in values.items()
        }
        times = converted["time_s"]
        if not math.isclose(times[0], 0.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("inlet velocity history must start at event time zero")
        if any(right <= left for left, right in zip(times, times[1:])):
            raise ValueError("inlet velocity history time must be strictly increasing")
        if any(value <= 0.0 for value in converted["u_m_s"]):
            raise ValueError(
                "inlet velocity history streamwise velocity must remain positive"
            )
        for name, channel in converted.items():
            object.__setattr__(self, name, channel)

    def velocity_at(self, time_s: float) -> tuple[float, float, float]:
        """Linearly interpolate the three components at one event time."""

        event_time = _finite(time_s, "inlet velocity event time")
        if event_time < self.time_s[0] or event_time > self.time_s[-1]:
            raise ValueError("inlet velocity time lies outside declared history")
        if self.interpolation == "piecewise_constant":
            index = min(
                int(np.searchsorted(self.time_s, event_time, side="right") - 1),
                len(self.time_s) - 1,
            )
            return tuple(
                float(channel[index])
                for channel in (self.u_m_s, self.v_m_s, self.w_m_s)
            )
        return tuple(
            float(np.interp(event_time, self.time_s, channel))
            for channel in (self.u_m_s, self.v_m_s, self.w_m_s)
        )


@dataclass(frozen=True)
class ResolvedSourceMomentumHistory3D:
    """Time-resolved source-jet velocity on the declared event clock.

    Unlike the inlet history, a source jet may point in any direction and may
    pass through zero.  The history is a closure for unresolved nozzle
    momentum; it does not claim that the nozzle diameter or discharge
    coefficient has been resolved by the grid.
    """

    time_s: tuple[float, ...]
    u_m_s: tuple[float, ...]
    v_m_s: tuple[float, ...]
    w_m_s: tuple[float, ...]
    history_id: str
    interpolation: str = "linear"

    def __post_init__(self) -> None:
        if not isinstance(self.history_id, str) or not self.history_id.strip():
            raise ValueError("source momentum history_id must be non-empty")
        if self.interpolation not in {"linear", "piecewise_constant"}:
            raise ValueError(
                "source momentum interpolation must be linear or piecewise_constant"
            )
        values = {
            "time_s": self.time_s,
            "u_m_s": self.u_m_s,
            "v_m_s": self.v_m_s,
            "w_m_s": self.w_m_s,
        }
        if any(not isinstance(item, tuple) for item in values.values()):
            raise TypeError("source momentum history channels must be tuples")
        lengths = {len(item) for item in values.values()}
        if len(lengths) != 1 or next(iter(lengths)) < 2:
            raise ValueError(
                "source momentum history requires equal channels of at least two values"
            )
        converted = {
            name: tuple(_finite(value, name) for value in channel)
            for name, channel in values.items()
        }
        times = converted["time_s"]
        if not math.isclose(times[0], 0.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError("source momentum history must start at event time zero")
        if any(right <= left for left, right in zip(times, times[1:])):
            raise ValueError("source momentum history time must be strictly increasing")
        for name, channel in converted.items():
            object.__setattr__(self, name, channel)

    def velocity_at(self, time_s: float) -> tuple[float, float, float]:
        """Interpolate the source-jet velocity at one event time."""

        event_time = _finite(time_s, "source momentum event time")
        if event_time < self.time_s[0] or event_time > self.time_s[-1]:
            raise ValueError("source momentum time lies outside declared history")
        if self.interpolation == "piecewise_constant":
            index = min(
                int(np.searchsorted(self.time_s, event_time, side="right") - 1),
                len(self.time_s) - 1,
            )
            return tuple(
                float(channel[index])
                for channel in (self.u_m_s, self.v_m_s, self.w_m_s)
            )
        return tuple(
            float(np.interp(event_time, self.time_s, channel))
            for channel in (self.u_m_s, self.v_m_s, self.w_m_s)
        )


@dataclass(frozen=True)
class ResolvedObstacleFlow3DConfig:
    domain_m: tuple[float, float, float]
    cells: tuple[int, int, int]
    wind_speed_m_s: float
    kinematic_viscosity_m2_s: float
    flow_steps: int
    lattice_inlet_speed: float
    duration_s: float
    scalar_time_step_s: float
    scalar_diffusivity_m2_s: float
    source_schedule: SourceRateSchedule
    source_position_m: tuple[float, float, float]
    source_sigma_m: float
    obstacles: tuple[SolidObstacle3D, ...]
    velocity_reference_scale_m_s: float | None = None
    source_total_mass_to_species_mass_ratio: float = 1.0
    circular_arc_representation: str = "finite_thickness_voxel"
    source_momentum_velocity_m_s: tuple[float, float, float] = (0.0, 0.0, 0.0)
    source_momentum_velocity_history: ResolvedSourceMomentumHistory3D | None = None
    source_sigma_xyz_m: tuple[float, float, float] | None = None
    source_kernel: str = "gaussian"
    receptor_vertical_sampling: str = "cell_center"
    scalar_advection_scheme: str = "muscl_van_leer"
    source_time_integration: str = "prestep"
    subgrid_model: str = "smagorinsky"
    flow_boundary_scheme: str = "equilibrium_copy"
    nonstreamwise_wall_boundary: str = (
        "ground_no_slip_lateral_top_free_slip"
    )
    inlet_profile_z_m: tuple[float, ...] | None = None
    inlet_profile_u_m_s: tuple[float, ...] | None = None
    inlet_temperature_profile_z_m: tuple[float, ...] | None = None
    inlet_temperature_profile_k: tuple[float, ...] | None = None
    obstacle_wall_temperature_k: float | None = None
    ground_temperature_k: float | None = None
    active_inlet_profile_y_m: tuple[float, ...] | None = None
    active_inlet_profile_u_m_s: tuple[float, ...] | None = None
    active_inlet_velocity_history: ResolvedInletVelocityHistory3D | None = None
    smagorinsky_constant: float = 0.0
    wale_constant: float = 0.5
    outlet_lattice_density: float = 1.0
    maximum_relaxation_time: float = 1.95
    maximum_lattice_mach_number: float = 0.15
    minimum_flow_steps: int = 500
    convergence_check_interval: int = 100
    flow_convergence_tolerance: float = 1.0e-6
    flow_average_last_steps: int = 0
    flow_statistical_convergence_tolerance: float = 0.01
    ambient_temperature_k: float = 293.15
    ambient_pressure_pa: float = 101_325.0
    ambient_air_density_kg_m3: float = 1.204
    thermal_buoyancy_coupling: bool = False
    source_temperature_k: float | None = None
    thermal_diffusivity_m2_s: float = 2.2e-5
    effective_schmidt_number: float | None = None
    effective_prandtl_number: float | None = None
    scalar_turbulent_diffusivity_model: str = (
        "legacy_total_kinematic_viscosity_ratio"
    )
    thermal_turbulent_diffusivity_model: str = (
        "legacy_total_kinematic_viscosity_ratio"
    )
    dry_air_heat_capacity_j_kg_k: float = 1005.0
    hydrogen_heat_capacity_j_kg_k: float = 14_300.0
    gravity_m_s2: float = 9.80665
    maximum_active_scalar_density_departure: float = 0.10
    low_mach_mass_flux_correction: bool = True
    active_momentum_solver: str = "d3q19_bgk"
    carrier_flow_initialization: str = "d3q19_spinup"
    active_provisional_face_projection: str = "incompressible_then_variable_density"
    pressure_correction_momentum_coupling: str = "equilibrium_replacement"
    variable_density_projection_linear_solver: str = "auto"
    thermodynamic_corrector_maximum_iterations: int = 8
    thermodynamic_corrector_minimum_iterations: int = 2
    thermodynamic_corrector_density_tolerance: float = 1.0e-6
    thermodynamic_corrector_velocity_tolerance: float = 1.0e-6
    require_thermodynamic_corrector_convergence: bool = False
    store_scalar_fields: bool = False
    transient_average_start_s: float | None = None
    initial_binary_cloud: RegularPolygonPrismInitialCloud | None = None
    binary_mixture_closure: BinaryIdealGasMixtureClosure | None = None
    initial_binary_cloud_momentum: str = "ambient_carried"
    background_surface_layer_mixing: NeutralSurfaceLayerMixing | None = None
    stratification_damping_alpha: float | None = None
    stratification_damping_exponent: float = 1.0
    stratification_shear_regularization_s_inv: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.domain_m, tuple) or len(self.domain_m) != 3:
            raise TypeError("domain_m must be a three-value tuple")
        domain = tuple(
            _finite(value, f"domain_m[{index}]")
            for index, value in enumerate(self.domain_m)
        )
        if any(value <= 0.0 for value in domain):
            raise ValueError("domain dimensions must be positive")
        object.__setattr__(self, "domain_m", domain)
        if (
            not isinstance(self.cells, tuple) or len(self.cells) != 3
            or any(isinstance(value, bool) or not isinstance(value, int) or value < 8
                   for value in self.cells)
        ):
            raise ValueError("cells must contain three integers of at least eight")
        nx, ny, nz = self.cells
        dx, dy, dz = domain[0] / nx, domain[1] / ny, domain[2] / nz
        if not (
            math.isclose(dx, dy, rel_tol=1.0e-9, abs_tol=1.0e-12)
            and math.isclose(dx, dz, rel_tol=1.0e-9, abs_tol=1.0e-12)
        ):
            raise ValueError("D3Q19 flow grid requires cubic cells")
        for name in (
            "kinematic_viscosity_m2_s", "lattice_inlet_speed",
            "duration_s", "scalar_time_step_s", "source_sigma_m",
            "maximum_relaxation_time", "flow_convergence_tolerance",
            "flow_statistical_convergence_tolerance",
            "ambient_temperature_k", "ambient_pressure_pa",
            "ambient_air_density_kg_m3",
        ):
            value = _finite(getattr(self, name), name)
            object.__setattr__(self, name, value)
            if value <= 0.0:
                raise ValueError(f"{name} must be positive")
        wind_speed = _finite(self.wind_speed_m_s, "wind_speed_m_s")
        if wind_speed < 0.0:
            raise ValueError("wind_speed_m_s must be non-negative")
        object.__setattr__(self, "wind_speed_m_s", wind_speed)
        velocity_scale = self.velocity_reference_scale_m_s
        if velocity_scale is None:
            if wind_speed <= 0.0:
                raise ValueError(
                    "zero-wind runs require velocity_reference_scale_m_s"
                )
            velocity_scale = wind_speed
        velocity_scale = _finite(
            velocity_scale, "velocity_reference_scale_m_s",
        )
        if velocity_scale <= 0.0:
            raise ValueError("velocity_reference_scale_m_s must be positive")
        object.__setattr__(self, "velocity_reference_scale_m_s", velocity_scale)
        if self.lattice_inlet_speed >= 0.08:
            raise ValueError("lattice_inlet_speed must remain below 0.08")
        scalar_diffusivity = _finite(
            self.scalar_diffusivity_m2_s, "scalar_diffusivity_m2_s"
        )
        if scalar_diffusivity < 0.0:
            raise ValueError("scalar_diffusivity_m2_s must be non-negative")
        object.__setattr__(self, "scalar_diffusivity_m2_s", scalar_diffusivity)
        for name in (
            "effective_schmidt_number", "effective_prandtl_number",
        ):
            value = getattr(self, name)
            if value is not None:
                value = _finite(value, name)
                if value <= 0.0:
                    raise ValueError(f"{name} must be positive when supplied")
                object.__setattr__(self, name, value)
        exponent = _finite(
            self.stratification_damping_exponent,
            "stratification_damping_exponent",
        )
        if exponent <= 0.0:
            raise ValueError("stratification_damping_exponent must be positive")
        object.__setattr__(self, "stratification_damping_exponent", exponent)
        if self.stratification_damping_alpha is not None:
            alpha = _finite(
                self.stratification_damping_alpha,
                "stratification_damping_alpha",
            )
            regularization = self.stratification_shear_regularization_s_inv
            if regularization is None:
                raise ValueError(
                    "Ri_g damping requires an explicit shear regularization"
                )
            regularization = _finite(
                regularization, "stratification_shear_regularization_s_inv",
            )
            if alpha <= 0.0 or regularization <= 0.0:
                raise ValueError("Ri_g damping parameters must be positive")
            if self.effective_schmidt_number is None:
                raise ValueError(
                    "Ri_g damping requires turbulent species diffusivity"
                )
            if not (self.thermal_buoyancy_coupling or self.initial_binary_cloud is not None):
                raise ValueError(
                    "Ri_g damping requires a density-resolving thermodynamic mode"
                )
            object.__setattr__(self, "stratification_damping_alpha", alpha)
            object.__setattr__(
                self, "stratification_shear_regularization_s_inv",
                regularization,
            )
        elif self.stratification_shear_regularization_s_inv is not None:
            raise ValueError(
                "shear regularization is only valid when Ri_g damping is enabled"
            )
        for name in (
            "scalar_turbulent_diffusivity_model",
            "thermal_turbulent_diffusivity_model",
        ):
            value = getattr(self, name)
            if value not in _TURBULENT_DIFFUSIVITY_MODELS:
                allowed = ", ".join(sorted(_TURBULENT_DIFFUSIVITY_MODELS))
                raise ValueError(f"{name} must be one of: {allowed}")
        for name in (
            "thermal_diffusivity_m2_s", "dry_air_heat_capacity_j_kg_k",
            "hydrogen_heat_capacity_j_kg_k",
            "maximum_active_scalar_density_departure",
            "maximum_lattice_mach_number",
            "thermodynamic_corrector_density_tolerance",
            "thermodynamic_corrector_velocity_tolerance",
        ):
            value = _finite(getattr(self, name), name)
            if value <= 0.0:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, value)
        gravity = _finite(self.gravity_m_s2, "gravity_m_s2")
        if gravity < 0.0:
            raise ValueError("gravity_m_s2 must be non-negative")
        object.__setattr__(self, "gravity_m_s2", gravity)
        smagorinsky = _finite(self.smagorinsky_constant, "smagorinsky_constant")
        if smagorinsky < 0.0 or smagorinsky > 0.30:
            raise ValueError("smagorinsky_constant must lie in [0, 0.30]")
        object.__setattr__(self, "smagorinsky_constant", smagorinsky)
        wale = _finite(self.wale_constant, "wale_constant")
        if wale < 0.0 or wale > 1.0:
            raise ValueError("wale_constant must lie in [0, 1]")
        object.__setattr__(self, "wale_constant", wale)
        if self.subgrid_model not in {"none", "smagorinsky", "wale"}:
            raise ValueError(
                "subgrid_model must be 'none', 'smagorinsky', or 'wale'"
            )
        if self.flow_boundary_scheme not in {
            "equilibrium_copy", "nonequilibrium_extrapolation",
            "characteristic_nonreflecting",
            "characteristic_multidimensional",
            "characteristic_multidimensional_extrapolated",
            "characteristic_multidimensional_vortical_limited",
            "characteristic_nscbc_convective",
        }:
            raise ValueError(
                "flow_boundary_scheme must be 'equilibrium_copy' or "
                "'nonequilibrium_extrapolation' or "
                "'characteristic_nonreflecting' or "
                "'characteristic_multidimensional' or "
                "'characteristic_multidimensional_extrapolated' or "
                "'characteristic_multidimensional_vortical_limited' or "
                "'characteristic_nscbc_convective'"
            )
        if self.nonstreamwise_wall_boundary not in {
            "ground_no_slip_lateral_top_free_slip",
            "all_no_slip_tunnel",
            "all_symmetry_free_slip",
            "all_pressure_open",
            "ground_free_slip_lateral_top_pressure_open",
        }:
            raise ValueError(
                "nonstreamwise_wall_boundary must be "
                "'ground_no_slip_lateral_top_free_slip' or "
                "'all_no_slip_tunnel', 'all_symmetry_free_slip', or "
                "'all_pressure_open', or 'ground_free_slip_lateral_top_pressure_open'"
            )
        outlet_density = _finite(
            self.outlet_lattice_density, "outlet_lattice_density"
        )
        if outlet_density <= 0.0:
            raise ValueError("outlet_lattice_density must be positive")
        object.__setattr__(self, "outlet_lattice_density", outlet_density)
        if self.maximum_relaxation_time <= 0.5:
            raise ValueError("maximum_relaxation_time must exceed 0.5")
        if self.maximum_lattice_mach_number > 0.30:
            raise ValueError("maximum_lattice_mach_number cannot exceed 0.30")
        for name, minimum in (
            ("flow_steps", 1), ("minimum_flow_steps", 1),
            ("convergence_check_interval", 1),
            ("thermodynamic_corrector_maximum_iterations", 1),
            ("thermodynamic_corrector_minimum_iterations", 1),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer of at least {minimum}")
        if self.minimum_flow_steps > self.flow_steps:
            raise ValueError("minimum_flow_steps cannot exceed flow_steps")
        if (
            isinstance(self.flow_average_last_steps, bool)
            or not isinstance(self.flow_average_last_steps, int)
            or self.flow_average_last_steps < 0
        ):
            raise ValueError("flow_average_last_steps must be a non-negative integer")
        if self.flow_average_last_steps == 1:
            raise ValueError("flow_average_last_steps must be zero or at least two")
        if self.flow_average_last_steps >= self.flow_steps:
            raise ValueError("flow_average_last_steps must be less than flow_steps")
        if (
            self.carrier_flow_initialization == "projected_inlet"
            and self.flow_average_last_steps != 0
        ):
            raise ValueError(
                "projected_inlet carrier-flow initialization does not have a "
                "D3Q19 averaging window"
            )
        if (
            self.flow_average_last_steps > 0
            and self.minimum_flow_steps != self.flow_steps
        ):
            raise ValueError(
                "time-averaged flow requires minimum_flow_steps to equal flow_steps"
            )
        if (
            self.thermodynamic_corrector_minimum_iterations
            > self.thermodynamic_corrector_maximum_iterations
        ):
            raise ValueError(
                "thermodynamic corrector minimum iterations cannot exceed maximum"
            )
        if self.scalar_advection_scheme not in {
            "upwind", "muscl_van_leer", "muscl_ssprk2",
        }:
            raise ValueError(
                "scalar_advection_scheme must be 'upwind', 'muscl_van_leer', "
                "or 'muscl_ssprk2'"
            )
        if self.receptor_vertical_sampling not in {
            "cell_center", "ground_linear",
        }:
            raise ValueError(
                "receptor_vertical_sampling must be 'cell_center' or "
                "'ground_linear'"
            )
        if self.source_time_integration not in {"prestep", "strang"}:
            raise ValueError(
                "source_time_integration must be 'prestep' or 'strang'"
            )
        if self.active_provisional_face_projection not in {
            "incompressible_then_variable_density", "variable_density_only",
        }:
            raise ValueError(
                "active_provisional_face_projection must be "
                "'incompressible_then_variable_density' or "
                "'variable_density_only'"
            )
        if self.active_momentum_solver not in {
            "d3q19_bgk", "finite_volume_projection",
        }:
            raise ValueError(
                "active_momentum_solver must be 'd3q19_bgk' or "
                "'finite_volume_projection'"
            )
        if self.circular_arc_representation not in {
            "finite_thickness_voxel", "zero_thickness_surface",
        }:
            raise ValueError(
                "circular_arc_representation must be "
                "'finite_thickness_voxel' or 'zero_thickness_surface'"
            )
        arc_walls = tuple(
            item for item in self.obstacles if isinstance(item, CircularArcWall)
        )
        if self.circular_arc_representation == "finite_thickness_voxel" and any(
            item.thickness_m == 0.0 for item in arc_walls
        ):
            raise ValueError(
                "finite_thickness_voxel requires positive CircularArcWall "
                "thickness_m; use zero_thickness_surface for an ideal wall"
            )
        if self.circular_arc_representation == "zero_thickness_surface" and any(
            item.thickness_m != 0.0 for item in arc_walls
        ):
            raise ValueError(
                "zero_thickness_surface requires CircularArcWall thickness_m=0"
            )
        if self.carrier_flow_initialization not in {
            "d3q19_spinup", "projected_inlet",
        }:
            raise ValueError(
                "carrier_flow_initialization must be 'd3q19_spinup' or "
                "'projected_inlet'"
            )
        if (
            self.carrier_flow_initialization == "projected_inlet"
            and self.active_momentum_solver != "finite_volume_projection"
        ):
            raise ValueError(
                "projected_inlet carrier-flow initialization requires the "
                "finite_volume_projection active momentum solver"
            )
        if self.pressure_correction_momentum_coupling not in {
            "equilibrium_replacement", "deferred_guo_force", "transport_only",
        }:
            raise ValueError(
                "pressure_correction_momentum_coupling must be "
                "'equilibrium_replacement', 'deferred_guo_force', or "
                "'transport_only'"
            )
        if self.variable_density_projection_linear_solver not in {
            "auto", "cached_direct", "cg", "amg_cg",
        }:
            raise ValueError(
                "variable_density_projection_linear_solver must be 'auto', "
                "'cached_direct', 'cg', or 'amg_cg'"
            )
        if not isinstance(self.source_schedule, SourceRateSchedule):
            raise TypeError("source_schedule must be a SourceRateSchedule")
        self.source_schedule.require_zero_endpoint("resolved 3-D obstacle source schedule")
        if self.source_schedule.duration_s > self.duration_s + 1.0e-12:
            raise ValueError("source schedule cannot outlast the transport run")
        source_mass_ratio = _finite(
            self.source_total_mass_to_species_mass_ratio,
            "source_total_mass_to_species_mass_ratio",
        )
        if source_mass_ratio < 1.0:
            raise ValueError(
                "source total-to-species mass ratio must be at least one"
            )
        if source_mass_ratio > 1.0 and self.initial_binary_cloud is None:
            raise ValueError(
                "diluted mixture source requires the binary-mixture closure"
            )
        object.__setattr__(
            self, "source_total_mass_to_species_mass_ratio", source_mass_ratio,
        )
        if not isinstance(self.source_position_m, tuple) or len(self.source_position_m) != 3:
            raise TypeError("source_position_m must be a three-value tuple")
        source = tuple(
            _finite(value, f"source_position_m[{index}]")
            for index, value in enumerate(self.source_position_m)
        )
        if any(value < 0.0 or value > bound for value, bound in zip(source, domain)):
            raise ValueError("source_position_m must lie inside the domain")
        object.__setattr__(self, "source_position_m", source)
        if (
            not isinstance(self.source_momentum_velocity_m_s, tuple)
            or len(self.source_momentum_velocity_m_s) != 3
        ):
            raise TypeError(
                "source_momentum_velocity_m_s must be a three-value tuple"
            )
        source_velocity = tuple(
            _finite(value, f"source_momentum_velocity_m_s[{index}]")
            for index, value in enumerate(self.source_momentum_velocity_m_s)
        )
        object.__setattr__(
            self, "source_momentum_velocity_m_s", source_velocity,
        )
        if self.source_sigma_xyz_m is not None:
            if (
                not isinstance(self.source_sigma_xyz_m, tuple)
                or len(self.source_sigma_xyz_m) != 3
            ):
                raise TypeError("source_sigma_xyz_m must be a three-value tuple")
            source_sigmas = tuple(
                _finite(value, f"source_sigma_xyz_m[{index}]")
                for index, value in enumerate(self.source_sigma_xyz_m)
            )
            if any(value <= 0.0 for value in source_sigmas):
                raise ValueError("source_sigma_xyz_m values must be positive")
            object.__setattr__(self, "source_sigma_xyz_m", source_sigmas)
        if self.source_kernel not in {
            "gaussian", "variance_matched_box", "cosine_tapered_box",
        }:
            raise ValueError(
                "source_kernel must be 'gaussian', 'variance_matched_box' "
                "or 'cosine_tapered_box'"
            )
        if (
            any(abs(value) > 0.0 for value in source_velocity)
            and self.active_momentum_solver != "finite_volume_projection"
        ):
            raise ValueError(
                "non-zero source_momentum_velocity_m_s requires the "
                "finite_volume_projection active momentum solver"
            )
        if (
            self.source_momentum_velocity_history is not None
            and any(
                abs(value) > 0.0
                for channel in (
                    self.source_momentum_velocity_history.u_m_s,
                    self.source_momentum_velocity_history.v_m_s,
                    self.source_momentum_velocity_history.w_m_s,
                )
                for value in channel
            )
            and self.active_momentum_solver != "finite_volume_projection"
        ):
            raise ValueError(
                "non-zero source_momentum_velocity_history requires the "
                "finite_volume_projection active momentum solver"
            )
        if not isinstance(self.obstacles, tuple):
            raise TypeError("obstacles must be a tuple")
        if any(
            not isinstance(item, (AxisAlignedCuboid, CircularArcWall, CrosswindCylinder))
            for item in self.obstacles
        ):
            raise TypeError(
                "obstacles must contain AxisAlignedCuboid or CircularArcWall values"
            )
        for obstacle in self.obstacles:
            bounds = solid_obstacle_bounds_m(obstacle)
            if not (
                2.0 * dx < bounds[0] < bounds[1] < domain[0] - 3.0 * dx
                and 2.0 * dy < bounds[2] < bounds[3] < domain[1] - 2.0 * dy
                and 0.0 <= bounds[4] < bounds[5] < domain[2] - 2.0 * dz
            ):
                raise ValueError("obstacles must leave resolved inlet, outlet, side and top fetch")
            if solid_obstacle_contains_point(obstacle, source):
                raise ValueError("source may not lie inside an obstacle")
        if (self.inlet_profile_z_m is None) != (self.inlet_profile_u_m_s is None):
            raise ValueError(
                "inlet_profile_z_m and inlet_profile_u_m_s must be supplied together"
            )
        if self.inlet_profile_z_m is not None:
            if not isinstance(self.inlet_profile_z_m, tuple) or not isinstance(
                self.inlet_profile_u_m_s, tuple
            ):
                raise TypeError("inlet profile coordinates and velocities must be tuples")
            if len(self.inlet_profile_z_m) < 2 or len(self.inlet_profile_z_m) != len(
                self.inlet_profile_u_m_s
            ):
                raise ValueError("inlet profile requires equally sized tuples")
            profile_z = tuple(_finite(value, "inlet_profile_z_m") for value in self.inlet_profile_z_m)
            profile_u = tuple(_finite(value, "inlet_profile_u_m_s") for value in self.inlet_profile_u_m_s)
            if any(right <= left for left, right in zip(profile_z, profile_z[1:])):
                raise ValueError("inlet_profile_z_m must be strictly increasing")
            if profile_z[0] < 0.0 or any(value < 0.0 for value in profile_u):
                raise ValueError("inlet profile heights and velocities must be non-negative")
            object.__setattr__(self, "inlet_profile_z_m", profile_z)
            object.__setattr__(self, "inlet_profile_u_m_s", profile_u)
        if (self.inlet_temperature_profile_z_m is None) != (
            self.inlet_temperature_profile_k is None
        ):
            raise ValueError(
                "inlet_temperature_profile_z_m and "
                "inlet_temperature_profile_k must be supplied together"
            )
        if self.inlet_temperature_profile_z_m is not None:
            if not isinstance(self.inlet_temperature_profile_z_m, tuple) or not isinstance(
                self.inlet_temperature_profile_k, tuple
            ):
                raise TypeError("inlet temperature profile coordinates and values must be tuples")
            if len(self.inlet_temperature_profile_z_m) < 2 or len(
                self.inlet_temperature_profile_z_m
            ) != len(self.inlet_temperature_profile_k):
                raise ValueError("inlet temperature profile requires equally sized tuples")
            profile_z_t = tuple(_finite(value, "inlet_temperature_profile_z_m") for value in self.inlet_temperature_profile_z_m)
            profile_t = tuple(_finite(value, "inlet_temperature_profile_k") for value in self.inlet_temperature_profile_k)
            if any(right <= left for left, right in zip(profile_z_t, profile_z_t[1:])):
                raise ValueError("inlet_temperature_profile_z_m must be strictly increasing")
            if profile_z_t[0] < 0.0 or any(value <= 0.0 for value in profile_t):
                raise ValueError("inlet temperature profile requires non-negative heights and positive kelvin")
            object.__setattr__(self, "inlet_temperature_profile_z_m", profile_z_t)
            object.__setattr__(self, "inlet_temperature_profile_k", profile_t)
            if not self.thermal_buoyancy_coupling:
                raise ValueError("inlet temperature profile requires thermal_buoyancy_coupling")
        for name in ("obstacle_wall_temperature_k", "ground_temperature_k"):
            value = getattr(self, name)
            if value is not None:
                value = _finite(value, name)
                if value <= 0.0:
                    raise ValueError(f"{name} must be positive")
                object.__setattr__(self, name, value)
                if not self.thermal_buoyancy_coupling:
                    raise ValueError(f"{name} requires thermal_buoyancy_coupling")
        if (self.active_inlet_profile_y_m is None) != (
            self.active_inlet_profile_u_m_s is None
        ):
            raise ValueError(
                "active_inlet_profile_y_m and active_inlet_profile_u_m_s "
                "must be supplied together"
            )
        if self.active_inlet_profile_y_m is not None:
            if not isinstance(self.active_inlet_profile_y_m, tuple) or not isinstance(
                self.active_inlet_profile_u_m_s, tuple
            ):
                raise TypeError("active inlet profile coordinates and velocities must be tuples")
            if len(self.active_inlet_profile_y_m) < 2 or len(
                self.active_inlet_profile_y_m
            ) != len(self.active_inlet_profile_u_m_s):
                raise ValueError("active inlet profile requires equally sized tuples")
            active_y = tuple(
                _finite(value, "active_inlet_profile_y_m")
                for value in self.active_inlet_profile_y_m
            )
            active_u = tuple(
                _finite(value, "active_inlet_profile_u_m_s")
                for value in self.active_inlet_profile_u_m_s
            )
            if any(right <= left for left, right in zip(active_y, active_y[1:])):
                raise ValueError("active_inlet_profile_y_m must be strictly increasing")
            if any(value <= 0.0 for value in active_u):
                raise ValueError("active_inlet_profile_u_m_s must be positive")
            object.__setattr__(self, "active_inlet_profile_y_m", active_y)
            object.__setattr__(self, "active_inlet_profile_u_m_s", active_u)
        if self.active_inlet_velocity_history is not None:
            if not isinstance(
                self.active_inlet_velocity_history,
                ResolvedInletVelocityHistory3D,
            ):
                raise TypeError(
                    "active_inlet_velocity_history must be "
                    "ResolvedInletVelocityHistory3D"
                )
            if self.active_inlet_velocity_history.time_s[-1] < self.duration_s:
                raise ValueError(
                    "active inlet velocity history must cover duration_s"
                )
        if self.source_momentum_velocity_history is not None:
            if not isinstance(
                self.source_momentum_velocity_history,
                ResolvedSourceMomentumHistory3D,
            ):
                raise TypeError(
                    "source_momentum_velocity_history must be "
                    "ResolvedSourceMomentumHistory3D"
                )
            if self.source_momentum_velocity_history.time_s[-1] < self.duration_s:
                raise ValueError(
                    "source momentum velocity history must cover duration_s"
                )
        for name in (
            "thermal_buoyancy_coupling", "low_mach_mass_flux_correction",
            "require_thermodynamic_corrector_convergence", "store_scalar_fields",
        ):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be boolean")
        if self.source_temperature_k is not None:
            source_temperature = _finite(
                self.source_temperature_k, "source_temperature_k"
            )
            if source_temperature < 14.0 or source_temperature > 1_000.0:
                raise ValueError("source_temperature_k must lie in [14, 1000] K")
            object.__setattr__(self, "source_temperature_k", source_temperature)
        if self.transient_average_start_s is not None:
            average_start = _finite(
                self.transient_average_start_s, "transient_average_start_s"
            )
            if not 0.0 <= average_start < self.duration_s:
                raise ValueError(
                    "transient_average_start_s must lie in [0, duration_s)"
                )
            object.__setattr__(self, "transient_average_start_s", average_start)
        if (
            self.thermal_buoyancy_coupling
            and self.source_temperature_k is None
            and (
                self.initial_binary_cloud is None
                or self.source_schedule.released_mass_kg > 0.0
            )
        ):
            raise ValueError(
                "source_temperature_k is required for thermal buoyancy coupling"
            )
        if (self.initial_binary_cloud is None) != (
            self.binary_mixture_closure is None
        ):
            raise ValueError(
                "initial_binary_cloud and binary_mixture_closure must be supplied together"
            )
        if self.initial_binary_cloud_momentum not in {
            "ambient_carried", "quiescent",
        }:
            raise ValueError(
                "initial_binary_cloud_momentum must be 'ambient_carried' or "
                "'quiescent'"
            )
        if (
            self.initial_binary_cloud is None
            and self.initial_binary_cloud_momentum != "ambient_carried"
        ):
            raise ValueError(
                "non-default initial_binary_cloud_momentum requires an initial "
                "binary cloud"
            )
        if self.initial_binary_cloud is not None:
            if not isinstance(
                self.initial_binary_cloud, RegularPolygonPrismInitialCloud,
            ):
                raise TypeError(
                    "initial_binary_cloud must be RegularPolygonPrismInitialCloud"
                )
            if not isinstance(
                self.binary_mixture_closure, BinaryIdealGasMixtureClosure,
            ):
                raise TypeError(
                    "binary_mixture_closure must be BinaryIdealGasMixtureClosure"
                )
            if not (
                self.thermal_buoyancy_coupling
                and self.low_mach_mass_flux_correction
                and self.active_momentum_solver == "finite_volume_projection"
            ):
                raise ValueError(
                    "initial binary cloud requires thermal coupling, variable-density "
                    "mass-flux correction and finite_volume_projection momentum"
                )
            binary_closure = self.binary_mixture_closure
            assert binary_closure is not None
            for configured, declared, label in (
                (self.ambient_temperature_k, binary_closure.ambient_temperature_k,
                 "ambient temperature"),
                (self.ambient_pressure_pa, binary_closure.ambient_pressure_pa,
                 "ambient pressure"),
                (self.ambient_air_density_kg_m3,
                 binary_closure.ambient_density_kg_m3, "ambient density"),
                (self.gravity_m_s2, binary_closure.gravity_m_s2, "gravity"),
            ):
                if not math.isclose(
                    configured, declared, rel_tol=1.0e-6, abs_tol=1.0e-10,
                ):
                    raise ValueError(
                        f"config and binary-mixture {label} must match"
                    )
        if (
            self.require_thermodynamic_corrector_convergence
            and not self.low_mach_mass_flux_correction
        ):
            raise ValueError(
                "corrector convergence cannot be required when mass-flux correction is disabled"
            )
        if self.active_momentum_solver == "finite_volume_projection" and not (
            self.thermal_buoyancy_coupling
            and self.low_mach_mass_flux_correction
        ):
            raise ValueError(
                "finite_volume_projection active momentum requires thermal "
                "coupling and variable-density mass-flux correction"
            )
        if self.wind_speed_m_s == 0.0:
            if not (
                self.active_momentum_solver == "finite_volume_projection"
                and self.carrier_flow_initialization == "projected_inlet"
            ):
                raise ValueError(
                    "zero ambient wind requires projected_inlet initialization "
                    "and finite_volume_projection momentum"
                )
            if self.inlet_profile_u_m_s is not None and any(
                value != 0.0 for value in self.inlet_profile_u_m_s
            ):
                raise ValueError(
                    "zero ambient wind requires a zero physical inlet profile"
                )
            if self.active_inlet_profile_y_m is not None:
                raise ValueError(
                    "zero ambient wind cannot use a forced-flow active inlet profile"
                )
            if self.active_inlet_velocity_history is not None:
                raise ValueError(
                    "zero ambient wind cannot use a forced-flow inlet history"
                )
        if (
            self.active_inlet_profile_y_m is not None
            and self.active_momentum_solver != "finite_volume_projection"
        ):
            raise ValueError(
                "active y-dependent inlet profiles require the "
                "finite_volume_projection active momentum solver"
            )
        if (
            self.active_inlet_velocity_history is not None
            and self.active_momentum_solver != "finite_volume_projection"
        ):
            raise ValueError(
                "time-dependent vector inlet requires the "
                "finite_volume_projection active momentum solver"
            )
        if self.background_surface_layer_mixing is not None:
            if not isinstance(
                self.background_surface_layer_mixing,
                NeutralSurfaceLayerMixing,
            ):
                raise TypeError(
                    "background_surface_layer_mixing must be "
                    "NeutralSurfaceLayerMixing"
                )
            if self.active_momentum_solver != "finite_volume_projection":
                raise ValueError(
                    "background surface-layer mixing requires the "
                    "finite_volume_projection active momentum solver"
                )
            if self.subgrid_model != "none":
                raise ValueError(
                    "the v1 background surface-layer closure cannot be combined "
                    "with a second SGS viscosity"
                )
            if (
                self.effective_schmidt_number is None
                or self.effective_prandtl_number is None
            ):
                raise ValueError(
                    "background surface-layer mixing requires declared turbulent "
                    "Schmidt and Prandtl numbers"
                )
            if (
                self.scalar_turbulent_diffusivity_model
                != "molecular_plus_subgrid_eddy_viscosity"
                or self.thermal_turbulent_diffusivity_model
                != "molecular_plus_subgrid_eddy_viscosity"
            ):
                raise ValueError(
                    "background surface-layer mixing requires molecular-plus-"
                    "eddy-viscosity scalar and thermal diffusion"
                )


@dataclass(frozen=True)
class ResolvedObstacleFlow3DRestartState:
    """Conserved and momentum state required to continue a binary FV run.

    The state intentionally excludes derived thermodynamic fields: they are
    reconstructed from species mass concentration and sensible enthalpy by
    the declared binary-mixture closure.  This prevents a checkpoint from
    carrying mutually inconsistent primitive and conserved variables.
    """

    schema: str
    configuration_sha256: str
    time_s: float
    concentration_kg_m3: np.ndarray
    sensible_enthalpy_density_j_m3: np.ndarray
    maximum_concentration_kg_m3: np.ndarray
    maximum_species_mole_fraction: np.ndarray
    velocity_u_m_s: np.ndarray
    velocity_v_m_s: np.ndarray
    velocity_w_m_s: np.ndarray
    u_face_m_s: np.ndarray
    v_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    projection_correction_potential_kg_m_s: np.ndarray | None

    def __post_init__(self) -> None:
        if self.schema != RESOLVED_OBSTACLE_FLOW_3D_RESTART_SCHEMA:
            raise ValueError("unsupported resolved 3-D restart schema")
        if (
            not isinstance(self.configuration_sha256, str)
            or len(self.configuration_sha256) != 64
            or any(
                character not in "0123456789abcdef"
                for character in self.configuration_sha256
            )
        ):
            raise ValueError(
                "restart configuration_sha256 must be a lowercase SHA-256"
            )
        time_value = _finite(self.time_s, "restart time_s")
        if time_value <= 0.0:
            raise ValueError("restart time_s must be positive")
        object.__setattr__(self, "time_s", time_value)
        cell_fields = (
            "concentration_kg_m3",
            "sensible_enthalpy_density_j_m3",
            "maximum_concentration_kg_m3",
            "maximum_species_mole_fraction",
            "velocity_u_m_s",
            "velocity_v_m_s",
            "velocity_w_m_s",
        )
        shapes = set()
        for name in cell_fields:
            value = np.asarray(getattr(self, name), dtype=float)
            if value.ndim != 3 or np.any(~np.isfinite(value)):
                raise ValueError(f"restart {name} must be a finite 3-D array")
            shapes.add(value.shape)
            object.__setattr__(self, name, value.copy())
        if len(shapes) != 1:
            raise ValueError("restart cell fields must share one shape")
        nz, ny, nx = next(iter(shapes))
        expected_faces = {
            "u_face_m_s": (nz, ny, nx + 1),
            "v_face_m_s": (nz, ny + 1, nx),
            "w_face_m_s": (nz + 1, ny, nx),
        }
        for name, expected in expected_faces.items():
            value = np.asarray(getattr(self, name), dtype=float)
            if value.shape != expected or np.any(~np.isfinite(value)):
                raise ValueError(
                    f"restart {name} must be finite with shape {expected}"
                )
            object.__setattr__(self, name, value.copy())
        potential = self.projection_correction_potential_kg_m_s
        if potential is not None:
            value = np.asarray(potential, dtype=float)
            if value.shape != (nz, ny, nx) or np.any(~np.isfinite(value)):
                raise ValueError(
                    "restart projection correction potential must be finite "
                    f"with shape {(nz, ny, nx)}"
                )
            object.__setattr__(
                self, "projection_correction_potential_kg_m_s", value.copy(),
            )


def _restart_configuration_sha256(
    config: ResolvedObstacleFlow3DConfig,
) -> str:
    """Hash the complete declared configuration for strict restart matching."""

    payload = json.dumps(
        {
            "solver_schema": RESOLVED_OBSTACLE_FLOW_3D_SCHEMA,
            "configuration": asdict(config),
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def rebind_restart_for_corrector_tolerance_audit(
    state: ResolvedObstacleFlow3DRestartState,
    *,
    source_config: ResolvedObstacleFlow3DConfig,
    target_config: ResolvedObstacleFlow3DConfig,
) -> ResolvedObstacleFlow3DRestartState:
    """Rebind one state only across the two nonlinear tolerance fields.

    This narrow bridge supports a matched strict/candidate postcheck starting
    from one identical conserved state.  It is deliberately not a general
    restart-hash override: the source hash must be authentic and every other
    declared configuration field must compare exactly.
    """

    if not isinstance(state, ResolvedObstacleFlow3DRestartState):
        raise TypeError("state must be ResolvedObstacleFlow3DRestartState")
    source_hash = _restart_configuration_sha256(source_config)
    if state.configuration_sha256 != source_hash:
        raise ValueError("restart does not match the declared source configuration")
    source = asdict(source_config)
    target = asdict(target_config)
    changed = {
        name for name in source
        if source[name] != target[name]
    }
    allowed = {
        "thermodynamic_corrector_density_tolerance",
        "thermodynamic_corrector_velocity_tolerance",
    }
    if not changed or not changed.issubset(allowed):
        unexpected = sorted(changed - allowed)
        raise ValueError(
            "corrector-tolerance audit rebind requires only the two tolerance "
            f"fields to change; unexpected changes: {unexpected}"
        )
    if (
        target_config.thermodynamic_corrector_density_tolerance
        != target_config.thermodynamic_corrector_velocity_tolerance
    ):
        raise ValueError("target corrector tolerances must be identical")
    return replace(
        state, configuration_sha256=_restart_configuration_sha256(target_config),
    )


def write_resolved_obstacle_flow_3d_restart(
    path: str | Path,
    state: ResolvedObstacleFlow3DRestartState,
) -> Path:
    """Atomically persist one restart without pickle-dependent payloads."""

    if not isinstance(state, ResolvedObstacleFlow3DRestartState):
        raise TypeError("state must be ResolvedObstacleFlow3DRestartState")
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    potential = state.projection_correction_potential_kg_m_s
    try:
        with temporary.open("xb") as stream:
            np.savez_compressed(
                stream,
                schema=np.asarray(state.schema),
                configuration_sha256=np.asarray(state.configuration_sha256),
                time_s=np.asarray(state.time_s),
                concentration_kg_m3=state.concentration_kg_m3,
                sensible_enthalpy_density_j_m3=(
                    state.sensible_enthalpy_density_j_m3
                ),
                maximum_concentration_kg_m3=(
                    state.maximum_concentration_kg_m3
                ),
                maximum_species_mole_fraction=(
                    state.maximum_species_mole_fraction
                ),
                velocity_u_m_s=state.velocity_u_m_s,
                velocity_v_m_s=state.velocity_v_m_s,
                velocity_w_m_s=state.velocity_w_m_s,
                u_face_m_s=state.u_face_m_s,
                v_face_m_s=state.v_face_m_s,
                w_face_m_s=state.w_face_m_s,
                projection_potential_present=np.asarray(potential is not None),
                projection_correction_potential_kg_m_s=(
                    np.empty((0,), dtype=float)
                    if potential is None else potential
                ),
            )
        os.replace(temporary, destination)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise
    return destination


def read_resolved_obstacle_flow_3d_restart(
    path: str | Path,
) -> ResolvedObstacleFlow3DRestartState:
    """Read and validate a restart written by the atomic NPZ contract."""

    source = Path(path).resolve()
    with np.load(source, allow_pickle=False) as payload:
        present = bool(payload["projection_potential_present"].item())
        potential = (
            np.asarray(
                payload["projection_correction_potential_kg_m_s"],
                dtype=float,
            )
            if present else None
        )
        return ResolvedObstacleFlow3DRestartState(
            schema=str(payload["schema"].item()),
            configuration_sha256=str(
                payload["configuration_sha256"].item()
            ),
            time_s=float(payload["time_s"].item()),
            concentration_kg_m3=payload["concentration_kg_m3"],
            sensible_enthalpy_density_j_m3=(
                payload["sensible_enthalpy_density_j_m3"]
            ),
            maximum_concentration_kg_m3=(
                payload["maximum_concentration_kg_m3"]
            ),
            maximum_species_mole_fraction=(
                payload["maximum_species_mole_fraction"]
            ),
            velocity_u_m_s=payload["velocity_u_m_s"],
            velocity_v_m_s=payload["velocity_v_m_s"],
            velocity_w_m_s=payload["velocity_w_m_s"],
            u_face_m_s=payload["u_face_m_s"],
            v_face_m_s=payload["v_face_m_s"],
            w_face_m_s=payload["w_face_m_s"],
            projection_correction_potential_kg_m_s=potential,
        )


@dataclass(frozen=True)
class ResolvedObstacleFlow3DResult:
    schema: str
    x_m: np.ndarray
    y_m: np.ndarray
    z_m: np.ndarray
    solid_mask: np.ndarray
    velocity_u_m_s: np.ndarray
    velocity_v_m_s: np.ndarray
    velocity_w_m_s: np.ndarray
    u_face_m_s: np.ndarray
    v_face_m_s: np.ndarray
    w_face_m_s: np.ndarray
    pressure_pa_gauge: np.ndarray
    pressure_coefficient: np.ndarray
    eddy_viscosity_m2_s: np.ndarray
    time_s: np.ndarray
    concentration_kg_m3: np.ndarray
    maximum_concentration_kg_m3: np.ndarray
    concentration_fields_kg_m3: np.ndarray | None
    sensible_enthalpy_density_j_m3: np.ndarray | None
    temperature_k: np.ndarray | None
    hydrogen_mole_fraction: np.ndarray | None
    maximum_hydrogen_mole_fraction: np.ndarray | None
    species_name: str
    species_mole_fraction: np.ndarray | None
    maximum_species_mole_fraction: np.ndarray | None
    mixture_density_kg_m3: np.ndarray | None
    buoyancy_acceleration_m_s2: np.ndarray | None
    pressure_correction_pa: np.ndarray | None
    continuity_residual_kg_m3_s: np.ndarray | None
    temperature_fields_k: np.ndarray | None
    receptor_traces_kg_m3: dict[str, np.ndarray]
    receptor_traces_vol_percent: dict[str, np.ndarray]
    receptor_temperature_traces_k: dict[str, np.ndarray]
    injected_mass_kg: np.ndarray
    domain_mass_kg: np.ndarray
    outflow_mass_kg: np.ndarray
    mass_residual_kg: np.ndarray
    injected_sensible_enthalpy_j: np.ndarray | None
    domain_sensible_enthalpy_j: np.ndarray | None
    outflow_sensible_enthalpy_j: np.ndarray | None
    energy_residual_j: np.ndarray | None
    flow_residual_history: np.ndarray
    restart_state: ResolvedObstacleFlow3DRestartState | None
    diagnostics: dict[str, object]
    ambient_temperature_k: float
    ambient_pressure_pa: float
    time_average_start_s: float | None = None
    time_average_sample_count: int = 0
    time_averaged_velocity_u_m_s: np.ndarray | None = None
    time_averaged_velocity_v_m_s: np.ndarray | None = None
    time_averaged_velocity_w_m_s: np.ndarray | None = None
    time_averaged_concentration_kg_m3: np.ndarray | None = None
    time_averaged_temperature_k: np.ndarray | None = None

    @property
    def concentration_vol_percent(self) -> np.ndarray:
        mole_fraction = (
            self.species_mole_fraction
            if self.species_mole_fraction is not None
            else self.hydrogen_mole_fraction
        )
        if mole_fraction is not None:
            return np.clip(mole_fraction * 100.0, 0.0, 100.0)
        molar_density = self.ambient_pressure_pa / (
            R_GAS_J_MOL_K * self.ambient_temperature_k
        )
        return np.clip(
            self.concentration_kg_m3 / HYDROGEN_MOLAR_MASS_KG_MOL
            / molar_density * 100.0,
            0.0, 100.0,
        )

    @property
    def maximum_concentration_vol_percent(self) -> np.ndarray:
        mole_fraction = (
            self.maximum_species_mole_fraction
            if self.maximum_species_mole_fraction is not None
            else self.maximum_hydrogen_mole_fraction
        )
        if mole_fraction is not None:
            return np.clip(
                mole_fraction * 100.0, 0.0, 100.0,
            )
        molar_density = self.ambient_pressure_pa / (
            R_GAS_J_MOL_K * self.ambient_temperature_k
        )
        return np.clip(
            self.maximum_concentration_kg_m3 / HYDROGEN_MOLAR_MASS_KG_MOL
            / molar_density * 100.0,
            0.0, 100.0,
        )

    @property
    def maximum_mass_residual_kg(self) -> float:
        return float(np.max(np.abs(self.mass_residual_kg)))

    @property
    def maximum_energy_residual_j(self) -> float | None:
        if self.energy_residual_j is None:
            return None
        return float(np.max(np.abs(self.energy_residual_j)))


def _equilibrium_d3q19(
    rho: np.ndarray, ux: np.ndarray, uy: np.ndarray, uz: np.ndarray,
) -> np.ndarray:
    cu = (
        _CX3[:, None, None, None] * ux[None]
        + _CY3[:, None, None, None] * uy[None]
        + _CZ3[:, None, None, None] * uz[None]
    )
    speed2 = ux**2 + uy**2 + uz**2
    return _W3[:, None, None, None] * rho[None] * (
        1.0 + 3.0 * cu + 4.5 * cu**2 - 1.5 * speed2[None]
    )


def _moments_d3q19(
    distribution: np.ndarray,
    solid: np.ndarray,
    force_lattice: tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray] = (
        0.0, 0.0, 0.0,
    ),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rho = np.maximum(np.sum(distribution, axis=0), 1.0e-14)
    acceleration = tuple(
        np.broadcast_to(np.asarray(component, dtype=float), rho.shape)
        for component in force_lattice
    )
    ux = np.tensordot(_CX3, distribution, axes=(0, 0)) / rho + 0.5 * acceleration[0]
    uy = np.tensordot(_CY3, distribution, axes=(0, 0)) / rho + 0.5 * acceleration[1]
    uz = np.tensordot(_CZ3, distribution, axes=(0, 0)) / rho + 0.5 * acceleration[2]
    ux[solid] = uy[solid] = uz[solid] = 0.0
    return rho, ux, uy, uz


def _guo_force_d3q19(
    rho: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    acceleration_lattice: tuple[
        float | np.ndarray, float | np.ndarray, float | np.ndarray
    ],
    relaxation_time: np.ndarray | float,
) -> np.ndarray:
    acceleration = np.stack([
        np.broadcast_to(np.asarray(component, dtype=float), rho.shape)
        for component in acceleration_lattice
    ])
    force = rho[None] * acceleration
    result = np.empty((19, *rho.shape), dtype=float)
    prefactor = 1.0 - 0.5 / relaxation_time
    for direction, (cx, cy, cz) in enumerate(_C3):
        cu = cx * ux + cy * uy + cz * uz
        projection = (
            (3.0 * (cx - ux) + 9.0 * cu * cx) * force[0]
            + (3.0 * (cy - uy) + 9.0 * cu * cy) * force[1]
            + (3.0 * (cz - uz) + 9.0 * cu * cz) * force[2]
        )
        result[direction] = _W3[direction] * prefactor * projection
    return result


def _solid_mask(config: ResolvedObstacleFlow3DConfig) -> tuple[np.ndarray, ...]:
    nx, ny, nz = config.cells
    dx = config.domain_m[0] / nx
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dx
    z = (np.arange(nz) + 0.5) * dx
    solid = np.zeros((nz, ny, nx), dtype=bool)
    for obstacle in config.obstacles:
        if (
            config.circular_arc_representation == "zero_thickness_surface"
            and isinstance(obstacle, CircularArcWall)
        ):
            continue
        solid |= solid_obstacle_voxel_mask(obstacle, x, y, z)
    return x, y, z, solid


def _circular_arc_blocked_faces(
    config: ResolvedObstacleFlow3DConfig,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build exact mid-surface crossings for the FV face-flux operators."""
    blocked_x = np.zeros((z.size, y.size, x.size + 1), dtype=bool)
    blocked_y = np.zeros((z.size, y.size + 1, x.size), dtype=bool)
    blocked_z = np.zeros((z.size + 1, y.size, x.size), dtype=bool)
    for obstacle in config.obstacles:
        if not isinstance(obstacle, CircularArcWall):
            continue
        faces = circular_arc_surface_face_masks(
            x, y, z,
            center_xy_m=obstacle.center_m,
            radius_m=obstacle.radius_m,
            z_min_m=obstacle.z_min_m,
            z_max_m=obstacle.z_max_m,
            arc_start_deg=obstacle.arc_start_deg,
            arc_end_deg=obstacle.arc_end_deg,
        )
        blocked_x |= faces[0]
        blocked_y |= faces[1]
    return blocked_x, blocked_y, blocked_z


def _circular_arc_blocked_links(
    config: ResolvedObstacleFlow3DConfig,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> np.ndarray | None:
    """Return first exact D3Q19 arc-hit fraction in surface mode, if any."""
    if config.circular_arc_representation != "zero_thickness_surface":
        return None
    fractions = np.full((len(_C3), z.size, y.size, x.size), np.nan)
    found_arc = False
    for obstacle in config.obstacles:
        if not isinstance(obstacle, CircularArcWall):
            continue
        found_arc = True
        wall_fractions = circular_arc_surface_link_fractions(
            x, y, z, _C3[1:],
            center_xy_m=obstacle.center_m,
            radius_m=obstacle.radius_m,
            z_min_m=obstacle.z_min_m,
            z_max_m=obstacle.z_max_m,
            arc_start_deg=obstacle.arc_start_deg,
            arc_end_deg=obstacle.arc_end_deg,
        )
        fractions[1:] = np.fmin(fractions[1:], wall_fractions)
    return fractions if found_arc else None


def _obstacle_voxelization_diagnostics(
    config: ResolvedObstacleFlow3DConfig,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> tuple[list[dict[str, object]], float, float]:
    """Return declared/effective obstacle bounds and voxelisation error."""
    spacing = config.domain_m[0] / config.cells[0]
    records: list[dict[str, object]] = []
    maximum_face_error = 0.0
    maximum_size_error = 0.0
    for obstacle in config.obstacles:
        if (
            config.circular_arc_representation == "zero_thickness_surface"
            and isinstance(obstacle, CircularArcWall)
        ):
            records.append({
                "label": obstacle.label,
                "geometry_type": type(obstacle).__name__,
                "geometry_representation": "zero_thickness_surface",
                "curved_surface_voxelization": False,
                "center_m": list(obstacle.center_m),
                "radius_m": obstacle.radius_m,
                "declared_radial_thickness_m": obstacle.thickness_m,
                "representation_note": "ideal surface; no voxel solid volume",
                "solid_cell_count": 0,
                "surface_link_and_face_masks_are_fractional": False,
            })
            continue
        obstacle_mask = solid_obstacle_voxel_mask(obstacle, x, y, z)
        selected = (
            np.flatnonzero(np.any(obstacle_mask, axis=(0, 1))),
            np.flatnonzero(np.any(obstacle_mask, axis=(0, 2))),
            np.flatnonzero(np.any(obstacle_mask, axis=(1, 2))),
        )
        if any(indices.size == 0 for indices in selected):
            raise ValueError(
                f"obstacle {obstacle.label!r} has no voxelized solid cells"
            )
        effective = tuple(
            value
            for coordinates, indices in zip((x, y, z), selected)
            for value in (
                float(coordinates[indices[0]] - 0.5 * spacing),
                float(coordinates[indices[-1]] + 0.5 * spacing),
            )
        )
        declared = solid_obstacle_bounds_m(obstacle)
        face_errors = tuple(
            effective_value - declared_value
            for effective_value, declared_value in zip(effective, declared)
        )
        declared_sizes = (
            declared[1] - declared[0],
            declared[3] - declared[2],
            declared[5] - declared[4],
        )
        effective_sizes = (
            effective[1] - effective[0],
            effective[3] - effective[2],
            effective[5] - effective[4],
        )
        size_errors = tuple(
            effective_value - declared_value
            for effective_value, declared_value in zip(
                effective_sizes, declared_sizes
            )
        )
        maximum_face_error = max(
            maximum_face_error, max(abs(value) for value in face_errors)
        )
        if isinstance(obstacle, CircularArcWall):
            # Axis-aligned bounds can coincide with cell faces while the
            # curved surface remains stair-stepped.  Never label that case as
            # exactly aligned; half a cell is the declared geometric floor.
            maximum_face_error = max(maximum_face_error, 0.5 * spacing)
        maximum_size_error = max(
            maximum_size_error, max(abs(value) for value in size_errors)
        )
        records.append({
            "label": obstacle.label,
            "geometry_type": type(obstacle).__name__,
            "curved_surface_voxelization": isinstance(obstacle, CircularArcWall),
            "declared_bounds_m": list(declared),
            "effective_voxel_bounds_m": list(effective),
            "signed_face_errors_m": list(face_errors),
            "signed_size_errors_m": list(size_errors),
            "solid_cell_counts_xyz": [int(indices.size) for indices in selected],
            "solid_cell_count": int(np.count_nonzero(obstacle_mask)),
        })
    return records, maximum_face_error, maximum_size_error


def _inlet_profile(
    config: ResolvedObstacleFlow3DConfig, z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if config.inlet_profile_z_m is None:
        physical = np.full_like(z, config.wind_speed_m_s)
    else:
        assert config.inlet_profile_u_m_s is not None
        physical = np.interp(
            z,
            np.asarray(config.inlet_profile_z_m),
            np.asarray(config.inlet_profile_u_m_s),
        )
    assert config.velocity_reference_scale_m_s is not None
    lattice = (
        physical * config.lattice_inlet_speed
        / config.velocity_reference_scale_m_s
    )
    if float(np.max(lattice)) >= 0.10:
        raise ValueError("mapped inlet profile exceeds lattice speed 0.10")
    return physical, lattice


def _active_inlet_plane(
    config: ResolvedObstacleFlow3DConfig,
    y: np.ndarray,
    z_profile_m_s: np.ndarray,
) -> np.ndarray:
    """Build the finite-volume z-y inlet plane from declared mean profiles."""
    plane = np.broadcast_to(
        np.asarray(z_profile_m_s, dtype=float)[:, None],
        (z_profile_m_s.size, y.size),
    ).copy()
    if config.active_inlet_profile_y_m is None:
        return plane
    assert config.active_inlet_profile_u_m_s is not None
    y_profile = np.interp(
        y,
        np.asarray(config.active_inlet_profile_y_m, dtype=float),
        np.asarray(config.active_inlet_profile_u_m_s, dtype=float),
    )
    assert config.velocity_reference_scale_m_s is not None
    plane *= y_profile[None, :] / config.velocity_reference_scale_m_s
    return plane


def _active_inlet_velocity_planes(
    config: ResolvedObstacleFlow3DConfig,
    reference_u_plane_m_s: np.ndarray,
    event_time_s: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map a vector history onto the declared spatial inlet shape.

    The history is a point/mean record, so one normalized spatial shape is
    applied to all three components.  This preserves the declared wall/height
    profile but is not a synthetic turbulence inflow generator.
    """

    reference = np.asarray(reference_u_plane_m_s, dtype=float)
    if config.active_inlet_velocity_history is None:
        return reference, np.zeros_like(reference), np.zeros_like(reference)
    u_value, v_value, w_value = (
        config.active_inlet_velocity_history.velocity_at(event_time_s)
    )
    assert config.velocity_reference_scale_m_s is not None
    shape = reference / config.velocity_reference_scale_m_s
    return shape * u_value, shape * v_value, shape * w_value


def _coerce_inlet_plane(
    inlet_physical: np.ndarray, nz: int, ny: int,
) -> np.ndarray:
    inlet = np.asarray(inlet_physical, dtype=float)
    if inlet.shape == (nz,):
        return np.broadcast_to(inlet[:, None], (nz, ny))
    if inlet.shape == (nz, ny):
        return inlet
    raise ValueError(
        f"inlet_physical must have shape {(nz,)} or {(nz, ny)}"
    )


def _multidimensional_characteristic_state(
    density: np.ndarray,
    outward_velocity: np.ndarray,
    tangential_velocity_a: np.ndarray,
    tangential_velocity_b: np.ndarray,
    *,
    reference_density: float,
    reference_outward_velocity: np.ndarray | float,
    reference_tangential_velocity_a: np.ndarray | float = 0.0,
    reference_tangential_velocity_b: np.ndarray | float = 0.0,
    acoustic_fraction: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an angle-aware isothermal acoustic outlet state.

    The one-dimensional characteristic assumes normal acoustic impedance,
    ``u_n' = c rho'/rho0``.  For a plane wave incident at angle theta from
    the normal, the exact linear relation is instead
    ``u_n' = c cos(theta) rho'/rho0``.  The local velocity perturbation gives
    ``cos(theta)`` without a fitted angle.  The outgoing invariant is retained
    and the incoming acoustic amplitude is set to zero under that impedance.

    An optional local ``acoustic_fraction`` blends this acoustic target with
    the extrapolated state.  The production vortical candidate derives that
    fraction from local dilatation/pressure-gradient versus curl magnitudes;
    it is not an empirical coefficient.
    """
    rho = np.maximum(np.asarray(density, dtype=float), 1.0e-14)
    normal = np.asarray(outward_velocity, dtype=float)
    tangent_a = np.asarray(tangential_velocity_a, dtype=float)
    tangent_b = np.asarray(tangential_velocity_b, dtype=float)
    reference_normal = np.broadcast_to(
        np.asarray(reference_outward_velocity, dtype=float), normal.shape,
    )
    reference_a = np.broadcast_to(
        np.asarray(reference_tangential_velocity_a, dtype=float), normal.shape,
    )
    reference_b = np.broadcast_to(
        np.asarray(reference_tangential_velocity_b, dtype=float), normal.shape,
    )
    sound_speed = math.sqrt(_CS2_3D)
    normal_perturbation = normal - reference_normal
    tangent_a_perturbation = tangent_a - reference_a
    tangent_b_perturbation = tangent_b - reference_b
    velocity_amplitude = np.sqrt(
        normal_perturbation**2
        + tangent_a_perturbation**2
        + tangent_b_perturbation**2
    )
    acoustic_log_density = np.log(rho / reference_density)
    signal_scale = np.maximum(
        velocity_amplitude, sound_speed * np.abs(acoustic_log_density),
    )
    resolved_signal = signal_scale > 1.0e-12
    incidence_cosine = np.ones_like(normal)
    np.divide(
        np.abs(normal_perturbation), velocity_amplitude,
        out=incidence_cosine,
        where=resolved_signal & (velocity_amplitude > 1.0e-14),
    )
    incidence_cosine = np.clip(incidence_cosine, 0.05, 1.0)
    outgoing_perturbation = (
        normal_perturbation + sound_speed * acoustic_log_density
    )
    target_log_density = outgoing_perturbation / (
        sound_speed * (1.0 + incidence_cosine)
    )
    target_normal_perturbation = (
        incidence_cosine * outgoing_perturbation
        / (1.0 + incidence_cosine)
    )
    target_outward_velocity = (
        reference_normal + target_normal_perturbation
    )
    if acoustic_fraction is not None:
        acoustic_weight = np.asarray(acoustic_fraction, dtype=float)
        if acoustic_weight.shape != normal.shape:
            acoustic_weight = np.broadcast_to(
                acoustic_weight, normal.shape,
            )
        if np.any(~np.isfinite(acoustic_weight)):
            raise ValueError("acoustic_fraction must be finite")
        if np.any((acoustic_weight < 0.0) | (acoustic_weight > 1.0)):
            raise ValueError("acoustic_fraction must be between zero and one")
        target_log_density = (
            acoustic_log_density
            + acoustic_weight * (
                target_log_density - acoustic_log_density
            )
        )
        target_outward_velocity = (
            normal
            + acoustic_weight * (
                target_outward_velocity - normal
            )
        )
    target_density = reference_density * np.exp(target_log_density)
    return target_density, target_outward_velocity, incidence_cosine


def _limited_boundary_extrapolation(
    first_interior: np.ndarray,
    second_interior: np.ndarray,
    third_interior: np.ndarray,
) -> np.ndarray:
    """TVD-extrapolate one plane without creating a new component extremum.

    A raw ``2*f1-f2`` boundary reconstruction is second order on a linear
    profile but can amplify a finite vortical disturbance.  The two adjacent
    one-sided slopes provide a coefficient-free minmod limiter: smooth
    monotone profiles retain second-order extrapolation, while a slope reversal
    falls back to the nearest interior state.
    """
    first = np.asarray(first_interior, dtype=float)
    second = np.asarray(second_interior, dtype=float)
    third = np.asarray(third_interior, dtype=float)
    if first.shape != second.shape or first.shape != third.shape:
        raise ValueError("boundary extrapolation planes must have equal shapes")
    near_slope = first - second
    inner_slope = second - third
    same_sign = near_slope * inner_slope > 0.0
    limited_slope = np.where(
        same_sign,
        np.sign(near_slope) * np.minimum(
            np.abs(near_slope), np.abs(inner_slope),
        ),
        0.0,
    )
    return first + limited_slope


def _boundary_acoustic_fraction(
    first_distribution: np.ndarray,
    second_distribution: np.ndarray,
    *,
    normal_spatial_axis: int,
    normal_coordinate_sign: float,
    energy_weighted: bool = False,
    include_pressure_gradient: bool = True,
) -> np.ndarray:
    """Separate locally irrotational/acoustic and vortical boundary content.

    The dimensionless blend is ``A / (A + |curl u|)``. ``A`` is velocity
    dilatation and can optionally include the isothermal pressure-gradient
    scale ``c_s |grad ln rho|``. Pure acoustic modes are curl-free and receive
    the characteristic correction; incompressible vortical modes are
    extrapolated. The NSCBC candidate disables the pressure term so a
    centrifugal vortex-pressure gradient is not misclassified as sound.
    """
    first = np.asarray(first_distribution, dtype=float)
    second = np.asarray(second_distribution, dtype=float)
    if first.shape != second.shape or first.ndim != 4:
        raise ValueError("boundary distribution planes must be equal 4-D arrays")
    if normal_spatial_axis not in {0, 1, 2}:
        raise ValueError("normal_spatial_axis must be zero, one or two")
    if normal_coordinate_sign not in {-1.0, 1.0}:
        raise ValueError("normal_coordinate_sign must be -1 or 1")

    def state(distribution: np.ndarray) -> tuple[np.ndarray, ...]:
        density = np.maximum(np.sum(distribution, axis=0), 1.0e-14)
        ux = np.tensordot(_CX3, distribution, axes=(0, 0)) / density
        uy = np.tensordot(_CY3, distribution, axes=(0, 0)) / density
        uz = np.tensordot(_CZ3, distribution, axes=(0, 0)) / density
        return density, ux, uy, uz

    first_state = state(first)
    second_state = state(second)

    def derivative(
        first_field: np.ndarray,
        second_field: np.ndarray,
        axis: int,
    ) -> np.ndarray:
        if axis == normal_spatial_axis:
            return normal_coordinate_sign * (first_field - second_field)
        if first_field.shape[axis] < 2:
            return np.zeros_like(first_field)
        return np.gradient(first_field, axis=axis)

    log_density_first = np.log(first_state[0])
    log_density_second = np.log(second_state[0])
    density_gradient = tuple(
        derivative(log_density_first, log_density_second, axis)
        for axis in range(3)
    )
    velocity_gradient = tuple(
        tuple(
            derivative(first_state[component], second_state[component], axis)
            for axis in range(3)
        )
        for component in range(1, 4)
    )
    du, dv, dw = velocity_gradient
    divergence = du[2] + dv[1] + dw[0]
    curl_x = dw[1] - dv[0]
    curl_y = du[0] - dw[2]
    curl_z = dv[2] - du[1]
    vorticity = np.sqrt(curl_x**2 + curl_y**2 + curl_z**2)
    acoustic = np.abs(divergence)
    if include_pressure_gradient:
        acoustic = np.sqrt(
            acoustic**2
            + _CS2_3D * sum(item**2 for item in density_gradient)
        )
    if energy_weighted:
        acoustic = acoustic**2
        vorticity = vorticity**2
    denominator = acoustic + vorticity
    fraction = np.zeros_like(acoustic)
    threshold = 1.0e-28 if energy_weighted else 1.0e-14
    np.divide(
        acoustic, denominator, out=fraction, where=denominator > threshold,
    )
    return np.clip(fraction, 0.0, 1.0)


def _convective_boundary_distribution(
    previous_boundary: np.ndarray,
    previous_interior: np.ndarray,
    *,
    normal_velocity_component: int,
    outward_sign: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Advance a boundary plane with a local one-way convective equation.

    In lattice units ``dt=dx=1``.  The local outward velocity is therefore the
    upwind Courant number. Backflow is not silently extrapolated from inside;
    its convective Courant number is zero and remains for the acoustic/inflow
    treatment to handle.
    """
    boundary = np.asarray(previous_boundary, dtype=float)
    interior = np.asarray(previous_interior, dtype=float)
    if boundary.shape != interior.shape or boundary.ndim != 4:
        raise ValueError("convective boundary planes must be equal 4-D arrays")
    if normal_velocity_component not in {0, 1, 2}:
        raise ValueError("normal_velocity_component must be zero, one or two")
    if outward_sign not in {-1.0, 1.0}:
        raise ValueError("outward_sign must be -1 or 1")
    density = np.maximum(np.sum(boundary, axis=0), 1.0e-14)
    directions = (_CX3, _CY3, _CZ3)
    normal_velocity = (
        np.tensordot(
            directions[normal_velocity_component], boundary, axes=(0, 0),
        ) / density
    )
    courant = np.clip(outward_sign * normal_velocity, 0.0, 1.0)
    advanced = boundary - courant[None] * (boundary - interior)
    return advanced, courant


def _nscbc_xmax_boundary_distribution(
    previous_boundary: np.ndarray,
    previous_interior: np.ndarray,
    *,
    reference_density: float,
    reference_velocity_x: np.ndarray | float,
) -> np.ndarray:
    """Advance an x-max plane with multidimensional isothermal NSCBC.

    ``R+ = u + c ln(rho)`` leaves through the subsonic x-max boundary and is
    upwinded from the interior. ``R-`` receives no normal incoming derivative,
    but both invariants retain their exact transverse Euler source terms. This
    allows pressure-balanced vortical motion to cross without being projected
    pointwise onto an acoustic impedance relation.
    """
    boundary = np.asarray(previous_boundary, dtype=float)
    interior = np.asarray(previous_interior, dtype=float)
    if boundary.shape != interior.shape or boundary.ndim != 4:
        raise ValueError("NSCBC boundary planes must be equal 4-D arrays")
    if boundary.shape[3] != 1:
        raise ValueError("NSCBC x-max planes must have singleton x dimension")
    reference_density = float(reference_density)
    if not math.isfinite(reference_density) or reference_density <= 0.0:
        raise ValueError("reference_density must be positive and finite")

    def state(distribution: np.ndarray) -> tuple[np.ndarray, ...]:
        density = np.maximum(np.sum(distribution, axis=0), 1.0e-14)
        ux = np.tensordot(_CX3, distribution, axes=(0, 0)) / density
        uy = np.tensordot(_CY3, distribution, axes=(0, 0)) / density
        uz = np.tensordot(_CZ3, distribution, axes=(0, 0)) / density
        return density, ux, uy, uz

    rho, ux, uy, uz = state(boundary)
    rho_inner, ux_inner, uy_inner, uz_inner = state(interior)
    sound_speed = math.sqrt(_CS2_3D)
    log_rho = np.log(rho / reference_density)
    log_rho_inner = np.log(rho_inner / reference_density)
    outgoing = ux + sound_speed * log_rho
    incoming = ux - sound_speed * log_rho
    outgoing_inner = ux_inner + sound_speed * log_rho_inner

    def tangent_gradient(field: np.ndarray, axis: int) -> np.ndarray:
        if field.shape[axis] < 2:
            return np.zeros_like(field)
        return np.gradient(field, axis=axis)

    outgoing_y = tangent_gradient(outgoing, 1)
    outgoing_z = tangent_gradient(outgoing, 0)
    incoming_y = tangent_gradient(incoming, 1)
    incoming_z = tangent_gradient(incoming, 0)
    uy_y = tangent_gradient(uy, 1)
    uy_z = tangent_gradient(uy, 0)
    uz_y = tangent_gradient(uz, 1)
    uz_z = tangent_gradient(uz, 0)
    log_rho_y = tangent_gradient(log_rho, 1)
    log_rho_z = tangent_gradient(log_rho, 0)
    transverse_divergence = uy_y + uz_z

    outgoing_new = (
        outgoing
        - (ux + sound_speed) * (outgoing - outgoing_inner)
        - uy * outgoing_y
        - uz * outgoing_z
        - sound_speed * transverse_divergence
    )
    incoming_new = (
        incoming
        - uy * incoming_y
        - uz * incoming_z
        + sound_speed * transverse_divergence
    )
    target_ux = 0.5 * (outgoing_new + incoming_new)
    target_log_rho = (outgoing_new - incoming_new) / (2.0 * sound_speed)
    target_rho = reference_density * np.exp(target_log_rho)

    # Tangential momentum is a convective/vorticity characteristic with the
    # transverse pressure and cross-advection terms retained.
    target_uy = (
        uy
        - ux * (uy - uy_inner)
        - uy * uy_y
        - uz * uy_z
        - _CS2_3D * log_rho_y
    )
    target_uz = (
        uz
        - ux * (uz - uz_inner)
        - uy * uz_y
        - uz * uz_z
        - _CS2_3D * log_rho_z
    )

    convected_distribution, _ = _convective_boundary_distribution(
        boundary,
        interior,
        normal_velocity_component=0,
        outward_sign=1.0,
    )
    convected_rho, convected_ux, convected_uy, convected_uz = state(
        convected_distribution
    )
    nonequilibrium = convected_distribution - _equilibrium_d3q19(
        convected_rho, convected_ux, convected_uy, convected_uz,
    )
    reference_velocity = np.broadcast_to(
        np.asarray(reference_velocity_x, dtype=float), target_ux.shape,
    )
    # Preserve an exactly uniform carrier even at roundoff.
    uniform = (
        np.abs(rho - reference_density) <= 1.0e-14
    ) & (
        np.abs(ux - reference_velocity) <= 1.0e-14
    ) & (
        np.abs(uy) <= 1.0e-14
    ) & (
        np.abs(uz) <= 1.0e-14
    )
    target_rho = np.where(uniform, reference_density, target_rho)
    target_ux = np.where(uniform, reference_velocity, target_ux)
    target_uy = np.where(uniform, 0.0, target_uy)
    target_uz = np.where(uniform, 0.0, target_uz)
    return _equilibrium_d3q19(
        target_rho, target_ux, target_uy, target_uz,
    ) + nonequilibrium


def _apply_outer_boundaries_d3q19(
    distribution: np.ndarray,
    inlet_u_lattice: np.ndarray,
    *,
    boundary_scheme: str,
    outlet_lattice_density: float,
    nonstreamwise_wall_boundary: str = (
        "ground_no_slip_lateral_top_free_slip"
    ),
    pressure_open_nonstream_faces: tuple[str, ...] = (
        "y_min", "y_max", "z_min", "z_max",
    ),
    apply_streamwise_boundaries: bool = True,
    previous_distribution: np.ndarray | None = None,
) -> None:
    """Apply declared streamwise and wall boundary operators.

    The non-equilibrium extrapolation option follows the decomposition of Guo,
    Zheng and Shi (2002): the prescribed macroscopic boundary equilibrium is
    combined with the non-equilibrium population from the adjacent fluid
    plane.  It is an opt-in validation candidate; the historical equilibrium/
    copy boundary remains the default until canonical and external gates pass.
    """
    if nonstreamwise_wall_boundary not in {
        "ground_no_slip_lateral_top_free_slip",
        "all_no_slip_tunnel",
        "all_symmetry_free_slip",
        "all_pressure_open",
        "ground_free_slip_lateral_top_pressure_open",
    }:
        raise ValueError("unsupported nonstreamwise wall boundary")
    if boundary_scheme not in {
        "equilibrium_copy", "nonequilibrium_extrapolation",
        "characteristic_nonreflecting",
        "characteristic_multidimensional",
        "characteristic_multidimensional_extrapolated",
        "characteristic_multidimensional_vortical_limited",
        "characteristic_nscbc_convective",
    }:
        raise ValueError(f"unsupported D3Q19 flow boundary scheme {boundary_scheme!r}")
    valid_nonstream_faces = {"y_min", "y_max", "z_min", "z_max"}
    if (
        not isinstance(pressure_open_nonstream_faces, tuple)
        or len(set(pressure_open_nonstream_faces))
        != len(pressure_open_nonstream_faces)
        or any(
            face not in valid_nonstream_faces
            for face in pressure_open_nonstream_faces
        )
    ):
        raise ValueError(
            "pressure_open_nonstream_faces must contain unique y/z face names"
        )
    if not isinstance(apply_streamwise_boundaries, bool):
        raise TypeError("apply_streamwise_boundaries must be boolean")
    if previous_distribution is not None:
        previous_distribution = np.asarray(previous_distribution, dtype=float)
        if previous_distribution.shape != distribution.shape:
            raise ValueError(
                "previous_distribution must match the current distribution"
            )

    def apply_nonstreamwise_walls() -> None:
        if nonstreamwise_wall_boundary in {"all_pressure_open", "ground_free_slip_lateral_top_pressure_open"}:
            # This is a fixed-pressure/zero-gradient candidate, not yet a
            # qualified non-reflecting boundary.  The adjacent plane supplies
            # velocity and non-equilibrium stress; ambient lattice density is
            # prescribed only for non-equilibrium extrapolation.
            def apply_plane(axis: int, boundary: int, neighbour: int) -> None:
                population_axis = axis + 1
                neighbour_distribution = np.take(
                    distribution, [neighbour], axis=population_axis,
                ).copy()
                if boundary == 0:
                    second_neighbour = 2
                else:
                    second_neighbour = -3
                acoustic_fraction = None
                if boundary_scheme in {
                    "characteristic_multidimensional_extrapolated",
                    "characteristic_multidimensional_vortical_limited",
                    "characteristic_nscbc_convective",
                }:
                    boundary_size = distribution.shape[population_axis]
                    second_distribution = np.take(
                        distribution, [second_neighbour],
                        axis=population_axis,
                    ).copy()
                    if (
                        boundary_scheme
                        in {
                            "characteristic_multidimensional_vortical_limited",
                            "characteristic_nscbc_convective",
                        }
                        and boundary_size >= 3
                    ):
                        acoustic_fraction = _boundary_acoustic_fraction(
                            neighbour_distribution,
                            second_distribution,
                            normal_spatial_axis=axis,
                            normal_coordinate_sign=(
                                -1.0 if boundary == 0 else 1.0
                            ),
                            energy_weighted=(
                                boundary_scheme
                                == "characteristic_nscbc_convective"
                            ),
                            include_pressure_gradient=(
                                boundary_scheme
                                != "characteristic_nscbc_convective"
                            ),
                        )
                    if boundary_size >= 4:
                        third_neighbour = 3 if boundary == 0 else -4
                        if (
                            boundary_scheme
                            in {
                                "characteristic_multidimensional_vortical_limited",
                                "characteristic_nscbc_convective",
                            }
                        ):
                            third_distribution = np.take(
                                distribution, [third_neighbour],
                                axis=population_axis,
                            ).copy()
                            sampled_distribution = _limited_boundary_extrapolation(
                                neighbour_distribution,
                                second_distribution,
                                third_distribution,
                            )
                        else:
                            sampled_distribution = (
                                2.0 * neighbour_distribution
                                - second_distribution
                            )
                    elif boundary_size >= 3:
                        sampled_distribution = (
                            2.0 * neighbour_distribution
                            - second_distribution
                            if boundary_scheme
                            == "characteristic_multidimensional_extrapolated"
                            else neighbour_distribution
                        )
                    else:
                        sampled_distribution = neighbour_distribution
                else:
                    sampled_distribution = neighbour_distribution
                if boundary_scheme == "equilibrium_copy":
                    boundary_distribution = neighbour_distribution
                else:
                    neighbour_rho = np.maximum(
                        np.sum(sampled_distribution, axis=0), 1.0e-14,
                    )
                    neighbour_ux = np.tensordot(
                        _CX3, sampled_distribution, axes=(0, 0),
                    ) / neighbour_rho
                    neighbour_uy = np.tensordot(
                        _CY3, sampled_distribution, axes=(0, 0),
                    ) / neighbour_rho
                    neighbour_uz = np.tensordot(
                        _CZ3, sampled_distribution, axes=(0, 0),
                    ) / neighbour_rho
                    nonequilibrium = sampled_distribution - _equilibrium_d3q19(
                        neighbour_rho, neighbour_ux, neighbour_uy, neighbour_uz,
                    )
                    if boundary_scheme in {
                        "characteristic_nonreflecting",
                        "characteristic_multidimensional",
                        "characteristic_multidimensional_extrapolated",
                        "characteristic_multidimensional_vortical_limited",
                        "characteristic_nscbc_convective",
                    }:
                        sound_speed = math.sqrt(_CS2_3D)
                        component_velocity = (
                            neighbour_uz if axis == 0 else neighbour_uy
                        )
                        outward_sign = -1.0 if boundary == 0 else 1.0
                        outward_velocity = outward_sign * component_velocity
                        if boundary_scheme in {
                            "characteristic_multidimensional",
                            "characteristic_multidimensional_extrapolated",
                            "characteristic_multidimensional_vortical_limited",
                            "characteristic_nscbc_convective",
                        }:
                            if axis == 0:
                                reference_ux = np.full_like(
                                    neighbour_ux,
                                    inlet_u_lattice[neighbour],
                                )
                                tangent_a = neighbour_ux
                                tangent_b = neighbour_uy
                            else:
                                reference_ux = np.broadcast_to(
                                    inlet_u_lattice[:, None, None],
                                    neighbour_ux.shape,
                                )
                                tangent_a = neighbour_ux
                                tangent_b = neighbour_uz
                            (
                                boundary_density,
                                target_outward_velocity,
                                _incidence_cosine,
                            ) = _multidimensional_characteristic_state(
                                neighbour_rho,
                                outward_velocity,
                                tangent_a,
                                tangent_b,
                                reference_density=outlet_lattice_density,
                                reference_outward_velocity=0.0,
                                reference_tangential_velocity_a=reference_ux,
                                reference_tangential_velocity_b=0.0,
                                acoustic_fraction=acoustic_fraction,
                            )
                        else:
                            outgoing_invariant = (
                                outward_velocity
                                + sound_speed * np.log(neighbour_rho)
                            )
                            incoming_invariant = (
                                -sound_speed * math.log(outlet_lattice_density)
                            )
                            target_outward_velocity = 0.5 * (
                                outgoing_invariant + incoming_invariant
                            )
                            boundary_density = np.exp(
                                (outgoing_invariant - incoming_invariant)
                                / (2.0 * sound_speed)
                            )
                        if axis == 0:
                            target_ux = neighbour_ux
                            target_uy = neighbour_uy
                            target_uz = (
                                outward_sign * target_outward_velocity
                            )
                        else:
                            target_ux = neighbour_ux
                            target_uy = (
                                outward_sign * target_outward_velocity
                            )
                            target_uz = neighbour_uz
                    else:
                        boundary_density = np.full_like(
                            neighbour_rho, outlet_lattice_density,
                        )
                        target_ux = neighbour_ux
                        target_uy = neighbour_uy
                        target_uz = neighbour_uz
                    boundary_distribution = _equilibrium_d3q19(
                        boundary_density,
                        target_ux, target_uy, target_uz,
                    ) + nonequilibrium
                boundary_size = distribution.shape[population_axis]
                boundary_index = boundary if boundary >= 0 else boundary_size - 1
                target = [slice(None)] * distribution.ndim
                target[population_axis] = slice(
                    boundary_index, boundary_index + 1,
                )
                distribution[tuple(target)] = boundary_distribution

            selected = set(pressure_open_nonstream_faces)
            if nonstreamwise_wall_boundary == "ground_free_slip_lateral_top_pressure_open":
                selected.discard("z_min")
            if "z_min" in selected:
                apply_plane(0, 0, 1)
            if "z_max" in selected:
                apply_plane(0, -1, -2)
            if "y_min" in selected:
                apply_plane(1, 0, 1)
            if "y_max" in selected:
                apply_plane(1, -1, -2)
            # Unselected faces are symmetry planes. This internal selection is
            # used by canonical tests to isolate one boundary orientation;
            # the production all_pressure_open mode selects all four faces.
            for direction, (cx, cy, cz) in enumerate(_C3):
                if "z_min" not in selected and cz > 0:
                    reflected = _DIRECTION_INDEX[(cx, cy, -cz)]
                    distribution[direction, 0, :, :] = distribution[
                        reflected, 0, :, :
                    ].copy()
                if "z_max" not in selected and cz < 0:
                    reflected = _DIRECTION_INDEX[(cx, cy, -cz)]
                    distribution[direction, -1, :, :] = distribution[
                        reflected, -1, :, :
                    ].copy()
                if "y_min" not in selected and cy > 0:
                    reflected = _DIRECTION_INDEX[(cx, -cy, cz)]
                    distribution[direction, :, 0, :] = distribution[
                        reflected, :, 0, :
                    ].copy()
                if "y_max" not in selected and cy < 0:
                    reflected = _DIRECTION_INDEX[(cx, -cy, cz)]
                    distribution[direction, :, -1, :] = distribution[
                        reflected, :, -1, :
                    ].copy()
            return
        if nonstreamwise_wall_boundary == "all_symmetry_free_slip":
            # Symmetry/free-slip is a finite-domain surrogate for distant open
            # faces. It reverses only the normal population component and is
            # deliberately not described as a non-reflecting open boundary.
            for direction, (cx, cy, cz) in enumerate(_C3):
                if cz > 0:
                    reflected = _DIRECTION_INDEX[(cx, cy, -cz)]
                    distribution[direction, 0, :, :] = distribution[
                        reflected, 0, :, :
                    ].copy()
                if cz < 0:
                    reflected = _DIRECTION_INDEX[(cx, cy, -cz)]
                    distribution[direction, -1, :, :] = distribution[
                        reflected, -1, :, :
                    ].copy()
                if cy > 0:
                    reflected = _DIRECTION_INDEX[(cx, -cy, cz)]
                    distribution[direction, :, 0, :] = distribution[
                        reflected, :, 0, :
                    ].copy()
                if cy < 0:
                    reflected = _DIRECTION_INDEX[(cx, -cy, cz)]
                    distribution[direction, :, -1, :] = distribution[
                        reflected, :, -1, :
                    ].copy()
            return
        # No-slip lower z wall: incoming populations reverse every component.
        for direction, (_cx, _cy, cz) in enumerate(_C3):
            if cz > 0:
                distribution[direction, 0, :, :] = distribution[
                    _OPPOSITE3[direction], 0, :, :
                ].copy()
        if nonstreamwise_wall_boundary == "all_no_slip_tunnel":
            # Full reversal at tunnel sidewalls and roof retains tangential
            # viscous stress instead of treating a confined tunnel as open air.
            for direction, (_cx, cy, cz) in enumerate(_C3):
                if cz < 0:
                    distribution[direction, -1, :, :] = distribution[
                        _OPPOSITE3[direction], -1, :, :
                    ].copy()
                if cy > 0:
                    distribution[direction, :, 0, :] = distribution[
                        _OPPOSITE3[direction], :, 0, :
                    ].copy()
                if cy < 0:
                    distribution[direction, :, -1, :] = distribution[
                        _OPPOSITE3[direction], :, -1, :
                    ].copy()
            return
        # Open-atmosphere free-slip top and sides reverse only the normal component.
        for direction, (cx, cy, cz) in enumerate(_C3):
            if cz < 0:
                reflected = _DIRECTION_INDEX[(cx, cy, -cz)]
                distribution[direction, -1, :, :] = distribution[
                    reflected, -1, :, :
                ].copy()
            if cy > 0:
                reflected = _DIRECTION_INDEX[(cx, -cy, cz)]
                distribution[direction, :, 0, :] = distribution[
                    reflected, :, 0, :
                ].copy()
            if cy < 0:
                reflected = _DIRECTION_INDEX[(cx, -cy, cz)]
                distribution[direction, :, -1, :] = distribution[
                    reflected, :, -1, :
                ].copy()

    nz, ny = distribution.shape[1:3]
    ux_in = np.broadcast_to(inlet_u_lattice[:, None, None], (nz, ny, 1))
    zero = np.zeros_like(ux_in)
    pressure_open = nonstreamwise_wall_boundary in {"all_pressure_open", "ground_free_slip_lateral_top_pressure_open"}
    if pressure_open:
        # Apply y/z first so the prescribed x inlet and pressure outlet retain
        # precedence at the twelve shared edge lines and eight corner nodes.
        apply_nonstreamwise_walls()
    if not apply_streamwise_boundaries:
        if not pressure_open:
            apply_nonstreamwise_walls()
        return
    if boundary_scheme == "equilibrium_copy":
        if not pressure_open:
            apply_nonstreamwise_walls()
        rho_in = np.ones((nz, ny, 1), dtype=float)
        distribution[:, :, :, 0] = _equilibrium_d3q19(
            rho_in, ux_in, zero, zero,
        )[:, :, :, 0]
        distribution[:, :, :, -1] = distribution[:, :, :, -2]
        return
    if boundary_scheme not in {
        "nonequilibrium_extrapolation", "characteristic_nonreflecting",
        "characteristic_multidimensional",
        "characteristic_multidimensional_extrapolated",
        "characteristic_multidimensional_vortical_limited",
        "characteristic_nscbc_convective",
    }:
        raise ValueError(f"unsupported D3Q19 flow boundary scheme {boundary_scheme!r}")

    inlet_neighbor = distribution[:, :, :, 1:2].copy()
    inlet_rho = np.maximum(np.sum(inlet_neighbor, axis=0), 1.0e-14)
    inlet_ux = np.tensordot(_CX3, inlet_neighbor, axes=(0, 0)) / inlet_rho
    inlet_uy = np.tensordot(_CY3, inlet_neighbor, axes=(0, 0)) / inlet_rho
    inlet_uz = np.tensordot(_CZ3, inlet_neighbor, axes=(0, 0)) / inlet_rho
    inlet_nonequilibrium = inlet_neighbor - _equilibrium_d3q19(
        inlet_rho, inlet_ux, inlet_uy, inlet_uz,
    )
    distribution[:, :, :, 0] = (
        _equilibrium_d3q19(inlet_rho, ux_in, zero, zero)
        + inlet_nonequilibrium
    )[:, :, :, 0]

    if boundary_scheme == "characteristic_nscbc_convective":
        history = (
            distribution
            if previous_distribution is None
            else previous_distribution
        )
        nscbc_distribution = _nscbc_xmax_boundary_distribution(
            history[:, :, :, -1:],
            history[:, :, :, -2:-1],
            reference_density=outlet_lattice_density,
            reference_velocity_x=ux_in,
        )
        distribution[:, :, :, -1] = nscbc_distribution[:, :, :, 0]
        if not pressure_open:
            apply_nonstreamwise_walls()
        return

    outlet_neighbor = distribution[:, :, :, -2:-1].copy()
    outlet_acoustic_fraction = None
    if boundary_scheme in {
        "characteristic_multidimensional_extrapolated",
        "characteristic_multidimensional_vortical_limited",
        "characteristic_nscbc_convective",
    }:
        outlet_second = distribution[:, :, :, -3:-2].copy()
        if (
            boundary_scheme in {
                "characteristic_multidimensional_vortical_limited",
                "characteristic_nscbc_convective",
            }
        ):
            outlet_acoustic_fraction = _boundary_acoustic_fraction(
                outlet_neighbor,
                outlet_second,
                normal_spatial_axis=2,
                normal_coordinate_sign=1.0,
                energy_weighted=(
                    boundary_scheme == "characteristic_nscbc_convective"
                ),
                include_pressure_gradient=(
                    boundary_scheme != "characteristic_nscbc_convective"
                ),
            )
        if distribution.shape[3] >= 4:
            if (
                boundary_scheme
                == "characteristic_multidimensional_vortical_limited"
            ):
                outlet_sample = _limited_boundary_extrapolation(
                    outlet_neighbor,
                    outlet_second,
                    distribution[:, :, :, -4:-3],
                )
            else:
                outlet_sample = 2.0 * outlet_neighbor - outlet_second
        else:
            outlet_sample = (
                2.0 * outlet_neighbor - outlet_second
                if distribution.shape[3] >= 3
                and boundary_scheme
                == "characteristic_multidimensional_extrapolated"
                else outlet_neighbor
            )
    else:
        outlet_sample = outlet_neighbor
    neighbor_rho = np.maximum(np.sum(outlet_sample, axis=0), 1.0e-14)
    neighbor_ux = np.tensordot(_CX3, outlet_sample, axes=(0, 0)) / neighbor_rho
    neighbor_uy = np.tensordot(_CY3, outlet_sample, axes=(0, 0)) / neighbor_rho
    neighbor_uz = np.tensordot(_CZ3, outlet_sample, axes=(0, 0)) / neighbor_rho
    outlet_nonequilibrium = outlet_sample - _equilibrium_d3q19(
        neighbor_rho, neighbor_ux, neighbor_uy, neighbor_uz,
    )
    convective_state = None
    if boundary_scheme == "characteristic_nscbc_convective":
        if previous_distribution is None:
            previous_distribution = distribution
        convective_distribution, _convective_courant = (
            _convective_boundary_distribution(
                previous_distribution[:, :, :, -1:],
                previous_distribution[:, :, :, -2:-1],
                normal_velocity_component=0,
                outward_sign=1.0,
            )
        )
        convective_rho = np.maximum(
            np.sum(convective_distribution, axis=0), 1.0e-14,
        )
        convective_ux = np.tensordot(
            _CX3, convective_distribution, axes=(0, 0),
        ) / convective_rho
        convective_uy = np.tensordot(
            _CY3, convective_distribution, axes=(0, 0),
        ) / convective_rho
        convective_uz = np.tensordot(
            _CZ3, convective_distribution, axes=(0, 0),
        ) / convective_rho
        convective_nonequilibrium = (
            convective_distribution
            - _equilibrium_d3q19(
                convective_rho,
                convective_ux,
                convective_uy,
                convective_uz,
            )
        )
        convective_state = (
            convective_rho,
            convective_ux,
            convective_uy,
            convective_uz,
            convective_nonequilibrium,
        )
    if boundary_scheme in {
        "characteristic_nonreflecting", "characteristic_multidimensional",
        "characteristic_multidimensional_extrapolated",
        "characteristic_multidimensional_vortical_limited",
        "characteristic_nscbc_convective",
    }:
        sound_speed = math.sqrt(_CS2_3D)
        # Linearise the incoming characteristic around the declared carrier
        # co-flow, not around stagnant ambient. This preserves a uniform mean
        # through-flow while still allowing its outgoing perturbation to leave.
        reference_outlet_ux = np.broadcast_to(
            inlet_u_lattice[:, None, None], neighbor_ux.shape,
        )
        if boundary_scheme in {
            "characteristic_multidimensional",
            "characteristic_multidimensional_extrapolated",
            "characteristic_multidimensional_vortical_limited",
            "characteristic_nscbc_convective",
        }:
            outlet_rho, target_ux, _incidence_cosine = (
                _multidimensional_characteristic_state(
                    neighbor_rho,
                    neighbor_ux,
                    neighbor_uy,
                    neighbor_uz,
                    reference_density=outlet_lattice_density,
                    reference_outward_velocity=reference_outlet_ux,
                    acoustic_fraction=(
                        None
                        if boundary_scheme == "characteristic_nscbc_convective"
                        else outlet_acoustic_fraction
                    ),
                )
            )
            if boundary_scheme == "characteristic_nscbc_convective":
                assert convective_state is not None
                assert outlet_acoustic_fraction is not None
                (
                    convective_rho,
                    convective_ux,
                    convective_uy,
                    convective_uz,
                    convective_nonequilibrium,
                ) = convective_state
                weight = outlet_acoustic_fraction
                outlet_rho = np.exp(
                    np.log(convective_rho)
                    + weight * (
                        np.log(np.maximum(outlet_rho, 1.0e-14))
                        - np.log(convective_rho)
                    )
                )
                target_ux = (
                    convective_ux + weight * (target_ux - convective_ux)
                )
                neighbor_uy = convective_uy
                neighbor_uz = convective_uz
                outlet_nonequilibrium = convective_nonequilibrium
        else:
            outgoing_invariant = (
                neighbor_ux + sound_speed * np.log(neighbor_rho)
            )
            incoming_invariant = (
                reference_outlet_ux
                - sound_speed * math.log(outlet_lattice_density)
            )
            target_ux = 0.5 * (outgoing_invariant + incoming_invariant)
            outlet_rho = np.exp(
                (outgoing_invariant - incoming_invariant)
                / (2.0 * sound_speed)
            )
    else:
        outlet_rho = np.full_like(neighbor_rho, outlet_lattice_density)
        target_ux = neighbor_ux
    distribution[:, :, :, -1] = (
        _equilibrium_d3q19(
            outlet_rho, target_ux, neighbor_uy, neighbor_uz,
        ) + outlet_nonequilibrium
    )[:, :, :, 0]

    # At edges and corners wall impermeability has precedence over the open
    # boundary prescription.  Interior inlet/outlet nodes retain the exact
    # prescribed velocity/density moments.
    if not pressure_open:
        apply_nonstreamwise_walls()


def _stream_d3q19(
    post_collision: np.ndarray,
    solid: np.ndarray,
    blocked_links: np.ndarray | None = None,
) -> np.ndarray:
    """Stream populations and apply halfway no-slip bounce-back on blocked links."""
    if post_collision.shape != (19, *solid.shape):
        raise ValueError("D3Q19 populations must have shape (19, nz, ny, nx)")
    if blocked_links is not None and blocked_links.shape != post_collision.shape:
        raise ValueError("D3Q19 link fractions must match population shape")
    streamed = np.empty_like(post_collision)
    for direction, (cx, cy, cz) in enumerate(_C3):
        streamed[direction] = np.roll(
            post_collision[direction], shift=(cz, cy, cx), axis=(0, 1, 2),
        )
    for direction in range(1, 19):
        cx, cy, cz = _C3[direction]
        destination_solid = np.roll(
            solid, shift=(-cz, -cy, -cx), axis=(0, 1, 2),
        )
        bounce = (~solid) & destination_solid
        streamed[_OPPOSITE3[direction], bounce] = post_collision[direction, bounce]
        if blocked_links is not None:
            arc_bounce = (
                (~solid) & (~destination_solid)
                & np.isfinite(blocked_links[direction])
            )
            q = blocked_links[direction]
            # Bouzidi-Firdaouss-Lallemand linear interpolation. If q<1/2
            # would sample behind a separate voxel solid, fall back to the
            # robust halfway rule for that link.
            upstream_solid = np.roll(
                solid, shift=(cz, cy, cx), axis=(0, 1, 2),
            )
            q_effective = np.where((q < 0.5) & upstream_solid, 0.5, q)
            near = arc_bounce & (q_effective < 0.5)
            far = arc_bounce & ~near
            if np.any(near):
                upstream_population = np.roll(
                    post_collision[direction],
                    shift=(cz, cy, cx), axis=(0, 1, 2),
                )
                streamed[_OPPOSITE3[direction], near] = (
                    2.0 * q_effective[near] * post_collision[direction, near]
                    + (1.0 - 2.0 * q_effective[near])
                    * upstream_population[near]
                )
            if np.any(far):
                q_far = q_effective[far]
                streamed[_OPPOSITE3[direction], far] = (
                    post_collision[direction, far] / (2.0 * q_far)
                    + ((2.0 * q_far - 1.0) / (2.0 * q_far))
                    * post_collision[_OPPOSITE3[direction], far]
                )
    return streamed


def _advance_d3q19(
    distribution: np.ndarray,
    *,
    solid: np.ndarray,
    inlet_u_lattice: np.ndarray,
    lattice_viscosity: float,
    flow_boundary_scheme: str,
    outlet_lattice_density: float,
    nonstreamwise_wall_boundary: str,
    subgrid_model: str,
    smagorinsky_constant: float,
    wale_constant: float,
    maximum_relaxation_time: float,
    acceleration_lattice: tuple[
        float | np.ndarray, float | np.ndarray, float | np.ndarray
    ] = (0.0, 0.0, 0.0),
    pressure_open_nonstream_faces: tuple[str, ...] = (
        "y_min", "y_max", "z_min", "z_max",
    ),
    apply_streamwise_boundaries: bool = True,
    blocked_links: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Advance one forced D3Q19 step and return distribution, tau and SGS nu."""
    rho, ux, uy, uz = _moments_d3q19(
        distribution, solid, acceleration_lattice,
    )
    equilibrium = _equilibrium_d3q19(rho, ux, uy, uz)
    tau, eddy_viscosity_lattice = subgrid_relaxation_time_lattice(
        ux, uy, uz, solid,
        base_lattice_viscosity=lattice_viscosity,
        model=subgrid_model,
        smagorinsky_constant=smagorinsky_constant,
        wale_constant=wale_constant,
        maximum_relaxation_time=maximum_relaxation_time,
    )
    forcing = _guo_force_d3q19(
        rho, ux, uy, uz, acceleration_lattice, tau,
    )
    post_collision = distribution - (distribution - equilibrium) / tau[None] + forcing
    streamed = _stream_d3q19(post_collision, solid, blocked_links)
    _apply_outer_boundaries_d3q19(
        streamed,
        inlet_u_lattice,
        boundary_scheme=flow_boundary_scheme,
        outlet_lattice_density=outlet_lattice_density,
        nonstreamwise_wall_boundary=nonstreamwise_wall_boundary,
        pressure_open_nonstream_faces=pressure_open_nonstream_faces,
        apply_streamwise_boundaries=apply_streamwise_boundaries,
        previous_distribution=distribution,
    )
    return streamed, tau, eddy_viscosity_lattice


def _projected_faces_from_cell_velocity(
    *,
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    solid: np.ndarray,
    inlet_physical: np.ndarray,
    spacing_m: float,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    nonstreamwise_wall_boundary: str = (
        "ground_no_slip_lateral_top_free_slip"
    ),
) -> tuple[np.ndarray, ...]:
    """Build impermeable faces and project them to cellwise mass consistency."""
    u_face, v_face, w_face, raw_divergence = _raw_faces_from_cell_velocity(
        u=u, v=v, w=w, solid=solid, inlet_physical=inlet_physical,
        spacing_m=spacing_m,
        nonstreamwise_wall_boundary=nonstreamwise_wall_boundary,
        blocked_faces=blocked_faces,
    )
    dx = spacing_m
    (
        u_face, v_face, w_face,
        divergence_before_projection, divergence_after_projection,
        projection_iterations,
    ) = _mass_consistent_face_projection_3d(
        u_face, v_face, w_face, solid, np.asarray((dx, dx, dx)),
        blocked_faces=blocked_faces,
    )
    return (
        u_face, v_face, w_face, raw_divergence,
        divergence_before_projection, divergence_after_projection,
        projection_iterations,
    )


def _raw_faces_from_cell_velocity(
    *,
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    solid: np.ndarray,
    inlet_physical: np.ndarray,
    spacing_m: float,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    nonstreamwise_wall_boundary: str = (
        "ground_no_slip_lateral_top_free_slip"
    ),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate momentum-predictor velocity without a pressure projection."""
    nz, ny, nx = solid.shape
    dx = spacing_m
    u_face = np.empty((nz, ny, nx + 1), dtype=float)
    v_face = np.zeros((nz, ny + 1, nx), dtype=float)
    w_face = np.zeros((nz + 1, ny, nx), dtype=float)
    u_face[:, :, 0] = _coerce_inlet_plane(inlet_physical, nz, ny)
    u_face[:, :, -1] = u[:, :, -1]
    u_face[:, :, 1:-1] = 0.5 * (u[:, :, :-1] + u[:, :, 1:])
    v_face[:, 1:-1, :] = 0.5 * (v[:, :-1, :] + v[:, 1:, :])
    w_face[1:-1, :, :] = 0.5 * (w[:-1, :, :] + w[1:, :, :])
    if nonstreamwise_wall_boundary in {"all_pressure_open", "ground_free_slip_lateral_top_pressure_open"}:
        v_face[:, 0, :] = v[:, 0, :]
        v_face[:, -1, :] = v[:, -1, :]
        if nonstreamwise_wall_boundary == "all_pressure_open":
            w_face[0, :, :] = w[0, :, :]
        w_face[-1, :, :] = w[-1, :, :]
    u_face[:, :, 1:-1][solid[:, :, :-1] | solid[:, :, 1:]] = 0.0
    v_face[:, 1:-1, :][solid[:, :-1, :] | solid[:, 1:, :]] = 0.0
    w_face[1:-1, :, :][solid[:-1, :, :] | solid[1:, :, :]] = 0.0
    if blocked_faces is not None:
        blocked_x, blocked_y, blocked_z = blocked_faces
        if (
            blocked_x.shape != u_face.shape
            or blocked_y.shape != v_face.shape
            or blocked_z.shape != w_face.shape
        ):
            raise ValueError("blocked face masks must match staggered flow faces")
        u_face[blocked_x] = 0.0
        v_face[blocked_y] = 0.0
        w_face[blocked_z] = 0.0
    v_face[:, 0, :][solid[:, 0, :]] = 0.0
    v_face[:, -1, :][solid[:, -1, :]] = 0.0
    w_face[0, :, :][solid[0, :, :]] = 0.0
    w_face[-1, :, :][solid[-1, :, :]] = 0.0
    raw_divergence = (
        (u_face[:, :, 1:] - u_face[:, :, :-1]) / dx
        + (v_face[:, 1:, :] - v_face[:, :-1, :]) / dx
        + (w_face[1:, :, :] - w_face[:-1, :, :]) / dx
    )
    raw_divergence[solid] = 0.0
    return u_face, v_face, w_face, raw_divergence


def _increment_projected_faces_from_cell_delta(
    *,
    previous_u_face: np.ndarray,
    previous_v_face: np.ndarray,
    previous_w_face: np.ndarray,
    delta_u: np.ndarray,
    delta_v: np.ndarray,
    delta_w: np.ndarray,
    solid: np.ndarray,
    spacing_m: float,
    inlet_physical: np.ndarray | None = None,
    blocked_faces: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    nonstreamwise_wall_boundary: str = (
        "ground_no_slip_lateral_top_free_slip"
    ),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Preserve projected face state and add only the momentum-predictor delta.

    Reconstructing absolute faces from cell averages is not an identity even at
    zero time step.  Repeating that map creates time-step-dependent filtering.
    The incremental map is exactly invariant when the momentum update is zero.
    """
    dx = spacing_m
    u_face = np.asarray(previous_u_face, dtype=float).copy()
    v_face = np.asarray(previous_v_face, dtype=float).copy()
    w_face = np.asarray(previous_w_face, dtype=float).copy()
    fluid_x = (~solid[:, :, :-1]) & (~solid[:, :, 1:])
    fluid_y = (~solid[:, :-1, :]) & (~solid[:, 1:, :])
    fluid_z = (~solid[:-1, :, :]) & (~solid[1:, :, :])
    if blocked_faces is not None:
        blocked_x, blocked_y, blocked_z = blocked_faces
        if (
            blocked_x.shape != u_face.shape
            or blocked_y.shape != v_face.shape
            or blocked_z.shape != w_face.shape
        ):
            raise ValueError("blocked face masks must match staggered flow faces")
        fluid_x &= ~blocked_x[:, :, 1:-1]
        fluid_y &= ~blocked_y[:, 1:-1, :]
        fluid_z &= ~blocked_z[1:-1, :, :]
    u_increment = 0.5 * (delta_u[:, :, :-1] + delta_u[:, :, 1:])
    v_increment = 0.5 * (delta_v[:, :-1, :] + delta_v[:, 1:, :])
    w_increment = 0.5 * (delta_w[:-1, :, :] + delta_w[1:, :, :])
    u_face[:, :, 1:-1][fluid_x] += u_increment[fluid_x]
    v_face[:, 1:-1, :][fluid_y] += v_increment[fluid_y]
    w_face[1:-1, :, :][fluid_z] += w_increment[fluid_z]
    outlet_fluid = ~solid[:, :, -1]
    u_face[:, :, -1][outlet_fluid] += delta_u[:, :, -1][outlet_fluid]
    if nonstreamwise_wall_boundary in {"all_pressure_open", "ground_free_slip_lateral_top_pressure_open"}:
        low_y_fluid = ~solid[:, 0, :]
        high_y_fluid = ~solid[:, -1, :]
        low_z_fluid = ~solid[0, :, :]
        high_z_fluid = ~solid[-1, :, :]
        v_face[:, 0, :][low_y_fluid] += delta_v[:, 0, :][low_y_fluid]
        v_face[:, -1, :][high_y_fluid] += delta_v[:, -1, :][high_y_fluid]
        if nonstreamwise_wall_boundary == "all_pressure_open":
            w_face[0, :, :][low_z_fluid] += delta_w[0, :, :][low_z_fluid]
        w_face[-1, :, :][high_z_fluid] += delta_w[-1, :, :][high_z_fluid]
    u_face[:, :, 1:-1][~fluid_x] = 0.0
    v_face[:, 1:-1, :][~fluid_y] = 0.0
    w_face[1:-1, :, :][~fluid_z] = 0.0
    u_face[:, :, 0][solid[:, :, 0]] = 0.0
    if inlet_physical is not None:
        u_face[:, :, 0] = _coerce_inlet_plane(
            inlet_physical, solid.shape[0], solid.shape[1],
        )
        u_face[:, :, 0][solid[:, :, 0]] = 0.0
    u_face[:, :, -1][solid[:, :, -1]] = 0.0
    if nonstreamwise_wall_boundary not in {"all_pressure_open", "ground_free_slip_lateral_top_pressure_open"}:
        v_face[:, 0, :] = 0.0
        v_face[:, -1, :] = 0.0
        w_face[0, :, :] = 0.0
        w_face[-1, :, :] = 0.0
    else:
        v_face[:, 0, :][solid[:, 0, :]] = 0.0
        v_face[:, -1, :][solid[:, -1, :]] = 0.0
        w_face[0, :, :][solid[0, :, :]] = 0.0
        w_face[-1, :, :][solid[-1, :, :]] = 0.0
        if nonstreamwise_wall_boundary == "ground_free_slip_lateral_top_pressure_open":
            w_face[0, :, :] = 0.0
    raw_divergence = (
        (u_face[:, :, 1:] - u_face[:, :, :-1]) / dx
        + (v_face[:, 1:, :] - v_face[:, :-1, :]) / dx
        + (w_face[1:, :, :] - w_face[:-1, :, :]) / dx
    )
    raw_divergence[solid] = 0.0
    if blocked_faces is not None:
        u_face[blocked_x] = 0.0
        v_face[blocked_y] = 0.0
        w_face[blocked_z] = 0.0
    return u_face, v_face, w_face, raw_divergence


def _initialize_projected_inlet_flow(
    config: ResolvedObstacleFlow3DConfig,
) -> tuple[np.ndarray, ...]:
    """Build a conservative full-site seed for the finite-volume solver.

    A molecular-viscosity atmospheric flow on a metre-scale grid cannot be
    represented by the BGK relaxation used by the optional D3Q19 spin-up: its
    relaxation time approaches the singular value 0.5.  The finite-volume
    momentum path does not need that lattice spin-up.  This initializer maps
    the declared vector inlet through the obstacle geometry, projects the face
    field to discrete continuity, and lets the physical-viscosity momentum
    equations perform all subsequent evolution.

    The returned field is an initial condition, not a converged wake and not a
    turbulence model.  Diagnostics intentionally report zero completed spin-up
    steps and no lattice relaxation time.
    """
    x, y, z, solid = _solid_mask(config)
    dx = config.domain_m[0] / config.cells[0]
    blocked_faces = _circular_arc_blocked_faces(config, x, y, z)
    inlet_physical, inlet_lattice = _inlet_profile(config, z)
    reference_plane = _active_inlet_plane(config, y, inlet_physical)
    inlet_u, inlet_v, inlet_w = _active_inlet_velocity_planes(
        config, reference_plane, 0.0,
    )
    nz, ny, nx = solid.shape
    u = np.broadcast_to(inlet_u[:, :, None], (nz, ny, nx)).copy()
    v = np.broadcast_to(inlet_v[:, :, None], (nz, ny, nx)).copy()
    w = np.broadcast_to(inlet_w[:, :, None], (nz, ny, nx)).copy()
    u[solid] = v[solid] = w[solid] = 0.0
    (
        u_face, v_face, w_face, raw_divergence,
        divergence_before_projection, divergence_after_projection,
        projection_iterations,
    ) = _projected_faces_from_cell_velocity(
        u=u,
        v=v,
        w=w,
        solid=solid,
        inlet_physical=inlet_u,
        spacing_m=dx,
        blocked_faces=blocked_faces,
        nonstreamwise_wall_boundary=config.nonstreamwise_wall_boundary,
    )
    # Keep cell momentum consistent with the projected transport faces.  Solid
    # cells remain excluded and are restored to exactly zero afterwards.
    u = 0.5 * (u_face[:, :, :-1] + u_face[:, :, 1:])
    v = 0.5 * (v_face[:, :-1, :] + v_face[:, 1:, :])
    w = 0.5 * (w_face[:-1, :, :] + w_face[1:, :, :])
    u[solid] = v[solid] = w[solid] = 0.0
    assert config.velocity_reference_scale_m_s is not None
    lattice_reference_dt = (
        config.lattice_inlet_speed * dx / config.velocity_reference_scale_m_s
    )
    velocity_to_lattice = lattice_reference_dt / dx
    lattice_density = np.ones_like(u)
    distribution = _equilibrium_d3q19(
        lattice_density,
        u * velocity_to_lattice,
        v * velocity_to_lattice,
        w * velocity_to_lattice,
    )
    pressure = np.zeros_like(u)
    pressure_coefficient = np.zeros_like(u)
    pressure[solid] = pressure_coefficient[solid] = np.nan
    eddy_viscosity = np.zeros_like(u)
    return (
        x, y, z, solid, u, v, w, u_face, v_face, w_face,
        pressure, pressure_coefficient, eddy_viscosity,
        np.empty((0, 2), dtype=float), lattice_density,
        lattice_reference_dt, None, None, 0, raw_divergence,
        divergence_before_projection, divergence_after_projection,
        projection_iterations, distribution, None, 0,
    )


def _solve_flow(
    config: ResolvedObstacleFlow3DConfig,
) -> tuple[np.ndarray, ...]:
    if config.carrier_flow_initialization == "projected_inlet":
        return _initialize_projected_inlet_flow(config)
    x, y, z, solid = _solid_mask(config)
    blocked_links = _circular_arc_blocked_links(config, x, y, z)
    dx = config.domain_m[0] / config.cells[0]
    assert config.velocity_reference_scale_m_s is not None
    physical_dt = (
        config.lattice_inlet_speed * dx
        / config.velocity_reference_scale_m_s
    )
    inlet_physical, inlet_lattice = _inlet_profile(config, z)
    lattice_viscosity = config.kinematic_viscosity_m2_s * physical_dt / dx**2
    base_tau = 0.5 + 3.0 * lattice_viscosity
    if base_tau <= 0.505:
        raise ValueError(
            "D3Q19 relaxation time is too close to 0.5; increase viscosity, "
            "lattice speed, or cell resolution"
        )
    nz, ny, nx = solid.shape
    rho = np.ones_like(solid, dtype=float)
    ux = np.broadcast_to(inlet_lattice[:, None, None], (nz, ny, nx)).copy()
    uy = np.zeros_like(ux)
    uz = np.zeros_like(ux)
    ux[solid] = 0.0
    distribution = _equilibrium_d3q19(rho, ux, uy, uz)
    previous = (ux.copy(), uy.copy(), uz.copy())
    residual_history: list[tuple[float, float]] = []
    actual_steps = config.flow_steps
    eddy_viscosity_lattice = np.zeros_like(ux)
    maximum_tau = base_tau
    average_steps = config.flow_average_last_steps
    averaging_start_step = config.flow_steps - average_steps + 1
    average_count = 0
    average_rho = np.zeros_like(rho)
    average_u = np.zeros_like(ux)
    average_v = np.zeros_like(uy)
    average_w = np.zeros_like(uz)
    average_eddy_viscosity = np.zeros_like(ux)
    first_half_u = np.zeros_like(ux)
    first_half_v = np.zeros_like(uy)
    first_half_w = np.zeros_like(uz)
    second_half_u = np.zeros_like(ux)
    second_half_v = np.zeros_like(uy)
    second_half_w = np.zeros_like(uz)
    first_half_count = 0
    second_half_count = 0

    for step in range(1, config.flow_steps + 1):
        distribution, tau, eddy_viscosity_lattice = _advance_d3q19(
            distribution,
            solid=solid,
            inlet_u_lattice=inlet_lattice,
            lattice_viscosity=lattice_viscosity,
            flow_boundary_scheme=config.flow_boundary_scheme,
            outlet_lattice_density=config.outlet_lattice_density,
            nonstreamwise_wall_boundary=config.nonstreamwise_wall_boundary,
            subgrid_model=config.subgrid_model,
            smagorinsky_constant=config.smagorinsky_constant,
            wale_constant=config.wale_constant,
            maximum_relaxation_time=config.maximum_relaxation_time,
            blocked_links=blocked_links,
        )
        maximum_tau = max(maximum_tau, float(np.max(tau[~solid])))
        current_moments = None
        if average_steps and step >= averaging_start_step:
            current_moments = _moments_d3q19(distribution, solid)
            current_rho, current_u, current_v, current_w = current_moments
            average_rho += current_rho
            average_u += current_u
            average_v += current_v
            average_w += current_w
            average_eddy_viscosity += eddy_viscosity_lattice
            average_count += 1
            if average_count <= average_steps // 2:
                first_half_u += current_u
                first_half_v += current_v
                first_half_w += current_w
                first_half_count += 1
            else:
                second_half_u += current_u
                second_half_v += current_v
                second_half_w += current_w
                second_half_count += 1
        if step % config.convergence_check_interval == 0 or step == config.flow_steps:
            if not np.all(np.isfinite(distribution)):
                raise RuntimeError(
                    "D3Q19 flow became non-finite; increase viscosity or reduce lattice speed"
                )
            if current_moments is None:
                current_moments = _moments_d3q19(distribution, solid)
            _, current_u, current_v, current_w = current_moments
            change = np.sqrt(
                (current_u - previous[0])**2
                + (current_v - previous[1])**2
                + (current_w - previous[2])**2
            )
            reference = max(
                float(np.max(np.sqrt(current_u**2 + current_v**2 + current_w**2))),
                float(np.max(inlet_lattice)),
                1.0e-15,
            )
            residual = float(np.max(change[~solid]) / reference)
            residual_history.append((float(step), residual))
            previous = (current_u.copy(), current_v.copy(), current_w.copy())
            if step >= config.minimum_flow_steps and residual <= config.flow_convergence_tolerance:
                actual_steps = step
                break

    rho, ux, uy, uz = _moments_d3q19(distribution, solid)
    statistical_residual = None
    if average_count:
        if average_count != average_steps:
            raise RuntimeError("D3Q19 flow averaging window was not completed")
        rho = average_rho / average_count
        ux = average_u / average_count
        uy = average_v / average_count
        uz = average_w / average_count
        eddy_viscosity_lattice = average_eddy_viscosity / average_count
        if first_half_count <= 0 or second_half_count <= 0:
            raise RuntimeError("D3Q19 flow averaging halves are empty")
        first_u = first_half_u / first_half_count
        first_v = first_half_v / first_half_count
        first_w = first_half_w / first_half_count
        second_u = second_half_u / second_half_count
        second_v = second_half_v / second_half_count
        second_w = second_half_w / second_half_count
        mean_change = np.sqrt(
            (second_u - first_u) ** 2
            + (second_v - first_v) ** 2
            + (second_w - first_w) ** 2
        )
        mean_reference = max(
            float(np.max(np.sqrt(ux**2 + uy**2 + uz**2))),
            float(np.max(inlet_lattice)),
            1.0e-15,
        )
        statistical_residual = float(np.max(mean_change[~solid]) / mean_reference)
    velocity_scale = dx / physical_dt
    u, v, w = ux * velocity_scale, uy * velocity_scale, uz * velocity_scale
    u[solid] = v[solid] = w[solid] = 0.0
    eddy_viscosity_physical = eddy_viscosity_lattice * dx**2 / physical_dt
    pressure_scale = config.ambient_air_density_kg_m3 * velocity_scale**2
    outlet_reference = float(np.mean(rho[:, :, -1]))
    pressure = _CS2_3D * (rho - outlet_reference) * pressure_scale
    assert config.velocity_reference_scale_m_s is not None
    dynamic_pressure = (
        0.5 * config.ambient_air_density_kg_m3
        * config.velocity_reference_scale_m_s**2
    )
    pressure_coefficient = pressure / dynamic_pressure
    pressure[solid] = pressure_coefficient[solid] = np.nan

    (
        u_face, v_face, w_face, raw_divergence,
        divergence_before_projection, divergence_after_projection,
        projection_iterations,
    ) = _projected_faces_from_cell_velocity(
        u=u, v=v, w=w, solid=solid, inlet_physical=inlet_physical,
        spacing_m=dx,
        nonstreamwise_wall_boundary=config.nonstreamwise_wall_boundary,
    )
    return (
        x, y, z, solid, u, v, w, u_face, v_face, w_face,
        pressure, pressure_coefficient, eddy_viscosity_physical,
        np.asarray(residual_history), rho, physical_dt, base_tau, maximum_tau,
        actual_steps, raw_divergence, divergence_before_projection,
        divergence_after_projection, projection_iterations, distribution,
        statistical_residual, average_count,
    )


def solve_resolved_obstacle_flow_3d(
    config: ResolvedObstacleFlow3DConfig,
    *,
    receptors: Sequence[FastObstacleReceptor3D] = (),
    handoff_request: BinaryMixtureHandoffPlane3DRequest | None = None,
    handoff_observer: Callable[[BinaryMixtureHandoffPlane3DFrame], None] | None = None,
    restart_state: ResolvedObstacleFlow3DRestartState | None = None,
    stop_time_s: float | None = None,
    face_flux_callback=None,
    enthalpy_flux_callback=None,
    interface_state_callback=None,
    interface_conserved_state_callback=None,
    interface_momentum_callback=None,
    interface_conserved_mass_callback=None,
) -> ResolvedObstacleFlow3DResult:
    if not isinstance(config, ResolvedObstacleFlow3DConfig):
        raise TypeError("config must be ResolvedObstacleFlow3DConfig")
    if face_flux_callback is not None and not callable(face_flux_callback):
        raise TypeError("face_flux_callback must be callable")
    if enthalpy_flux_callback is not None and not callable(enthalpy_flux_callback):
        raise TypeError("enthalpy_flux_callback must be callable")
    if interface_state_callback is not None and not callable(interface_state_callback):
        raise TypeError("interface_state_callback must be callable")
    if (
        interface_conserved_state_callback is not None
        and not callable(interface_conserved_state_callback)
    ):
        raise TypeError("interface_conserved_state_callback must be callable")
    if interface_momentum_callback is not None and not callable(interface_momentum_callback):
        raise TypeError("interface_momentum_callback must be callable")
    if interface_conserved_mass_callback is not None and not callable(interface_conserved_mass_callback):
        raise TypeError("interface_conserved_mass_callback must be callable")
    checked_receptors = tuple(receptors)
    if any(not isinstance(item, FastObstacleReceptor3D) for item in checked_receptors):
        raise TypeError("receptors must contain FastObstacleReceptor3D values")
    if len({item.label for item in checked_receptors}) != len(checked_receptors):
        raise ValueError("receptor labels must be unique")
    if (handoff_request is None) != (handoff_observer is None):
        raise ValueError(
            "handoff_request and handoff_observer must be supplied together"
        )
    if handoff_request is not None and not isinstance(
        handoff_request, BinaryMixtureHandoffPlane3DRequest,
    ):
        raise TypeError(
            "handoff_request must be BinaryMixtureHandoffPlane3DRequest"
        )
    if handoff_observer is not None and not callable(handoff_observer):
        raise TypeError("handoff_observer must be callable")
    if restart_state is not None:
        if not isinstance(restart_state, ResolvedObstacleFlow3DRestartState):
            raise TypeError(
                "restart_state must be ResolvedObstacleFlow3DRestartState"
            )
        if not (
            config.thermal_buoyancy_coupling
            and config.initial_binary_cloud is not None
            and config.binary_mixture_closure is not None
            and config.active_momentum_solver == "finite_volume_projection"
            and config.low_mach_mass_flux_correction
        ):
            raise ValueError(
                "restart currently requires the conservative binary-mixture "
                "finite_volume_projection path"
            )
        if restart_state.configuration_sha256 != _restart_configuration_sha256(
            config
        ):
            raise ValueError(
                "restart configuration SHA-256 does not match the current "
                "solver configuration"
            )
    integration_start_time_s = (
        0.0 if restart_state is None else restart_state.time_s
    )
    integration_stop_time_s = (
        config.duration_s
        if stop_time_s is None
        else _finite(stop_time_s, "stop_time_s")
    )
    if not (
        integration_start_time_s < integration_stop_time_s
        <= config.duration_s + 1.0e-12
    ):
        raise ValueError(
            "stop_time_s must be after the restart time and no later than "
            "config.duration_s"
        )
    integration_duration_s = (
        integration_stop_time_s - integration_start_time_s
    )
    wall_started = time.perf_counter()
    flow_started = time.perf_counter()
    (
        x, y, z, solid, u, v, w, u_face, v_face, w_face,
        pressure, pressure_coefficient, eddy_viscosity,
        residual_history, lattice_density, flow_dt, base_tau, maximum_tau,
        actual_flow_steps, raw_divergence, divergence_before_projection,
        divergence_after_projection, projection_iterations, distribution,
        statistical_residual, average_count,
    ) = _solve_flow(config)
    flow_runtime = time.perf_counter() - flow_started
    dx = config.domain_m[0] / config.cells[0]
    inlet_physical, inlet_lattice = _inlet_profile(config, z)
    active_inlet_physical = _active_inlet_plane(
        config, y, inlet_physical,
    )
    (
        current_inlet_u_physical,
        current_inlet_v_physical,
        current_inlet_w_physical,
    ) = _active_inlet_velocity_planes(
        config, active_inlet_physical, 0.0,
    )
    obstacle_blocked_faces = _circular_arc_blocked_faces(config, x, y, z)
    obstacle_blocked_links = _circular_arc_blocked_links(config, x, y, z)
    zero_thickness_arc_present = (
        config.circular_arc_representation == "zero_thickness_surface"
        and any(isinstance(item, CircularArcWall) for item in config.obstacles)
    )
    embedded_arc_surface_area_m2 = 0.0
    declared_arc_surface_area_m2 = 0.0
    embedded_arc_surface_cell_count = 0
    embedded_arc_minimum_wall_distance_m = math.inf
    embedded_arc_area_by_cell = (
        np.zeros_like(solid, dtype=float) if zero_thickness_arc_present else None
    )
    embedded_arc_area_normal_x_by_cell = (
        np.zeros_like(solid, dtype=float) if zero_thickness_arc_present else None
    )
    embedded_arc_area_normal_y_by_cell = (
        np.zeros_like(solid, dtype=float) if zero_thickness_arc_present else None
    )
    embedded_arc_wall_distance_by_cell = (
        np.full(solid.shape, np.inf, dtype=float)
        if zero_thickness_arc_present else None
    )
    if zero_thickness_arc_present:
        for obstacle in config.obstacles:
            if not isinstance(obstacle, CircularArcWall):
                continue
            cell_area, _area_normal_x, _area_normal_y, _wall_distance = (
                circular_arc_surface_cell_metrics(
                    x, y, z,
                    center_xy_m=obstacle.center_m,
                    radius_m=obstacle.radius_m,
                    z_min_m=obstacle.z_min_m,
                    z_max_m=obstacle.z_max_m,
                    arc_start_deg=obstacle.arc_start_deg,
                    arc_end_deg=obstacle.arc_end_deg,
                )
            )
            embedded_arc_surface_area_m2 += float(np.sum(cell_area))
            assert embedded_arc_area_by_cell is not None
            assert embedded_arc_area_normal_x_by_cell is not None
            assert embedded_arc_area_normal_y_by_cell is not None
            assert embedded_arc_wall_distance_by_cell is not None
            embedded_arc_area_by_cell += cell_area
            embedded_arc_area_normal_x_by_cell += _area_normal_x
            embedded_arc_area_normal_y_by_cell += _area_normal_y
            cut_surface_cells = cell_area > 0.0
            embedded_arc_surface_cell_count += int(np.count_nonzero(cut_surface_cells))
            embedded_arc_minimum_wall_distance_m = min(
                embedded_arc_minimum_wall_distance_m,
                float(np.min(_wall_distance[cut_surface_cells])),
            )
            embedded_arc_wall_distance_by_cell[cut_surface_cells] = np.minimum(
                embedded_arc_wall_distance_by_cell[cut_surface_cells],
                _wall_distance[cut_surface_cells],
            )
            declared_arc_surface_area_m2 += (
                obstacle.radius_m
                * math.radians(obstacle.angular_span_deg)
                * (obstacle.z_max_m - obstacle.z_min_m)
            )
    if config.active_momentum_solver == "finite_volume_projection":
        # The complete z-y mean inlet is the authoritative boundary for
        # momentum, pressure correction and scalar transport, independent of
        # whether the initial field came from D3Q19 or a face projection.
        u_face[:, :, 0] = current_inlet_u_physical
        (
            u_face, v_face, w_face, raw_divergence,
            divergence_before_projection, divergence_after_projection,
            projection_iterations,
        ) = _projected_faces_from_cell_velocity(
            u=u, v=v, w=w, solid=solid,
            inlet_physical=current_inlet_u_physical,
            spacing_m=dx,
            blocked_faces=obstacle_blocked_faces,
            nonstreamwise_wall_boundary=config.nonstreamwise_wall_boundary,
        )
        u = 0.5 * (u_face[:, :, :-1] + u_face[:, :, 1:])
        v = 0.5 * (v_face[:, :-1, :] + v_face[:, 1:, :])
        w = 0.5 * (w_face[:-1, :, :] + w_face[1:, :, :])
        u[solid] = v[solid] = w[solid] = 0.0
    projection_open_faces = (
        ("x_max", "y_min", "y_max", "z_min", "z_max")
        if config.nonstreamwise_wall_boundary == "all_pressure_open"
        else ("x_max", "y_min", "y_max", "z_max")
        if config.nonstreamwise_wall_boundary == "ground_free_slip_lateral_top_pressure_open"
        else ("x_max",)
    )
    background_eddy_viscosity = np.zeros_like(eddy_viscosity)
    if config.background_surface_layer_mixing is not None:
        background_eddy_viscosity = (
            config.background_surface_layer_mixing.eddy_viscosity_m2_s(
                z, solid,
            )
        )
        eddy_viscosity = eddy_viscosity + background_eddy_viscosity
    applied_scalar_diffusivity = _effective_transport_diffusivity(
        molecular_diffusivity_m2_s=config.scalar_diffusivity_m2_s,
        effective_dimensionless_number=config.effective_schmidt_number,
        base_kinematic_viscosity_m2_s=config.kinematic_viscosity_m2_s,
        subgrid_eddy_viscosity_m2_s=eddy_viscosity,
        model=config.scalar_turbulent_diffusivity_model,
    )
    applied_thermal_diffusivity = _effective_transport_diffusivity(
        molecular_diffusivity_m2_s=config.thermal_diffusivity_m2_s,
        effective_dimensionless_number=config.effective_prandtl_number,
        base_kinematic_viscosity_m2_s=config.kinematic_viscosity_m2_s,
        subgrid_eddy_viscosity_m2_s=eddy_viscosity,
        model=config.thermal_turbulent_diffusivity_model,
    )
    (
        obstacle_voxelization,
        maximum_obstacle_face_error,
        maximum_obstacle_size_error,
    ) = _obstacle_voxelization_diagnostics(config, x, y, z)
    cell_volume = dx**3
    maximum_lattice_mach = float(np.max(
        np.sqrt(u[~solid]**2 + v[~solid]**2 + w[~solid]**2)
        * flow_dt / dx / math.sqrt(_CS2_3D)
    ))
    if maximum_lattice_mach > config.maximum_lattice_mach_number:
        raise RuntimeError(
            "resolved 3-D flow exceeded the declared lattice Mach gate: "
            f"{maximum_lattice_mach:.6g} > {config.maximum_lattice_mach_number:.6g}"
        )
    if config.thermal_buoyancy_coupling:
        active_dt = (
            flow_dt
            if config.active_momentum_solver == "d3q19_bgk"
            else config.scalar_time_step_s
        )
        step_ratio = integration_duration_s / active_dt
        steps = int(round(step_ratio))
        if steps < 1 or not math.isclose(
            step_ratio, steps, rel_tol=1.0e-8, abs_tol=1.0e-10,
        ):
            raise ValueError(
                "thermal buoyancy coupling requires duration_s to be an integer "
                f"multiple of the active momentum step {active_dt:.12g} s"
            )
        if config.active_momentum_solver == "d3q19_bgk" and not math.isclose(
            config.scalar_time_step_s, flow_dt, rel_tol=1.0e-8, abs_tol=1.0e-12,
        ):
            raise ValueError(
                "thermal buoyancy coupling requires scalar_time_step_s to equal "
                f"the D3Q19 physical step {flow_dt:.12g} s"
            )
        dt = active_dt
    else:
        steps = int(math.ceil(
            integration_duration_s / config.scalar_time_step_s
        ))
        dt = integration_duration_s / steps
    advective_number = dt * (
        float(np.max(np.abs(u_face))) / dx
        + float(np.max(np.abs(v_face))) / dx
        + float(np.max(np.abs(w_face))) / dx
    )
    diffusion_number = 6.0 * float(
        np.max(np.asarray(applied_scalar_diffusivity))
    ) * dt / dx**2
    thermal_diffusion_number = (
        6.0 * float(np.max(np.asarray(applied_thermal_diffusivity))) * dt / dx**2
        if config.thermal_buoyancy_coupling else 0.0
    )
    if advective_number + max(diffusion_number, thermal_diffusion_number) > 1.0 + 1.0e-12:
        raise ValueError(
            "resolved 3-D scalar stability limit exceeded: "
            f"advection={advective_number:.6g}, diffusion={diffusion_number:.6g}, "
            f"thermal_diffusion={thermal_diffusion_number:.6g}"
        )
    source_sigma_xyz = (
        (config.source_sigma_m, config.source_sigma_m, config.source_sigma_m)
        if config.source_sigma_xyz_m is None
        else config.source_sigma_xyz_m
    )
    if config.source_kernel == "gaussian":
        weights = (
            _source_weights(
                x, y, z, solid, config.source_position_m, config.source_sigma_m,
            )
            if config.source_sigma_xyz_m is None
            else _source_weights_xyz(
                x, y, z, solid, config.source_position_m, source_sigma_xyz,
            )
        )
    elif config.source_kernel == "variance_matched_box":
        weights = _source_weights_variance_matched_box(
            x, y, z, solid, config.source_position_m, source_sigma_xyz,
        )
    else:
        weights = _source_weights_cosine_tapered_box(
            x, y, z, solid, config.source_position_m, source_sigma_xyz,
        )
    stencils = [
        _receptor_stencil(
            x, y, z, solid, receptor,
            vertical_sampling=config.receptor_vertical_sampling,
        )
        for receptor in checked_receptors
    ]
    binary_initial_cloud = config.initial_binary_cloud is not None
    if handoff_request is not None and not binary_initial_cloud:
        raise ValueError(
            "conservative binary handoff recording requires an initial binary cloud"
        )
    concentration = np.zeros_like(solid, dtype=float)
    maximum_concentration = concentration.copy()
    closure = None
    sensible_enthalpy = None
    maximum_mole_fraction = None
    thermodynamic_state = None
    temperature_fields = None
    injected_energy = domain_energy = outflow_energy = energy_residual = None
    cumulative_injected_energy = 0.0
    cumulative_wall_heat_energy_j = 0.0
    maximum_wall_heat_rate_w = 0.0
    cumulative_outflow_energy = 0.0
    source_specific_enthalpy = None
    initial_cloud_record = None
    initial_cloud_maximum_speed_m_s = None
    initial_cloud_maximum_momentum_density_kg_m2_s = None
    if config.thermal_buoyancy_coupling:
        if binary_initial_cloud:
            assert config.binary_mixture_closure is not None
            assert config.initial_binary_cloud is not None
            closure = config.binary_mixture_closure
            initial_cloud_record = initialise_regular_polygon_prism_cloud(
                x, y, z,
                cloud=config.initial_binary_cloud,
                closure=closure,
                solid_mask=solid,
            )
            concentration = (
                initial_cloud_record.species_mass_concentration_kg_m3.copy()
            )
            sensible_enthalpy = (
                initial_cloud_record.sensible_enthalpy_density_j_m3.copy()
            )
            thermodynamic_state = initial_cloud_record.thermodynamic_state
            if (
                config.inlet_temperature_profile_z_m is not None
                and not np.any(concentration > 0.0)
            ):
                initial_temperature = np.broadcast_to(
                    np.interp(
                        z,
                        config.inlet_temperature_profile_z_m,
                        config.inlet_temperature_profile_k,
                    )[:, None, None],
                    concentration.shape,
                ).copy()
                initial_temperature[solid] = config.ambient_temperature_k
                thermodynamic_state = closure.state_from_mole_fraction_and_temperature(
                    np.zeros_like(concentration), initial_temperature,
                )
                sensible_enthalpy = (
                    thermodynamic_state.sensible_enthalpy_density_j_m3.copy()
                )
                sensible_enthalpy[solid] = 0.0
            if config.source_schedule.released_mass_kg > 0.0:
                assert config.source_temperature_k is not None
                species_mass_fraction = (
                    1.0 / config.source_total_mass_to_species_mass_ratio
                )
                source_mixture_heat_capacity = (
                    species_mass_fraction
                    * closure.species_heat_capacity_j_kg_k
                    + (1.0 - species_mass_fraction)
                    * closure.carrier_heat_capacity_j_kg_k
                )
                source_specific_enthalpy = (
                    config.source_total_mass_to_species_mass_ratio
                    * source_mixture_heat_capacity
                    * (config.source_temperature_k - config.ambient_temperature_k)
                )
            maximum_mole_fraction = (
                thermodynamic_state.species_mole_fraction.copy()
            )
            if config.initial_binary_cloud_momentum == "quiescent":
                # The finite cloud replaces ambient gas in proportion to its
                # exact cell-volume fraction.  Its own initial momentum is
                # zero; only the ambient share of a cut cell retains the
                # carrier-flow momentum.  Dividing that retained momentum by
                # the reconstructed mixture density gives a conservative
                # cell velocity without a concentration-derived tuning term.
                cloud_fraction = initial_cloud_record.volume_fraction
                mixture_density = thermodynamic_state.mixture_density_kg_m3
                retained_ambient_momentum_fraction = (
                    (1.0 - cloud_fraction)
                    * closure.ambient_density_kg_m3
                    / mixture_density
                )
                u *= retained_ambient_momentum_fraction
                v *= retained_ambient_momentum_fraction
                w *= retained_ambient_momentum_fraction
                u[solid] = v[solid] = w[solid] = 0.0
                (
                    u_face, v_face, w_face, raw_divergence,
                ) = _raw_faces_from_cell_velocity(
                    u=u, v=v, w=w, solid=solid,
                    inlet_physical=current_inlet_u_physical,
                    spacing_m=dx,
                    blocked_faces=obstacle_blocked_faces,
                    nonstreamwise_wall_boundary=(
                        config.nonstreamwise_wall_boundary
                    ),
                )
                cloud_cells = cloud_fraction > 1.0e-12
                initial_cloud_maximum_speed_m_s = float(np.max(np.sqrt(
                    u[cloud_cells] ** 2
                    + v[cloud_cells] ** 2
                    + w[cloud_cells] ** 2
                )))
                initial_cloud_maximum_momentum_density_kg_m2_s = float(
                    np.max(
                        mixture_density[cloud_cells]
                        * np.sqrt(
                            u[cloud_cells] ** 2
                            + v[cloud_cells] ** 2
                            + w[cloud_cells] ** 2
                        )
                    )
                )
            else:
                cloud_cells = initial_cloud_record.volume_fraction > 1.0e-12
                initial_cloud_maximum_speed_m_s = float(np.max(np.sqrt(
                    u[cloud_cells] ** 2
                    + v[cloud_cells] ** 2
                    + w[cloud_cells] ** 2
                )))
                initial_cloud_maximum_momentum_density_kg_m2_s = float(
                    np.max(
                        thermodynamic_state.mixture_density_kg_m3[cloud_cells]
                        * np.sqrt(
                            u[cloud_cells] ** 2
                            + v[cloud_cells] ** 2
                            + w[cloud_cells] ** 2
                        )
                    )
                )
        else:
            closure = HydrogenAirActiveScalarClosure(
                ambient_temperature_k=config.ambient_temperature_k,
                ambient_pressure_pa=config.ambient_pressure_pa,
                ambient_air_density_kg_m3=config.ambient_air_density_kg_m3,
                dry_air_heat_capacity_j_kg_k=config.dry_air_heat_capacity_j_kg_k,
                hydrogen_heat_capacity_j_kg_k=config.hydrogen_heat_capacity_j_kg_k,
                gravity_m_s2=config.gravity_m_s2,
                maximum_relative_density_departure=(
                    config.maximum_active_scalar_density_departure
                ),
            )
            assert config.source_temperature_k is not None
            source_specific_enthalpy = (
                closure.source_specific_sensible_enthalpy_j_kg(
                    config.source_temperature_k
                )
            )
            sensible_enthalpy = np.zeros_like(concentration)
            if config.inlet_temperature_profile_z_m is not None:
                profile_temperature = np.interp(
                    z,
                    config.inlet_temperature_profile_z_m,
                    config.inlet_temperature_profile_k,
                )
                sensible_enthalpy = np.broadcast_to(
                    (
                        config.ambient_air_density_kg_m3
                        * config.dry_air_heat_capacity_j_kg_k
                        * (profile_temperature - config.ambient_temperature_k)
                    )[:, None, None],
                    concentration.shape,
                ).copy()
                sensible_enthalpy[solid] = 0.0
            thermodynamic_state = closure.state_from_conserved(
                concentration, sensible_enthalpy,
            )
            maximum_mole_fraction = (
                thermodynamic_state.hydrogen_mole_fraction.copy()
            )
        inlet_sensible_enthalpy_profile = None
        if config.inlet_temperature_profile_z_m is not None:
            inlet_temperature = np.interp(
                z,
                config.inlet_temperature_profile_z_m,
                config.inlet_temperature_profile_k,
            )
            inlet_sensible_enthalpy_profile = np.broadcast_to(
                (
                    config.ambient_air_density_kg_m3
                    * config.dry_air_heat_capacity_j_kg_k
                    * (inlet_temperature - config.ambient_temperature_k)
                )[:, None],
                (len(z), len(y)),
            ).copy()
        injected_energy = [0.0]
        domain_energy = [float(np.sum(sensible_enthalpy) * cell_volume)]
        outflow_energy = [0.0]
        energy_residual = [0.0]
        if config.store_scalar_fields:
            temperature_fields = [thermodynamic_state.temperature_k.copy()]
    restart_time_s = 0.0
    restart_projection_potential = None
    if restart_state is not None:
        expected_shape = solid.shape
        if restart_state.concentration_kg_m3.shape != expected_shape:
            raise ValueError(
                "restart cell shape does not match the configured grid"
            )
        assert isinstance(closure, BinaryIdealGasMixtureClosure)
        restart_time_s = restart_state.time_s
        concentration = restart_state.concentration_kg_m3.copy()
        sensible_enthalpy = (
            restart_state.sensible_enthalpy_density_j_m3.copy()
        )
        maximum_concentration = (
            restart_state.maximum_concentration_kg_m3.copy()
        )
        maximum_mole_fraction = (
            restart_state.maximum_species_mole_fraction.copy()
        )
        u = restart_state.velocity_u_m_s.copy()
        v = restart_state.velocity_v_m_s.copy()
        w = restart_state.velocity_w_m_s.copy()
        u_face = restart_state.u_face_m_s.copy()
        v_face = restart_state.v_face_m_s.copy()
        w_face = restart_state.w_face_m_s.copy()
        restart_projection_potential = (
            None
            if restart_state.projection_correction_potential_kg_m_s is None
            else restart_state.projection_correction_potential_kg_m_s.copy()
        )
        thermodynamic_state = (
            closure
            .state_from_species_mass_concentration_and_sensible_enthalpy(
                concentration, sensible_enthalpy,
            )
        )
        (
            current_inlet_u_physical,
            current_inlet_v_physical,
            current_inlet_w_physical,
        ) = _active_inlet_velocity_planes(
            config, active_inlet_physical, restart_time_s,
        )
        u_face[:, :, 0] = current_inlet_u_physical
        eddy_viscosity = subgrid_eddy_viscosity_physical(
            u, v, w, solid,
            grid_spacing_m=dx,
            model=config.subgrid_model,
            smagorinsky_constant=config.smagorinsky_constant,
            wale_constant=config.wale_constant,
        ) + background_eddy_viscosity
    else:
        # A finite initial cloud is part of the maximum-field history at t=0.
        maximum_concentration = concentration.copy()
    maximum_gradient_richardson_number = 0.0
    if config.stratification_damping_alpha is not None:
        if thermodynamic_state is None:
            raise RuntimeError("Ri_g damping requires an initialized density field")
        (
            applied_scalar_diffusivity,
            applied_thermal_diffusivity,
            maximum_gradient_richardson_number,
        ) = _apply_vertical_stratification_sensitivity(
            applied_scalar_diffusivity,
            applied_thermal_diffusivity,
            molecular_scalar_diffusivity_m2_s=config.scalar_diffusivity_m2_s,
            molecular_thermal_diffusivity_m2_s=config.thermal_diffusivity_m2_s,
            density_kg_m3=thermodynamic_state.mixture_density_kg_m3,
            u_m_s=u,
            v_m_s=v,
            solid_mask=solid,
            grid_spacing_m=dx,
            config=config,
        )
    times = [restart_time_s]
    initial_species_mass = float(np.sum(concentration) * cell_volume)
    initial_sensible_enthalpy = (
        0.0 if sensible_enthalpy is None
        else float(np.sum(sensible_enthalpy) * cell_volume)
    )
    injected = [0.0]
    domain_mass = [initial_species_mass]
    outflow = [0.0]
    residual = [0.0]
    traces = {
        receptor.label: [max(0.0, float(sum(
            weight * concentration[k, j, i]
            for k, j, i, weight in stencil
        )))]
        for receptor, stencil in zip(checked_receptors, stencils)
    }
    initial_mole_fraction = None
    if thermodynamic_state is not None:
        initial_mole_fraction = (
            thermodynamic_state.species_mole_fraction
            if binary_initial_cloud
            else thermodynamic_state.hydrogen_mole_fraction
        )
    vol_percent_traces = {
        receptor.label: [
            0.0 if initial_mole_fraction is None else 100.0 * float(np.clip(
                sum(
                    weight * initial_mole_fraction[k, j, i]
                    for k, j, i, weight in stencil
                ), 0.0, 1.0,
            ))
        ]
        for receptor, stencil in zip(checked_receptors, stencils)
    }
    temperature_traces = {
        receptor.label: [
            config.ambient_temperature_k
            if thermodynamic_state is None else float(sum(
                weight * thermodynamic_state.temperature_k[k, j, i]
                for k, j, i, weight in stencil
            ))
        ]
        for receptor, stencil in zip(checked_receptors, stencils)
    }
    fields = [concentration.copy()] if config.store_scalar_fields else None
    transient_average_count = 0
    averaged_u = averaged_v = averaged_w = None
    averaged_concentration = averaged_temperature = None
    cumulative_injected = 0.0
    cumulative_outflow = 0.0
    maximum_advective_number = advective_number
    maximum_diffusion_number = diffusion_number
    maximum_thermal_diffusion_number = thermal_diffusion_number
    maximum_buoyancy_acceleration = 0.0
    maximum_embedded_arc_wall_shear_stress = 0.0
    maximum_embedded_arc_wall_friction_velocity = 0.0
    maximum_embedded_arc_wall_tangential_speed = 0.0
    cumulative_embedded_arc_wall_force_n = np.zeros(3, dtype=float)
    maximum_density_departure = 0.0
    active_scalar_out_of_range_cell_steps = 0
    corrector_iterations_per_step: list[int] = []
    projection_linear_iterations_per_step: list[int] = []
    corrector_nonconverged_steps = 0
    maximum_corrector_density_change = 0.0
    maximum_corrector_velocity_change = 0.0
    maximum_variable_density_continuity_residual = 0.0
    maximum_variable_density_global_residual = 0.0
    maximum_momentum_balance_residual = 0.0
    maximum_momentum_advective_courant = 0.0
    maximum_momentum_viscous_stability = 0.0
    maximum_eddy_viscosity_over_transient = float(
        np.max(eddy_viscosity[~solid])
    )
    minimum_transient_inlet_u = float(np.min(current_inlet_u_physical))
    maximum_transient_inlet_u = float(np.max(current_inlet_u_physical))
    maximum_transient_inlet_speed = float(np.max(np.sqrt(
        current_inlet_u_physical**2
        + current_inlet_v_physical**2
        + current_inlet_w_physical**2
    )))
    last_variable_density_solver_method = None
    last_variable_density_preconditioner = None
    last_variable_density_preconditioner_levels = None
    last_variable_density_preconditioner_operator_complexity = None
    last_variable_density_pyamg_version = None
    last_pressure_correction = None
    last_projection_potential = restart_projection_potential
    last_continuity_residual = None
    lattice_viscosity = config.kinematic_viscosity_m2_s * flow_dt / dx**2
    pending_pressure_acceleration_lattice = tuple(
        np.zeros_like(concentration) for _ in range(3)
    )
    handoff_x_face_index = None
    handoff_sample_interval_steps = None
    handoff_frame_count = 0
    if handoff_request is not None:
        x_index_float = handoff_request.x_face_m / dx
        if not math.isclose(
            x_index_float, round(x_index_float), rel_tol=0.0, abs_tol=1.0e-10,
        ):
            raise ValueError(
                "handoff x_face_m must coincide with a solver face"
            )
        handoff_x_face_index = int(round(x_index_float))
        if not 0 < handoff_x_face_index < solid.shape[2]:
            raise ValueError("handoff plane must be an interior solver x face")
        interval_float = handoff_request.sample_interval_s / dt
        if not math.isclose(
            interval_float, round(interval_float), rel_tol=0.0, abs_tol=1.0e-10,
        ):
            raise ValueError(
                "handoff sample_interval_s must be an integer multiple of the active step"
            )
        handoff_sample_interval_steps = int(round(interval_float))
        if handoff_sample_interval_steps < 1:
            raise ValueError("handoff sample interval must contain at least one step")

    def emit_handoff_frame(sample_time_s: float) -> None:
        nonlocal handoff_frame_count
        if handoff_request is None or handoff_observer is None:
            return
        assert handoff_x_face_index is not None
        assert isinstance(closure, BinaryIdealGasMixtureClosure)
        assert thermodynamic_state is not None
        assert sensible_enthalpy is not None
        frame = extract_binary_mixture_handoff_plane_3d(
            thermodynamic_state,
            sensible_enthalpy,
            u_face_m_s=u_face,
            cell_v_m_s=v,
            cell_w_m_s=w,
            solid_mask=solid,
            x_face_index=handoff_x_face_index,
            spacing_m=(dx, dx, dx),
            time_s=sample_time_s,
            time_step_s=dt,
            species_diffusivity_m2_s=applied_scalar_diffusivity,
            thermal_diffusivity_m2_s=applied_thermal_diffusivity,
            ambient_total_mass_flux_kg_m2_s=(
                closure.ambient_density_kg_m3 * current_inlet_u_physical
            ),
            ambient_v_m_s=current_inlet_v_physical,
            ambient_w_m_s=current_inlet_w_physical,
            advection_scheme=config.scalar_advection_scheme,
        )
        handoff_observer(frame)
        handoff_frame_count += 1

    if handoff_request is not None and handoff_request.include_initial:
        emit_handoff_frame(restart_time_s)
    transport_started = time.perf_counter()
    for step in range(steps):
        start = restart_time_s + step * dt
        end = restart_time_s + (step + 1) * dt
        (
            current_inlet_u_physical,
            current_inlet_v_physical,
            current_inlet_w_physical,
        ) = _active_inlet_velocity_planes(
            config, active_inlet_physical, 0.5 * (start + end),
        )
        minimum_transient_inlet_u = min(
            minimum_transient_inlet_u,
            float(np.min(current_inlet_u_physical)),
        )
        maximum_transient_inlet_u = max(
            maximum_transient_inlet_u,
            float(np.max(current_inlet_u_physical)),
        )
        maximum_transient_inlet_speed = max(
            maximum_transient_inlet_speed,
            float(np.max(np.sqrt(
                current_inlet_u_physical**2
                + current_inlet_v_physical**2
                + current_inlet_w_physical**2
            ))),
        )
        source_mass = config.source_schedule.mass_between(start, end)
        source_hydrogen_increment = source_mass * weights / cell_volume
        source_total_mass_increment = (
            source_hydrogen_increment
            * config.source_total_mass_to_species_mass_ratio
        )
        source_energy_increment = np.zeros_like(concentration)
        if sensible_enthalpy is not None and (
            not binary_initial_cloud or source_mass > 0.0
        ):
            assert source_specific_enthalpy is not None
            source_energy_increment = (
                source_hydrogen_increment * source_specific_enthalpy
            )
        corrector_transport_done = False
        source_pre_fraction = (
            0.5 if config.source_time_integration == "strang" else 1.0
        )
        source_post_fraction = 1.0 - source_pre_fraction
        if not (
            sensible_enthalpy is not None
            and config.low_mach_mass_flux_correction
        ):
            concentration += source_pre_fraction * source_hydrogen_increment
            if sensible_enthalpy is not None:
                sensible_enthalpy += source_pre_fraction * source_energy_increment
        cumulative_injected += source_mass
        if interface_conserved_state_callback is not None:
            updated_state = interface_conserved_state_callback(
                concentration.copy(),
                None if sensible_enthalpy is None else sensible_enthalpy.copy(),
                start,
                dt,
            )
            if not isinstance(updated_state, tuple) or len(updated_state) != 2:
                raise TypeError(
                    "interface_conserved_state_callback must return "
                    "(concentration, sensible_enthalpy)"
                )
            updated_concentration = np.asarray(updated_state[0], dtype=float)
            if (
                updated_concentration.shape != concentration.shape
                or not np.all(np.isfinite(updated_concentration))
                or np.any(updated_concentration < 0.0)
            ):
                raise ValueError(
                    "interface conserved concentration must be finite, "
                    "non-negative, and solver-shaped"
                )
            concentration = updated_concentration
            if sensible_enthalpy is None:
                if updated_state[1] is not None:
                    raise ValueError(
                        "non-thermal resolved runs must return None enthalpy"
                    )
            else:
                updated_enthalpy = np.asarray(updated_state[1], dtype=float)
                if (
                    updated_enthalpy.shape != sensible_enthalpy.shape
                    or not np.all(np.isfinite(updated_enthalpy))
                ):
                    raise ValueError(
                        "interface sensible enthalpy must be finite and "
                        "solver-shaped"
                    )
                sensible_enthalpy = updated_enthalpy
        if interface_conserved_mass_callback is not None:
            callback_mass, callback_enthalpy = interface_conserved_mass_callback(start, dt)
            callback_mass = float(callback_mass)
            callback_enthalpy = float(callback_enthalpy)
            if (
                not math.isfinite(callback_mass) or callback_mass < 0.0
                or not math.isfinite(callback_enthalpy)
            ):
                raise ValueError(
                    "interface conserved mass callback must return finite mass "
                    "and enthalpy"
                )
            cumulative_injected += callback_mass
            cumulative_injected_energy += callback_enthalpy
        if sensible_enthalpy is not None:
            assert closure is not None
            if binary_initial_cloud:
                forcing_state = thermodynamic_state
                assert forcing_state is not None
                if source_specific_enthalpy is not None:
                    cumulative_injected_energy += (
                        source_mass * source_specific_enthalpy
                    )
            else:
                assert source_specific_enthalpy is not None
                cumulative_injected_energy += source_mass * source_specific_enthalpy
                if config.low_mach_mass_flux_correction:
                    forcing_state = closure.state_from_conserved(
                        concentration
                        + source_pre_fraction * source_hydrogen_increment,
                        sensible_enthalpy
                        + source_pre_fraction * source_energy_increment,
                    )
                else:
                    forcing_state = closure.state_from_conserved(
                        concentration, sensible_enthalpy,
                    )
            buoyancy_physical = forcing_state.buoyancy_acceleration_m_s2.copy()
            buoyancy_physical[solid] = 0.0
            maximum_buoyancy_acceleration = max(
                maximum_buoyancy_acceleration,
                float(np.max(np.abs(buoyancy_physical[~solid]))),
            )
            velocity_scale = dx / flow_dt
            acceleration_lattice = tuple(
                np.zeros_like(concentration) for _ in range(3)
            )
            if config.active_momentum_solver == "finite_volume_projection":
                assert thermodynamic_state is not None
                old_u, old_v, old_w = u.copy(), v.copy(), w.copy()
                old_density = thermodynamic_state.mixture_density_kg_m3.copy()
                old_density[solid] = config.ambient_air_density_kg_m3
                rho_face_x, rho_face_y, rho_face_z = _face_density(
                    old_density, solid,
                )
                wall_force_x = np.zeros_like(old_density)
                wall_force_y = np.zeros_like(old_density)
                wall_force_z = np.zeros_like(old_density)
                if (
                    embedded_arc_area_by_cell is not None
                    and embedded_arc_area_normal_x_by_cell is not None
                    and embedded_arc_area_normal_y_by_cell is not None
                    and embedded_arc_wall_distance_by_cell is not None
                ):
                    wall_force, wall_diagnostics = (
                        smooth_wall_log_law_force_density_3d(
                            u, v, w, old_density,
                            surface_area_m2=embedded_arc_area_by_cell,
                            area_normal_x_m2=embedded_arc_area_normal_x_by_cell,
                            area_normal_y_m2=embedded_arc_area_normal_y_by_cell,
                            wall_distance_m=embedded_arc_wall_distance_by_cell,
                            spacing_m=(dx, dx, dx),
                            molecular_kinematic_viscosity_m2_s=(
                                config.kinematic_viscosity_m2_s
                            ),
                        )
                    )
                    wall_force_x, wall_force_y, wall_force_z = wall_force
                    maximum_embedded_arc_wall_shear_stress = max(
                        maximum_embedded_arc_wall_shear_stress,
                        float(wall_diagnostics["maximum_wall_shear_stress_pa"]),
                    )
                    maximum_embedded_arc_wall_friction_velocity = max(
                        maximum_embedded_arc_wall_friction_velocity,
                        float(wall_diagnostics["maximum_friction_velocity_m_s"]),
                    )
                    maximum_embedded_arc_wall_tangential_speed = max(
                        maximum_embedded_arc_wall_tangential_speed,
                        float(wall_diagnostics["maximum_tangential_speed_m_s"]),
                    )
                    cumulative_embedded_arc_wall_force_n += np.asarray((
                        wall_diagnostics["net_wall_force_x_n"],
                        wall_diagnostics["net_wall_force_y_n"],
                        wall_diagnostics["net_wall_force_z_n"],
                    )) * dt
                source_momentum_velocity = config.source_momentum_velocity_m_s
                if config.source_momentum_velocity_history is not None:
                    source_momentum_velocity = (
                        config.source_momentum_velocity_history.velocity_at(
                            0.5 * (start + end)
                        )
                    )
                momentum = advance_variable_density_momentum_3d(
                    old_density, u, v, w,
                    mass_flux_x_kg_m2_s=rho_face_x * u_face,
                    mass_flux_y_kg_m2_s=rho_face_y * v_face,
                    mass_flux_z_kg_m2_s=rho_face_z * w_face,
                    solid_mask=solid,
                    spacing_m=(dx, dx, dx),
                    time_step_s=dt,
                    kinematic_viscosity_m2_s=(
                        config.kinematic_viscosity_m2_s + eddy_viscosity
                    ),
                    acceleration_x_m_s2=wall_force_x / old_density,
                    acceleration_y_m_s2=wall_force_y / old_density,
                    acceleration_z_m_s2=(
                        buoyancy_physical + wall_force_z / old_density
                    ),
                    mass_source_rate_kg_m3_s=(
                        source_total_mass_increment / dt
                    ),
                    source_velocity_m_s=source_momentum_velocity,
                    inlet_u_m_s=current_inlet_u_physical,
                    inlet_v_m_s=current_inlet_v_physical,
                    inlet_w_m_s=current_inlet_w_physical,
                    nonstreamwise_wall_boundary=(
                        config.nonstreamwise_wall_boundary
                    ),
                )
                u = momentum.provisional_u_m_s
                v = momentum.provisional_v_m_s
                w = momentum.provisional_w_m_s
                u_face, v_face, w_face, raw_divergence = (
                    _increment_projected_faces_from_cell_delta(
                        previous_u_face=u_face,
                        previous_v_face=v_face,
                        previous_w_face=w_face,
                        delta_u=u - old_u,
                        delta_v=v - old_v,
                        delta_w=w - old_w,
                        solid=solid,
                        spacing_m=dx,
                        inlet_physical=current_inlet_u_physical,
                        blocked_faces=obstacle_blocked_faces,
                        nonstreamwise_wall_boundary=(
                            config.nonstreamwise_wall_boundary
                        ),
                    )
                )
                maximum_momentum_balance_residual = max(
                    maximum_momentum_balance_residual,
                    float(np.max(np.abs(
                        momentum.momentum_balance_residual_kg_m_s
                    ))),
                )
                maximum_momentum_advective_courant = max(
                    maximum_momentum_advective_courant,
                    momentum.advective_courant_number,
                )
                maximum_momentum_viscous_stability = max(
                    maximum_momentum_viscous_stability,
                    momentum.viscous_stability_number,
                )
                ux, uy, uz = (
                    u / velocity_scale,
                    v / velocity_scale,
                    w / velocity_scale,
                )
                eddy_viscosity = subgrid_eddy_viscosity_physical(
                    u, v, w, solid,
                    grid_spacing_m=dx,
                    model=config.subgrid_model,
                    smagorinsky_constant=config.smagorinsky_constant,
                    wale_constant=config.wale_constant,
                )
                eddy_viscosity = (
                    eddy_viscosity + background_eddy_viscosity
                )
            else:
                acceleration_lattice = (
                    pending_pressure_acceleration_lattice[0],
                    pending_pressure_acceleration_lattice[1],
                    buoyancy_physical * flow_dt**2 / dx
                    + pending_pressure_acceleration_lattice[2],
                )
                pending_pressure_acceleration_lattice = tuple(
                    np.zeros_like(concentration) for _ in range(3)
                )
                distribution, tau, eddy_viscosity_lattice = _advance_d3q19(
                    distribution,
                    solid=solid,
                    inlet_u_lattice=inlet_lattice,
                    lattice_viscosity=lattice_viscosity,
                    flow_boundary_scheme=config.flow_boundary_scheme,
                    outlet_lattice_density=config.outlet_lattice_density,
                    nonstreamwise_wall_boundary=(
                        config.nonstreamwise_wall_boundary
                    ),
                    subgrid_model=config.subgrid_model,
                    smagorinsky_constant=config.smagorinsky_constant,
                    wale_constant=config.wale_constant,
                    maximum_relaxation_time=config.maximum_relaxation_time,
                    acceleration_lattice=acceleration_lattice,
                    blocked_links=obstacle_blocked_links,
                )
                maximum_tau = max(maximum_tau, float(np.max(tau[~solid])))
                lattice_density, ux, uy, uz = _moments_d3q19(
                    distribution, solid,
                    acceleration_lattice,
                )
                u, v, w = (
                    ux * velocity_scale,
                    uy * velocity_scale,
                    uz * velocity_scale,
                )
                u[solid] = v[solid] = w[solid] = 0.0
                eddy_viscosity = eddy_viscosity_lattice * dx**2 / flow_dt
            maximum_eddy_viscosity_over_transient = max(
                maximum_eddy_viscosity_over_transient,
                float(np.max(eddy_viscosity[~solid])),
            )
            applied_scalar_diffusivity = _effective_transport_diffusivity(
                molecular_diffusivity_m2_s=config.scalar_diffusivity_m2_s,
                effective_dimensionless_number=config.effective_schmidt_number,
                base_kinematic_viscosity_m2_s=config.kinematic_viscosity_m2_s,
                subgrid_eddy_viscosity_m2_s=eddy_viscosity,
                model=config.scalar_turbulent_diffusivity_model,
            )
            applied_thermal_diffusivity = _effective_transport_diffusivity(
                molecular_diffusivity_m2_s=config.thermal_diffusivity_m2_s,
                effective_dimensionless_number=config.effective_prandtl_number,
                base_kinematic_viscosity_m2_s=config.kinematic_viscosity_m2_s,
                subgrid_eddy_viscosity_m2_s=eddy_viscosity,
                model=config.thermal_turbulent_diffusivity_model,
            )
            if config.stratification_damping_alpha is not None:
                if thermodynamic_state is None:
                    raise RuntimeError(
                        "Ri_g damping requires an initialized density field"
                    )
                (
                    applied_scalar_diffusivity,
                    applied_thermal_diffusivity,
                    step_gradient_richardson_maximum,
                ) = _apply_vertical_stratification_sensitivity(
                    applied_scalar_diffusivity,
                    applied_thermal_diffusivity,
                    molecular_scalar_diffusivity_m2_s=(
                        config.scalar_diffusivity_m2_s
                    ),
                    molecular_thermal_diffusivity_m2_s=(
                        config.thermal_diffusivity_m2_s
                    ),
                    density_kg_m3=thermodynamic_state.mixture_density_kg_m3,
                    u_m_s=u,
                    v_m_s=v,
                    solid_mask=solid,
                    grid_spacing_m=dx,
                    config=config,
                )
                maximum_gradient_richardson_number = max(
                    maximum_gradient_richardson_number,
                    step_gradient_richardson_maximum,
                )
            current_diffusion_number = 6.0 * float(np.max(
                np.asarray(applied_scalar_diffusivity)
            )) * dt / dx**2
            current_thermal_diffusion_number = 6.0 * float(np.max(
                np.asarray(applied_thermal_diffusivity)
            )) * dt / dx**2
            maximum_diffusion_number = max(
                maximum_diffusion_number, current_diffusion_number,
            )
            maximum_thermal_diffusion_number = max(
                maximum_thermal_diffusion_number,
                current_thermal_diffusion_number,
            )
            if (
                config.low_mach_mass_flux_correction
                and config.active_momentum_solver
                == "finite_volume_projection"
            ):
                raw_maximum = float(np.max(np.abs(raw_divergence[~solid])))
                divergence_before_projection = raw_maximum
                divergence_after_projection = raw_maximum
                projection_iterations = 0
            elif (
                config.low_mach_mass_flux_correction
                and config.active_provisional_face_projection
                == "variable_density_only"
            ):
                u_face, v_face, w_face, raw_divergence = (
                    _raw_faces_from_cell_velocity(
                        u=u, v=v, w=w, solid=solid,
                        inlet_physical=current_inlet_u_physical, spacing_m=dx,
                        blocked_faces=obstacle_blocked_faces,
                        nonstreamwise_wall_boundary=(
                            config.nonstreamwise_wall_boundary
                        ),
                    )
                )
                raw_maximum = float(np.max(np.abs(raw_divergence[~solid])))
                divergence_before_projection = raw_maximum
                divergence_after_projection = raw_maximum
                projection_iterations = 0
            else:
                (
                    u_face, v_face, w_face, raw_divergence,
                    divergence_before_projection, divergence_after_projection,
                    projection_iterations,
                ) = _projected_faces_from_cell_velocity(
                    u=u, v=v, w=w, solid=solid,
                    inlet_physical=current_inlet_u_physical, spacing_m=dx,
                    blocked_faces=obstacle_blocked_faces,
                    nonstreamwise_wall_boundary=(
                        config.nonstreamwise_wall_boundary
                    ),
                )
            wall_heat_increment_j = 0.0
            if (
                config.obstacle_wall_temperature_k is not None
                or config.ground_temperature_k is not None
            ):
                assert thermodynamic_state is not None
                density_for_heat = thermodynamic_state.mixture_density_kg_m3
                if binary_initial_cloud:
                    heat_capacity_for_heat = (
                        thermodynamic_state.mixture_heat_capacity_j_kg_k
                    )
                else:
                    volumetric_heat_capacity = (
                        config.ambient_air_density_kg_m3
                        * config.dry_air_heat_capacity_j_kg_k
                        + concentration * config.hydrogen_heat_capacity_j_kg_k
                    )
                    heat_capacity_for_heat = (
                        volumetric_heat_capacity
                        / np.maximum(density_for_heat, 1.0e-12)
                    )
                wall_heat = isothermal_wall_heat_source_3d(
                    thermodynamic_state.temperature_k,
                    density_for_heat,
                    heat_capacity_for_heat,
                    applied_thermal_diffusivity,
                    solid,
                    spacing_m=(dx, dx, dx),
                    obstacle_wall_temperature_k=config.obstacle_wall_temperature_k,
                    ground_temperature_k=config.ground_temperature_k,
                )
                source_energy_increment += (
                    wall_heat.sensible_enthalpy_source_rate_w_m3 * dt
                )
                wall_heat_increment_j = wall_heat.wall_heat_rate_w * dt
                cumulative_wall_heat_energy_j += wall_heat_increment_j
                maximum_wall_heat_rate_w = max(
                    maximum_wall_heat_rate_w, abs(wall_heat.wall_heat_rate_w),
                )
            if config.low_mach_mass_flux_correction:
                if binary_initial_cloud:
                    assert isinstance(closure, BinaryIdealGasMixtureClosure)
                    binary_inlet_fluxes = {}
                    if config.inlet_temperature_profile_z_m is not None:
                        inlet_temperature = np.interp(
                            z,
                            config.inlet_temperature_profile_z_m,
                            config.inlet_temperature_profile_k,
                        )
                        inlet_state = closure.state_from_mole_fraction_and_temperature(
                            np.zeros_like(inlet_temperature), inlet_temperature,
                        )
                        binary_inlet_fluxes = {
                            "inlet_total_mass_flux_kg_m2_s": (
                                inlet_state.mixture_density_kg_m3[:, None]
                                * current_inlet_u_physical
                            ),
                            "inlet_species_mass_flux_kg_m2_s": np.zeros_like(
                                current_inlet_u_physical
                            ),
                            "inlet_sensible_enthalpy_flux_w_m2": (
                                inlet_state.sensible_enthalpy_density_j_m3[:, None]
                                * current_inlet_u_physical
                            ),
                        }
                    corrector = solve_binary_mixture_initial_cloud_step_3d(
                        concentration,
                        sensible_enthalpy,
                        provisional_u_face_m_s=u_face,
                        provisional_v_face_m_s=v_face,
                        provisional_w_face_m_s=w_face,
                        solid_mask=solid,
                        spacing_m=(dx, dx, dx),
                        time_step_s=dt,
                        species_diffusivity_m2_s=applied_scalar_diffusivity,
                        thermal_diffusivity_m2_s=applied_thermal_diffusivity,
                        closure=closure,
                        advection_scheme=config.scalar_advection_scheme,
                        source_time_integration=config.source_time_integration,
                        maximum_iterations=(
                            config.thermodynamic_corrector_maximum_iterations
                        ),
                        minimum_iterations=(
                            config.thermodynamic_corrector_minimum_iterations
                        ),
                        relative_density_tolerance=(
                            config.thermodynamic_corrector_density_tolerance
                        ),
                        relative_mass_flux_tolerance=(
                            config.thermodynamic_corrector_velocity_tolerance
                        ),
                        open_boundary_faces=projection_open_faces,
                        projection_linear_solver=(
                            config.variable_density_projection_linear_solver
                        ),
                        initial_correction_potential_kg_m_s=(
                            last_projection_potential
                        ),
                        species_source_increment_kg_m3=(
                            source_hydrogen_increment
                        ),
                        total_mass_source_increment_kg_m3=(
                            source_total_mass_increment
                        ),
                        sensible_enthalpy_source_increment_j_m3=(
                            source_energy_increment
                        ),
                        **binary_inlet_fluxes,
                        blocked_faces=obstacle_blocked_faces,
                    )
                else:
                    assert isinstance(closure, HydrogenAirActiveScalarClosure)
                    corrector = solve_low_mach_thermal_corrector_3d(
                        concentration,
                        sensible_enthalpy,
                        hydrogen_source_increment_kg_m3=source_hydrogen_increment,
                        sensible_enthalpy_source_increment_j_m3=(
                            source_energy_increment
                        ),
                        provisional_u_face_m_s=u_face,
                        provisional_v_face_m_s=v_face,
                        provisional_w_face_m_s=w_face,
                        solid_mask=solid,
                        spacing_m=(dx, dx, dx),
                        time_step_s=dt,
                        species_diffusivity_m2_s=applied_scalar_diffusivity,
                        thermal_diffusivity_m2_s=applied_thermal_diffusivity,
                        closure=closure,
                        advection_scheme=config.scalar_advection_scheme,
                        source_time_integration=config.source_time_integration,
                        inlet_sensible_enthalpy_density_j_m3=(
                            inlet_sensible_enthalpy_profile
                        ),
                        maximum_iterations=(
                            config.thermodynamic_corrector_maximum_iterations
                        ),
                        minimum_iterations=(
                            config.thermodynamic_corrector_minimum_iterations
                        ),
                        relative_density_tolerance=(
                            config.thermodynamic_corrector_density_tolerance
                        ),
                        relative_velocity_tolerance=(
                            config.thermodynamic_corrector_velocity_tolerance
                        ),
                        open_boundary_faces=projection_open_faces,
                        initial_correction_potential_kg_m_s=(
                            last_projection_potential
                        ),
                        projection_linear_solver=(
                            config.variable_density_projection_linear_solver
                        ),
                        blocked_faces=obstacle_blocked_faces,
                        face_flux_callback=face_flux_callback,
                        enthalpy_flux_callback=enthalpy_flux_callback,
                    )
                if (
                    config.require_thermodynamic_corrector_convergence
                    and not corrector.converged
                ):
                    raise RuntimeError(
                        "thermodynamic mass-flux corrector did not converge in "
                        f"{corrector.iterations} iterations"
                    )
                corrector_iterations_per_step.append(corrector.iterations)
                projection_linear_iterations_per_step.append(int(
                    np.sum(corrector.projection_iterations_history)
                ))
                if not corrector.converged:
                    corrector_nonconverged_steps += 1
                final_corrector_change = corrector.iteration_history[-1]
                maximum_corrector_density_change = max(
                    maximum_corrector_density_change,
                    float(final_corrector_change[0]),
                )
                maximum_corrector_velocity_change = max(
                    maximum_corrector_velocity_change,
                    float(final_corrector_change[1]),
                )
                maximum_variable_density_continuity_residual = max(
                    maximum_variable_density_continuity_residual,
                    float(final_corrector_change[2]),
                )
                maximum_variable_density_global_residual = max(
                    maximum_variable_density_global_residual,
                    abs(
                        corrector.mass_projection
                        .global_mass_balance_residual_kg_s
                    ),
                )
                last_variable_density_solver_method = (
                    corrector.mass_projection.solver_method
                )
                last_variable_density_preconditioner = (
                    corrector.mass_projection.preconditioner
                )
                last_variable_density_preconditioner_levels = (
                    corrector.mass_projection.preconditioner_level_count
                )
                last_variable_density_preconditioner_operator_complexity = (
                    corrector.mass_projection
                    .preconditioner_operator_complexity
                )
                last_variable_density_pyamg_version = (
                    corrector.mass_projection.pyamg_version
                )
                concentration = (
                    corrector.species_mass_concentration_kg_m3
                    if binary_initial_cloud
                    else corrector.hydrogen_mass_concentration_kg_m3
                )
                sensible_enthalpy = (
                    corrector.sensible_enthalpy_density_j_m3
                )
                thermodynamic_state = corrector.thermodynamic_state
                u_face = corrector.u_face_m_s
                v_face = corrector.v_face_m_s
                w_face = corrector.w_face_m_s
                step_outflow = (
                    corrector.species_boundary_outflow_mass_kg
                    if binary_initial_cloud
                    else corrector.hydrogen_outflow_mass_kg
                )
                step_energy_outflow = (
                    corrector.sensible_enthalpy_boundary_outflow_j
                    if binary_initial_cloud
                    else corrector.sensible_enthalpy_outflow_j
                )
                last_pressure_correction = (
                    corrector.mass_projection.pressure_correction_pa.copy()
                    if binary_initial_cloud
                    else corrector.pressure_correction_pa.copy()
                )
                last_projection_potential = (
                    corrector.mass_projection
                    .correction_potential_kg_m_s.copy()
                )
                last_continuity_residual = (
                    corrector.mass_projection
                    .continuity_residual_after_kg_m3_s.copy()
                )

                # Couple the pressure projection back to D3Q19.  The historical
                # equilibrium replacement is retained as the default candidate.
                # A deferred Guo impulse provides a force-consistent alternative;
                # transport_only is a diagnostic that deliberately omits momentum
                # feedback while retaining conservative variable-density fluxes.
                corrected_u = 0.5 * (u_face[:, :, :-1] + u_face[:, :, 1:])
                corrected_v = 0.5 * (v_face[:, :-1, :] + v_face[:, 1:, :])
                corrected_w = 0.5 * (w_face[:-1, :, :] + w_face[1:, :, :])
                corrected_u[solid] = corrected_v[solid] = corrected_w[solid] = 0.0
                coupling = config.pressure_correction_momentum_coupling
                if (
                    config.active_momentum_solver == "d3q19_bgk"
                    and coupling == "equilibrium_replacement"
                ):
                    _, raw_ux, raw_uy, raw_uz = _moments_d3q19(
                        distribution, solid,
                    )
                    target_raw_ux = (
                        corrected_u / velocity_scale
                        - 0.5 * acceleration_lattice[0]
                    )
                    target_raw_uy = (
                        corrected_v / velocity_scale
                        - 0.5 * acceleration_lattice[1]
                    )
                    target_raw_uz = (
                        corrected_w / velocity_scale
                        - 0.5 * acceleration_lattice[2]
                    )
                    target_raw_ux[solid] = 0.0
                    target_raw_uy[solid] = 0.0
                    target_raw_uz[solid] = 0.0
                    distribution += _equilibrium_d3q19(
                        lattice_density,
                        target_raw_ux, target_raw_uy, target_raw_uz,
                    ) - _equilibrium_d3q19(
                        lattice_density, raw_ux, raw_uy, raw_uz,
                    )
                    if not np.all(np.isfinite(distribution)):
                        raise RuntimeError(
                            "pressure-corrected D3Q19 population became non-finite"
                        )
                elif (
                    config.active_momentum_solver == "d3q19_bgk"
                    and coupling == "deferred_guo_force"
                ):
                    pending_pressure_acceleration_lattice = (
                        (corrected_u - u) * flow_dt / dx,
                        (corrected_v - v) * flow_dt / dx,
                        (corrected_w - w) * flow_dt / dx,
                    )
                u, v, w = corrected_u, corrected_v, corrected_w
                ux, uy, uz = (
                    u / velocity_scale,
                    v / velocity_scale,
                    w / velocity_scale,
                )
                corrector_transport_done = True

        velocity_scale = dx / flow_dt
        if interface_momentum_callback is not None:
            updated_momentum = interface_momentum_callback(
                u.copy(), v.copy(), w.copy(),
                u_face.copy(), v_face.copy(), w_face.copy(),
                start,
                dt,
            )
            if not isinstance(updated_momentum, tuple) or len(updated_momentum) != 6:
                raise TypeError(
                    "interface_momentum_callback must return six velocity arrays"
                )
            arrays = tuple(np.asarray(value, dtype=float) for value in updated_momentum)
            expected = (
                u.shape, v.shape, w.shape,
                u_face.shape, v_face.shape, w_face.shape,
            )
            if any(array.shape != shape for array, shape in zip(arrays, expected)):
                raise ValueError("interface momentum arrays must match solver shapes")
            if any(not np.all(np.isfinite(array)) for array in arrays):
                raise ValueError("interface momentum arrays must be finite")
            u, v, w, u_face, v_face, w_face = arrays
            u[solid] = v[solid] = w[solid] = 0.0
            ux, uy, uz = (
                u / velocity_scale,
                v / velocity_scale,
                w / velocity_scale,
            )

        current_diffusion_number = 6.0 * float(
            np.max(np.asarray(applied_scalar_diffusivity))
        ) * dt / dx**2
        current_thermal_diffusion_number = (
            6.0 * float(np.max(np.asarray(applied_thermal_diffusivity)))
            * dt / dx**2
            if sensible_enthalpy is not None else 0.0
        )
        # Non-thermal runs do not enter the corrector branch, so derive the
        # lattice velocity here as well before applying the common stability
        # gates and optional interface momentum callback.
        ux, uy, uz = (
            u / velocity_scale,
            v / velocity_scale,
            w / velocity_scale,
        )
        current_lattice_mach = float(np.max(
            np.sqrt(ux[~solid]**2 + uy[~solid]**2 + uz[~solid]**2)
            / math.sqrt(_CS2_3D)
        ))
        maximum_lattice_mach = max(
            maximum_lattice_mach, current_lattice_mach,
        )
        if current_lattice_mach > config.maximum_lattice_mach_number:
            raise RuntimeError(
                "active buoyancy/corrector exceeded the declared lattice Mach gate: "
                f"{current_lattice_mach:.6g} > "
                f"{config.maximum_lattice_mach_number:.6g}"
            )
        current_advective_number = dt * (
            float(np.max(np.abs(u_face))) / dx
            + float(np.max(np.abs(v_face))) / dx
            + float(np.max(np.abs(w_face))) / dx
        )
        maximum_advective_number = max(
            maximum_advective_number, current_advective_number,
        )
        if current_advective_number + max(
            current_diffusion_number, current_thermal_diffusion_number,
        ) > 1.0 + 1.0e-12:
            raise RuntimeError(
                "active thermal scalar exceeded its explicit stability limit: "
                f"advection={current_advective_number:.6g}, "
                f"diffusion={current_diffusion_number:.6g}, "
                "thermal_diffusion="
                f"{current_thermal_diffusion_number:.6g}"
            )
        if not corrector_transport_done:
            if interface_state_callback is not None:
                updated = interface_state_callback(concentration.copy(), start, dt)
                updated = np.asarray(updated, dtype=float)
                if updated.shape != concentration.shape or not np.all(np.isfinite(updated)) or np.any(updated < 0.0):
                    raise ValueError("interface_state_callback must return a finite non-negative concentration field")
                concentration = updated
            concentration, step_outflow = _transport_step(
                concentration,
                u_face=u_face, v_face=v_face, w_face=w_face,
                solid=solid, diffusivity=applied_scalar_diffusivity,
                dt=dt, spacing=np.asarray((dx, dx, dx)),
                advection_scheme=config.scalar_advection_scheme,
                open_boundary_faces=projection_open_faces,
                blocked_faces=obstacle_blocked_faces,
                face_flux_callback=face_flux_callback,
            )
            if source_post_fraction:
                concentration += source_post_fraction * source_hydrogen_increment
        minimum = float(np.min(concentration))
        if minimum < -1.0e-12:
            raise RuntimeError(
                f"resolved 3-D concentration became negative: {minimum:.6g}"
            )
        concentration = np.maximum(concentration, 0.0)
        np.maximum(maximum_concentration, concentration, out=maximum_concentration)
        if not corrector_transport_done:
            step_energy_outflow = 0.0
        if sensible_enthalpy is not None and not corrector_transport_done:
            sensible_enthalpy, step_energy_outflow = _transport_step(
                sensible_enthalpy,
                u_face=u_face, v_face=v_face, w_face=w_face,
                solid=solid, diffusivity=applied_thermal_diffusivity,
                dt=dt, spacing=np.asarray((dx, dx, dx)),
                signed_boundary_flux=True,
                advection_scheme=config.scalar_advection_scheme,
                open_boundary_faces=projection_open_faces,
                blocked_faces=obstacle_blocked_faces,
                face_flux_callback=enthalpy_flux_callback,
            )
            if source_post_fraction:
                sensible_enthalpy += source_post_fraction * source_energy_increment
            if not np.all(np.isfinite(sensible_enthalpy)):
                raise RuntimeError("resolved 3-D sensible enthalpy became non-finite")
            assert closure is not None
            if binary_initial_cloud:
                thermodynamic_state = (
                    closure
                    .state_from_species_mass_concentration_and_sensible_enthalpy(
                        concentration, sensible_enthalpy,
                    )
                )
            else:
                thermodynamic_state = closure.state_from_conserved(
                    concentration, sensible_enthalpy,
                )
        if sensible_enthalpy is not None:
            assert thermodynamic_state is not None
            maximum_density_departure = max(
                maximum_density_departure,
                float(np.max(np.abs(
                    thermodynamic_state.relative_density_departure[~solid]
                ))),
            )
            if not binary_initial_cloud:
                active_scalar_out_of_range_cell_steps += int(np.count_nonzero(
                    (~thermodynamic_state.active_scalar_applicability_mask)
                    & (~solid)
                ))
            assert maximum_mole_fraction is not None
            current_mole_fraction = (
                thermodynamic_state.species_mole_fraction
                if binary_initial_cloud
                else thermodynamic_state.hydrogen_mole_fraction
            )
            np.maximum(
                maximum_mole_fraction,
                current_mole_fraction,
                out=maximum_mole_fraction,
            )
        if (
            handoff_request is not None
            and handoff_sample_interval_steps is not None
            and (
                (step + 1) % handoff_sample_interval_steps == 0
                or step == steps - 1
            )
        ):
            emit_handoff_frame(end)
        if (
            config.transient_average_start_s is not None
            and end + 1.0e-12 >= config.transient_average_start_s
        ):
            if averaged_u is None:
                averaged_u = np.zeros_like(u)
                averaged_v = np.zeros_like(v)
                averaged_w = np.zeros_like(w)
                averaged_concentration = np.zeros_like(concentration)
                if thermodynamic_state is not None:
                    averaged_temperature = np.zeros_like(
                        thermodynamic_state.temperature_k
                    )
            averaged_u += u
            averaged_v += v
            averaged_w += w
            averaged_concentration += concentration
            if averaged_temperature is not None and thermodynamic_state is not None:
                averaged_temperature += thermodynamic_state.temperature_k
            transient_average_count += 1
        cumulative_outflow += step_outflow
        cumulative_outflow_energy += step_energy_outflow
        mass = float(np.sum(concentration) * cell_volume)
        times.append(end)
        injected.append(cumulative_injected)
        domain_mass.append(mass)
        outflow.append(cumulative_outflow)
        residual.append(
            mass + cumulative_outflow
            - initial_species_mass - cumulative_injected
        )
        if sensible_enthalpy is not None:
            assert injected_energy is not None
            assert domain_energy is not None
            assert outflow_energy is not None
            assert energy_residual is not None
            energy = float(np.sum(sensible_enthalpy) * cell_volume)
            injected_energy.append(
                cumulative_injected_energy + cumulative_wall_heat_energy_j
            )
            domain_energy.append(energy)
            outflow_energy.append(cumulative_outflow_energy)
            energy_residual.append(
                energy + cumulative_outflow_energy
                - initial_sensible_enthalpy - cumulative_injected_energy
                - cumulative_wall_heat_energy_j
            )
        for receptor, stencil in zip(checked_receptors, stencils):
            traces[receptor.label].append(max(0.0, float(sum(
                weight * concentration[k, j, i]
                for k, j, i, weight in stencil
            ))))
            if thermodynamic_state is not None:
                current_mole_fraction = (
                    thermodynamic_state.species_mole_fraction
                    if binary_initial_cloud
                    else thermodynamic_state.hydrogen_mole_fraction
                )
                vol_percent_traces[receptor.label].append(100.0 * float(np.clip(
                    sum(
                        weight * current_mole_fraction[k, j, i]
                        for k, j, i, weight in stencil
                    ), 0.0, 1.0,
                )))
                temperature_traces[receptor.label].append(float(sum(
                    weight * thermodynamic_state.temperature_k[k, j, i]
                    for k, j, i, weight in stencil
                )))
            else:
                molar_density = config.ambient_pressure_pa / (
                    R_GAS_J_MOL_K * config.ambient_temperature_k
                )
                receptor_concentration = traces[receptor.label][-1]
                vol_percent_traces[receptor.label].append(float(np.clip(
                    receptor_concentration / HYDROGEN_MOLAR_MASS_KG_MOL
                    / molar_density * 100.0, 0.0, 100.0,
                )))
                temperature_traces[receptor.label].append(
                    config.ambient_temperature_k
                )
        if fields is not None:
            fields.append(concentration.copy())
        if temperature_fields is not None:
            assert thermodynamic_state is not None
            temperature_fields.append(thermodynamic_state.temperature_k.copy())
    transport_runtime = time.perf_counter() - transport_started
    if config.thermal_buoyancy_coupling:
        velocity_scale = dx / flow_dt
        pressure_scale = config.ambient_air_density_kg_m3 * velocity_scale**2
        outlet_reference = float(np.mean(lattice_density[:, :, -1]))
        pressure = _CS2_3D * (
            lattice_density - outlet_reference
        ) * pressure_scale
        if last_pressure_correction is not None:
            pressure += last_pressure_correction
        assert config.velocity_reference_scale_m_s is not None
        dynamic_pressure = (
            0.5 * config.ambient_air_density_kg_m3
            * config.velocity_reference_scale_m_s**2
        )
        pressure_coefficient = pressure / dynamic_pressure
        pressure[solid] = pressure_coefficient[solid] = np.nan
    final_concentration = concentration.copy()
    final_concentration[solid] = np.nan
    final_maximum = maximum_concentration.copy()
    final_maximum[solid] = np.nan
    final_enthalpy = final_temperature = final_mole_fraction = None
    final_maximum_mole_fraction = final_mixture_density = final_buoyancy = None
    final_pressure_correction = final_continuity_residual = None
    if last_pressure_correction is not None:
        final_pressure_correction = last_pressure_correction.copy()
        final_pressure_correction[solid] = np.nan
    if last_continuity_residual is not None:
        final_continuity_residual = last_continuity_residual.copy()
        final_continuity_residual[solid] = np.nan
    if sensible_enthalpy is not None:
        assert thermodynamic_state is not None
        assert maximum_mole_fraction is not None
        final_enthalpy = sensible_enthalpy.copy()
        final_temperature = thermodynamic_state.temperature_k.copy()
        final_mole_fraction = (
            thermodynamic_state.species_mole_fraction.copy()
            if binary_initial_cloud
            else thermodynamic_state.hydrogen_mole_fraction.copy()
        )
        final_maximum_mole_fraction = maximum_mole_fraction.copy()
        final_mixture_density = thermodynamic_state.mixture_density_kg_m3.copy()
        final_buoyancy = thermodynamic_state.buoyancy_acceleration_m_s2.copy()
        for array in (
            final_enthalpy, final_temperature, final_mole_fraction,
            final_maximum_mole_fraction, final_mixture_density, final_buoyancy,
        ):
            array[solid] = np.nan
    fluid = ~solid
    (
        background_eddy_viscosity_minimum,
        background_eddy_viscosity_maximum,
        background_eddy_viscosity_mean,
    ) = _diffusivity_statistics(background_eddy_viscosity, fluid)
    scalar_diffusivity_minimum, scalar_diffusivity_maximum, scalar_diffusivity_mean = (
        _diffusivity_statistics(applied_scalar_diffusivity, fluid)
    )
    thermal_diffusivity_minimum, thermal_diffusivity_maximum, thermal_diffusivity_mean = (
        _diffusivity_statistics(applied_thermal_diffusivity, fluid)
    )
    density_variation = float(
        (np.max(lattice_density[fluid]) - np.min(lattice_density[fluid]))
        / np.mean(lattice_density[fluid])
    )
    if config.obstacles:
        obstacle = config.obstacles[0]
        (
            obstacle_x_min, obstacle_x_max,
            obstacle_y_min, obstacle_y_max,
            obstacle_z_min, obstacle_z_max,
        ) = solid_obstacle_bounds_m(obstacle)
        obstacle_height = obstacle_z_max - obstacle_z_min
        lee = (
            (x[None, None, :] > obstacle_x_max)
            & (
                x[None, None, :]
                <= obstacle_x_max + 4.0 * obstacle_height
            )
            & (y[None, :, None] >= obstacle_y_min)
            & (y[None, :, None] <= obstacle_y_max)
            & (z[:, None, None] <= 1.5 * obstacle_z_max)
            & fluid
        )
    else:
        lee = np.zeros_like(solid)
    warnings = [
        "quantitative facility use is withheld pending canonical and external validation",
    ]
    source_sigma_xyz = (
        config.source_sigma_xyz_m
        if config.source_sigma_xyz_m is not None
        else (config.source_sigma_m,) * 3
    )
    source_fwhm_cells = _source_kernel_fwhm_cells(source_sigma_xyz, dx)
    source_kernel_underresolved = min(source_fwhm_cells) < 2.0
    if source_kernel_underresolved:
        warnings.append(
            "source Gaussian FWHM is below two cells on at least one axis; "
            "source shape is under-resolved and concentration peaks are "
            "regularization-sensitive"
        )
    if config.carrier_flow_initialization == "projected_inlet":
        warnings.append(
            "the projected-inlet carrier field is a conservative initial "
            "condition for finite-volume momentum, not a converged wake or "
            "an independently qualified turbulence closure"
        )
    if config.flow_boundary_scheme == "nonequilibrium_extrapolation":
        warnings.append(
            "non-equilibrium extrapolation inlet/pressure outlet is an opt-in "
            "candidate pending open-boundary canonical and AIJ verification"
        )
    elif config.flow_boundary_scheme in {
        "characteristic_nonreflecting", "characteristic_multidimensional",
    }:
        warnings.append(
            "characteristic pressure-open treatment is an opt-in component "
            "candidate and is not a qualified facility boundary"
        )
    else:
        warnings.append(
            "equilibrium inlet and copy outlet require independent boundary verification"
        )
    if config.active_inlet_profile_y_m is not None:
        warnings.append(
            "the active y-dependent inlet is a declared mean-velocity profile; "
            "it does not supply fan turbulence RMS, integral length scale, "
            "spectrum or dissipation"
        )
    if config.active_inlet_velocity_history is not None:
        warnings.append(
            "the time-dependent three-component inlet applies one measured "
            "vector history over the complete inlet-plane shape; it preserves "
            "temporal variation but does not reconstruct spatial turbulence"
        )
    if config.background_surface_layer_mixing is not None:
        warnings.append(
            "neutral Monin-Obukhov surface-layer mixing uses declared u_star "
            "and z0; no stability correction is inferred from the reported "
            "D/E class because a Monin-Obukhov length was not measured"
        )
        warnings.append(
            "the background surface-layer closure supplies vertical ambient "
            "eddy viscosity but no obstacle-generated turbulence production; "
            "the declared Schmidt/Prandtl envelope requires independent scoring"
        )
    if average_count:
        warnings.append(
            "reported carrier velocity is a finite-window time mean; statistical "
            "convergence compares the first and second halves of that window"
        )
    if maximum_obstacle_face_error > 1.0e-9 * dx:
        warnings.append(
            "one or more declared obstacle faces are not aligned with cell "
            "faces; voxelization error is reported and quantitative grid "
            "comparison requires aligned geometry"
        )
    if config.subgrid_model == "wale" and config.wale_constant > 0.0:
        warnings.append(
            "WALE is an opt-in candidate with invariant checks only; D3Q19 "
            "channel and wall-mounted-cube validation remains pending"
        )
    elif config.subgrid_model == "smagorinsky" and config.smagorinsky_constant > 0.0:
        warnings.append(
            "Smagorinsky closure has no wall damping and is not calibrated or validated"
        )
    elif config.background_surface_layer_mixing is not None:
        warnings.append(
            "no obstacle-scale subgrid closure is active; the atmospheric "
            "background closure does not replace wake-turbulence validation"
        )
    else:
        warnings.append("no subgrid turbulence closure is active")
    if config.thermal_buoyancy_coupling:
        if config.low_mach_mass_flux_correction:
            if binary_initial_cloud:
                warnings.append(
                    "binary initial-cloud transport conserves species mass, uses "
                    "full-mixture EOS density and projects total mass flux; its "
                    "constant-property ideal-gas closure remains research-only"
                )
            else:
                warnings.append(
                    "mass-flux continuity is enforced, but carrier thermodynamics "
                    "remains an active-scalar approximation rather than a full "
                    "variable-property mixture solve"
                )
            if config.active_momentum_solver == "finite_volume_projection":
                warnings.append(
                    "finite-volume variable-density momentum is an opt-in "
                    "candidate pending canonical, grid/time and external validation"
                )
            elif config.pressure_correction_momentum_coupling == "transport_only":
                warnings.append(
                    "transport_only pressure correction is a causal diagnostic "
                    "that omits pressure-to-momentum feedback and is not a "
                    "promotable physical model"
                )
            else:
                warnings.append(
                    "the selected pressure-to-D3Q19 feedback has not passed the "
                    "active three-level time gate and remains withheld"
                )
            if (
                config.active_provisional_face_projection
                == "variable_density_only"
            ):
                warnings.append(
                    "variable_density_only provisional projection is a rejected "
                    "diagnostic candidate and did not cure active time sensitivity"
                )
        else:
            warnings.append(
                "active-scalar carrier inventory is not a variable-density "
                "continuity/EOS solve"
            )
        warnings.append(
            "constant heat capacities omit real-gas and phase-change source thermodynamics"
        )
        if config.stratification_damping_alpha is not None:
            warnings.append(
                "opt-in Richardson-number vertical turbulent-diffusion damping "
                "is a cross-domain sensitivity candidate; its coefficients and "
                "dense-gas predictions have not been independently validated"
            )
        if corrector_nonconverged_steps:
            warnings.append(
                "thermodynamic mass-flux corrector did not meet its fixed-point "
                "tolerance on one or more steps; affected results are withheld"
            )
        if active_scalar_out_of_range_cell_steps and not binary_initial_cloud:
            warnings.append(
                "active-scalar density-departure limit was exceeded; affected results are withheld"
            )
    else:
        warnings.append(
            "isothermal one-way coupling omits cryogenic energy and buoyancy feedback"
        )
    if config.effective_schmidt_number is not None:
        if (
            config.scalar_turbulent_diffusivity_model
            == "legacy_total_kinematic_viscosity_ratio"
        ):
            warnings.append(
                "legacy uniform effective Schmidt nu/Sc scalar diffusion is "
                "retained only for v14 reproduction; it is a numerical "
                "surrogate, not an SGS closure"
            )
        else:
            warnings.append(
                "molecular plus local SGS eddy-viscosity/Sc scalar transport is "
                "a gradient-diffusion candidate pending external validation"
            )
    if config.effective_prandtl_number is not None:
        if (
            config.thermal_turbulent_diffusivity_model
            == "legacy_total_kinematic_viscosity_ratio"
        ):
            warnings.append(
                "legacy uniform effective Prandtl nu/Pr thermal diffusion is "
                "retained only for v14 reproduction; it is a numerical "
                "surrogate, not an SGS closure"
            )
        else:
            warnings.append(
                "molecular plus local SGS eddy-viscosity/Pr thermal transport is "
                "a gradient-diffusion candidate pending external validation"
            )
    completed_restart_state = None
    if binary_initial_cloud and config.active_momentum_solver == "finite_volume_projection":
        assert sensible_enthalpy is not None
        assert maximum_mole_fraction is not None
        completed_restart_state = ResolvedObstacleFlow3DRestartState(
            schema=RESOLVED_OBSTACLE_FLOW_3D_RESTART_SCHEMA,
            configuration_sha256=_restart_configuration_sha256(config),
            time_s=float(times[-1]),
            concentration_kg_m3=concentration,
            sensible_enthalpy_density_j_m3=sensible_enthalpy,
            maximum_concentration_kg_m3=maximum_concentration,
            maximum_species_mole_fraction=maximum_mole_fraction,
            velocity_u_m_s=u,
            velocity_v_m_s=v,
            velocity_w_m_s=w,
            u_face_m_s=u_face,
            v_face_m_s=v_face,
            w_face_m_s=w_face,
            projection_correction_potential_kg_m_s=last_projection_potential,
        )
    return ResolvedObstacleFlow3DResult(
        schema=RESOLVED_OBSTACLE_FLOW_3D_SCHEMA,
        x_m=x, y_m=y, z_m=z, solid_mask=solid,
        velocity_u_m_s=u, velocity_v_m_s=v, velocity_w_m_s=w,
        u_face_m_s=u_face, v_face_m_s=v_face, w_face_m_s=w_face,
        pressure_pa_gauge=pressure, pressure_coefficient=pressure_coefficient,
        eddy_viscosity_m2_s=eddy_viscosity,
        time_s=np.asarray(times),
        concentration_kg_m3=final_concentration,
        maximum_concentration_kg_m3=final_maximum,
        concentration_fields_kg_m3=(None if fields is None else np.asarray(fields)),
        sensible_enthalpy_density_j_m3=final_enthalpy,
        temperature_k=final_temperature,
        hydrogen_mole_fraction=(
            None if binary_initial_cloud else final_mole_fraction
        ),
        maximum_hydrogen_mole_fraction=(
            None if binary_initial_cloud else final_maximum_mole_fraction
        ),
        species_name=(
            config.binary_mixture_closure.species_name
            if binary_initial_cloud and config.binary_mixture_closure is not None
            else "hydrogen"
        ),
        species_mole_fraction=final_mole_fraction,
        maximum_species_mole_fraction=final_maximum_mole_fraction,
        mixture_density_kg_m3=final_mixture_density,
        buoyancy_acceleration_m_s2=final_buoyancy,
        pressure_correction_pa=final_pressure_correction,
        continuity_residual_kg_m3_s=final_continuity_residual,
        temperature_fields_k=(
            None if temperature_fields is None else np.asarray(temperature_fields)
        ),
        receptor_traces_kg_m3={label: np.asarray(values) for label, values in traces.items()},
        receptor_traces_vol_percent={
            label: np.asarray(values) for label, values in vol_percent_traces.items()
        },
        receptor_temperature_traces_k={
            label: np.asarray(values) for label, values in temperature_traces.items()
        },
        injected_mass_kg=np.asarray(injected),
        domain_mass_kg=np.asarray(domain_mass),
        outflow_mass_kg=np.asarray(outflow),
        mass_residual_kg=np.asarray(residual),
        injected_sensible_enthalpy_j=(
            None if injected_energy is None else np.asarray(injected_energy)
        ),
        domain_sensible_enthalpy_j=(
            None if domain_energy is None else np.asarray(domain_energy)
        ),
        outflow_sensible_enthalpy_j=(
            None if outflow_energy is None else np.asarray(outflow_energy)
        ),
        energy_residual_j=(
            None if energy_residual is None else np.asarray(energy_residual)
        ),
        flow_residual_history=residual_history,
        restart_state=completed_restart_state,
        diagnostics={
            "fidelity_tier": "tier_3_resolved_viscous_3d_prototype",
            "restart_applied": restart_state is not None,
            "restart_time_s": restart_time_s,
            "integration_stop_time_s": integration_stop_time_s,
            "applicability": "conditional",
            "quantitative_prediction_allowed": False,
            "nonlocal_scalar_routing": False,
            "local_conservative_species_transport": True,
            "local_conservative_sensible_enthalpy_transport": bool(
                config.thermal_buoyancy_coupling
            ),
            "interface_conserved_state_callback_applied": (
                interface_conserved_state_callback is not None
            ),
            "interface_momentum_callback_applied": (
                interface_momentum_callback is not None
            ),
            "interface_conserved_mass_callback_applied": (
                interface_conserved_mass_callback is not None
            ),
            "ambient_wind_speed_m_s": config.wind_speed_m_s,
            "velocity_reference_scale_m_s": (
                config.velocity_reference_scale_m_s
            ),
            "reference_dynamic_pressure_pa": (
                0.5 * config.ambient_air_density_kg_m3
                * config.velocity_reference_scale_m_s**2
            ),
            "zero_ambient_wind_mode": config.wind_speed_m_s == 0.0,
            "flow_model": (
                "conservative finite-volume variable-density momentum predictor"
                if config.thermal_buoyancy_coupling
                and config.active_momentum_solver == "finite_volume_projection"
                else (
                "D3Q19 BGK with Guo active-scalar buoyancy forcing"
                if config.thermal_buoyancy_coupling
                else "isothermal low-Mach D3Q19 BGK lattice-Boltzmann"
                )
            ),
            "pressure_model": (
                "isothermal lattice pressure plus variable-density mass-flux correction"
                if config.thermal_buoyancy_coupling
                and config.low_mach_mass_flux_correction
                else "isothermal lattice equation of state p=rho*c_s^2"
            ),
            "turbulence_model": (
                "neutral Monin-Obukhov background surface-layer eddy viscosity"
                if config.background_surface_layer_mixing is not None else (
                    "explicit Smagorinsky-Lilly strain viscosity"
                    if config.subgrid_model == "smagorinsky"
                    and config.smagorinsky_constant > 0.0
                    else (
                        "wall-adapting local eddy-viscosity (WALE)"
                        if config.subgrid_model == "wale"
                        and config.wale_constant > 0.0
                        else "none"
                    )
                )
            ),
            "background_surface_layer_mixing": (
                None if config.background_surface_layer_mixing is None
                else config.background_surface_layer_mixing.record()
            ),
            "background_eddy_viscosity_minimum_m2_s": (
                background_eddy_viscosity_minimum
            ),
            "background_eddy_viscosity_maximum_m2_s": (
                background_eddy_viscosity_maximum
            ),
            "background_eddy_viscosity_mean_m2_s": (
                background_eddy_viscosity_mean
            ),
            "subgrid_model": config.subgrid_model,
            "smagorinsky_constant": config.smagorinsky_constant,
            "wale_constant": config.wale_constant,
            "wall_damping": False,
            "wall_adapting_subgrid_model": bool(
                config.subgrid_model == "wale" and config.wale_constant > 0.0
            ),
            "scalar_model": "local conservative finite-volume advection-diffusion",
            "source_spatial_discretization": (
                "cell-integrated variance-matched compact box over Cartesian cells"
                if config.source_kernel == "variance_matched_box"
                else (
                    "cell-integrated cosine-tapered compact box over Cartesian cells"
                    if config.source_kernel == "cosine_tapered_box"
                    else (
                    "cell-integrated separable anisotropic Gaussian over Cartesian cells"
                    if config.source_sigma_xyz_m is not None
                    else "cell-integrated separable isotropic Gaussian over Cartesian cells"
                    )
                )
            ),
            "source_kernel": config.source_kernel,
            "source_schedule_id": config.source_schedule.source_id,
            "source_total_mass_to_species_mass_ratio": (
                config.source_total_mass_to_species_mass_ratio
            ),
            "source_schedule_mass_basis": (
                "tracked_species_mass; total mixture source mass is ratio-scaled"
                if config.source_total_mass_to_species_mass_ratio > 1.0
                else "tracked_species_mass_equals_total_source_mass"
            ),
            "source_momentum_velocity_m_s": list(
                config.source_momentum_velocity_m_s
            ),
            "source_momentum_velocity_history_id": (
                None if config.source_momentum_velocity_history is None
                else config.source_momentum_velocity_history.history_id
            ),
            "source_momentum_velocity_history_sample_count": (
                0 if config.source_momentum_velocity_history is None
                else len(config.source_momentum_velocity_history.time_s)
            ),
            "source_momentum_velocity_history_component_ranges_m_s": (
                None if config.source_momentum_velocity_history is None
                else [
                    [
                        float(min(channel)),
                        float(max(channel)),
                    ] for channel in (
                        config.source_momentum_velocity_history.u_m_s,
                        config.source_momentum_velocity_history.v_m_s,
                        config.source_momentum_velocity_history.w_m_s,
                    )
                ]
            ),
            "source_sigma_xyz_m": (
                None if config.source_sigma_xyz_m is None
                else list(config.source_sigma_xyz_m)
            ),
            "source_kernel_fwhm_cell_count_xyz": list(source_fwhm_cells),
            "source_kernel_underresolved": source_kernel_underresolved,
            "source_kernel_resolution_criterion": (
                "FWHM at least 2 grid cells in every axis; diagnostic threshold only"
            ),
            "scalar_advection_scheme": config.scalar_advection_scheme,
            "source_time_integration": config.source_time_integration,
            "active_momentum_solver": config.active_momentum_solver,
            "active_inlet_profile_model": (
                "declared_y_profile_times_z_profile"
                if config.active_inlet_profile_y_m is not None
                else "uniform_in_y_times_z_profile"
            ),
            "active_inlet_profile_y_m": (
                None if config.active_inlet_profile_y_m is None
                else list(config.active_inlet_profile_y_m)
            ),
            "active_inlet_profile_u_m_s": (
                None if config.active_inlet_profile_u_m_s is None
                else list(config.active_inlet_profile_u_m_s)
            ),
            "inlet_temperature_boundary_model": (
                "conservative_x_min_sensible_enthalpy_inflow; linearly interpolated in height with endpoint hold"
                if config.inlet_temperature_profile_z_m is not None else None
            ),
            "inlet_temperature_profile_z_m": (
                None if config.inlet_temperature_profile_z_m is None
                else list(config.inlet_temperature_profile_z_m)
            ),
            "inlet_temperature_profile_k": (
                None if config.inlet_temperature_profile_k is None
                else list(config.inlet_temperature_profile_k)
            ),
            "initial_temperature_profile_applied": bool(
                config.inlet_temperature_profile_z_m is not None
            ),
            "active_inlet_plane_minimum_m_s": float(
                np.min(active_inlet_physical)
            ),
            "active_inlet_plane_maximum_m_s": float(
                np.max(active_inlet_physical)
            ),
            "active_inlet_plane_mean_m_s": float(
                np.mean(active_inlet_physical)
            ),
            "active_inlet_velocity_history_id": (
                None if config.active_inlet_velocity_history is None
                else config.active_inlet_velocity_history.history_id
            ),
            "active_inlet_velocity_history_sample_count": (
                0 if config.active_inlet_velocity_history is None
                else len(config.active_inlet_velocity_history.time_s)
            ),
            "active_inlet_velocity_history_time_range_s": (
                None if config.active_inlet_velocity_history is None
                else [
                    config.active_inlet_velocity_history.time_s[0],
                    config.active_inlet_velocity_history.time_s[-1],
                ]
            ),
            "active_inlet_velocity_history_component_ranges_m_s": (
                None if config.active_inlet_velocity_history is None
                else {
                    "u": [
                        min(config.active_inlet_velocity_history.u_m_s),
                        max(config.active_inlet_velocity_history.u_m_s),
                    ],
                    "v": [
                        min(config.active_inlet_velocity_history.v_m_s),
                        max(config.active_inlet_velocity_history.v_m_s),
                    ],
                    "w": [
                        min(config.active_inlet_velocity_history.w_m_s),
                        max(config.active_inlet_velocity_history.w_m_s),
                    ],
                }
            ),
            "active_inlet_temporal_mapping": (
                (
                    "linear_interpolation"
                    if config.active_inlet_velocity_history.interpolation == "linear"
                    else "piecewise_constant"
                ) + "_single_vector_times_declared_spatial_shape"
                if config.active_inlet_velocity_history is not None
                else "steady"
            ),
            "minimum_applied_transient_inlet_u_m_s": minimum_transient_inlet_u,
            "maximum_applied_transient_inlet_u_m_s": maximum_transient_inlet_u,
            "maximum_applied_transient_inlet_speed_m_s": (
                maximum_transient_inlet_speed
            ),
            "source_momentum_velocity_m_s": list(
                config.source_momentum_velocity_m_s
            ),
            "active_provisional_face_projection": (
                config.active_provisional_face_projection
            ),
            "pressure_correction_momentum_coupling": (
                config.pressure_correction_momentum_coupling
            ),
            "scalar_diffusivity_model": (
                config.scalar_turbulent_diffusivity_model
                if config.effective_schmidt_number is not None
                else "declared constant diffusivity"
            ),
            "declared_scalar_diffusivity_m2_s": (
                config.scalar_diffusivity_m2_s
            ),
            "applied_scalar_diffusivity_m2_s": (
                float(applied_scalar_diffusivity)
                if np.ndim(applied_scalar_diffusivity) == 0 else None
            ),
            "applied_scalar_diffusivity_minimum_m2_s": (
                scalar_diffusivity_minimum
            ),
            "applied_scalar_diffusivity_maximum_m2_s": (
                scalar_diffusivity_maximum
            ),
            "applied_scalar_diffusivity_mean_m2_s": scalar_diffusivity_mean,
            "scalar_diffusivity_anisotropy": (
                "none"
                if config.stratification_damping_alpha is None
                else "horizontal molecular+neutral turbulent; vertical molecular+Ri_g-damped turbulent"
            ),
            "stratification_damping_scope": (
                "disabled"
                if config.stratification_damping_alpha is None
                else "species and sensible-enthalpy turbulent diffusion only; momentum eddy viscosity is unchanged"
            ),
            "stratification_damping_update": (
                "disabled"
                if config.stratification_damping_alpha is None
                else "once per active step from pre-transport mixture density and current cell-centred horizontal velocity (lagged Ri_g)"
            ),
            "stratification_damping_alpha": config.stratification_damping_alpha,
            "stratification_damping_exponent": (
                config.stratification_damping_exponent
            ),
            "stratification_shear_regularization_s_inv": (
                config.stratification_shear_regularization_s_inv
            ),
            "maximum_gradient_richardson_number": (
                maximum_gradient_richardson_number
            ),
            "effective_schmidt_number": config.effective_schmidt_number,
            "thermal_scalar_model": (
                "local conservative signed sensible-enthalpy advection-diffusion"
                if config.thermal_buoyancy_coupling else "disabled"
            ),
            "thermal_diffusivity_model": (
                config.thermal_turbulent_diffusivity_model
                if config.effective_prandtl_number is not None
                else "declared constant diffusivity"
            ),
            "declared_thermal_diffusivity_m2_s": (
                config.thermal_diffusivity_m2_s
            ),
            "applied_thermal_diffusivity_m2_s": (
                float(applied_thermal_diffusivity)
                if np.ndim(applied_thermal_diffusivity) == 0 else None
            ),
            "applied_thermal_diffusivity_minimum_m2_s": (
                thermal_diffusivity_minimum
            ),
            "applied_thermal_diffusivity_maximum_m2_s": (
                thermal_diffusivity_maximum
            ),
            "applied_thermal_diffusivity_mean_m2_s": thermal_diffusivity_mean,
            "effective_prandtl_number": config.effective_prandtl_number,
            "reported_cell_velocity": (
                "projected finite-volume variable-density momentum velocity"
                if config.thermal_buoyancy_coupling
                and config.active_momentum_solver == "finite_volume_projection"
                else (
                "face-mass-flux-corrected velocity reconciled with D3Q19 momentum"
                if config.thermal_buoyancy_coupling
                and config.low_mach_mass_flux_correction
                else "unprojected D3Q19 momentum moments"
                )
            ),
            "scalar_face_velocity": (
                "variable-density mass-flux pressure projection"
                if config.thermal_buoyancy_coupling
                and config.low_mach_mass_flux_correction
                else "mass-consistent Helmholtz projection"
            ),
            "energy_equation": config.thermal_buoyancy_coupling,
            "density_buoyancy_feedback": config.thermal_buoyancy_coupling,
            "variable_density_continuity_equation": bool(
                config.thermal_buoyancy_coupling
                and config.low_mach_mass_flux_correction
            ),
            "thermodynamic_corrector": (
                "species/enthalpy/EOS/mass-flux fixed point"
                if config.thermal_buoyancy_coupling
                and config.low_mach_mass_flux_correction
                else "disabled"
            ),
            "pressure_correction_to_momentum": bool(
                config.thermal_buoyancy_coupling
                and config.low_mach_mass_flux_correction
                and (
                    config.active_momentum_solver == "finite_volume_projection"
                    or config.pressure_correction_momentum_coupling
                    != "transport_only"
                )
            ),
            "thermodynamic_closure": (
                "constant-pressure conservative binary ideal-gas mixture"
                if binary_initial_cloud else (
                    "constant-cp ideal-gas H2/air active-scalar carrier approximation"
                    if config.thermal_buoyancy_coupling else "none"
                )
            ),
            "initial_cloud_model": (
                None if not binary_initial_cloud
                else "exact-cell-intersection regular-polygon prism"
            ),
            "initial_cloud_label": (
                None if config.initial_binary_cloud is None
                else config.initial_binary_cloud.label
            ),
            "initial_cloud_exact_geometry_volume_m3": (
                None if initial_cloud_record is None
                else initial_cloud_record.exact_geometry_volume_m3
            ),
            "initial_cloud_represented_geometry_volume_m3": (
                None if initial_cloud_record is None
                else initial_cloud_record.represented_geometry_volume_m3
            ),
            "initial_cloud_relative_geometry_volume_error": (
                None if initial_cloud_record is None
                else initial_cloud_record.relative_geometry_volume_error
            ),
            "initial_cloud_species_mass_kg": initial_species_mass,
            "initial_binary_cloud_momentum": (
                None if initial_cloud_record is None
                else config.initial_binary_cloud_momentum
            ),
            "initial_cloud_maximum_speed_m_s": (
                initial_cloud_maximum_speed_m_s
            ),
            "initial_cloud_maximum_momentum_density_kg_m2_s": (
                initial_cloud_maximum_momentum_density_kg_m2_s
            ),
            "binary_handoff_recording_enabled": handoff_request is not None,
            "binary_handoff_x_face_m": (
                None if handoff_request is None else handoff_request.x_face_m
            ),
            "binary_handoff_x_face_index": handoff_x_face_index,
            "binary_handoff_sample_interval_s": (
                None if handoff_request is None
                else handoff_request.sample_interval_s
            ),
            "binary_handoff_frame_count": handoff_frame_count,
            "obstacle_count": len(config.obstacles),
            "circular_arc_representation": config.circular_arc_representation,
            "d3q19_zero_thickness_arc_link_bounce_back_applied": bool(
                obstacle_blocked_links is not None
                and (
                    config.active_momentum_solver == "d3q19_bgk"
                    or config.carrier_flow_initialization == "d3q19_spinup"
                )
            ),
            "d3q19_zero_thickness_arc_link_rule": (
                "Bouzidi-Firdaouss-Lallemand linear interpolation; "
                "halfway fallback if upstream interpolation node is solid"
                if obstacle_blocked_links is not None
                else None
            ),
            "d3q19_zero_thickness_arc_link_count": int(
                0 if obstacle_blocked_links is None
                else np.count_nonzero(np.isfinite(obstacle_blocked_links))
            ),
            "embedded_arc_blocked_face_count": int(sum(
                np.count_nonzero(mask) for mask in obstacle_blocked_faces
            )),
            "embedded_arc_surface_area_m2": embedded_arc_surface_area_m2,
            "embedded_arc_declared_surface_area_m2": declared_arc_surface_area_m2,
            "embedded_arc_surface_cell_count": embedded_arc_surface_cell_count,
            "embedded_arc_minimum_wall_distance_m": (
                embedded_arc_minimum_wall_distance_m
                if zero_thickness_arc_present else None
            ),
            "embedded_arc_domain_area_coverage_fraction": (
                embedded_arc_surface_area_m2 / declared_arc_surface_area_m2
                if declared_arc_surface_area_m2 > 0.0 else None
            ),
            "embedded_arc_face_masks_applied_to_scalar_fluxes": True,
            "embedded_arc_fractional_cut_cell_areas": False,
            "embedded_arc_tangential_no_slip_stress_applied": bool(
                not zero_thickness_arc_present
                or config.active_momentum_solver in {
                    "d3q19_bgk", "finite_volume_projection",
                }
            ),
            "embedded_arc_momentum_boundary_model": (
                "smooth-wall Reichardt immersed shear stress plus blocked "
                "normal faces; roughness sensitivity still required"
                if zero_thickness_arc_present and (
                    config.active_momentum_solver == "finite_volume_projection"
                ) else (
                    "Bouzidi-Firdaouss-Lallemand interpolated no-slip links"
                    if zero_thickness_arc_present else (
                        "finite-thickness voxel no-slip boundary"
                        if any(isinstance(item, CircularArcWall) for item in config.obstacles)
                        else "not applicable: no circular-arc obstacle"
                    )
                )
            ),
            "finite_lateral_bypass": bool(config.obstacles),
            "solid_velocity_boundary": "link bounce-back no-slip",
            "ground_boundary": (
                "fixed-pressure/zero-gradient open candidate"
                if config.nonstreamwise_wall_boundary == "all_pressure_open"
                else (
                    "specular symmetry/free-slip surrogate"
                    if config.nonstreamwise_wall_boundary
                    in {"all_symmetry_free_slip", "ground_free_slip_lateral_top_pressure_open"}
                    else "on-node bounce-back no-slip"
                )
            ),
            "nonstreamwise_wall_boundary": config.nonstreamwise_wall_boundary,
            "lateral_and_top_boundary": (
                "full bounce-back no-slip tunnel walls"
                if config.nonstreamwise_wall_boundary == "all_no_slip_tunnel"
                else (
                    "fixed-pressure/zero-gradient open candidate"
                    if config.nonstreamwise_wall_boundary
                    in {"all_pressure_open", "ground_free_slip_lateral_top_pressure_open"}
                    else (
                        "all-face specular symmetry/free-slip surrogate"
                        if config.nonstreamwise_wall_boundary
                        == "all_symmetry_free_slip"
                        else "specular free slip"
                    )
                )
            ),
            "pressure_projection_open_boundary_faces": list(
                projection_open_faces
            ),
            "pressure_open_boundary_qualification": (
                (
                    "multidimensional_characteristic_candidate_"
                    "angle_sweep_corner_and_domain_validation_required"
                    if config.flow_boundary_scheme in {
                        "characteristic_multidimensional",
                        "characteristic_multidimensional_extrapolated",
                        "characteristic_multidimensional_vortical_limited",
                        "characteristic_nscbc_convective",
                    }
                    else (
                        "normal_and_30deg_linear_acoustic_gates_passed_"
                        "angle_sweep_failed_corner_and_domain_pending"
                        if config.flow_boundary_scheme
                        == "characteristic_nonreflecting"
                        else "fixed_pressure_candidate_reflection_gate_failed"
                    )
                )
                if config.nonstreamwise_wall_boundary in {"all_pressure_open", "ground_free_slip_lateral_top_pressure_open"}
                else "not_applicable"
            ),
            "flow_boundary_scheme": config.flow_boundary_scheme,
            "inlet_boundary": (
                "non-equilibrium extrapolation velocity profile"
                if config.flow_boundary_scheme in {
                    "nonequilibrium_extrapolation",
                    "characteristic_nonreflecting",
                    "characteristic_multidimensional",
                    "characteristic_multidimensional_extrapolated",
                    "characteristic_multidimensional_vortical_limited",
                    "characteristic_nscbc_convective",
                }
                else "equilibrium velocity profile"
            ),
            "outlet_boundary": (
                (
                    (
                        "time-dependent convective vortical mode with "
                        "incoming-acoustic characteristic suppression"
                        if config.flow_boundary_scheme
                        == "characteristic_nscbc_convective"
                        else (
                            "TVD-limited, acoustic/vortical-separated angle-aware "
                            "isothermal characteristic incoming-wave suppression"
                            if config.flow_boundary_scheme
                            == "characteristic_multidimensional_vortical_limited"
                            else (
                                "second-order-extrapolated angle-aware isothermal "
                                "characteristic incoming-wave suppression"
                                if config.flow_boundary_scheme
                                == "characteristic_multidimensional_extrapolated"
                                else "angle-aware isothermal characteristic incoming-wave suppression"
                            )
                        )
                    )
                )
                if config.flow_boundary_scheme in {
                    "characteristic_multidimensional",
                    "characteristic_multidimensional_extrapolated",
                    "characteristic_multidimensional_vortical_limited",
                    "characteristic_nscbc_convective",
                }
                else (
                    "isothermal characteristic incoming-wave suppression"
                    if config.flow_boundary_scheme
                    == "characteristic_nonreflecting"
                    else (
                        "non-equilibrium extrapolation fixed lattice density"
                        if config.flow_boundary_scheme
                        == "nonequilibrium_extrapolation"
                        else "zero-gradient population copy"
                    )
                )
            ),
            "outlet_lattice_density": config.outlet_lattice_density,
            "cell_size_m": dx,
            "obstacle_voxelization": obstacle_voxelization,
            "maximum_obstacle_face_error_m": maximum_obstacle_face_error,
            "maximum_obstacle_face_error_cells": (
                maximum_obstacle_face_error / dx
            ),
            "maximum_obstacle_size_error_m": maximum_obstacle_size_error,
            "maximum_obstacle_size_error_cells": (
                maximum_obstacle_size_error / dx
            ),
            "obstacle_grid_alignment_exact": bool(
                maximum_obstacle_face_error <= 1.0e-9 * dx
            ) if config.circular_arc_representation == "finite_thickness_voxel"
            else None,
            "obstacle_grid_alignment_metric_applicable": (
                config.circular_arc_representation == "finite_thickness_voxel"
            ),
            "ambient_air_density_kg_m3": config.ambient_air_density_kg_m3,
            "carrier_flow_initialization": config.carrier_flow_initialization,
            "flow_physical_time_step_s": flow_dt,
            "scalar_time_step_s": dt,
            "flow_steps_requested": config.flow_steps,
            "flow_steps_completed": actual_flow_steps,
            "flow_field_statistic": (
                "time_mean" if average_count else "final_instantaneous"
            ),
            "flow_average_last_steps": config.flow_average_last_steps,
            "flow_average_sample_count": average_count,
            "flow_statistical_convergence_tolerance": (
                config.flow_statistical_convergence_tolerance
            ),
            "flow_statistical_residual": statistical_residual,
            "instantaneous_flow_converged": bool(
                residual_history.size
                and residual_history[-1, 1] <= config.flow_convergence_tolerance
            ),
            "flow_converged": bool(
                (
                    statistical_residual
                    <= config.flow_statistical_convergence_tolerance
                )
                if statistical_residual is not None
                else (
                    residual_history.size
                    and residual_history[-1, 1]
                    <= config.flow_convergence_tolerance
                )
            ),
            "final_flow_residual": (
                float(statistical_residual)
                if statistical_residual is not None
                else (
                    None
                    if not residual_history.size
                    else float(residual_history[-1, 1])
                )
            ),
            "final_instantaneous_flow_residual": (
                None if not residual_history.size else float(residual_history[-1, 1])
            ),
            "flow_convergence_criterion": (
                "first-half versus second-half time-mean velocity"
                if statistical_residual is not None
                else "velocity change over convergence-check interval"
            ),
            "base_lattice_relaxation_time": base_tau,
            "maximum_lattice_relaxation_time": maximum_tau,
            "maximum_lattice_mach_number": (
                maximum_lattice_mach
            ),
            "maximum_allowed_lattice_mach_number": (
                config.maximum_lattice_mach_number
            ),
            "maximum_relative_lattice_density_variation": density_variation,
            "reynolds_number_obstacle_height": (
                None
                if not config.obstacles
                else (
                    config.wind_speed_m_s
                    * obstacle_height
                    / config.kinematic_viscosity_m2_s
                )
            ),
            "maximum_eddy_viscosity_m2_s": float(np.max(eddy_viscosity[fluid])),
            "maximum_eddy_viscosity_over_transient_m2_s": (
                maximum_eddy_viscosity_over_transient
            ),
            "maximum_absolute_raw_face_divergence_s_1": float(
                np.max(np.abs(raw_divergence[fluid]))
            ),
            "maximum_absolute_divergence_before_projection_s_1": (
                divergence_before_projection
            ),
            "maximum_absolute_divergence_s_1": divergence_after_projection,
            "face_projection_iterations": projection_iterations,
            "recirculating_lee_cell_count": int(np.count_nonzero(lee & (u < 0.0))),
            "minimum_lee_streamwise_velocity_m_s": (
                None if not np.any(lee) else float(np.min(u[lee]))
            ),
            "maximum_mass_residual_kg": float(np.max(np.abs(residual))),
            "maximum_relative_mass_residual": float(
                np.max(np.abs(residual))
                / max(abs(initial_species_mass + cumulative_injected), 1.0e-30)
            ),
            "maximum_energy_residual_j": (
                None if energy_residual is None
                else float(np.max(np.abs(energy_residual)))
            ),
            "maximum_relative_energy_residual": (
                None if energy_residual is None or injected_energy is None
                else float(
                    np.max(np.abs(energy_residual))
                    / max(
                        abs(initial_sensible_enthalpy + injected_energy[-1]),
                        1.0e-30,
                    )
                )
            ),
            "advective_number": maximum_advective_number,
            "diffusion_number": maximum_diffusion_number,
            "thermal_diffusion_number": maximum_thermal_diffusion_number,
            "source_temperature_k": config.source_temperature_k,
            "source_specific_sensible_enthalpy_j_kg": source_specific_enthalpy,
            "isothermal_wall_heat_boundary_model": (
                "first_order_half_cell_conduction_to_fluid_cells"
                if config.obstacle_wall_temperature_k is not None
                or config.ground_temperature_k is not None else None
            ),
            "obstacle_wall_temperature_k": config.obstacle_wall_temperature_k,
            "ground_temperature_k": config.ground_temperature_k,
            "cumulative_wall_heat_added_j": cumulative_wall_heat_energy_j,
            "maximum_absolute_wall_heat_rate_w": maximum_wall_heat_rate_w,
            "maximum_absolute_buoyancy_acceleration_m_s2": (
                maximum_buoyancy_acceleration
            ),
            "embedded_arc_wall_stress_closure": (
                "smooth_wall_reichardt"
                if zero_thickness_arc_present
                and config.active_momentum_solver == "finite_volume_projection"
                else (
                    "d3q19_bfl_interpolated_bounceback"
                    if zero_thickness_arc_present else None
                )
            ),
            "embedded_arc_wall_roughness_assumption": (
                "hydraulically smooth; roughness sensitivity not yet validated"
                if zero_thickness_arc_present
                and config.active_momentum_solver == "finite_volume_projection"
                else None
            ),
            "maximum_embedded_arc_wall_shear_stress_pa": (
                maximum_embedded_arc_wall_shear_stress
            ),
            "maximum_embedded_arc_wall_friction_velocity_m_s": (
                maximum_embedded_arc_wall_friction_velocity
            ),
            "maximum_embedded_arc_wall_tangential_speed_m_s": (
                maximum_embedded_arc_wall_tangential_speed
            ),
            "cumulative_embedded_arc_wall_impulse_x_n_s": float(
                cumulative_embedded_arc_wall_force_n[0]
            ),
            "cumulative_embedded_arc_wall_impulse_y_n_s": float(
                cumulative_embedded_arc_wall_force_n[1]
            ),
            "cumulative_embedded_arc_wall_impulse_z_n_s": float(
                cumulative_embedded_arc_wall_force_n[2]
            ),
            "maximum_relative_mixture_density_departure": (
                maximum_density_departure
            ),
            "maximum_active_scalar_density_departure": (
                config.maximum_active_scalar_density_departure
            ),
            "active_scalar_out_of_range_cell_steps": (
                active_scalar_out_of_range_cell_steps
            ),
            "active_scalar_applicability_passed": bool(
                not config.thermal_buoyancy_coupling
                or (
                    active_scalar_out_of_range_cell_steps == 0
                    and corrector_nonconverged_steps == 0
                )
            ),
            "thermodynamic_corrector_maximum_iterations_used": (
                None if not corrector_iterations_per_step
                else max(corrector_iterations_per_step)
            ),
            "thermodynamic_corrector_mean_iterations_used": (
                None if not corrector_iterations_per_step
                else float(np.mean(corrector_iterations_per_step))
            ),
            "thermodynamic_corrector_iterations_per_step": (
                list(corrector_iterations_per_step)
            ),
            "thermodynamic_corrector_nonconverged_steps": (
                corrector_nonconverged_steps
            ),
            "maximum_final_corrector_relative_density_change": (
                maximum_corrector_density_change
            ),
            "maximum_final_corrector_relative_velocity_change": (
                maximum_corrector_velocity_change
            ),
            "maximum_variable_density_continuity_residual_kg_m3_s": (
                maximum_variable_density_continuity_residual
            ),
            "maximum_variable_density_global_residual_kg_s": (
                maximum_variable_density_global_residual
            ),
            "maximum_momentum_balance_residual_kg_m_s": (
                maximum_momentum_balance_residual
            ),
            "maximum_momentum_advective_courant_number": (
                maximum_momentum_advective_courant
            ),
            "maximum_momentum_viscous_stability_number": (
                maximum_momentum_viscous_stability
            ),
            "variable_density_linear_solver": (
                last_variable_density_solver_method
            ),
            "variable_density_requested_linear_solver": (
                config.variable_density_projection_linear_solver
            ),
            "variable_density_preconditioner": (
                last_variable_density_preconditioner
            ),
            "variable_density_preconditioner_level_count": (
                last_variable_density_preconditioner_levels
            ),
            "variable_density_preconditioner_operator_complexity": (
                last_variable_density_preconditioner_operator_complexity
            ),
            "variable_density_pyamg_version": (
                last_variable_density_pyamg_version
            ),
            "variable_density_projection_warm_start": bool(
                config.low_mach_mass_flux_correction
            ),
            "variable_density_projection_maximum_linear_iterations_per_step": (
                None if not projection_linear_iterations_per_step
                else max(projection_linear_iterations_per_step)
            ),
            "variable_density_projection_mean_linear_iterations_per_step": (
                None if not projection_linear_iterations_per_step
                else float(np.mean(projection_linear_iterations_per_step))
            ),
            "variable_density_projection_linear_iterations_per_step": (
                list(projection_linear_iterations_per_step)
            ),
            "coupled_transient_steps": (
                steps if config.thermal_buoyancy_coupling else 0
            ),
            "flow_solve_runtime_s": flow_runtime,
            "transport_runtime_s": transport_runtime,
            "total_runtime_s": time.perf_counter() - wall_started,
            "warnings": warnings,
        },
        ambient_temperature_k=config.ambient_temperature_k,
        ambient_pressure_pa=config.ambient_pressure_pa,
        time_average_start_s=(
            config.transient_average_start_s
            if transient_average_count > 0 else None
        ),
        time_average_sample_count=transient_average_count,
        time_averaged_velocity_u_m_s=(
            None if averaged_u is None else averaged_u / transient_average_count
        ),
        time_averaged_velocity_v_m_s=(
            None if averaged_v is None else averaged_v / transient_average_count
        ),
        time_averaged_velocity_w_m_s=(
            None if averaged_w is None else averaged_w / transient_average_count
        ),
        time_averaged_concentration_kg_m3=(
            None if averaged_concentration is None
            else averaged_concentration / transient_average_count
        ),
        time_averaged_temperature_k=(
            None if averaged_temperature is None
            else averaged_temperature / transient_average_count
        ),
    )


__all__ = [
    "RESOLVED_OBSTACLE_FLOW_3D_SCHEMA",
    "RESOLVED_OBSTACLE_FLOW_3D_RESTART_SCHEMA",
    "ResolvedInletVelocityHistory3D",
    "ResolvedSourceMomentumHistory3D",
    "ResolvedObstacleFlow3DConfig",
    "ResolvedObstacleFlow3DRestartState",
    "ResolvedObstacleFlow3DResult",
    "read_resolved_obstacle_flow_3d_restart",
    "rebind_restart_for_corrector_tolerance_audit",
    "solve_resolved_obstacle_flow_3d",
    "write_resolved_obstacle_flow_3d_restart",
]
