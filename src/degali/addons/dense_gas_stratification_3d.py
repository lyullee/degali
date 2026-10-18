"""Observation-blind dense-layer stratification diagnostics.

This module diagnoses local stable density stratification without changing
momentum, species, or energy fluxes.  It intentionally contains no turbulent
damping law and no coefficient inferred from concentration observations.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np


DENSE_GAS_STRATIFICATION_3D_SCHEMA = "degali.dense-gas-stratification-3d.v1"
GRADIENT_RICHARDSON_SUMMARY_SCHEMA = (
    "degali.gradient-richardson-summary-3d.v1"
)


def _readonly(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=float).copy()
    result.setflags(write=False)
    return result


def finite_volume_vertical_gradient(
    field: np.ndarray,
    solid_mask: np.ndarray,
    *,
    grid_spacing_m: float,
) -> np.ndarray:
    """Return a cell-centred vertical gradient using FV face geometry.

    Interior fluid-fluid face values are arithmetic averages.  Domain and
    fluid-solid faces copy the adjacent fluid-cell value, which is the same
    zero-normal-gradient convention used by the viscous/diffusive finite-
    volume operators.  The array order is ``(z, y, x)``.
    """

    values = np.asarray(field, dtype=float).copy()
    solid = np.asarray(solid_mask)
    spacing = float(grid_spacing_m)
    if values.ndim != 3:
        raise ValueError("field must be a three-dimensional array")
    if solid.dtype != np.bool_ or solid.shape != values.shape:
        raise ValueError("solid_mask must be boolean with the field shape")
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("grid_spacing_m must be positive and finite")

    # Resolved obstacle fields intentionally carry NaN (or another sentinel)
    # in solid cells because those cells have no thermodynamic state.  The FV
    # face policy below never samples a solid value for a fluid-solid face, so
    # sanitize solid entries only to keep domain-face arithmetic finite.  A
    # non-finite fluid entry remains a hard input error.
    fluid = ~solid
    if not np.all(np.isfinite(values[fluid])):
        raise ValueError("fluid field entries must be finite")
    values[solid] = 0.0

    faces = np.zeros(
        (values.shape[0] + 1, values.shape[1], values.shape[2]), dtype=float,
    )
    faces[0] = values[0]
    faces[-1] = values[-1]
    lower = values[:-1]
    upper = values[1:]
    lower_fluid = fluid[:-1]
    upper_fluid = fluid[1:]
    both = lower_fluid & upper_fluid
    only_lower = lower_fluid & ~upper_fluid
    only_upper = ~lower_fluid & upper_fluid
    interior = faces[1:-1]
    interior[both] = 0.5 * (lower[both] + upper[both])
    interior[only_lower] = lower[only_lower]
    interior[only_upper] = upper[only_upper]
    gradient = (faces[1:] - faces[:-1]) / spacing
    gradient[solid] = 0.0
    return gradient


@dataclass(frozen=True)
class GradientRichardsonDiagnostic3D:
    density_vertical_gradient_kg_m4: np.ndarray
    horizontal_shear_squared_s_inv2: np.ndarray
    stable_buoyancy_frequency_squared_s_inv2: np.ndarray
    gradient_richardson_number: np.ndarray
    fluid_mask: np.ndarray
    grid_spacing_m: float
    shear_regularization_s_inv: float
    gravitational_acceleration_m_s2: float

    def record(self) -> dict[str, object]:
        return {
            "schema": DENSE_GAS_STRATIFICATION_3D_SCHEMA,
            "mode": "diagnostic_only_no_flux_modification",
            "formula": (
                "Ri_g=max(0,-(g/rho)*(drho/dz))/"
                "((du/dz)^2+(dv/dz)^2+epsilon_s^2)"
            ),
            "grid_spacing_m": self.grid_spacing_m,
            "shear_regularization_s_inv": self.shear_regularization_s_inv,
            "gravitational_acceleration_m_s2": (
                self.gravitational_acceleration_m_s2
            ),
            "solid_cell_value_policy": "zero_and_excluded_from_summaries",
            "fluxes_modified": False,
            "observation_derived_coefficients": [],
        }


def diagnose_gradient_richardson_3d(
    density_kg_m3: np.ndarray,
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    solid_mask: np.ndarray,
    *,
    grid_spacing_m: float,
    shear_regularization_s_inv: float,
    gravitational_acceleration_m_s2: float = 9.80665,
) -> GradientRichardsonDiagnostic3D:
    """Diagnose local stable-density Richardson number on a cubic FV grid."""

    density = np.asarray(density_kg_m3, dtype=float)
    u = np.asarray(u_m_s, dtype=float)
    v = np.asarray(v_m_s, dtype=float)
    solid = np.asarray(solid_mask)
    if density.ndim != 3 or u.shape != density.shape or v.shape != density.shape:
        raise ValueError("density and horizontal velocity must share one 3-D shape")
    if solid.dtype != np.bool_ or solid.shape != density.shape:
        raise ValueError("solid_mask must be boolean with the field shape")
    fluid = ~solid
    fluid = ~solid
    if (
        not np.all(np.isfinite(density[fluid]))
        or not np.all(np.isfinite(u[fluid]))
        or not np.all(np.isfinite(v[fluid]))
        or np.any(density[fluid] <= 0.0)
    ):
        raise ValueError("fluid density must be positive and all fields finite")
    epsilon = float(shear_regularization_s_inv)
    gravity = float(gravitational_acceleration_m_s2)
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("shear_regularization_s_inv must be positive and finite")
    if not math.isfinite(gravity) or gravity <= 0.0:
        raise ValueError(
            "gravitational_acceleration_m_s2 must be positive and finite"
        )

    density_gradient = finite_volume_vertical_gradient(
        density, solid, grid_spacing_m=grid_spacing_m,
    )
    du_dz = finite_volume_vertical_gradient(
        u, solid, grid_spacing_m=grid_spacing_m,
    )
    dv_dz = finite_volume_vertical_gradient(
        v, solid, grid_spacing_m=grid_spacing_m,
    )
    shear_squared = du_dz**2 + dv_dz**2
    buoyancy_squared = np.zeros_like(density)
    buoyancy_squared[fluid] = np.maximum(
        0.0,
        -gravity * density_gradient[fluid] / density[fluid],
    )
    richardson = np.zeros_like(density)
    richardson[fluid] = buoyancy_squared[fluid] / (
        shear_squared[fluid] + epsilon**2
    )
    if not np.all(np.isfinite(richardson)) or np.any(richardson < 0.0):
        raise RuntimeError("gradient Richardson diagnostic became invalid")
    return GradientRichardsonDiagnostic3D(
        density_vertical_gradient_kg_m4=_readonly(density_gradient),
        horizontal_shear_squared_s_inv2=_readonly(shear_squared),
        stable_buoyancy_frequency_squared_s_inv2=_readonly(buoyancy_squared),
        gradient_richardson_number=_readonly(richardson),
        fluid_mask=np.asarray(fluid, dtype=bool),
        grid_spacing_m=float(grid_spacing_m),
        shear_regularization_s_inv=epsilon,
        gravitational_acceleration_m_s2=gravity,
    )


def summarize_gradient_richardson_regions(
    diagnostic: GradientRichardsonDiagnostic3D,
    region_masks: Mapping[str, np.ndarray],
    *,
    neutral_scalar_eddy_diffusivity_m2_s: np.ndarray,
    cell_volume_m3: float,
) -> dict[str, object]:
    """Summarize ``Ri_g`` and undamped neutral mixing in fixed regions.

    The neutral diffusivity is reported beside the diagnostic; it is never
    multiplied by a damping function here.
    """

    diffusivity = np.asarray(neutral_scalar_eddy_diffusivity_m2_s, dtype=float)
    shape = diagnostic.gradient_richardson_number.shape
    volume = float(cell_volume_m3)
    if diffusivity.shape != shape or not np.all(np.isfinite(diffusivity)):
        raise ValueError("neutral diffusivity must be finite with diagnostic shape")
    if np.any(diffusivity < 0.0):
        raise ValueError("neutral diffusivity must be non-negative")
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("cell_volume_m3 must be positive and finite")
    regions: dict[str, object] = {}
    for label, supplied_mask in region_masks.items():
        if not isinstance(label, str) or not label:
            raise ValueError("region labels must be non-empty strings")
        mask = np.asarray(supplied_mask)
        if mask.dtype != np.bool_ or mask.shape != shape:
            raise ValueError("region masks must be boolean with diagnostic shape")
        selected = mask & diagnostic.fluid_mask
        if not np.any(selected):
            raise ValueError(f"region {label} contains no fluid cells")
        ri = diagnostic.gradient_richardson_number[selected]
        k_neutral = diffusivity[selected]
        regions[label] = {
            "fluid_cell_count": int(ri.size),
            "sampling_volume_m3": float(ri.size * volume),
            "gradient_richardson_mean": float(np.mean(ri)),
            "gradient_richardson_p50": float(np.quantile(ri, 0.50)),
            "gradient_richardson_p90": float(np.quantile(ri, 0.90)),
            "gradient_richardson_p95": float(np.quantile(ri, 0.95)),
            "gradient_richardson_maximum": float(np.max(ri)),
            "fraction_at_or_above_0p25": float(np.mean(ri >= 0.25)),
            "neutral_scalar_eddy_diffusivity_mean_m2_s": float(
                np.mean(k_neutral)
            ),
            "neutral_scalar_eddy_diffusivity_p95_m2_s": float(
                np.quantile(k_neutral, 0.95)
            ),
            "stable_neutral_diffusivity_volume_integral_m5_s": float(
                np.sum(k_neutral[ri > 0.0]) * volume
            ),
        }
    return {
        "schema": GRADIENT_RICHARDSON_SUMMARY_SCHEMA,
        "diagnostic": diagnostic.record(),
        "regions": regions,
        "damping_applied": False,
    }


def stratification_damped_diffusivity_tensor_3d(
    molecular_diffusivity_m2_s: float | np.ndarray,
    neutral_turbulent_diffusivity_m2_s: np.ndarray,
    gradient_richardson_number: np.ndarray,
    fluid_mask: np.ndarray,
    *,
    alpha: float,
    exponent: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a sensitivity-only anisotropic diffusivity tensor.

    Horizontal turbulent diffusivity is retained; only the *vertical
    turbulent* contribution is smoothly reduced as
    ``Kz_turb = Kturb / (1 + alpha*max(Ri_g, 0))**exponent``. Molecular
    diffusivity is added unchanged in all three directions. This Richardson-
    dependent shape is inspired by vertical mixing closures used in other
    stratified-flow domains (notably Pacanowski--Philander ocean models), but
    its parameters are not established for dense-gas releases. Callers must
    therefore treat it as an opt-in sensitivity candidate, never a validated
    default. There is no hard turbulence cutoff.

    The return order is ``(Dx, Dy, Dz)`` to match Cartesian face directions;
    solid-cell diffusivities are zero and the input arrays are not mutated.
    """
    turbulent = np.asarray(neutral_turbulent_diffusivity_m2_s, dtype=float)
    richardson = np.asarray(gradient_richardson_number, dtype=float)
    fluid = np.asarray(fluid_mask)
    if turbulent.ndim != 3 or richardson.shape != turbulent.shape:
        raise ValueError("turbulent diffusivity and Ri_g must share a 3-D shape")
    if fluid.dtype != np.bool_ or fluid.shape != turbulent.shape:
        raise ValueError("fluid_mask must be boolean with the field shape")
    if (
        not np.all(np.isfinite(turbulent[fluid]))
        or np.any(turbulent[fluid] < 0.0)
        or not np.all(np.isfinite(richardson[fluid]))
        or np.any(richardson[fluid] < 0.0)
    ):
        raise ValueError("fluid turbulent diffusivity and Ri_g must be finite and non-negative")
    alpha_value = float(alpha)
    exponent_value = float(exponent)
    if not math.isfinite(alpha_value) or alpha_value <= 0.0:
        raise ValueError("alpha must be positive and finite")
    if not math.isfinite(exponent_value) or exponent_value <= 0.0:
        raise ValueError("exponent must be positive and finite")
    molecular = np.asarray(molecular_diffusivity_m2_s, dtype=float)
    if molecular.ndim == 0:
        molecular = np.full(turbulent.shape, float(molecular), dtype=float)
    elif molecular.shape != turbulent.shape:
        raise ValueError("molecular diffusivity must be scalar or match the 3-D field")
    if not np.all(np.isfinite(molecular[fluid])) or np.any(molecular[fluid] < 0.0):
        raise ValueError("fluid molecular diffusivity must be finite and non-negative")

    damping = np.zeros_like(turbulent)
    damping[fluid] = np.power(
        1.0 + alpha_value * richardson[fluid], -exponent_value,
    )
    horizontal = np.zeros_like(turbulent)
    vertical = np.zeros_like(turbulent)
    horizontal[fluid] = molecular[fluid] + turbulent[fluid]
    vertical[fluid] = molecular[fluid] + damping[fluid] * turbulent[fluid]
    return horizontal, horizontal.copy(), vertical


__all__ = [
    "DENSE_GAS_STRATIFICATION_3D_SCHEMA",
    "GRADIENT_RICHARDSON_SUMMARY_SCHEMA",
    "GradientRichardsonDiagnostic3D",
    "finite_volume_vertical_gradient",
    "diagnose_gradient_richardson_3d",
    "summarize_gradient_richardson_regions",
    "stratification_damped_diffusivity_tensor_3d",
]
