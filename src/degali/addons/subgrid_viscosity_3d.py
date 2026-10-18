"""Local algebraic subgrid viscosities for the resolved D3Q19 tier.

The WALE invariant follows Nicoud and Ducros (1999), with unit lattice filter
width.  It is a research closure until the channel and wall-mounted-cube gates
are passed; implementing the invariant is not itself validation.
"""

from __future__ import annotations

import math

import numpy as np


SUBGRID_VISCOSITY_3D_SCHEMA = "degali.subgrid-viscosity-3d.v1"


def _checked_velocity(
    ux: np.ndarray, uy: np.ndarray, uz: np.ndarray, solid_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    arrays = tuple(np.asarray(value, dtype=float) for value in (ux, uy, uz))
    if arrays[0].ndim != 3 or any(value.shape != arrays[0].shape for value in arrays):
        raise ValueError("velocity components must share one three-dimensional shape")
    if any(not np.all(np.isfinite(value)) for value in arrays):
        raise ValueError("velocity components must be finite")
    solid = np.asarray(solid_mask)
    if solid.dtype != np.bool_ or solid.shape != arrays[0].shape:
        raise ValueError("solid_mask must be boolean with the velocity shape")
    return arrays[0], arrays[1], arrays[2], solid


def velocity_gradient_tensor_lattice(
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    solid_mask: np.ndarray,
) -> np.ndarray:
    """Return ``g[..., i, j] = d u_i / d x_j`` on a unit lattice."""
    ux, uy, uz, _solid = _checked_velocity(ux, uy, uz, solid_mask)
    du_dz, du_dy, du_dx = np.gradient(ux)
    dv_dz, dv_dy, dv_dx = np.gradient(uy)
    dw_dz, dw_dy, dw_dx = np.gradient(uz)
    gradient = np.empty((*ux.shape, 3, 3), dtype=float)
    gradient[..., 0, 0] = du_dx
    gradient[..., 0, 1] = du_dy
    gradient[..., 0, 2] = du_dz
    gradient[..., 1, 0] = dv_dx
    gradient[..., 1, 1] = dv_dy
    gradient[..., 1, 2] = dv_dz
    gradient[..., 2, 0] = dw_dx
    gradient[..., 2, 1] = dw_dy
    gradient[..., 2, 2] = dw_dz
    return gradient


def smagorinsky_eddy_viscosity_lattice(
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    solid_mask: np.ndarray,
    *,
    constant: float,
) -> np.ndarray:
    """Return the standard strain-only Smagorinsky eddy viscosity."""
    ux, uy, uz, solid = _checked_velocity(ux, uy, uz, solid_mask)
    coefficient = float(constant)
    if not math.isfinite(coefficient) or coefficient < 0.0:
        raise ValueError("constant must be finite and non-negative")
    if coefficient == 0.0:
        return np.zeros_like(ux)
    gradient = velocity_gradient_tensor_lattice(ux, uy, uz, solid)
    strain = 0.5 * (gradient + np.swapaxes(gradient, -1, -2))
    strain_magnitude = np.sqrt(2.0 * np.sum(strain**2, axis=(-2, -1)))
    viscosity = coefficient**2 * strain_magnitude
    viscosity[solid] = 0.0
    return viscosity


def wale_eddy_viscosity_lattice(
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    solid_mask: np.ndarray,
    *,
    constant: float = 0.5,
) -> np.ndarray:
    """Return the WALE eddy viscosity with unit lattice filter width."""
    ux, uy, uz, solid = _checked_velocity(ux, uy, uz, solid_mask)
    coefficient = float(constant)
    if not math.isfinite(coefficient) or coefficient < 0.0:
        raise ValueError("constant must be finite and non-negative")
    if coefficient == 0.0:
        return np.zeros_like(ux)
    gradient = velocity_gradient_tensor_lattice(ux, uy, uz, solid)
    strain = 0.5 * (gradient + np.swapaxes(gradient, -1, -2))
    gradient_squared = np.einsum(
        "...ik,...kj->...ij", gradient, gradient,
    )
    symmetric_squared = 0.5 * (
        gradient_squared + np.swapaxes(gradient_squared, -1, -2)
    )
    trace = np.trace(gradient_squared, axis1=-2, axis2=-1)
    identity = np.eye(3, dtype=float)
    traceless_squared = symmetric_squared - (
        trace[..., None, None] * identity / 3.0
    )
    strain_squared = np.sum(strain**2, axis=(-2, -1))
    wale_squared = np.sum(traceless_squared**2, axis=(-2, -1))
    numerator = wale_squared**1.5
    denominator = strain_squared**2.5 + wale_squared**1.25
    viscosity = np.zeros_like(ux)
    positive = denominator > np.finfo(float).tiny
    viscosity[positive] = (
        coefficient**2 * numerator[positive] / denominator[positive]
    )
    viscosity[solid] = 0.0
    if not np.all(np.isfinite(viscosity)) or np.any(viscosity < 0.0):
        raise RuntimeError("WALE eddy viscosity became invalid")
    return viscosity


def subgrid_eddy_viscosity_physical(
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    w_m_s: np.ndarray,
    solid_mask: np.ndarray,
    *,
    grid_spacing_m: float,
    model: str,
    smagorinsky_constant: float,
    wale_constant: float,
) -> np.ndarray:
    """Return a local SGS viscosity in physical units on a cubic grid.

    The lattice helpers differentiate with respect to the cell index and use a
    unit filter width.  For a physical velocity field, multiplying their
    result by the cubic filter width recovers ``(C Delta)^2 |grad u|`` in
    ``m2/s``.
    """
    spacing = float(grid_spacing_m)
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("grid_spacing_m must be positive and finite")
    if model not in {"none", "smagorinsky", "wale"}:
        raise ValueError("model must be 'none', 'smagorinsky', or 'wale'")
    if model == "none":
        ux, _uy, _uz, _solid = _checked_velocity(
            u_m_s, v_m_s, w_m_s, solid_mask,
        )
        return np.zeros_like(ux)
    if model == "smagorinsky":
        index_viscosity = smagorinsky_eddy_viscosity_lattice(
            u_m_s, v_m_s, w_m_s, solid_mask,
            constant=smagorinsky_constant,
        )
    else:
        index_viscosity = wale_eddy_viscosity_lattice(
            u_m_s, v_m_s, w_m_s, solid_mask,
            constant=wale_constant,
        )
    return index_viscosity * spacing


def subgrid_relaxation_time_lattice(
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    solid_mask: np.ndarray,
    *,
    base_lattice_viscosity: float,
    model: str,
    smagorinsky_constant: float,
    wale_constant: float,
    maximum_relaxation_time: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return BGK relaxation time and selected local SGS viscosity."""
    base = float(base_lattice_viscosity)
    maximum_tau = float(maximum_relaxation_time)
    if not math.isfinite(base) or base < 0.0:
        raise ValueError("base_lattice_viscosity must be finite and non-negative")
    if not math.isfinite(maximum_tau) or maximum_tau <= 0.5:
        raise ValueError("maximum_relaxation_time must exceed 0.5")
    if model not in {"none", "smagorinsky", "wale"}:
        raise ValueError("model must be 'none', 'smagorinsky', or 'wale'")
    if model == "none":
        eddy = np.zeros_like(np.asarray(ux, dtype=float))
    elif model == "smagorinsky":
        eddy = smagorinsky_eddy_viscosity_lattice(
            ux, uy, uz, solid_mask, constant=smagorinsky_constant,
        )
    else:
        eddy = wale_eddy_viscosity_lattice(
            ux, uy, uz, solid_mask, constant=wale_constant,
        )
    tau = np.minimum(0.5 + 3.0 * (base + eddy), maximum_tau)
    return tau, eddy


__all__ = [
    "SUBGRID_VISCOSITY_3D_SCHEMA",
    "velocity_gradient_tensor_lattice",
    "smagorinsky_eddy_viscosity_lattice",
    "wale_eddy_viscosity_lattice",
    "subgrid_eddy_viscosity_physical",
    "subgrid_relaxation_time_lattice",
]
