"""Conservative post-impact wall-jet and vertical-profile operators.

These operators cover the gap between a downward/impinging source and a
wind-aligned scalar transport plane.  They redistribute an already resolved
H2 field; they do not infer a source rate, fit a concentration multiplier, or
claim to be a three-dimensional wake solver.  Every angular sector, mixing
length, and vertical profile is therefore an explicit transport hypothesis
that must be validated separately from the source ledger.
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter1d


POST_IMPACT_WALLJET_SCHEMA = "degali.post-impact-walljet-operator.v1"


def _grid(values: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    field = np.asarray(values, dtype=float)
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    if field.shape != (y_values.size, x_values.size):
        raise ValueError("field shape must be (len(y), len(x))")
    if x_values.size < 2 or y_values.size < 2:
        raise ValueError("x and y must each contain at least two nodes")
    if np.any(~np.isfinite(field)) or np.any(field < 0.0) or np.any(field > 1.0):
        raise ValueError("field must contain finite values in [0, 1]")
    if np.any(~np.isfinite(x_values)) or np.any(~np.isfinite(y_values)):
        raise ValueError("x and y must be finite")
    if np.any(np.diff(x_values) <= 0.0) or np.any(np.diff(y_values) <= 0.0):
        raise ValueError("x and y must be strictly increasing")
    return field, x_values, y_values


def _integral(field: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    return float(np.trapezoid(np.trapezoid(field, x, axis=1), y, axis=0))


def _renormalize(field: np.ndarray, reference: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    before = _integral(reference, x, y)
    after = _integral(field, x, y)
    if before <= 1.0e-300:
        return np.zeros_like(field)
    if after <= 1.0e-300:
        raise ValueError("wall-jet rotation moved all scalar inventory outside the grid")
    result = np.maximum(field, 0.0) * (before / after)
    if float(np.max(result)) > 1.0 + 1.0e-8:
        raise ValueError("wall-jet grid is too narrow to conserve a bounded scalar field")
    return np.clip(result, 0.0, 1.0)


def _rotated_average(
    base_field: np.ndarray,
    *,
    x: np.ndarray,
    y: np.ndarray,
    angles_rad: np.ndarray,
    weights: np.ndarray,
    origin_x: float,
    origin_y: float,
) -> np.ndarray:
    values, x_values, y_values = _grid(base_field, x, y)
    if not math.isfinite(origin_x) or not math.isfinite(origin_y):
        raise ValueError("origin_x and origin_y must be finite")
    interpolator = RegularGridInterpolator(
        (y_values, x_values), values, bounds_error=False, fill_value=0.0,
    )
    xx, yy = np.meshgrid(x_values, y_values)
    local_x = xx - float(origin_x)
    local_y = yy - float(origin_y)
    output = np.zeros_like(values)
    for angle, weight in zip(angles_rad, weights, strict=True):
        rotated_x = local_x * math.cos(float(angle)) + local_y * math.sin(float(angle))
        rotated_y = -local_x * math.sin(float(angle)) + local_y * math.cos(float(angle))
        output += float(weight) * interpolator(
            np.column_stack(((rotated_y + origin_y).ravel(), (rotated_x + origin_x).ravel()))
        ).reshape(values.shape)
    return _renormalize(output, values, x_values, y_values)


def _sector_nodes(half_angle_deg: float, nodes: int) -> tuple[np.ndarray, np.ndarray]:
    if not math.isfinite(float(half_angle_deg)) or not 0.0 < half_angle_deg < 90.0:
        raise ValueError("half_angle_deg must lie in (0, 90)")
    if isinstance(nodes, bool) or not isinstance(nodes, int) or nodes < 3 or nodes % 2 == 0:
        raise ValueError("nodes must be an odd integer >= 3")
    legendre, weights = np.polynomial.legendre.leggauss(nodes)
    return math.radians(float(half_angle_deg)) * legendre, 0.5 * weights


def radial_walljet_field(
    base_field: np.ndarray,
    *,
    x: np.ndarray,
    y: np.ndarray,
    half_angle_deg: float = 60.0,
    nodes: int = 9,
    origin_x: float = 0.0,
    origin_y: float = 0.0,
) -> np.ndarray:
    """Average a wind-aligned field over a declared radial wall-jet sector."""

    angles, weights = _sector_nodes(half_angle_deg, nodes)
    return _rotated_average(
        base_field, x=x, y=y, angles_rad=angles, weights=weights,
        origin_x=origin_x, origin_y=origin_y,
    )


def biased_radial_walljet_field(
    base_field: np.ndarray,
    *,
    x: np.ndarray,
    y: np.ndarray,
    half_angle_deg: float = 60.0,
    nodes: int = 17,
    bias_kappa: float = 0.0,
    bias_center_deg: float = 0.0,
    origin_x: float = 0.0,
    origin_y: float = 0.0,
) -> np.ndarray:
    """Radial wall-jet with a source/vector-momentum von-Mises bias."""

    if not math.isfinite(float(bias_kappa)) or bias_kappa < 0.0:
        raise ValueError("bias_kappa must be finite and non-negative")
    if not math.isfinite(float(bias_center_deg)):
        raise ValueError("bias_center_deg must be finite")
    angles, weights = _sector_nodes(half_angle_deg, nodes)
    weights = weights * np.exp(float(bias_kappa) * np.cos(angles - math.radians(bias_center_deg)))
    weights /= float(np.sum(weights))
    return _rotated_average(
        base_field, x=x, y=y, angles_rad=angles, weights=weights,
        origin_x=origin_x, origin_y=origin_y,
    )


def bifurcated_walljet_field(
    base_field: np.ndarray,
    *,
    x: np.ndarray,
    y: np.ndarray,
    half_angle_deg: float = 60.0,
    split_fraction: float = 0.5,
    spread_fraction: float = 0.25,
    bias_kappa: float = 0.0,
    bias_center_deg: float = 0.0,
    lobe_center_offset_deg: float = 0.0,
    nodes: int = 33,
    origin_x: float = 0.0,
    origin_y: float = 0.0,
) -> np.ndarray:
    """Redistribute a field into two symmetric, normalized impact lobes."""

    if not 0.0 < split_fraction < 1.0:
        raise ValueError("split_fraction must lie in (0, 1)")
    if not 0.0 < spread_fraction <= 1.0:
        raise ValueError("spread_fraction must lie in (0, 1]")
    if not math.isfinite(float(lobe_center_offset_deg)):
        raise ValueError("lobe_center_offset_deg must be finite")
    if not math.isfinite(float(bias_kappa)) or bias_kappa < 0.0:
        raise ValueError("bias_kappa must be finite and non-negative")
    if not math.isfinite(float(bias_center_deg)):
        raise ValueError("bias_center_deg must be finite")
    angles, quadrature = _sector_nodes(half_angle_deg, nodes)
    sector = math.radians(float(half_angle_deg))
    centres = math.radians(float(lobe_center_offset_deg)) + np.asarray(
        (-split_fraction * sector, split_fraction * sector), dtype=float
    )
    sigma = spread_fraction * sector
    density = np.exp(-0.5 * ((angles[:, None] - centres[None, :]) / sigma) ** 2).sum(axis=1)
    if bias_kappa > 0.0:
        density *= np.exp(float(bias_kappa) * np.cos(angles - math.radians(bias_center_deg)))
    weights = quadrature * density
    weights /= float(np.sum(weights))
    return _rotated_average(
        base_field, x=x, y=y, angles_rad=angles, weights=weights,
        origin_x=origin_x, origin_y=origin_y,
    )


def _stability_factor(stability: str) -> float:
    value = str(stability).strip().upper()
    if value in {"A", "B"}:
        return 0.85
    if value in {"C", "D"}:
        return 1.0
    if value in {"E", "F"}:
        return 1.15
    raise ValueError("stability must be one of A, B, C, D, E or F")


def radial_walljet_mixing(
    field: np.ndarray,
    *,
    x: np.ndarray,
    y: np.ndarray,
    height_m: float,
    wind_speed_m_s: float,
    jet_speed_m_s: float,
    roughness_m: float = 0.03,
    source_height_m: float = 0.05,
    reference_height_m: float = 10.0,
    stability: str = "D",
    start_m: float = 0.0,
    cutoff_m: float | None = None,
    core_aware: bool = False,
) -> np.ndarray:
    """Apply source-derived lateral wall-jet mixing with column closure."""

    values, x_values, y_values = _grid(field, x, y)
    for name, value in (
        ("height_m", height_m), ("wind_speed_m_s", wind_speed_m_s),
        ("jet_speed_m_s", jet_speed_m_s), ("roughness_m", roughness_m),
        ("source_height_m", source_height_m),
        ("reference_height_m", reference_height_m), ("start_m", start_m),
    ):
        if not math.isfinite(float(value)):
            raise ValueError(f"{name} must be finite")
    if wind_speed_m_s <= 0.0 or jet_speed_m_s <= 0.0:
        raise ValueError("wind and jet speeds must be positive")
    if roughness_m <= 0.0 or source_height_m < 0.0 or reference_height_m <= roughness_m:
        raise ValueError("roughness/source/reference heights are inconsistent")
    if start_m < 0.0:
        raise ValueError("start_m must be non-negative")
    if cutoff_m is not None and (not math.isfinite(float(cutoff_m)) or cutoff_m < start_m):
        raise ValueError("cutoff_m must be finite and not below start_m")
    if y_values.size < 2:
        raise ValueError("y must contain at least two nodes")
    dy = float(np.median(np.diff(y_values)))
    u_star = 0.4 * float(wind_speed_m_s) / math.log(
        (float(reference_height_m) + roughness_m) / roughness_m
    )
    vertical_mix_length = max(float(height_m) - source_height_m + roughness_m, roughness_m)
    phi_m = _stability_factor(stability)
    output = values.copy()
    core_half_width = max(float(height_m), source_height_m + roughness_m) if core_aware else None
    for index, station in enumerate(x_values):
        distance = max(float(station) - float(start_m), 0.0)
        if distance <= 0.0 or (cutoff_m is not None and station > cutoff_m):
            continue
        mixing_length = math.sqrt(max(float(station) + roughness_m, roughness_m) * vertical_mix_length)
        sigma_y = math.sqrt(max(2.0 * u_star * mixing_length * distance / (phi_m * jet_speed_m_s), 0.0))
        if sigma_y <= 1.0e-12:
            continue
        incoming = np.maximum(output[:, index], 0.0)
        if core_half_width is None:
            mixed = gaussian_filter1d(incoming, sigma=sigma_y / dy, mode="constant", cval=0.0)
            before = float(np.trapezoid(incoming, y_values))
            after = float(np.trapezoid(np.maximum(mixed, 0.0), y_values))
            if before > 1.0e-300 and after > 1.0e-300:
                output[:, index] = np.maximum(mixed, 0.0) * before / after
        else:
            core_mask = np.abs(y_values) <= core_half_width
            core = np.where(core_mask, incoming, 0.0)
            wing = np.where(core_mask, 0.0, incoming)
            before = float(np.trapezoid(wing, y_values))
            mixed = gaussian_filter1d(wing, sigma=sigma_y / dy, mode="constant", cval=0.0)
            mixed = np.where(core_mask, 0.0, np.maximum(mixed, 0.0))
            after = float(np.trapezoid(mixed, y_values))
            if before > 1.0e-300 and after > 1.0e-300:
                output[:, index] = core + mixed * before / after
    return np.clip(output, 0.0, 1.0)


def _layer_thickness(heights_m: Sequence[float]) -> np.ndarray:
    heights = np.asarray(heights_m, dtype=float)
    if heights.ndim != 1 or heights.size < 2 or np.any(~np.isfinite(heights)):
        raise ValueError("heights_m must contain at least two finite nodes")
    if np.any(np.diff(heights) <= 0.0):
        raise ValueError("heights_m must be strictly increasing")
    edges = np.empty(heights.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (heights[:-1] + heights[1:])
    edges[0] = heights[0] - 0.5 * (heights[1] - heights[0])
    edges[-1] = heights[-1] + 0.5 * (heights[-1] - heights[-2])
    return np.diff(edges)


def source_state_profile_weights(
    heights_m: Sequence[float],
    x_m: Sequence[float],
    source_x_m: Sequence[float],
    source_depth_m: Sequence[float],
    source_centre_m: Sequence[float],
    *,
    ground_reflection: bool = True,
) -> np.ndarray:
    """Return normalized ground-reflected vertical finite-volume weights."""

    heights = np.asarray(heights_m, dtype=float)
    x_values = np.asarray(x_m, dtype=float)
    source_x = np.asarray(source_x_m, dtype=float)
    source_depth = np.asarray(source_depth_m, dtype=float)
    source_centre = np.asarray(source_centre_m, dtype=float)
    if x_values.ndim != 1 or np.any(~np.isfinite(x_values)):
        raise ValueError("x_m must be a finite one-dimensional array")
    if np.any(np.diff(x_values) <= 0.0):
        raise ValueError("x_m must be strictly increasing")
    if source_x.size < 2 or any(array.size != source_x.size for array in (source_depth, source_centre)):
        raise ValueError("source state arrays must have equal length >= 2")
    if any(np.any(~np.isfinite(array)) for array in (source_x, source_depth, source_centre)):
        raise ValueError("source state arrays must be finite")
    if np.any(np.diff(source_x) <= 0.0) or np.any(source_depth <= 0.0) or np.any(source_centre < 0.0):
        raise ValueError("source x must increase, depth must be positive, and centre must be non-negative")
    thickness = _layer_thickness(heights)
    depth = np.interp(x_values, source_x, source_depth)
    centre = np.minimum(np.interp(x_values, source_x, source_centre), 1.5 * depth)
    sigma = np.where(
        centre <= 0.5 * depth,
        np.maximum(depth - centre, 1.0e-12) / math.sqrt(3.0),
        0.5 * depth / math.sqrt(3.0),
    )
    z = heights[:, None]
    log_profile = -0.5 * ((z - centre[None, :]) / sigma[None, :]) ** 2
    if ground_reflection:
        log_profile = np.logaddexp(
            log_profile,
            -0.5 * ((z + centre[None, :]) / sigma[None, :]) ** 2,
        )
    profile = np.exp(log_profile - np.max(log_profile, axis=0, keepdims=True))
    normal = np.sum(thickness[:, None] * profile, axis=0)
    return profile / np.maximum(normal[None, :], 1.0e-300)


def redistribute_vertical_profile(
    fields: Mapping[float, np.ndarray],
    heights_m: Sequence[float],
    *,
    x_m: Sequence[float],
    source_x_m: Sequence[float],
    source_depth_m: Sequence[float],
    source_centre_m: Sequence[float],
    ground_reflection: bool = True,
) -> dict[float, np.ndarray]:
    """Redistribute existing multi-height fields without changing inventory."""

    heights = tuple(float(value) for value in heights_m)
    weights = source_state_profile_weights(
        heights, x_m, source_x_m, source_depth_m, source_centre_m,
        ground_reflection=ground_reflection,
    )
    arrays = [np.asarray(fields[height], dtype=float) for height in heights]
    if any(array.ndim != 2 for array in arrays) or any(np.any(~np.isfinite(array)) or np.any(array < 0.0) for array in arrays):
        raise ValueError("height fields must be finite, non-negative two-dimensional arrays")
    if any(array.shape != arrays[0].shape for array in arrays):
        raise ValueError("height fields must share a grid shape")
    x_values = np.asarray(x_m, dtype=float)
    if arrays[0].shape[1] != x_values.size:
        raise ValueError("field x dimension must match x_m")
    thickness = _layer_thickness(heights)
    inventory = sum(thickness[index] * arrays[index] for index in range(len(heights)))
    output = {
        height: inventory * weights[index][None, :]
        for index, height in enumerate(heights)
    }
    after = sum(thickness[index] * output[height] for index, height in enumerate(heights))
    if not np.allclose(inventory, after, rtol=2.0e-12, atol=2.0e-12):
        raise RuntimeError("vertical profile redistribution failed finite-volume closure")
    return output


__all__ = [
    "POST_IMPACT_WALLJET_SCHEMA", "radial_walljet_field",
    "biased_radial_walljet_field", "bifurcated_walljet_field",
    "radial_walljet_mixing", "source_state_profile_weights",
    "redistribute_vertical_profile",
]
