"""Finite-volume initial-cloud contracts for dense binary-gas releases."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .binary_mixture_thermodynamics import (
    BinaryIdealGasMixtureClosure,
    BinaryIdealGasMixtureState,
)


@dataclass(frozen=True)
class RegularPolygonPrismInitialCloud:
    """Uniform binary mixture in a vertical regular-polygon prism.

    ``flat_to_flat_m`` is the distance between opposite faces.  This convention
    matches the approximately 2,050 m3 Thorney Island gas holder when the
    twelve-sided holder is 14 m across and filled to 13 m.
    """

    center_xy_m: tuple[float, float]
    flat_to_flat_m: float
    z_min_m: float
    z_max_m: float
    species_mole_fraction: float
    temperature_k: float
    sides: int = 12
    rotation_deg: float = 0.0
    label: str = "initial_cloud"

    def __post_init__(self) -> None:
        if not isinstance(self.center_xy_m, tuple) or len(self.center_xy_m) != 2:
            raise TypeError("center_xy_m must be a two-value tuple")
        center = tuple(float(value) for value in self.center_xy_m)
        if any(not math.isfinite(value) for value in center):
            raise ValueError("center_xy_m must be finite")
        object.__setattr__(self, "center_xy_m", center)
        for name in (
            "flat_to_flat_m", "z_min_m", "z_max_m",
            "species_mole_fraction", "temperature_k", "rotation_deg",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        if self.flat_to_flat_m <= 0.0:
            raise ValueError("flat_to_flat_m must be positive")
        if self.z_min_m < 0.0 or self.z_max_m <= self.z_min_m:
            raise ValueError("cloud vertical bounds must satisfy 0 <= z_min < z_max")
        if not 0.0 <= self.species_mole_fraction <= 1.0:
            raise ValueError("species_mole_fraction must lie in [0, 1]")
        if self.temperature_k <= 0.0:
            raise ValueError("temperature_k must be positive")
        if isinstance(self.sides, bool) or not isinstance(self.sides, int):
            raise TypeError("sides must be an integer")
        if self.sides < 4 or self.sides % 2:
            raise ValueError("sides must be an even integer of at least four")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("label must be a non-empty string")

    @property
    def apothem_m(self) -> float:
        return 0.5 * self.flat_to_flat_m

    @property
    def circumradius_m(self) -> float:
        return self.apothem_m / math.cos(math.pi / self.sides)

    @property
    def horizontal_area_m2(self) -> float:
        return (
            self.sides * self.apothem_m**2 * math.tan(math.pi / self.sides)
        )

    @property
    def exact_volume_m3(self) -> float:
        return self.horizontal_area_m2 * (self.z_max_m - self.z_min_m)

    def vertices_xy_m(self) -> np.ndarray:
        # Half a sector offset makes rotation_deg the outward normal of one face.
        angles = (
            math.radians(self.rotation_deg)
            + math.pi / self.sides
            + 2.0 * math.pi * np.arange(self.sides) / self.sides
        )
        return np.column_stack((
            self.center_xy_m[0] + self.circumradius_m * np.cos(angles),
            self.center_xy_m[1] + self.circumradius_m * np.sin(angles),
        ))


@dataclass(frozen=True)
class FiniteVolumeInitialCloud3D:
    """Cell-averaged conserved initial state and its geometry audit."""

    volume_fraction: np.ndarray
    species_mass_concentration_kg_m3: np.ndarray
    sensible_enthalpy_density_j_m3: np.ndarray
    thermodynamic_state: BinaryIdealGasMixtureState
    exact_geometry_volume_m3: float
    represented_geometry_volume_m3: float
    relative_geometry_volume_error: float
    integrated_species_mass_kg: float
    exact_species_mass_kg: float
    relative_species_mass_error: float


def _clip_polygon_axis(
    polygon: list[tuple[float, float]],
    *,
    axis: int,
    bound: float,
    keep_greater: bool,
) -> list[tuple[float, float]]:
    if not polygon:
        return []
    result: list[tuple[float, float]] = []

    def inside(point: tuple[float, float]) -> bool:
        return point[axis] >= bound if keep_greater else point[axis] <= bound

    previous = polygon[-1]
    previous_inside = inside(previous)
    for current in polygon:
        current_inside = inside(current)
        if current_inside != previous_inside:
            denominator = current[axis] - previous[axis]
            if abs(denominator) > 1.0e-15:
                fraction = (bound - previous[axis]) / denominator
                intersection = (
                    previous[0] + fraction * (current[0] - previous[0]),
                    previous[1] + fraction * (current[1] - previous[1]),
                )
                result.append(intersection)
        if current_inside:
            result.append(current)
        previous = current
        previous_inside = current_inside
    return result


def _polygon_rectangle_intersection_area(
    vertices: np.ndarray,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
) -> float:
    polygon = [tuple(map(float, row)) for row in vertices]
    polygon = _clip_polygon_axis(
        polygon, axis=0, bound=x_min, keep_greater=True,
    )
    polygon = _clip_polygon_axis(
        polygon, axis=0, bound=x_max, keep_greater=False,
    )
    polygon = _clip_polygon_axis(
        polygon, axis=1, bound=y_min, keep_greater=True,
    )
    polygon = _clip_polygon_axis(
        polygon, axis=1, bound=y_max, keep_greater=False,
    )
    if len(polygon) < 3:
        return 0.0
    x = np.asarray([point[0] for point in polygon])
    y = np.asarray([point[1] for point in polygon])
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _uniform_spacing(coordinates: np.ndarray, name: str) -> float:
    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 1 or len(values) < 2 or np.any(~np.isfinite(values)):
        raise ValueError(f"{name} must be a finite one-dimensional grid")
    spacing = float(values[1] - values[0])
    if spacing <= 0.0 or not np.allclose(
        np.diff(values), spacing, rtol=1.0e-12, atol=1.0e-14,
    ):
        raise ValueError(f"{name} must be strictly increasing and uniform")
    return spacing


def initialise_regular_polygon_prism_cloud(
    x_m: np.ndarray,
    y_m: np.ndarray,
    z_m: np.ndarray,
    *,
    cloud: RegularPolygonPrismInitialCloud,
    closure: BinaryIdealGasMixtureClosure,
    solid_mask: np.ndarray | None = None,
) -> FiniteVolumeInitialCloud3D:
    """Create conservative cell averages using exact polygon/box intersections."""
    if not isinstance(cloud, RegularPolygonPrismInitialCloud):
        raise TypeError("cloud must be RegularPolygonPrismInitialCloud")
    if not isinstance(closure, BinaryIdealGasMixtureClosure):
        raise TypeError("closure must be BinaryIdealGasMixtureClosure")
    x = np.asarray(x_m, dtype=float)
    y = np.asarray(y_m, dtype=float)
    z = np.asarray(z_m, dtype=float)
    dx = _uniform_spacing(x, "x_m")
    dy = _uniform_spacing(y, "y_m")
    dz = _uniform_spacing(z, "z_m")
    shape = (len(z), len(y), len(x))
    if solid_mask is None:
        solid = np.zeros(shape, dtype=bool)
    else:
        solid = np.asarray(solid_mask)
        if solid.shape != shape or solid.dtype != np.bool_:
            raise ValueError("solid_mask must be boolean with shape (z, y, x)")

    vertices = cloud.vertices_xy_m()
    horizontal_fraction = np.zeros((len(y), len(x)), dtype=float)
    polygon_x_min, polygon_y_min = np.min(vertices, axis=0)
    polygon_x_max, polygon_y_max = np.max(vertices, axis=0)
    candidate_i = np.flatnonzero(
        (x + 0.5 * dx > polygon_x_min) & (x - 0.5 * dx < polygon_x_max)
    )
    candidate_j = np.flatnonzero(
        (y + 0.5 * dy > polygon_y_min) & (y - 0.5 * dy < polygon_y_max)
    )
    for j in candidate_j:
        for i in candidate_i:
            area = _polygon_rectangle_intersection_area(
                vertices,
                x[i] - 0.5 * dx,
                x[i] + 0.5 * dx,
                y[j] - 0.5 * dy,
                y[j] + 0.5 * dy,
            )
            horizontal_fraction[j, i] = min(max(area / (dx * dy), 0.0), 1.0)
    z_overlap = np.maximum(
        0.0,
        np.minimum(z + 0.5 * dz, cloud.z_max_m)
        - np.maximum(z - 0.5 * dz, cloud.z_min_m),
    )
    volume_fraction = (
        z_overlap[:, None, None] / dz * horizontal_fraction[None, :, :]
    )
    if np.any((volume_fraction > 1.0e-14) & solid):
        raise ValueError("initial cloud intersects a solid obstacle cell")

    cloud_state = closure.state_from_mole_fraction_and_temperature(
        cloud.species_mole_fraction, cloud.temperature_k,
    )
    cloud_species_concentration = float(
        cloud_state.species_mass_concentration_kg_m3
    )
    cloud_enthalpy = float(cloud_state.sensible_enthalpy_density_j_m3)
    species_concentration = volume_fraction * cloud_species_concentration
    sensible_enthalpy = volume_fraction * cloud_enthalpy
    state = closure.state_from_species_mass_concentration_and_sensible_enthalpy(
        species_concentration, sensible_enthalpy,
    )
    cell_volume = dx * dy * dz
    represented_volume = float(np.sum(volume_fraction) * cell_volume)
    exact_volume = cloud.exact_volume_m3
    represented_species_mass = float(np.sum(species_concentration) * cell_volume)
    exact_species_mass = cloud_species_concentration * exact_volume
    return FiniteVolumeInitialCloud3D(
        volume_fraction=volume_fraction,
        species_mass_concentration_kg_m3=species_concentration,
        sensible_enthalpy_density_j_m3=sensible_enthalpy,
        thermodynamic_state=state,
        exact_geometry_volume_m3=exact_volume,
        represented_geometry_volume_m3=represented_volume,
        relative_geometry_volume_error=(represented_volume / exact_volume - 1.0),
        integrated_species_mass_kg=represented_species_mass,
        exact_species_mass_kg=exact_species_mass,
        relative_species_mass_error=(
            represented_species_mass / exact_species_mass - 1.0
            if exact_species_mass > 0.0 else 0.0
        ),
    )


__all__ = [
    "FiniteVolumeInitialCloud3D",
    "RegularPolygonPrismInitialCloud",
    "initialise_regular_polygon_prism_cloud",
]
