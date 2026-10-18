"""Local obstacle-patch geometry and conservative coarse/patch remapping."""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class ObstaclePatchBounds:
    """Axis-aligned local patch around an obstacle in metres."""

    x_min: float
    x_max: float
    y_min: float
    y_max: float
    z_min: float
    z_max: float

    @property
    def volume_m3(self) -> float:
        return (self.x_max - self.x_min) * (self.y_max - self.y_min) * (self.z_max - self.z_min)

def cartesian_face_coordinates(domain_m: tuple[float, float, float],
                               cells: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
    """Return x/y/z face coordinates using the obstacle solver convention."""
    if len(domain_m) != 3 or len(cells) != 3 or any(int(n) < 1 for n in cells):
        raise ValueError("domain_m and cells must contain three positive values")
    return tuple(np.linspace(0.0, float(length), int(count) + 1)
                 for length, count in zip(domain_m, cells))

@dataclass
class PatchFluxAccumulator:
    """Callable face-flux collector for one local patch."""
    x_face_m: np.ndarray
    y_face_m: np.ndarray
    z_face_m: np.ndarray
    patch: ObstaclePatchBounds

    def __post_init__(self) -> None:
        self.net_outward_kg_s: list[float] = []
        self.face_histories: list[dict[str, float]] = []

    def __call__(self, flux_x: np.ndarray, flux_y: np.ndarray, flux_z: np.ndarray) -> None:
        record = extract_patch_boundary_flux(
            flux_x=flux_x, flux_y=flux_y, flux_z=flux_z,
            x_face_m=self.x_face_m, y_face_m=self.y_face_m,
            z_face_m=self.z_face_m, patch=self.patch)
        self.face_histories.append(record)
        self.net_outward_kg_s.append(record["net_outward_kg_s"])

    @property
    def rates_kg_s(self) -> np.ndarray:
        """Return the accumulated net outward rates as a numeric history."""
        return np.asarray(self.net_outward_kg_s, dtype=float)


def make_obstacle_patch_bounds(
    bounds: tuple[float, float, float, float, float, float],
    *,
    upstream_m: float,
    downstream_m: float,
    lateral_m: float,
    vertical_m: float,
    domain_m: tuple[float, float, float],
) -> ObstaclePatchBounds:
    """Expand obstacle bounds and clip to a domain without empty dimensions."""
    values = tuple(float(value) for value in bounds)
    if len(values) != 6 or any(not math.isfinite(value) for value in values):
        raise ValueError("bounds must contain six finite values")
    x0, x1, y0, y1, z0, z1 = values
    if not (x1 > x0 and y1 > y0 and z1 > z0):
        raise ValueError("obstacle bounds must be strictly increasing")
    up, down, lat, vert = map(float, (upstream_m, downstream_m, lateral_m, vertical_m))
    if min(up, down, lat, vert) < 0.0 or not all(math.isfinite(v) for v in (up, down, lat, vert)):
        raise ValueError("patch buffers must be finite and non-negative")
    Lx, Ly, Lz = map(float, domain_m)
    if min(Lx, Ly, Lz) <= 0.0:
        raise ValueError("domain dimensions must be positive")
    return ObstaclePatchBounds(
        max(0.0, x0 - up), min(Lx, x1 + down),
        max(0.0, y0 - lat), min(Ly, y1 + lat),
        max(0.0, z0 - vert), min(Lz, z1 + vert),
    )


def conservative_patch_remap(
    coarse_mass_kg: np.ndarray,
    *,
    coarse_cell_volume_m3: float,
    patch_shape: tuple[int, int, int],
) -> np.ndarray:
    """Create a uniform patch field preserving total coarse mass."""
    coarse = np.asarray(coarse_mass_kg, dtype=float)
    if np.any(~np.isfinite(coarse)) or np.any(coarse < 0.0):
        raise ValueError("coarse mass must be finite and non-negative")
    volume = float(coarse_cell_volume_m3)
    if not math.isfinite(volume) or volume <= 0.0:
        raise ValueError("coarse_cell_volume_m3 must be positive")
    if len(patch_shape) != 3 or any(int(n) < 1 for n in patch_shape):
        raise ValueError("patch_shape must contain three positive sizes")
    total = float(np.sum(coarse))
    return np.full(tuple(int(n) for n in patch_shape), total / np.prod(patch_shape), dtype=float)

def overlap_mass_remap_to_patch(
    coarse_mass_kg: np.ndarray, *, x_m: np.ndarray, y_m: np.ndarray, z_m: np.ndarray,
    patch: ObstaclePatchBounds, patch_shape: tuple[int, int, int],
) -> tuple[np.ndarray, float]:
    """Remap only coarse mass at cell centres inside the local patch."""
    mass = np.asarray(coarse_mass_kg, dtype=float)
    x, y, z = (np.asarray(v, dtype=float) for v in (x_m, y_m, z_m))
    if mass.shape != (len(z), len(y), len(x)):
        raise ValueError("coarse mass shape must be (z, y, x) matching coordinate axes")
    if not np.all(np.isfinite(mass)) or np.any(mass < 0):
        raise ValueError("coarse mass must be finite and non-negative")
    mask = ((x[None, None, :] >= patch.x_min) & (x[None, None, :] <= patch.x_max)
            & (y[None, :, None] >= patch.y_min) & (y[None, :, None] <= patch.y_max)
            & (z[:, None, None] >= patch.z_min) & (z[:, None, None] <= patch.z_max))
    selected = float(np.sum(mass[mask]))
    shape = tuple(int(n) for n in patch_shape)
    if len(shape) != 3 or any(n < 1 for n in shape):
        raise ValueError("patch_shape must contain three positive sizes")
    return np.full(shape, selected / np.prod(shape), dtype=float), selected

def extract_patch_boundary_flux(*, flux_x: np.ndarray, flux_y: np.ndarray,
                                flux_z: np.ndarray, x_face_m: np.ndarray,
                                y_face_m: np.ndarray, z_face_m: np.ndarray,
                                patch: ObstaclePatchBounds,
                                face_tolerance_m: float = 1e-10) -> dict[str, float]:
    """Integrate scalar face fluxes on the six local-patch boundary planes."""
    fx, fy, fz = (np.asarray(v, dtype=float) for v in (flux_x, flux_y, flux_z))
    xf, yf, zf = (np.asarray(v, dtype=float) for v in (x_face_m, y_face_m, z_face_m))
    if fx.ndim != 3 or fy.ndim != 3 or fz.ndim != 3:
        raise ValueError("face flux arrays must be three-dimensional")
    tol = float(face_tolerance_m)
    if not math.isfinite(tol) or tol < 0:
        raise ValueError("face_tolerance_m must be finite and non-negative")
    def plane(values, target):
        indices = np.flatnonzero(np.isclose(values, target, rtol=0, atol=tol))
        if len(indices) != 1:
            raise ValueError("patch boundary does not coincide with exactly one face plane")
        return int(indices[0])
    ix0, ix1 = plane(xf, patch.x_min), plane(xf, patch.x_max)
    iy0, iy1 = plane(yf, patch.y_min), plane(yf, patch.y_max)
    iz0, iz1 = plane(zf, patch.z_min), plane(zf, patch.z_max)
    # Flux arrays are mass/area/time. Face areas are inferred from coordinate spacing.
    dy = float(np.median(np.diff(yf))); dz = float(np.median(np.diff(zf)));
    dx = float(np.median(np.diff(xf)))
    xm = (xf[:-1] >= patch.x_min) & (xf[:-1] < patch.x_max)
    ym = (yf[:-1] >= patch.y_min) & (yf[:-1] < patch.y_max)
    zm = (zf[:-1] >= patch.z_min) & (zf[:-1] < patch.z_max)
    result = {
        "x_min_kg_s": float(np.sum(fx[np.ix_(zm, ym, [ix0])]) * dy * dz),
        "x_max_kg_s": float(np.sum(fx[np.ix_(zm, ym, [ix1])]) * dy * dz),
        "y_min_kg_s": float(np.sum(fy[np.ix_(zm, [iy0], xm)]) * dx * dz),
        "y_max_kg_s": float(np.sum(fy[np.ix_(zm, [iy1], xm)]) * dx * dz),
        "z_min_kg_s": float(np.sum(fz[np.ix_([iz0], ym, xm)]) * dx * dy),
        "z_max_kg_s": float(np.sum(fz[np.ix_([iz1], ym, xm)]) * dx * dy),
    }
    result["net_outward_kg_s"] = (result["x_max_kg_s"] - result["x_min_kg_s"]
        + result["y_max_kg_s"] - result["y_min_kg_s"]
        + result["z_max_kg_s"] - result["z_min_kg_s"])
    return result


__all__ = ["ObstaclePatchBounds", "PatchFluxAccumulator", "cartesian_face_coordinates",
           "make_obstacle_patch_bounds",
           "conservative_patch_remap", "overlap_mass_remap_to_patch"]
__all__.append("extract_patch_boundary_flux")
