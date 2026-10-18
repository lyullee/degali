"""Calculate and render a pressure-projected continuous H2 wall scenario.

The authoritative calculation uses DEGALI's variable-density momentum,
pressure-projection and thermal/species corrector path.  The barrier therefore
changes the carrier velocity and pressure field; it is not merely a mask on
hydrogen concentration.  Spatial upsampling and Gaussian smoothing remain
output-only display operations and never feed back into the solver.

The result is a coarse-grid, physics-bearing research demonstration.  It is
not external experimental validation and not LES/RANS CFD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.ndimage import gaussian_filter, zoom

from degali.addons.continuous_barrier_resolved_3d import (
    audit_continuous_barrier_physics,
    build_continuous_barrier_resolved_config,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedObstacleFlow3DConfig,
    solve_resolved_obstacle_flow_3d,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs" / "continuous-barrier-3d-animation"
MOLAR_MASS_H2_KG_MOL = 2.01588e-3
MOLAR_MASS_AIR_KG_MOL = 28.96546e-3

CAMERA_PRESETS: dict[str, dict[str, float | str]] = {
    "upwind-oblique": {
        "elev": 27.0,
        "azim": -66.0,
        "label": "upwind oblique",
    },
    "side-profile": {
        "elev": 16.0,
        "azim": -90.0,
        "label": "crosswind side profile",
    },
    "elevated-leeward": {
        "elev": 38.0,
        "azim": 116.0,
        "label": "elevated leeward",
    },
    "wall-face-elevated": {
        "elev": 34.0,
        "azim": -48.0,
        "label": "elevated windward wall face",
    },
}
DEFAULT_CAMERA_PRESETS = (
    "upwind-oblique",
    "side-profile",
    "wall-face-elevated",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _barrier_segments(barrier: dict[str, object]) -> list[dict[str, float]]:
    segments = barrier.get("segments")
    if segments is None:
        return [{
            "y_min_m": float(barrier["y_min_m"]),
            "y_max_m": float(barrier["y_max_m"]),
        }]
    return [
        {
            "y_min_m": float(segment["y_min_m"]),
            "y_max_m": float(segment["y_max_m"]),
        }
        for segment in segments
    ]


def _obstacle_blocks(barrier: dict[str, object]) -> list[dict[str, float | str]]:
    buildings = barrier.get("buildings")
    if buildings is not None:
        return [
            {
                "id": str(block.get("id", f"building-{index + 1}")),
                "x_m": float(block["x_m"]),
                "thickness_m": float(block["thickness_m"]),
                "y_min_m": float(block["y_min_m"]),
                "y_max_m": float(block["y_max_m"]),
                "height_m": float(block["height_m"]),
            }
            for index, block in enumerate(buildings)
        ]
    return [
        {
            "id": f"wall-{index + 1}",
            "x_m": float(barrier["x_m"]),
            "thickness_m": float(barrier["thickness_m"]),
            "y_min_m": segment["y_min_m"],
            "y_max_m": segment["y_max_m"],
            "height_m": float(barrier["height_m"]),
        }
        for index, segment in enumerate(_barrier_segments(barrier))
    ]


def _barrier_mask(
    domain_m: tuple[float, float, float],
    cells: tuple[int, int, int],
    *,
    x_m: float,
    thickness_m: float,
    y_min_m: float,
    y_max_m: float,
    height_m: float,
    segments: list[dict[str, float]] | None = None,
    buildings: list[dict[str, float | str]] | None = None,
    **_metadata: object,
) -> np.ndarray:
    nx, ny, nz = cells
    lx, ly, lz = domain_m
    x = (np.arange(nx) + 0.5) * lx / nx
    y = (np.arange(ny) + 0.5) * ly / ny
    z = (np.arange(nz) + 0.5) * lz / nz
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    wall = np.zeros_like(xx, dtype=bool)
    blocks = buildings or [{
        "x_m": x_m,
        "thickness_m": thickness_m,
        "y_min_m": y_min_m,
        "y_max_m": y_max_m,
        "height_m": height_m,
    } for _segment in (segments or [{"y_min_m": y_min_m, "y_max_m": y_max_m}])]
    if buildings is None and segments is not None:
        blocks = [{
            "x_m": x_m,
            "thickness_m": thickness_m,
            "y_min_m": segment["y_min_m"],
            "y_max_m": segment["y_max_m"],
            "height_m": height_m,
        } for segment in segments]
    for block in blocks:
        wall |= (
            (xx >= float(block["x_m"]))
            & (xx <= float(block["x_m"]) + float(block["thickness_m"]))
            & (yy >= float(block["y_min_m"]))
            & (yy <= float(block["y_max_m"]))
            & (zz <= float(block["height_m"]))
        )
    return wall


def _config(
    args: argparse.Namespace,
) -> tuple[ResolvedObstacleFlow3DConfig, dict[str, object]]:
    return build_continuous_barrier_resolved_config(
        duration_s=args.duration_s,
        release_duration_s=args.release_duration_s,
        release_rate_kg_s=args.release_rate_kg_s,
        wind_speed_m_s=args.wind_speed_m_s,
        barrier_layout=getattr(args, "barrier_layout", "single"),
        store_scalar_fields=True,
    )


def _mass_concentration_to_vol_percent(
    h2_mass_concentration_kg_m3: np.ndarray,
    config: ResolvedObstacleFlow3DConfig,
) -> np.ndarray:
    # This exactly matches HydrogenAirActiveScalarClosure: the active-scalar
    # carrier inventory is the declared ambient dry-air density per cell.
    concentration = np.maximum(
        np.asarray(h2_mass_concentration_kg_m3, dtype=float), 0.0,
    )
    h2_moles = concentration / MOLAR_MASS_H2_KG_MOL
    air_moles = config.ambient_air_density_kg_m3 / MOLAR_MASS_AIR_KG_MOL
    return 100.0 * np.divide(
        h2_moles,
        h2_moles + air_moles,
        out=np.zeros_like(h2_moles),
        where=(h2_moles + air_moles) > 0.0,
    )


def _display_field(
    field_kg_m3: np.ndarray,
    config: ResolvedObstacleFlow3DConfig,
    *,
    upscale: float,
    smooth_sigma_cells: float,
    barrier: dict[str, object] | None = None,
) -> np.ndarray:
    concentration = _mass_concentration_to_vol_percent(field_kg_m3, config)
    if upscale != 1.0:
        concentration = zoom(concentration, upscale, order=1, mode="nearest")
    linear_concentration = concentration.copy()
    solid = None
    wall_shadow_guard = None
    if barrier is not None:
        nz, ny, nx = concentration.shape
        solid = _barrier_mask(
            config.domain_m,
            (nx, ny, nz),
            **barrier,
        )
        # Geometry-aware display filtering starts from an exactly empty solid.
        # With the truncated kernel below, a two-display-cell-thick wall blocks
        # any direct upstream-to-downstream smoothing, while the central gap
        # remains a connected fluid corridor.
        concentration[solid] = 0.0
        linear_concentration[solid] = 0.0
        display_dx = config.domain_m[0] / nx
        guard_depth_m = max(config.domain_m[0] / config.cells[0], 2.0 * display_dx)
        x = (np.arange(nx) + 0.5) * display_dx
        y = (np.arange(ny) + 0.5) * config.domain_m[1] / ny
        z = (np.arange(nz) + 0.5) * config.domain_m[2] / nz
        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        wall_shadow_guard = np.zeros_like(solid)
        for block in _obstacle_blocks(barrier):
            wall_shadow_guard |= (
                (xx > float(block["x_m"]) + float(block["thickness_m"]))
                & (
                    xx
                    <= float(block["x_m"])
                    + float(block["thickness_m"])
                    + guard_depth_m
                )
                & (yy >= float(block["y_min_m"]))
                & (yy <= float(block["y_max_m"]))
                & (zz <= float(block["height_m"]))
            )
    if smooth_sigma_cells > 0.0:
        concentration = gaussian_filter(
            concentration,
            smooth_sigma_cells,
            mode="nearest",
            truncate=2.0,
        )
    concentration = np.clip(concentration, 0.0, 100.0).astype(np.float32)
    if wall_shadow_guard is not None:
        # A display blur may otherwise raise a real but tiny downstream value
        # by borrowing colour from the wall's upstream face.  In the first
        # physical cell behind each solid segment, never show more than the
        # linearly interpolated raw solver field.  True gap/overtop transport
        # remains visible because it is already present in that raw field.
        concentration[wall_shadow_guard] = np.minimum(
            concentration[wall_shadow_guard],
            linear_concentration[wall_shadow_guard],
        )
    if solid is not None:
        # Interpolation/smoothing is a display operation and must not paint
        # concentration inside an impermeable solid.  Reapply the exact
        # physical obstacle after all visual filtering.
        concentration[solid] = 0.0
    return concentration


def _temporal_upsample(
    snapshots: list[np.ndarray],
    snapshot_time_s: np.ndarray,
    factor: int,
) -> tuple[list[np.ndarray], np.ndarray]:
    """Linearly add display frames without altering solver snapshots."""
    if factor < 1:
        raise ValueError("temporal interpolation factor must be at least one")
    if factor == 1:
        return snapshots, np.asarray(snapshot_time_s, dtype=float)
    expanded_fields: list[np.ndarray] = []
    expanded_times: list[float] = []
    for index in range(len(snapshots) - 1):
        left = snapshots[index]
        right = snapshots[index + 1]
        for substep in range(factor):
            fraction = substep / factor
            expanded_fields.append(
                ((1.0 - fraction) * left + fraction * right).astype(np.float32)
            )
            expanded_times.append(
                float(
                    (1.0 - fraction) * snapshot_time_s[index]
                    + fraction * snapshot_time_s[index + 1]
                )
            )
    expanded_fields.append(snapshots[-1])
    expanded_times.append(float(snapshot_time_s[-1]))
    return expanded_fields, np.asarray(expanded_times, dtype=float)


def _cloud_isosurface_triangles(
    field_vol_percent: np.ndarray,
    domain_m: tuple[float, float, float],
    *,
    level_vol_percent: float,
    vertical_exaggeration: float,
) -> np.ndarray:
    """Return a closed, smooth display mesh in matplotlib x-y-z order."""
    from skimage.measure import marching_cubes

    field = np.asarray(field_vol_percent, dtype=float)
    if field.ndim != 3 or min(field.shape) < 2:
        raise ValueError("isosurface field must be a three-dimensional grid")
    level = float(level_vol_percent)
    if not np.isfinite(level) or level <= 0.0:
        raise ValueError("isosurface level must be positive and finite")
    if float(np.max(field)) < level:
        return np.empty((0, 3, 3), dtype=float)
    nz, ny, nx = field.shape
    lx, ly, lz = (float(value) for value in domain_m)
    dx, dy, dz = lx / nx, ly / ny, lz / nz
    # The zero padding closes clouds that meet an open domain boundary.  Its
    # origin represents a ghost-cell centre at -0.5 cell, hence the offset.
    padded = np.pad(field, 1, mode="constant", constant_values=0.0)
    vertices_zyx, faces, _normals, _values = marching_cubes(
        padded,
        level=level,
        spacing=(dz, dy, dx),
        allow_degenerate=False,
    )
    vertices = np.column_stack((
        vertices_zyx[:, 2] - 0.5 * dx,
        vertices_zyx[:, 1] - 0.5 * dy,
        (vertices_zyx[:, 0] - 0.5 * dz) * vertical_exaggeration,
    ))
    return vertices[faces]


def _split_gap_diagnostics(
    snapshots: list[np.ndarray],
    config: ResolvedObstacleFlow3DConfig,
    barrier: dict[str, object],
) -> dict[str, object]:
    """Report raw-solver passage at the first cell plane behind the opening."""
    if barrier["layout"] != "split-gap":
        return {}
    nx, ny, nz = config.cells
    dx = config.domain_m[0] / nx
    dy = config.domain_m[1] / ny
    dz = config.domain_m[2] / nz
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    z = (np.arange(nz) + 0.5) * dz
    downstream_indices = np.flatnonzero(
        x > float(barrier["x_m"]) + float(barrier["thickness_m"])
    )
    ix = int(downstream_indices[0])
    gap_y = (y > float(barrier["gap_y_min_m"])) & (
        y < float(barrier["gap_y_max_m"])
    )
    below_top = z <= float(barrier["height_m"])
    raw = np.stack(snapshots, axis=0)
    gap_plane = raw[:, below_top, :, ix][:, :, gap_y]
    inventory = np.sum(gap_plane, axis=(1, 2)) * dy * dz
    return {
        "first_downstream_plane_x_m": float(x[ix]),
        "gap_solver_cells_y": int(np.count_nonzero(gap_y)),
        "gap_below_wall_top_peak_concentration_kg_m3": float(np.max(gap_plane)),
        "gap_below_wall_top_peak_plane_inventory_kg_m": float(np.max(inventory)),
        "gap_below_wall_top_final_plane_inventory_kg_m": float(inventory[-1]),
        "interpretation": (
            "Non-zero raw-solver concentration behind the opening demonstrates "
            "transport through the fluid gap; this is not inferred from display interpolation."
        ),
    }


def _render(
    snapshots: list[np.ndarray],
    snapshot_time_s: np.ndarray,
    config: ResolvedObstacleFlow3DConfig,
    barrier: dict[str, object],
    args: argparse.Namespace,
) -> tuple[Path, Path, dict[str, float]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation, colors
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    first = _display_field(
        snapshots[0], config, upscale=args.render_upscale,
        smooth_sigma_cells=args.smooth_sigma_cells,
        barrier=barrier,
    )
    nz, ny, nx = first.shape
    x = np.linspace(0.5 * config.domain_m[0] / nx, config.domain_m[0] - 0.5 * config.domain_m[0] / nx, nx)
    y = np.linspace(0.5 * config.domain_m[1] / ny, config.domain_m[1] - 0.5 * config.domain_m[1] / ny, ny)
    z = np.linspace(0.5 * config.domain_m[2] / nz, config.domain_m[2] - 0.5 * config.domain_m[2] / nz, nz)
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    px, py, pz = xx.ravel(), yy.ravel(), zz.ravel()

    maxima = []
    display_counts = []
    wall_face_counts = []
    windward_ix = int(np.flatnonzero(x < float(barrier["x_m"]))[-1])
    wall_height_indices = np.flatnonzero(z <= float(barrier["height_m"]))
    for raw in snapshots:
        display = _display_field(
            raw, config, upscale=args.render_upscale,
            smooth_sigma_cells=args.smooth_sigma_cells,
            barrier=barrier,
        )
        maxima.append(float(display.max()))
        display_counts.append(int(np.count_nonzero(display >= args.display_threshold_vol_percent)))
        wall_face_count = 0
        for segment in _barrier_segments(barrier):
            segment_y = np.flatnonzero(
                (y >= segment["y_min_m"]) & (y <= segment["y_max_m"])
            )
            wall_face_count += int(
                np.count_nonzero(
                    display[np.ix_(wall_height_indices, segment_y, [windward_ix])]
                    >= args.display_threshold_vol_percent
                )
            )
        wall_face_counts.append(wall_face_count)
    positive_maxima = [value for value in maxima if value > 0.0]
    if getattr(args, "cloud_style", "points") == "isosurface":
        vmax = max(
            args.display_threshold_vol_percent * 1.01,
            float(np.percentile(positive_maxima or [1.0], 95)),
        )
    else:
        vmax = max(4.0, float(np.percentile(positive_maxima or [4.0], 95)))
    norm = colors.LogNorm(vmin=args.display_threshold_vol_percent, vmax=vmax)
    cmap = plt.get_cmap("turbo")
    surface_levels = (
        args.display_threshold_vol_percent,
        math.sqrt(args.display_threshold_vol_percent * vmax),
        max(args.display_threshold_vol_percent * 1.02, 0.62 * vmax),
    )
    surface_alphas = (0.34, 0.52, 0.78)

    fig = plt.figure(figsize=(11.2, 7.0), facecolor="#07111f")
    # Matplotlib's automatic 3-D collection sorting can draw transparent
    # upwind points on top of an opaque wall in a leeward view.  Explicit
    # ordering is required here because that artefact looks exactly like gas
    # penetrating the solid.
    axis = fig.add_subplot(111, projection="3d", computed_zorder=False)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(scalar, ax=axis, fraction=0.025, pad=0.06)
    colorbar.set_label("Display-interpolated H$_2$ concentration (vol%)", color="white")
    colorbar.ax.tick_params(colors="white")
    colorbar.outline.set_edgecolor("#94a3b8")

    def draw(index: int) -> list[object]:
        axis.clear()
        # A filled translucent ground surface suffers from the same collection
        # sorting problem and can appear above low gas.  A z=0 reference grid
        # communicates the floor without hiding the plume.
        for ground_x in np.linspace(0.0, config.domain_m[0], 8):
            axis.plot(
                [ground_x, ground_x], [0.0, config.domain_m[1]], [0.0, 0.0],
                color="#64748b", alpha=0.28, linewidth=0.55, zorder=0,
            )
        for ground_y in np.linspace(0.0, config.domain_m[1], 6):
            axis.plot(
                [0.0, config.domain_m[0]], [ground_y, ground_y], [0.0, 0.0],
                color="#64748b", alpha=0.28, linewidth=0.55, zorder=0,
            )
        field = _display_field(
            snapshots[index], config, upscale=args.render_upscale,
            smooth_sigma_cells=args.smooth_sigma_cells,
            barrier=barrier,
        )
        flat = field.ravel()
        selected = np.flatnonzero(flat >= args.display_threshold_vol_percent)
        if selected.size > args.max_points:
            stride = int(math.ceil(selected.size / args.max_points))
            selected = selected[::stride]
        wall_x0 = float(barrier["x_m"])
        wall_x1 = wall_x0 + float(barrier["thickness_m"])

        def cluster_depth_zorder(x_position_m: float) -> float:
            if args.camera_preset == "elevated-leeward":
                return 100.0 + x_position_m
            if args.camera_preset == "upwind-oblique":
                return 300.0 - x_position_m
            # Side profile is an x-z projection; keep all opaque buildings in
            # front of the integrated crosswind point cloud.
            return 100.0

        def scatter_cells(indices: np.ndarray, *, zorder: float) -> None:
            if not indices.size:
                return
            concentration = flat[indices]
            rgba = cmap(norm(concentration))
            alpha_min = float(getattr(args, "gas_alpha_min", 0.32))
            alpha_max = float(getattr(args, "gas_alpha_max", 0.96))
            rgba[:, 3] = np.clip(
                alpha_min + (alpha_max - alpha_min) * norm(concentration),
                alpha_min,
                alpha_max,
            )
            axis.scatter(
                px[indices], py[indices],
                pz[indices] * args.vertical_exaggeration,
                c=rgba,
                s=float(getattr(args, "gas_point_size", 9.0)),
                depthshade=False, linewidths=0.0,
                zorder=zorder,
            )

        near_lee = np.asarray([], dtype=int)
        if getattr(args, "cloud_style", "points") == "isosurface":
            for surface_index, (level, alpha) in enumerate(zip(
                surface_levels, surface_alphas,
            )):
                triangles = _cloud_isosurface_triangles(
                    field,
                    config.domain_m,
                    level_vol_percent=level,
                    vertical_exaggeration=args.vertical_exaggeration,
                )
                if not triangles.size:
                    continue
                cloud = Poly3DCollection(
                    triangles,
                    facecolor=cmap(norm(level)),
                    edgecolor="none",
                    linewidth=0.0,
                    alpha=alpha,
                    antialiased=True,
                    zorder=2.0 + 0.2 * surface_index,
                )
                cloud.set_zsort("average")
                axis.add_collection3d(cloud)
        elif (
            selected.size
            and args.camera_preset == "elevated-leeward"
            and barrier["layout"] != "building-cluster"
        ):
            # Camera is in the lee (+x): draw the far/upwind cloud first, then
            # the wall, then only the genuinely downstream cloud.  This makes
            # the wall an opaque occluder instead of allowing an upstream
            # transparent point cloud to show through it.
            far_upwind = selected[px[selected] < wall_x1]
            near_lee = selected[px[selected] >= wall_x1]
            scatter_cells(far_upwind, zorder=2.0)
        elif selected.size and barrier["layout"] == "building-cluster":
            # One collection for the whole plume cannot be depth-sorted
            # reliably against several opaque Poly3D buildings.  Thin x slabs
            # allow deterministic front-to-back layering for the two oblique
            # views, so gas behind a building cannot be painted through it.
            for x_value in np.unique(px[selected]):
                slab = selected[np.isclose(px[selected], x_value)]
                scatter_cells(
                    slab,
                    zorder=(
                        50.0
                        if args.camera_preset == "side-profile"
                        else cluster_depth_zorder(float(x_value))
                    ),
                )
        elif selected.size:
            scatter_cells(selected, zorder=4.0)
        for block in _obstacle_blocks(barrier):
            block_x0 = float(block["x_m"])
            block_x1 = block_x0 + float(block["thickness_m"])
            wall_y0 = float(block["y_min_m"])
            wall_y1 = float(block["y_max_m"])
            block_zorder = (
                cluster_depth_zorder(0.5 * (block_x0 + block_x1))
                if barrier["layout"] == "building-cluster"
                else 5.0
            )
            axis.bar3d(
                [block_x0], [wall_y0], [0.0],
                [block["thickness_m"]],
                [wall_y1 - wall_y0],
                [float(block["height_m"]) * args.vertical_exaggeration],
                color=(
                    "#475569"
                    if barrier["layout"] == "building-cluster"
                    else "#4b5563"
                ),
                alpha=1.0, shade=True,
                edgecolor="#f8fafc", linewidth=0.55,
                zorder=block_zorder,
            )
            axis.plot(
                [block_x0, block_x1, block_x1, block_x0, block_x0],
                [wall_y0, wall_y0, wall_y1, wall_y1, wall_y0],
                [0.0] * 5,
                color="#f8fafc",
                linewidth=1.4,
                zorder=block_zorder + 0.1,
            )
            if (
                args.camera_preset == "wall-face-elevated"
                and barrier["layout"] != "building-cluster"
            ):
                segment_y = np.flatnonzero(
                    (y >= wall_y0) & (y <= wall_y1)
                )
                face_value = field[
                    np.ix_(wall_height_indices, segment_y, [windward_ix])
                ][:, :, 0]
                face_y, face_z = np.meshgrid(
                    y[segment_y], z[wall_height_indices], indexing="xy"
                )
                face_x = np.full_like(face_y, float(barrier["x_m"]) - 1.0e-3)
                face_rgba = cmap(norm(np.maximum(
                    face_value, args.display_threshold_vol_percent
                )))
                face_rgba[..., 3] = np.where(
                    face_value >= args.display_threshold_vol_percent, 0.98, 0.0
                )
                axis.plot_surface(
                    face_x,
                    face_y,
                    face_z * args.vertical_exaggeration,
                    facecolors=face_rgba,
                    shade=False,
                    linewidth=0.0,
                    antialiased=False,
                    zorder=6.0,
                )
        if near_lee.size:
            scatter_cells(near_lee, zorder=8.0)
        sx, sy, sz = config.source_position_m
        if barrier["layout"] == "building-cluster":
            source_zorder = cluster_depth_zorder(float(sx)) + 0.2
        else:
            source_zorder = (
                3.0 if args.camera_preset == "elevated-leeward" else 8.0
            )
        axis.scatter(
            [sx], [sy], [sz * args.vertical_exaggeration], marker="*", s=110,
            color="white", edgecolor="#111827", linewidth=0.7,
            zorder=source_zorder,
        )
        axis.plot(
            [sx, sx], [sy, sy], [0.0, sz * args.vertical_exaggeration],
            color="#f8fafc", linewidth=0.8, alpha=0.9,
            zorder=source_zorder,
        )
        axis.plot(
            [sx, barrier["x_m"]], [sy, sy], [0.0, 0.0],
            color="#fbbf24", linewidth=1.4, linestyle="--",
            zorder=source_zorder,
        )
        if args.camera_preset == "wall-face-elevated":
            axis.set_xlim(15.0, 145.0)
        else:
            x_limit = 175.0 if args.view_mode == "barrier-detail" else config.domain_m[0]
            axis.set_xlim(0.0, x_limit)
        axis.set_ylim(0.0, config.domain_m[1])
        axis.set_zlim(0.0, config.domain_m[2] * args.vertical_exaggeration)
        axis.set_xlabel("Downwind x (m)", color="white", labelpad=8)
        axis.set_ylabel("Crosswind y (m)", color="white", labelpad=8)
        axis.set_zlabel("Height z (m)", color="white", labelpad=8)
        if args.camera_preset == "side-profile":
            axis.set_yticks([])
            axis.set_ylabel("")
        z_ticks = np.linspace(0.0, config.domain_m[2], 5)
        axis.set_zticks(z_ticks * args.vertical_exaggeration)
        axis.set_zticklabels([f"{value:g}" for value in z_ticks])
        axis.tick_params(colors="#cbd5e1")
        camera = CAMERA_PRESETS[args.camera_preset]
        axis.view_init(elev=float(camera["elev"]), azim=float(camera["azim"]))
        axis.set_box_aspect(
            (1.24, 0.92, 0.62)
            if args.view_mode == "barrier-detail"
            else (1.45, 0.82, 0.55)
        )
        axis.set_facecolor("#07111f")
        for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
            pane.set_facecolor((0.035, 0.075, 0.13, 1.0))
            pane.set_edgecolor("#334155")
        axis.grid(False)
        if barrier["layout"] == "split-gap":
            title_geometry = "two protective barriers with central gap"
        elif barrier["layout"] == "blocking-wall":
            title_geometry = "wide protective wall with open side bypass"
        elif barrier["layout"] == "building-cluster":
            title_geometry = "staggered four-building cluster"
        else:
            title_geometry = "straight protective barrier"
        axis.set_title(
            f"DEGALI continuous H$_2$ release — {title_geometry}",
            color="white", pad=16, fontsize=15, weight="bold",
        )
        release_state = "continuous release" if snapshot_time_s[index] <= args.release_duration_s else "release stopped"
        axis.text2D(
            0.02, 0.94,
            f"physical time {snapshot_time_s[index]:6.1f} s   |   {release_state}   |   max {float(field.max()):6.2f} vol%",
            transform=axis.transAxes, color="#e2e8f0", fontsize=10,
        )
        if barrier["layout"] == "split-gap":
            geometry_note = (
                f"two {barrier['height_m']:g} m-high wall segments   |   "
                f"central opening {barrier['gap_width_m']:g} m"
            )
        elif barrier["layout"] == "blocking-wall":
            geometry_note = (
                f"impermeable wall {barrier['height_m']:g} m high × "
                f"{barrier['y_max_m'] - barrier['y_min_m']:g} m wide   |   "
                "resolved top and side bypass"
            )
        elif barrier["layout"] == "building-cluster":
            heights = [
                float(block["height_m"])
                for block in _obstacle_blocks(barrier)
            ]
            geometry_note = (
                f"{len(heights)} staggered buildings   |   "
                f"heights {min(heights):g}–{max(heights):g} m   |   "
                "open street canyons"
            )
        else:
            geometry_note = (
                f"wall {barrier['height_m']:g} m high × "
                f"{barrier['y_max_m'] - barrier['y_min_m']:g} m wide"
            )
        axis.text2D(
            0.02, 0.055,
            f"source-to-wall distance {barrier['x_m'] - sx:g} m   |   "
            f"{geometry_note}   |   view: {camera['label']}",
            transform=axis.transAxes, color="#cbd5e1", fontsize=8.5,
        )
        if (
            args.camera_preset == "wall-face-elevated"
            and barrier["layout"] != "building-cluster"
        ):
            axis.text2D(
                0.02,
                0.09,
                "coloured wall face = concentration in the adjacent upstream solver cell",
                transform=axis.transAxes,
                color="#e2e8f0",
                fontsize=8.5,
            )
        if args.camera_preset == "elevated-leeward":
            routes = getattr(args, "transport_routes", None)
            route_note = ""
            if routes:
                total = max(float(routes["total_downstream_kg_s"]), 1.0e-30)
                top_fraction = 100.0 * float(routes["over_wall_top_kg_s"]) / total
                side_fraction = (
                    100.0
                    * float(routes["around_sides_or_open_gap_kg_s"])
                    / total
                )
                route_note = (
                    f"   |   final x={float(routes['plane_x_m']):g} m flux: "
                    f"direct 0, over top {top_fraction:.1f}%, sides {side_fraction:.1f}%"
                )
            depth_note = (
                "rear depth layers: upwind cloud → opaque buildings → lee cloud; "
                "shown flux split is across the first building"
                if barrier["layout"] == "building-cluster"
                else (
                    "rear depth order: upwind cloud → opaque wall → lee cloud; "
                    "visible lee gas must arrive over the top or around an open side"
                )
            )
            axis.text2D(
                0.02,
                0.09,
                depth_note + route_note,
                transform=axis.transAxes,
                color="#e2e8f0",
                fontsize=8.0,
            )
        axis.text2D(
            0.02, 0.02,
            "VARIABLE-DENSITY MOMENTUM + PRESSURE PROJECTION — coarse-grid research result, not validated LES/RANS CFD",
            transform=axis.transAxes, color="#fbbf24", fontsize=8.7, weight="bold",
        )
        return []

    movie = animation.FuncAnimation(
        fig, draw, frames=len(snapshots), interval=1000.0 / args.fps, blit=False
    )
    geometry_slug = {
        "split-gap": "split_barrier_gap",
        "blocking-wall": "blocking_wall",
        "building-cluster": "building_cluster",
    }.get(str(barrier["layout"]), "straight_barrier")
    gif = output_dir / (
        f"degali_continuous_release_{geometry_slug}_{args.camera_preset}_3d.gif"
    )
    movie.save(gif, writer=animation.PillowWriter(fps=args.fps), dpi=args.dpi)
    poster_index = int(np.argmax(
        wall_face_counts
        if args.camera_preset == "wall-face-elevated"
        else display_counts
    ))
    draw(poster_index)
    poster = output_dir / (
        f"degali_continuous_release_{geometry_slug}_{args.camera_preset}_3d_poster.png"
    )
    fig.savefig(poster, dpi=args.dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return gif, poster, {
        "maximum_display_concentration_vol_percent": float(max(maxima)),
        "poster_time_s": float(snapshot_time_s[poster_index]),
        "render_grid_cells": [int(nx), int(ny), int(nz)],
        **(
            {
                "cloud_surface_levels_vol_percent": [
                    float(value) for value in surface_levels
                ],
                "cloud_surface_alphas": list(surface_alphas),
            }
            if getattr(args, "cloud_style", "points") == "isosurface"
            else {}
        ),
    }


def build(args: argparse.Namespace) -> dict[str, object]:
    config, barrier = _config(args)
    started = time.perf_counter()
    result = solve_resolved_obstacle_flow_3d(config)
    solve_runtime = time.perf_counter() - started
    if result.concentration_fields_kg_m3 is None:
        raise RuntimeError("resolved barrier run did not retain scalar fields")
    target_times = np.arange(
        0.0,
        config.duration_s + 0.5 * args.snapshot_interval_s,
        args.snapshot_interval_s,
        dtype=float,
    )
    if target_times.size == 0 or target_times[-1] < config.duration_s - 1.0e-12:
        target_times = np.append(target_times, config.duration_s)
    target_times = np.minimum(target_times, config.duration_s)
    snapshot_indices = np.unique(np.clip(
        np.searchsorted(result.time_s, target_times, side="left"),
        0,
        result.time_s.size - 1,
    ))
    snapshot_times = np.asarray(result.time_s[snapshot_indices], dtype=float)
    snapshots = [
        np.asarray(result.concentration_fields_kg_m3[index], dtype=np.float32)
        for index in snapshot_indices
    ]
    raw_snapshot_count = len(snapshots)
    physics_audit = audit_continuous_barrier_physics(result, barrier)
    if not physics_audit["all_internal_physics_checks_passed"]:
        failed = [
            name for name, passed in physics_audit["checks"].items() if not passed
        ]
        raise RuntimeError(
            "resolved barrier physics audit failed: " + ", ".join(failed)
        )
    split_gap_diagnostics = _split_gap_diagnostics(snapshots, config, barrier)
    snapshot_cache_path = args.output_dir.resolve() / "solver_snapshots_float32.npz"
    args.output_dir.resolve().mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        snapshot_cache_path,
        time_s=np.asarray(snapshot_times, dtype=float),
        fields_kg_m3=np.stack(snapshots, axis=0),
    )
    snapshots, display_times = _temporal_upsample(
        snapshots,
        snapshot_times,
        args.temporal_interpolation_factor,
    )
    if args.camera_preset == "all":
        camera_presets = (
            ("upwind-oblique", "side-profile", "elevated-leeward")
            if barrier["layout"] == "building-cluster"
            else DEFAULT_CAMERA_PRESETS
        )
    else:
        camera_presets = (args.camera_preset,)
    render_started = time.perf_counter()
    view_outputs: dict[str, dict[str, object]] = {}
    for camera_preset in camera_presets:
        render_args = argparse.Namespace(**vars(args))
        render_args.camera_preset = camera_preset
        render_args.transport_routes = physics_audit[
            "final_downstream_transport_routes"
        ]
        view_started = time.perf_counter()
        gif, poster, visual = _render(
            snapshots, display_times, config, barrier, render_args
        )
        view_outputs[camera_preset] = {
            "animation_gif": str(gif),
            "poster_png": str(poster),
            "runtime_s": time.perf_counter() - view_started,
            **visual,
        }
    render_runtime = time.perf_counter() - render_started

    source_to_wall = float(barrier["x_m"]) - config.source_position_m[0]
    metadata: dict[str, object] = {
        "schema": "degali.continuous-barrier-3d-visualization.v2",
        "classification": (
            "hypothetical_coarse_grid_variable_density_projection_scenario"
        ),
        "resolved_3d_finite_volume": True,
        "les_or_rans_cfd": False,
        "experimental_validation": False,
        "display_interpolation_changes_solver_state": False,
        "scenario": {
            "domain_m": list(config.domain_m),
            "solver_cells": list(config.cells),
            "duration_s": config.duration_s,
            "time_step_s": float(result.diagnostics["scalar_time_step_s"]),
            "source_position_m": list(config.source_position_m),
            "release_rate_kg_s": args.release_rate_kg_s,
            "release_duration_s": args.release_duration_s,
            "wind_speed_m_s": args.wind_speed_m_s,
            "barrier": barrier,
            "source_to_barrier_distance_m": source_to_wall,
        },
        "solver_audit": {
            "runtime_s": solve_runtime,
            "solver_snapshot_count": raw_snapshot_count,
            "injected_h2_mass_kg": float(result.injected_mass_kg[-1]),
            "outflow_h2_mass_kg": float(result.outflow_mass_kg[-1]),
            "remaining_h2_mass_kg": float(result.domain_mass_kg[-1]),
            "maximum_mass_residual_kg": float(result.maximum_mass_residual_kg),
            "maximum_momentum_advective_courant": float(
                result.diagnostics["maximum_momentum_advective_courant_number"]
            ),
            "maximum_scalar_advective_number": float(
                result.diagnostics["advective_number"]
            ),
            "maximum_scalar_diffusion_number": float(
                result.diagnostics["diffusion_number"]
            ),
            "solid_cell_count": int(np.count_nonzero(result.solid_mask)),
            "physics": physics_audit,
            **({"split_gap_transport": split_gap_diagnostics} if split_gap_diagnostics else {}),
        },
        "render": {
            "runtime_s": render_runtime,
            "fps": args.fps,
            "display_frame_count": len(snapshots),
            "temporal_interpolation_factor": args.temporal_interpolation_factor,
            "render_upscale": args.render_upscale,
            "smooth_sigma_cells": args.smooth_sigma_cells,
            "geometry_aware_smoothing": (
                "solid-mask-before-and-after-truncated-gaussian-plus-"
                "first-downstream-cell-raw-field-cap"
            ),
            "ground_reference_grid_rendered_at_z_m": 0.0,
            "elevated_leeward_depth_order": (
                "upwind_cloud_then_opaque_wall_then_downwind_cloud"
            ),
            **(
                {
                    "gap_display_cells_y": int(
                        round(
                            float(barrier["gap_width_m"])
                            / (config.domain_m[1] / (config.cells[1] * args.render_upscale))
                        )
                    )
                }
                if barrier["layout"] == "split-gap"
                else {}
            ),
            "display_threshold_vol_percent": args.display_threshold_vol_percent,
            "cloud_style": args.cloud_style,
            "gas_alpha_min": args.gas_alpha_min,
            "gas_alpha_max": args.gas_alpha_max,
            "gas_point_size": args.gas_point_size,
            "vertical_exaggeration": args.vertical_exaggeration,
            "view_mode": args.view_mode,
            "camera_presets": list(camera_presets),
            "views": view_outputs,
        },
        "implementation": {
            "solver_sha256": _sha256(
                ROOT / "src/degali/addons/resolved_obstacle_flow_3d.py"
            ),
            "scenario_sha256": _sha256(
                ROOT / "src/degali/addons/continuous_barrier_resolved_3d.py"
            ),
            "renderer_sha256": _sha256(Path(__file__).resolve()),
        },
        "outputs": {
            "animation_gif": {
                key: value["animation_gif"] for key, value in view_outputs.items()
            },
            "poster_png": {
                key: value["poster_png"] for key, value in view_outputs.items()
            },
            "solver_snapshot_cache_npz": str(snapshot_cache_path),
        },
    }
    metadata_path = args.output_dir.resolve() / "visualization_metadata.json"
    metadata["outputs"]["metadata_json"] = str(metadata_path)
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--duration-s", type=float, default=600.0)
    parser.add_argument("--release-duration-s", type=float, default=600.0)
    parser.add_argument("--release-rate-kg-s", type=float, default=0.10)
    parser.add_argument("--wind-speed-m-s", type=float, default=2.2)
    parser.add_argument(
        "--barrier-layout",
        choices=("single", "split-gap", "blocking-wall", "building-cluster"),
        default="single",
    )
    parser.add_argument("--snapshot-interval-s", type=float, default=5.0)
    parser.add_argument("--temporal-interpolation-factor", type=int, default=1)
    parser.add_argument("--render-upscale", type=float, default=1.5)
    parser.add_argument("--smooth-sigma-cells", type=float, default=0.55)
    parser.add_argument("--display-threshold-vol-percent", type=float, default=0.04)
    parser.add_argument(
        "--cloud-style", choices=("points", "isosurface"), default="points"
    )
    parser.add_argument("--gas-alpha-min", type=float, default=0.32)
    parser.add_argument("--gas-alpha-max", type=float, default=0.96)
    parser.add_argument("--gas-point-size", type=float, default=9.0)
    parser.add_argument("--vertical-exaggeration", type=float, default=5.0)
    parser.add_argument(
        "--view-mode", choices=("full", "barrier-detail"), default="full"
    )
    parser.add_argument(
        "--camera-preset",
        choices=(*CAMERA_PRESETS, "all"),
        default="upwind-oblique",
    )
    parser.add_argument("--max-points", type=int, default=10000)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=105)
    args = parser.parse_args()
    if not 0.0 < args.release_duration_s <= args.duration_s:
        parser.error("release-duration-s must lie inside or end with duration-s")
    if min(
        args.release_rate_kg_s,
        args.wind_speed_m_s,
        args.snapshot_interval_s,
        args.render_upscale,
        args.display_threshold_vol_percent,
        args.vertical_exaggeration,
    ) <= 0.0:
        parser.error("positive physical and rendering controls are required")
    if args.temporal_interpolation_factor < 1:
        parser.error("temporal-interpolation-factor must be at least one")
    if not 0.0 <= args.gas_alpha_min <= args.gas_alpha_max <= 1.0:
        parser.error("gas alpha must satisfy 0 <= min <= max <= 1")
    if args.gas_point_size <= 0.0:
        parser.error("gas-point-size must be positive")
    metadata = build(args)
    print(json.dumps(metadata, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
