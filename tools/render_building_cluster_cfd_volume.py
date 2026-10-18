"""Render the DEGALI building-cluster transient in a conventional CFD style.

This renderer deliberately separates calculation from presentation.  It reads
the saved, mass-conservative DEGALI concentration snapshots, converts them to
H2 volume percent, and applies display-only interpolation.  VTK then performs
composite ray-cast volume rendering with a fixed colour/opacity transfer
function.  A centreline vertical slice is included as the familiar CFD
contour view.  Neither interpolation nor rendering feeds back into the solver.

The result is CFD-style post-processing of a coarse DEGALI resolved field; it
is not a claim that the underlying calculation is LES or RANS CFD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import imageio.v2 as imageio
import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter, zoom

from degali.addons.continuous_barrier_resolved_3d import (
    build_continuous_barrier_resolved_config,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = (
    ROOT
    / "outputs"
    / "continuous-building-cluster-cloud-volume-2026-10-18"
    / "solver_snapshots_float32.npz"
)
DEFAULT_OUTPUT = (
    ROOT / "outputs" / "continuous-building-cluster-cfd-volume-2026-10-18"
)
MOLAR_MASS_H2_KG_MOL = 2.01588e-3
MOLAR_MASS_AIR_KG_MOL = 28.96546e-3


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _mass_concentration_to_vol_percent(
    concentration_kg_m3: np.ndarray,
    ambient_air_density_kg_m3: float,
) -> np.ndarray:
    concentration = np.maximum(np.asarray(concentration_kg_m3, dtype=float), 0.0)
    h2_moles = concentration / MOLAR_MASS_H2_KG_MOL
    air_moles = ambient_air_density_kg_m3 / MOLAR_MASS_AIR_KG_MOL
    return 100.0 * np.divide(
        h2_moles,
        h2_moles + air_moles,
        out=np.zeros_like(h2_moles),
        where=(h2_moles + air_moles) > 0.0,
    )


def _building_mask(
    shape_zyx: tuple[int, int, int],
    domain_m: tuple[float, float, float],
    buildings: tuple[dict[str, object], ...],
) -> np.ndarray:
    nz, ny, nx = shape_zyx
    lx, ly, lz = domain_m
    x = (np.arange(nx) + 0.5) * lx / nx
    y = (np.arange(ny) + 0.5) * ly / ny
    z = (np.arange(nz) + 0.5) * lz / nz
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    solid = np.zeros(shape_zyx, dtype=bool)
    for block in buildings:
        x0 = float(block["x_m"])
        x1 = x0 + float(block["thickness_m"])
        solid |= (
            (xx >= x0)
            & (xx <= x1)
            & (yy >= float(block["y_min_m"]))
            & (yy <= float(block["y_max_m"]))
            & (zz <= float(block["height_m"]))
        )
    return solid


def _display_field(
    field_kg_m3: np.ndarray,
    *,
    ambient_air_density_kg_m3: float,
    domain_m: tuple[float, float, float],
    buildings: tuple[dict[str, object], ...],
    upscale: int,
    smooth_sigma: float,
) -> np.ndarray:
    vol_percent = _mass_concentration_to_vol_percent(
        field_kg_m3,
        ambient_air_density_kg_m3,
    )
    if upscale > 1:
        vol_percent = zoom(vol_percent, upscale, order=1, mode="nearest")
    solid = _building_mask(vol_percent.shape, domain_m, buildings)
    vol_percent[solid] = 0.0
    if smooth_sigma > 0.0:
        # Normalized convolution retains near-wall gas while excluding solid
        # cells.  The buildings are several display cells thick, so this
        # compact filter does not bridge one side of a building to the other.
        fluid = (~solid).astype(float)
        numerator = gaussian_filter(
            vol_percent * fluid,
            smooth_sigma,
            mode="nearest",
            truncate=2.0,
        )
        denominator = gaussian_filter(
            fluid,
            smooth_sigma,
            mode="nearest",
            truncate=2.0,
        )
        vol_percent = np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator),
            where=denominator > 1.0e-8,
        )
        vol_percent[solid] = 0.0
    return np.asarray(vol_percent, dtype=np.float32)


def _opacity_transfer(
    *, threshold_vol_percent: float, maximum_vol_percent: float,
) -> np.ndarray:
    samples = np.linspace(0.0, maximum_vol_percent, 256)
    # Low concentrations remain visible but transparent; the compact core is
    # sufficiently opaque to show impingement and accumulation at buildings.
    anchors_x = np.asarray([
        0.0,
        threshold_vol_percent * 0.999,
        threshold_vol_percent,
        0.04,
        0.08,
        0.16,
        0.28,
        maximum_vol_percent,
    ])
    anchors_a = np.asarray([0.0, 0.0, 0.018, 0.035, 0.075, 0.16, 0.31, 0.50])
    # PyVista/VTK expects a 256-entry lookup transfer in the uint8 range.
    return np.rint(255.0 * np.interp(samples, anchors_x, anchors_a)).astype(np.uint8)


def _image_grid(
    field: np.ndarray,
    domain_m: tuple[float, float, float],
) -> pv.ImageData:
    nz, ny, nx = field.shape
    lx, ly, lz = domain_m
    spacing = (lx / nx, ly / ny, lz / nz)
    grid = pv.ImageData(
        dimensions=(nx, ny, nz),
        spacing=spacing,
        origin=(0.5 * spacing[0], 0.5 * spacing[1], 0.5 * spacing[2]),
    )
    grid.point_data["H2_vol_percent"] = field.ravel(order="C")
    return grid


def _interpolated_frames(
    times_s: np.ndarray,
    fields_kg_m3: np.ndarray,
    temporal_factor: int,
):
    if temporal_factor <= 1:
        for time_s, field in zip(times_s, fields_kg_m3, strict=True):
            yield float(time_s), np.asarray(field)
        return
    for index in range(len(times_s) - 1):
        for step in range(temporal_factor):
            weight = step / temporal_factor
            yield (
                float((1.0 - weight) * times_s[index] + weight * times_s[index + 1]),
                (1.0 - weight) * fields_kg_m3[index] + weight * fields_kg_m3[index + 1],
            )
    yield float(times_s[-1]), np.asarray(fields_kg_m3[-1])


def _add_scene_geometry(
    plotter: pv.Plotter,
    *,
    domain_m: tuple[float, float, float],
    buildings: tuple[dict[str, object], ...],
    source_position_m: tuple[float, float, float],
    wind_speed_m_s: float,
) -> None:
    lx, ly, _lz = domain_m
    ground = pv.Plane(
        center=(0.5 * lx, 0.5 * ly, 0.0),
        direction=(0.0, 0.0, 1.0),
        i_size=lx,
        j_size=ly,
        i_resolution=16,
        j_resolution=12,
    )
    plotter.add_mesh(
        ground,
        color="#202936",
        opacity=0.80,
        show_edges=True,
        edge_color="#3c4a59",
        line_width=0.5,
        lighting=False,
        name="ground",
    )
    for block in buildings:
        x0 = float(block["x_m"])
        x1 = x0 + float(block["thickness_m"])
        box = pv.Box(bounds=(
            x0,
            x1,
            float(block["y_min_m"]),
            float(block["y_max_m"]),
            0.0,
            float(block["height_m"]),
        ))
        plotter.add_mesh(
            box,
            color="#aeb8c4",
            opacity=1.0,
            show_edges=True,
            edge_color="#111820",
            line_width=1.5,
            ambient=0.30,
            diffuse=0.72,
            specular=0.20,
            specular_power=18.0,
            name=str(block["id"]),
            force_opaque=True,
        )
    source = pv.Sphere(radius=2.2, center=source_position_m)
    plotter.add_mesh(
        source,
        color="#fff3a6",
        emissive=True,
        name="continuous source",
    )
    wind_arrow = pv.Arrow(
        start=(7.0, 10.0, 43.0),
        direction=(1.0, 0.0, 0.0),
        scale=24.0,
        tip_length=0.18,
        tip_radius=0.08,
        shaft_radius=0.025,
    )
    plotter.add_mesh(
        wind_arrow,
        color="white",
        emissive=True,
        name="wind arrow",
    )
    plotter.add_point_labels(
        [(19.0, 10.0, 46.0)],
        [f"wind {wind_speed_m_s:.1f} m/s"],
        text_color="white",
        font_size=12,
        point_size=0,
        shape=None,
        always_visible=True,
        name="wind label",
    )


def render(args: argparse.Namespace) -> dict[str, object]:
    started = time.perf_counter()
    config, barrier = build_continuous_barrier_resolved_config(
        duration_s=120.0,
        release_duration_s=120.0,
        release_rate_kg_s=0.10,
        wind_speed_m_s=2.2,
        barrier_layout="building-cluster",
        store_scalar_fields=True,
    )
    buildings = tuple(barrier["buildings"])
    with np.load(args.snapshot_cache) as cache:
        times_s = np.asarray(cache["time_s"], dtype=float)
        fields_kg_m3 = np.asarray(cache["fields_kg_m3"], dtype=np.float32)
    expected_shape = (len(times_s), config.cells[2], config.cells[1], config.cells[0])
    if fields_kg_m3.shape != expected_shape:
        raise ValueError(
            f"snapshot shape {fields_kg_m3.shape} does not match {expected_shape}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    mp4_path = args.output_dir / "degali_building_cluster_cfd_volume.mp4"
    gif_path = args.output_dir / "degali_building_cluster_cfd_volume.gif"
    poster_path = args.output_dir / "degali_building_cluster_cfd_volume_poster.png"
    metadata_path = args.output_dir / "visualization_metadata.json"

    first_display = _display_field(
        fields_kg_m3[0],
        ambient_air_density_kg_m3=config.ambient_air_density_kg_m3,
        domain_m=config.domain_m,
        buildings=buildings,
        upscale=args.upscale,
        smooth_sigma=args.smooth_sigma,
    )
    grid = _image_grid(first_display, config.domain_m)
    opacity = _opacity_transfer(
        threshold_vol_percent=args.threshold_vol_percent,
        maximum_vol_percent=args.maximum_vol_percent,
    )

    plotter = pv.Plotter(
        off_screen=True,
        window_size=(args.width, args.height),
        lighting="three lights",
    )
    plotter.set_background("#07111d", top="#14283b")
    plotter.enable_depth_peeling(number_of_peels=12, occlusion_ratio=0.0)
    volume_actor = plotter.add_volume(
        grid,
        scalars="H2_vol_percent",
        cmap="turbo",
        clim=(0.0, args.maximum_vol_percent),
        opacity=opacity,
        n_colors=256,
        blending="composite",
        mapper="smart",
        opacity_unit_distance=2.5,
        shade=True,
        ambient=0.22,
        diffuse=0.70,
        specular=0.22,
        specular_power=18.0,
        scalar_bar_args={
            "title": "H2 concentration [vol%]",
            "vertical": True,
            "position_x": 0.86,
            "position_y": 0.19,
            "width": 0.055,
            "height": 0.62,
            "n_labels": 5,
            "fmt": "%.2f",
            "color": "white",
            "title_font_size": 15,
            "label_font_size": 12,
        },
        name="H2 volume",
    )
    _add_scene_geometry(
        plotter,
        domain_m=config.domain_m,
        buildings=buildings,
        source_position_m=config.source_position_m,
        wind_speed_m_s=config.wind_speed_m_s,
    )
    plotter.show_bounds(
        bounds=(0.0, config.domain_m[0], 0.0, config.domain_m[1], 0.0, config.domain_m[2]),
        grid="back",
        location="outer",
        xtitle="x [m]",
        ytitle="y [m]",
        ztitle="z [m]",
        color="#cad4df",
        font_size=10,
        all_edges=True,
    )
    plotter.add_text(
        "DEGALI transient H2 dispersion | volume + centreline contour",
        position="upper_left",
        font_size=14,
        color="white",
        name="title",
    )
    time_actor = plotter.add_text(
        "t = 0.0 s",
        position=(18, 55),
        font_size=18,
        color="#fff3a6",
        name="time",
    )
    display_dx_m = config.domain_m[0] / (config.cells[0] * args.upscale)
    plotter.add_text(
        (
            "Composite ray-cast display | solver grid 5 m | "
            f"display interpolation {display_dx_m:g} m"
        ),
        position="lower_left",
        font_size=10,
        color="#cad4df",
        name="method",
    )
    plotter.camera_position = [
        (245.0, -175.0, 118.0),
        (83.0, 60.0, 14.0),
        (0.0, 0.0, 1.0),
    ]
    plotter.camera.zoom(1.08)

    frames = list(_interpolated_frames(times_s, fields_kg_m3, args.temporal_factor))
    mp4_writer = imageio.get_writer(
        mp4_path,
        fps=args.fps,
        codec="libx264",
        quality=8,
        macro_block_size=1,
    )
    gif_writer = imageio.get_writer(
        gif_path,
        mode="I",
        duration=1000.0 / args.fps,
        loop=0,
        subrectangles=True,
    )
    maximum_display = 0.0
    for frame_index, (time_s, field_kg_m3) in enumerate(frames):
        display = _display_field(
            field_kg_m3,
            ambient_air_density_kg_m3=config.ambient_air_density_kg_m3,
            domain_m=config.domain_m,
            buildings=buildings,
            upscale=args.upscale,
            smooth_sigma=args.smooth_sigma,
        )
        maximum_display = max(maximum_display, float(np.max(display)))
        flat_display = display.ravel(order="C")
        grid.point_data["H2_vol_percent"] = flat_display
        grid.Modified()
        # ``add_volume`` may wrap the input in a mapper-owned image object.
        # Updating that object explicitly guarantees the time-varying scalar
        # field reaches VTK instead of leaving the initial zero frame resident.
        volume_grid = volume_actor.mapper.dataset
        volume_grid.point_data["H2_vol_percent"] = flat_display
        volume_grid.Modified()

        # Rebuild the section from the current point data.  It is deliberately
        # semi-transparent so the volume and the opaque buildings remain legible.
        section = grid.slice(
            normal=(0.0, 1.0, 0.0),
            origin=(0.5 * config.domain_m[0], 60.0, 0.5 * config.domain_m[2]),
        )
        plotter.add_mesh(
            section,
            scalars="H2_vol_percent",
            cmap="turbo",
            clim=(0.0, args.maximum_vol_percent),
            # Use the scalar transfer function here as well.  A constant
            # opacity would draw the whole zero-valued plane as a coloured
            # rectangle and hide the surrounding volume.
            opacity=np.clip(opacity.astype(float) * 1.35, 0.0, 220.0).astype(np.uint8),
            lighting=False,
            show_scalar_bar=False,
            name="centreline contour",
            remove_existing_actor=True,
        )
        time_actor.SetInput(f"t = {time_s:5.1f} s")
        plotter.render()
        image = plotter.screenshot(return_img=True)
        mp4_writer.append_data(image)
        gif_writer.append_data(image)
        if frame_index == len(frames) - 1:
            imageio.imwrite(poster_path, image)
        print(
            f"rendered frame {frame_index + 1}/{len(frames)} at t={time_s:.1f} s",
            flush=True,
        )
    mp4_writer.close()
    gif_writer.close()
    plotter.close()

    metadata = {
        "schema": "degali.cfd-style-volume-visualization.v1",
        "classification": (
            "CFD-style post-processing of a coarse DEGALI resolved field; "
            "not LES/RANS CFD and not external validation"
        ),
        "visual_reference_basis": [
            {
                "software": "ANSYS Fluent",
                "feature": "solution animations of contours, vectors, pathlines, particle tracks and scenes",
                "url": "https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_sec_solve_animate.html",
            },
            {
                "software": "ParaView",
                "feature": "composite volume rendering with colour and opacity transfer functions",
                "url": "https://docs.paraview.org/en/latest/UsersGuide/displayingData.html",
            },
        ],
        "source_snapshot_cache": str(args.snapshot_cache.resolve()),
        "source_snapshot_sha256": _sha256(args.snapshot_cache),
        "scenario": {
            "domain_m": list(config.domain_m),
            "solver_cells_xyz": list(config.cells),
            "solver_cell_size_m": [
                config.domain_m[i] / config.cells[i] for i in range(3)
            ],
            "release_rate_kg_s": 0.10,
            "release_duration_s": 120.0,
            "wind_speed_m_s": config.wind_speed_m_s,
            "source_position_m": list(config.source_position_m),
            "buildings": list(buildings),
        },
        "render": {
            "method": "VTK composite ray-cast volume plus vertical slice",
            "display_only_spatial_interpolation": True,
            "display_only_temporal_interpolation": True,
            "upscale": args.upscale,
            "display_grid_shape_zyx": list(first_display.shape),
            "smooth_sigma_display_cells": args.smooth_sigma,
            "threshold_vol_percent": args.threshold_vol_percent,
            "fixed_colour_range_vol_percent": [0.0, args.maximum_vol_percent],
            "frame_count": len(frames),
            "fps": args.fps,
            "video_duration_s": len(frames) / args.fps,
            "maximum_display_concentration_vol_percent": maximum_display,
            "window_size_px": [args.width, args.height],
            "runtime_s": time.perf_counter() - started,
        },
        "outputs": {
            "mp4": str(mp4_path.resolve()),
            "gif": str(gif_path.resolve()),
            "poster": str(poster_path.resolve()),
            "metadata": str(metadata_path.resolve()),
        },
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--width", type=int, default=1280)
    parser.add_argument("--height", type=int, default=720)
    parser.add_argument("--fps", type=int, default=6)
    parser.add_argument("--temporal-factor", type=int, default=2)
    parser.add_argument("--upscale", type=int, default=2)
    parser.add_argument("--smooth-sigma", type=float, default=0.65)
    parser.add_argument("--threshold-vol-percent", type=float, default=0.02)
    parser.add_argument("--maximum-vol-percent", type=float, default=0.40)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metadata = render(args)
    print(json.dumps(metadata["outputs"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
