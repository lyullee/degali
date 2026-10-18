"""Render an optional smooth 3-D visual reconstruction of a TI021 replay.

The input is the receptor time-series NPZ written by
``score_ti021_dense_layer_hybrid.py``.  The renderer performs anisotropic
inverse-distance interpolation between predicted receptor values, optional
Gaussian display smoothing, and temporal interpolation.  It never feeds the
interpolated values back into DEGALI and therefore does not change any score,
mass ledger, or solver result.

This is presentation-only output.  It is not a resolved 3-D flow solution,
CFD, or validation evidence.  The classification is embedded in every frame
and in the accompanying JSON metadata.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import NamedTuple

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "outputs"
    / "ti021-sequential-strict-f020-c1p7s0p5-independent-600s-2026-10-17.npz"
)
DEFAULT_WORKBOOK = ROOT / "ti21_input_ascii.xls"
DEFAULT_OUTPUT = ROOT / "outputs" / "ti021-interpolated-3d-animation"


class InterpolationPlan(NamedTuple):
    indices: np.ndarray
    weights: np.ndarray
    support: np.ndarray


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _frame_times(
    source_time_s: np.ndarray,
    *,
    frames: int,
    end_time_s: float | None,
    time_warp: str,
) -> np.ndarray:
    if frames < 2:
        raise ValueError("frames must be at least 2")
    start = float(source_time_s[0])
    end = float(source_time_s[-1] if end_time_s is None else end_time_s)
    if not start < end <= float(source_time_s[-1]):
        raise ValueError("end_time_s must lie inside the input time range")
    unit = np.linspace(0.0, 1.0, frames)
    if time_warp == "linear":
        warped = unit
    elif time_warp == "early_detail":
        # Half the rendered frames cover the first quarter of physical time.
        warped = unit * unit
    else:
        raise ValueError("time_warp must be 'linear' or 'early_detail'")
    return start + (end - start) * warped


def _interpolate_time(
    source_time_s: np.ndarray,
    traces_sensor_time: np.ndarray,
    frame_time_s: np.ndarray,
) -> np.ndarray:
    source_time_s = np.asarray(source_time_s, dtype=float)
    traces = np.asarray(traces_sensor_time, dtype=float)
    if source_time_s.ndim != 1 or traces.ndim != 2:
        raise ValueError("time must be 1-D and concentration traces must be 2-D")
    if traces.shape[1] != source_time_s.size:
        raise ValueError("trace time dimension does not match time_s")
    if np.any(np.diff(source_time_s) <= 0.0):
        raise ValueError("time_s must be strictly increasing")
    return np.stack(
        [np.interp(frame_time_s, source_time_s, np.maximum(trace, 0.0)) for trace in traces],
        axis=1,
    )


def _inverse_distance_plan(
    sensor_xyz_m: np.ndarray,
    target_xyz_m: np.ndarray,
    *,
    neighbours: int,
    power: float,
    vertical_scale: float,
    support_radius_m: float,
) -> InterpolationPlan:
    sensors = np.asarray(sensor_xyz_m, dtype=float)
    targets = np.asarray(target_xyz_m, dtype=float)
    if sensors.ndim != 2 or targets.ndim != 2 or sensors.shape[1] != 3 or targets.shape[1] != 3:
        raise ValueError("sensor and target coordinates must have shape (n, 3)")
    if neighbours < 1 or power <= 0.0 or vertical_scale <= 0.0 or support_radius_m <= 0.0:
        raise ValueError("interpolation controls must be positive")
    scale = np.asarray((1.0, 1.0, vertical_scale), dtype=float)
    tree = cKDTree(sensors * scale)
    k = min(int(neighbours), sensors.shape[0])
    distance, indices = tree.query(targets * scale, k=k)
    if k == 1:
        distance = distance[:, None]
        indices = indices[:, None]
    exact = distance[:, 0] <= 1.0e-12
    safe = np.maximum(distance, 1.0e-12)
    weights = safe ** (-float(power))
    weights /= weights.sum(axis=1, keepdims=True)
    if np.any(exact):
        weights[exact] = 0.0
        weights[exact, 0] = 1.0
    nearest_horizontal = np.linalg.norm(
        targets[:, None, :2] - sensors[indices[:, :1], :2], axis=2
    )[:, 0]
    # Smoothly fade the reconstruction outside the measured/predicted receptor
    # footprint instead of inventing a constant far-field extrapolation.
    ratio = np.clip(nearest_horizontal / float(support_radius_m), 0.0, 1.0)
    support = np.cos(0.5 * np.pi * ratio) ** 2
    support[nearest_horizontal >= support_radius_m] = 0.0
    return InterpolationPlan(indices=indices, weights=weights, support=support)


def _apply_plan(plan: InterpolationPlan, sensor_values: np.ndarray) -> np.ndarray:
    values = np.asarray(sensor_values, dtype=float)
    reconstructed = np.sum(values[plan.indices] * plan.weights, axis=1)
    return np.maximum(reconstructed * plan.support, 0.0)


def _sensor_positions(workbook: Path) -> np.ndarray:
    dataset = read_thorney_ti021_workbook(workbook)
    return np.asarray(
        [
            (60.0 + sensor.local_y_m, 320.0 - sensor.local_x_m, sensor.height_m)
            for sensor in dataset.dose_sensors
        ],
        dtype=float,
    )


def _grid(
    sensors: np.ndarray,
    *,
    nx: int,
    ny: int,
    nz: int,
    z_max_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if min(nx, ny, nz) < 2 or z_max_m <= 0.0:
        raise ValueError("grid dimensions must be >=2 and z_max_m must be positive")
    x = np.linspace(40.0, max(580.0, float(sensors[:, 0].max()) + 20.0), nx)
    y = np.linspace(
        min(0.0, float(sensors[:, 1].min()) - 20.0),
        max(640.0, float(sensors[:, 1].max()) + 20.0),
        ny,
    )
    z = np.linspace(0.0, z_max_m, nz)
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    return x, y, z, points


def reconstruct_fields(
    sensor_frame_values: np.ndarray,
    plan: InterpolationPlan,
    *,
    shape: tuple[int, int, int],
    smooth_sigma_cells: float,
) -> np.ndarray:
    values = np.asarray(sensor_frame_values, dtype=float)
    fields = np.empty((values.shape[0], *shape), dtype=np.float32)
    for frame, sensor_values in enumerate(values):
        field = _apply_plan(plan, sensor_values).reshape(shape)
        if smooth_sigma_cells > 0.0:
            field = gaussian_filter(field, smooth_sigma_cells, mode="nearest")
        fields[frame] = np.clip(field, 0.0, 100.0).astype(np.float32)
    return fields


def _draw_fence(axis, *, vertical_exaggeration: float) -> None:
    theta = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 72)
    height = np.linspace(0.0, 5.0, 6)
    tt, zz = np.meshgrid(theta, height)
    xx = 60.0 + 50.0 * np.cos(tt)
    yy = 320.0 + 50.0 * np.sin(tt)
    axis.plot_surface(
        xx,
        yy,
        zz * vertical_exaggeration,
        color="#4b5563",
        alpha=0.72,
        linewidth=0.0,
        shade=True,
    )


def _render(
    *,
    fields: np.ndarray,
    frame_time_s: np.ndarray,
    x_m: np.ndarray,
    y_m: np.ndarray,
    z_m: np.ndarray,
    output_dir: Path,
    fps: int,
    display_threshold_vol_percent: float,
    vertical_exaggeration: float,
    max_points: int,
    dpi: int,
) -> tuple[Path, Path, int]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation, colors

    output_dir.mkdir(parents=True, exist_ok=True)
    xx, yy, zz = np.meshgrid(x_m, y_m, z_m, indexing="xy")
    # meshgrid above is (ny, nx, nz); field storage is (nz, ny, nx).
    px = np.transpose(xx, (2, 0, 1)).ravel()
    py = np.transpose(yy, (2, 0, 1)).ravel()
    pz = np.transpose(zz, (2, 0, 1)).ravel()

    positive = fields[fields >= display_threshold_vol_percent]
    vmax = max(4.0, float(np.percentile(positive, 99.5)) if positive.size else 4.0)
    norm = colors.LogNorm(vmin=display_threshold_vol_percent, vmax=vmax)
    cmap = plt.get_cmap("turbo")

    fig = plt.figure(figsize=(11.2, 7.0), facecolor="#07111f")
    axis = fig.add_subplot(111, projection="3d")
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(scalar, ax=axis, fraction=0.025, pad=0.06)
    colorbar.set_label("Reconstructed released-mixture concentration (vol%)", color="white")
    colorbar.ax.tick_params(colors="white")
    colorbar.outline.set_edgecolor("#94a3b8")

    def draw(index: int) -> list[object]:
        axis.clear()
        field = fields[index]
        flat = field.ravel()
        selected = np.flatnonzero(flat >= display_threshold_vol_percent)
        if selected.size > max_points:
            stride = int(math.ceil(selected.size / max_points))
            selected = selected[::stride]
        if selected.size:
            concentration = flat[selected]
            rgba = cmap(norm(concentration))
            rgba[:, 3] = np.clip(0.12 + 0.65 * norm(concentration), 0.12, 0.77)
            axis.scatter(
                px[selected],
                py[selected],
                pz[selected] * vertical_exaggeration,
                c=rgba,
                s=7.0,
                depthshade=False,
                linewidths=0.0,
            )
        _draw_fence(axis, vertical_exaggeration=vertical_exaggeration)
        axis.scatter(
            [60.0], [320.0], [0.4 * vertical_exaggeration],
            marker="*", s=100, color="#ffffff", edgecolor="#111827", linewidth=0.7,
        )
        axis.set_xlim(float(x_m[0]), float(x_m[-1]))
        axis.set_ylim(float(y_m[0]), float(y_m[-1]))
        axis.set_zlim(0.0, float(z_m[-1] * vertical_exaggeration))
        axis.set_xlabel("Site x (m)", color="white", labelpad=8)
        axis.set_ylabel("Site y (m)", color="white", labelpad=8)
        axis.set_zlabel("Height z (m)", color="white", labelpad=8)
        z_ticks = np.linspace(0.0, float(z_m[-1]), 4)
        axis.set_zticks(z_ticks * vertical_exaggeration)
        axis.set_zticklabels([f"{value:g}" for value in z_ticks])
        axis.tick_params(colors="#cbd5e1")
        axis.view_init(elev=26.0, azim=-62.0)
        axis.set_box_aspect((1.25, 1.0, 0.48))
        axis.set_facecolor("#07111f")
        for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
            pane.set_facecolor((0.035, 0.075, 0.13, 1.0))
            pane.set_edgecolor("#334155")
        axis.grid(False)
        axis.set_title(
            "DEGALI TI021 — 3D visual reconstruction",
            color="white",
            pad=16,
            fontsize=15,
            weight="bold",
        )
        axis.text2D(
            0.02,
            0.94,
            f"physical time {frame_time_s[index]:6.1f} s   |   max {float(field.max()):6.2f} vol%",
            transform=axis.transAxes,
            color="#e2e8f0",
            fontsize=10,
        )
        axis.text2D(
            0.02,
            0.02,
            "VISUAL RECONSTRUCTION — interpolated receptor field, not resolved 3D CFD",
            transform=axis.transAxes,
            color="#fbbf24",
            fontsize=9,
            weight="bold",
        )
        axis.text2D(
            0.02,
            0.055,
            f"vertical exaggeration {vertical_exaggeration:g}×   |   display threshold {display_threshold_vol_percent:g} vol%",
            transform=axis.transAxes,
            color="#94a3b8",
            fontsize=8.5,
        )
        return []

    animation_object = animation.FuncAnimation(
        fig, draw, frames=fields.shape[0], interval=1000.0 / fps, blit=False
    )
    gif_path = output_dir / "degali_ti021_interpolated_3d.gif"
    animation_object.save(gif_path, writer=animation.PillowWriter(fps=fps), dpi=dpi)

    # Poster selects the frame with the largest displayed plume volume.
    volume_count = np.count_nonzero(fields >= display_threshold_vol_percent, axis=(1, 2, 3))
    poster_index = int(np.argmax(volume_count))
    draw(poster_index)
    poster_path = output_dir / "degali_ti021_interpolated_3d_poster.png"
    fig.savefig(poster_path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return gif_path, poster_path, poster_index


def build(args: argparse.Namespace) -> dict[str, object]:
    input_path = args.input_npz.resolve()
    workbook = args.workbook.resolve()
    with np.load(input_path) as archive:
        time_s = np.asarray(archive["time_s"], dtype=float)
        traces = np.asarray(archive["concentration_vol_percent"], dtype=float)
    sensors = _sensor_positions(workbook)
    if traces.shape[0] != sensors.shape[0]:
        raise ValueError("input trace count does not match TI021 sensor geometry")

    frame_time_s = _frame_times(
        time_s,
        frames=args.frames,
        end_time_s=args.end_time_s,
        time_warp=args.time_warp,
    )
    frame_values = _interpolate_time(time_s, traces, frame_time_s)
    x_m, y_m, z_m, target = _grid(
        sensors, nx=args.nx, ny=args.ny, nz=args.nz, z_max_m=args.z_max_m
    )
    plan = _inverse_distance_plan(
        sensors,
        target,
        neighbours=args.neighbours,
        power=args.idw_power,
        vertical_scale=args.vertical_distance_scale,
        support_radius_m=args.support_radius_m,
    )
    fields = reconstruct_fields(
        frame_values,
        plan,
        shape=(args.nz, args.ny, args.nx),
        smooth_sigma_cells=args.smooth_sigma_cells,
    )
    gif, poster, poster_index = _render(
        fields=fields,
        frame_time_s=frame_time_s,
        x_m=x_m,
        y_m=y_m,
        z_m=z_m,
        output_dir=args.output_dir.resolve(),
        fps=args.fps,
        display_threshold_vol_percent=args.display_threshold_vol_percent,
        vertical_exaggeration=args.vertical_exaggeration,
        max_points=args.max_points,
        dpi=args.dpi,
    )
    metadata: dict[str, object] = {
        "schema": "degali.ti021-interpolated-3d-visualization.v1",
        "classification": "presentation_only_interpolated_visual_reconstruction",
        "resolved_3d_physics": False,
        "cfd_claim_allowed": False,
        "validation_evidence_allowed": False,
        "interpolation_changes_solver_outputs": False,
        "source": {
            "npz": str(input_path),
            "npz_sha256": _sha256(input_path),
            "workbook": str(workbook),
            "workbook_sha256": _sha256(workbook),
            "sensor_count": int(sensors.shape[0]),
        },
        "render": {
            "frames": int(args.frames),
            "fps": int(args.fps),
            "physical_time_range_s": [float(frame_time_s[0]), float(frame_time_s[-1])],
            "time_warp": args.time_warp,
            "grid_cells": [int(args.nx), int(args.ny), int(args.nz)],
            "neighbours": int(args.neighbours),
            "idw_power": float(args.idw_power),
            "vertical_distance_scale": float(args.vertical_distance_scale),
            "support_radius_m": float(args.support_radius_m),
            "smooth_sigma_cells": float(args.smooth_sigma_cells),
            "display_threshold_vol_percent": float(args.display_threshold_vol_percent),
            "vertical_exaggeration": float(args.vertical_exaggeration),
            "poster_time_s": float(frame_time_s[poster_index]),
        },
        "outputs": {"animation_gif": str(gif), "poster_png": str(poster)},
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
    parser.add_argument("--input-npz", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--frames", type=int, default=72)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--end-time-s", type=float, default=600.0)
    parser.add_argument("--time-warp", choices=("linear", "early_detail"), default="early_detail")
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--ny", type=int, default=72)
    parser.add_argument("--nz", type=int, default=18)
    parser.add_argument("--z-max-m", type=float, default=9.0)
    parser.add_argument("--neighbours", type=int, default=12)
    parser.add_argument("--idw-power", type=float, default=2.0)
    parser.add_argument("--vertical-distance-scale", type=float, default=24.0)
    parser.add_argument("--support-radius-m", type=float, default=130.0)
    parser.add_argument("--smooth-sigma-cells", type=float, default=0.65)
    parser.add_argument("--display-threshold-vol-percent", type=float, default=0.08)
    parser.add_argument("--vertical-exaggeration", type=float, default=15.0)
    parser.add_argument("--max-points", type=int, default=9000)
    parser.add_argument("--dpi", type=int, default=105)
    args = parser.parse_args()
    if args.fps < 1 or args.max_points < 100 or args.dpi < 60:
        parser.error("fps, max-points, and dpi are outside their supported ranges")
    metadata = build(args)
    print(json.dumps(metadata, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
