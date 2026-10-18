"""Plot the completed TI021 obstacle-resolved maximum concentration field."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROPERTY_AUDIT = (
    ROOT / "outputs" / "thorney-ti021-r12-property-audit-2026-10-11-v1.json"
)

try:
    from tools.audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )
except ModuleNotFoundError:  # Direct ``python tools/<script>.py`` execution.
    from audit_thorney_ti021_handoff_streaming_preflight import (
        build_configuration,
    )

from degali.addons.site_geometry import CircularArcWall


def _projection_max_ignoring_nonfinite(
    field: np.ndarray, *, axis: int,
) -> np.ndarray:
    """Project a 3-D field without letting masked solid cells poison maxima."""
    values = np.asarray(field, dtype=float)
    if values.ndim != 3:
        raise ValueError("field must have [z, y, x] dimensions")
    finite = np.isfinite(values)
    maximum = np.max(np.where(finite, values, -np.inf), axis=axis)
    return np.where(np.any(finite, axis=axis), maximum, np.nan)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--field", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--field-view", choices=("endpoint", "time-maximum"),
        default="time-maximum",
    )
    parser.add_argument("--species-label", default="tracked species")
    parser.add_argument("--contour-label", default="reference contour")
    parser.add_argument(
        "--field-frame", choices=("inner-obstacle", "outer-posthandoff"),
        default="inner-obstacle",
    )
    parser.add_argument("--property-audit", type=Path, default=DEFAULT_PROPERTY_AUDIT)
    args = parser.parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle

    with np.load(args.field) as data:
        x = data["x_m"]
        y = data["y_m"]
        z = data["z_m"]
        crosswind_q = (
            data["crosswind_q_m"]
            if args.field_frame == "outer-posthandoff" and "crosswind_q_m" in data.files
            else y
        )
        field_key = (
            "concentration_vol_percent"
            if args.field_view == "endpoint"
            else "maximum_concentration_vol_percent"
        )
        field = data[field_key]
        solid = data["solid_mask"]

    arc = None
    if args.field_frame == "inner-obstacle":
        properties = json.loads(args.property_audit.read_text(encoding="utf-8"))
        spacing_m = float(np.median(np.diff(x)))
        config, _, _ = build_configuration(
            properties,
            cell_size_m=spacing_m,
            time_step_s=0.02,
            active_steps=1,
            corrector_tolerance=1.0e-4,
            circular_arc_representation="zero_thickness_surface",
        )
        arc = next(
            (item for item in config.obstacles if isinstance(item, CircularArcWall)),
            None,
        )

    # The saved field is indexed [z, y, x].  Use physically meaningful slices.
    # Solid obstacle cells are intentionally NaN in the stored fields.  Ignore
    # them while reducing, preserving all-solid projected cells as NaN.
    top = _projection_max_ignoring_nonfinite(field, axis=0)
    front = _projection_max_ignoring_nonfinite(field, axis=1)
    side = _projection_max_ignoring_nonfinite(field, axis=2)
    top_solid = np.any(solid, axis=0)
    front_solid = np.any(solid, axis=1)
    side_solid = np.any(solid, axis=2)
    levels = np.asarray([0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32])
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_under("#f8fafc")
    norm = LogNorm(vmin=0.1, vmax=32)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    if args.field_frame == "outer-posthandoff":
        panels = (
            (top, top_solid, x, crosswind_q, "Downwind/crosswind view (max over z)", "downwind x (m)", "crosswind q (m)"),
            (front, front_solid, x, z, "Downwind/vertical envelope (max over q)", "downwind x (m)", "z (m)"),
            (side, side_solid, crosswind_q, z, "Crosswind/vertical envelope (max over x)", "crosswind q (m)", "z (m)"),
        )
    else:
        panels = ((top, top_solid, x, y, "Plan view (max over z)", "x (m)", "y (m)"),
                  (front, front_solid, x, z, "Front view (max over y)", "x (m)", "z (m)"),
                  (side, side_solid, y, z, "Side view (max over x)", "y (m)", "z (m)"))
    image = None
    for axis, (values, mask, xx, yy, title, xlabel, ylabel) in zip(axes, panels):
        image = axis.contourf(xx, yy, values, levels=levels, cmap=cmap, norm=norm, extend="max")
        if float(np.nanmax(values)) >= 4.0:
            axis.contour(xx, yy, values, levels=[4.0], colors="white", linewidths=3.0)
            axis.contour(xx, yy, values, levels=[4.0], colors="#111827", linewidths=1.1)
        # Hatch the obstacle footprint so high concentration is not confused with solid cells.
        if np.any(mask):
            axis.contourf(xx, yy, mask.astype(float), levels=[0.5, 1.5], colors="none", hatches=["///"], alpha=0.0)
            axis.contour(xx, yy, mask.astype(float), levels=[0.5], colors="#111827", linewidths=1.2)
        if arc is not None:
            theta = np.deg2rad(np.linspace(
                arc.arc_start_deg,
                arc.arc_start_deg + arc.angular_span_deg,
                361,
            ))
            arc_x = arc.center_m[0] + arc.radius_m * np.cos(theta)
            arc_y = arc.center_m[1] + arc.radius_m * np.sin(theta)
            if title.startswith("Plan"):
                axis.plot(arc_x, arc_y, color="#111827", linewidth=2.0,
                          linestyle=(0, (5, 2)), label="idealized arc surface")
            elif title.startswith("Front"):
                axis.add_patch(Rectangle(
                    (float(np.min(arc_x)), arc.z_min_m),
                    float(np.ptp(arc_x)), arc.z_max_m - arc.z_min_m,
                    fill=False, edgecolor="#111827", linewidth=1.6,
                    linestyle=(0, (5, 2)),
                ))
            else:
                axis.add_patch(Rectangle(
                    (float(np.min(arc_y)), arc.z_min_m),
                    float(np.ptp(arc_y)), arc.z_max_m - arc.z_min_m,
                    fill=False, edgecolor="#111827", linewidth=1.6,
                    linestyle=(0, (5, 2)),
                ))
        axis.set_title(title, weight="bold")
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.15)
    field_label = (
        f"Endpoint {args.species_label} concentration (vol%)"
        if args.field_view == "endpoint"
        else f"Time-maximum {args.species_label} concentration (vol%)"
    )
    fig.suptitle(
        f"DEGALI TI021 {args.field_frame.replace('-', ' ')} transient field\n{field_label}",
        weight="bold", fontsize=15,
    )
    cbar = fig.colorbar(image, ax=axes, shrink=0.86, pad=0.02)
    cbar.set_label(f"{args.species_label} concentration (vol%)")
    overlay_note = (
        "dashed = idealized zero-thickness arc and its projections"
        if arc is not None else
        "post-handoff outer field only; obstacle is outside this domain"
        if args.field_frame == "outer-posthandoff" else
        "obstacle footprint shown by hatching"
    )
    fig.text(
        0.5, 0.01,
        f"White/black contour = 4 vol% {args.contour_label}; {overlay_note}",
        ha="center", fontsize=9,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    fig.savefig(args.output.with_suffix(".svg"))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
