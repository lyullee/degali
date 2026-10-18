"""Render the maximum H2 vol% field from the fast obstacle screen."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def render(field_path: Path, output_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from matplotlib.patches import Rectangle

    data = np.load(field_path)
    x = data["x_m"]
    y = data["y_m"]
    z = data["z_m"]
    field = np.asarray(data["maximum_h2_vol_percent"], dtype=float)
    mid = int(np.argmin(np.abs(y - 20.0)))
    slice_xz = field[:, mid, :]
    levels = np.asarray([0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 10.0])
    norm = BoundaryNorm(levels, ncolors=256, clip=False)
    figure, axis = plt.subplots(figsize=(10.5, 4.8), constrained_layout=True)
    image = axis.contourf(x, z, slice_xz, levels=levels, norm=norm, cmap="magma", extend="max")
    if np.nanmax(slice_xz) >= 4.0:
        axis.contour(x, z, slice_xz, levels=[4.0], colors="cyan", linewidths=2.0)
        axis.text(31.0, 11.0, "4 vol%", color="cyan", weight="bold")
    axis.add_patch(Rectangle((30.0, 0.0), 10.0, 10.0, facecolor="#263238", edgecolor="white", hatch="////", zorder=5))
    axis.set(xlabel="Downwind x (m)", ylabel="Height z (m)", title="DEGALI fast H2 obstacle screen — maximum vol% at y=20 m")
    axis.set_xlim(float(x[0]), float(x[-1])); axis.set_ylim(float(z[0]), float(z[-1]))
    figure.colorbar(image, ax=axis, label="H2 concentration (vol%)")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--field", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.field, args.output)


if __name__ == "__main__":
    main()
