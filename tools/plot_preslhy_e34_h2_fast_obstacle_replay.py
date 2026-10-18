"""Render a conditional PRESLHY H2 obstacle replay field."""

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
    x, y, z = data["x_m"], data["y_m"], data["z_m"]
    field = np.asarray(data["maximum_h2_vol_percent"], dtype=float)
    mid = int(np.argmin(np.abs(y - 0.8)))
    section = field[:, mid, :]
    levels = np.asarray([0.0, 1.0, 4.0, 10.0, 25.0, 50.0, 100.0])
    norm = BoundaryNorm(levels, ncolors=256, clip=False)
    fig, ax = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    image = ax.contourf(x, z, section, levels=levels, norm=norm, cmap="magma", extend="max")
    if np.nanmax(section) >= 4.0:
        ax.contour(x, z, section, levels=[4.0], colors="cyan", linewidths=2.0)
        ax.text(1.01, 0.82, "4 vol%", color="cyan", weight="bold")
    ax.add_patch(Rectangle((1.0, 0.0), 0.4, 0.5, facecolor="#263238", edgecolor="white", hatch="////", zorder=5))
    ax.set(xlabel="Downwind x (m)", ylabel="Height z (m)", title="PRESLHY E3.4 conditional H2 replay — maximum vol% at y=0.8 m")
    ax.set_xlim(float(x[0]), float(x[-1])); ax.set_ylim(float(z[0]), float(z[-1]))
    fig.colorbar(image, ax=ax, label="H2 concentration (vol%)")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--field", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.field, args.output)


if __name__ == "__main__":
    main()
