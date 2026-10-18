"""Plot the resolved FFI maximum H2 field and 4 vol% contour."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lateral-m", type=float, default=0.45)
    args = parser.parse_args()
    with np.load(args.input) as data:
        x = np.asarray(data["x_m"], dtype=float)
        y = np.asarray(data["y_m"], dtype=float)
        z = np.asarray(data["z_m"], dtype=float)
        field = np.asarray(data["maximum_concentration_vol_percent"], dtype=float)
        solid = np.asarray(data["solid_mask"], dtype=bool)
    lateral_index = int(np.argmin(np.abs(y - args.lateral_m)))
    values = field[:, lateral_index, :].copy()
    obstacle = solid[:, lateral_index, :]
    values[obstacle] = np.nan
    xx, zz = np.meshgrid(x, z)
    fig, axis = plt.subplots(figsize=(11, 4.8))
    image = axis.pcolormesh(xx, zz, values, shading="nearest", cmap="magma", vmin=0.0)
    finite = np.isfinite(values)
    if np.any(finite & (values >= 4.0)) and np.nanmin(values) < 4.0 < np.nanmax(values):
        axis.contour(xx, zz, values, levels=[4.0], colors="cyan", linewidths=1.6)
    axis.contourf(xx, zz, obstacle, levels=[0.5, 1.5], colors=["#4b5563"], alpha=0.9)
    axis.set_xlabel("streamwise x (m)")
    axis.set_ylabel("height z (m)")
    axis.set_title(f"FFI O1 resolved buoyant maximum H₂ field at y={y[lateral_index]:.3f} m")
    axis.grid(alpha=0.2)
    fig.colorbar(image, ax=axis, label="H₂ (vol%)")
    fig.text(
        0.5, 0.01,
        "Conditional sensitivity geometry; cyan contour = 4 vol%; not a facility validation score.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    print({"output": str(args.output.resolve()), "lateral_index": lateral_index})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
