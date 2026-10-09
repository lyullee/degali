"""Render a compact convergence/conservation figure from the 3-D audit JSON.

The figure is a numerical-verification diagnostic only.  It deliberately does
not label the error as experimental accuracy or sensor uncertainty.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _load(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "degali.transient-dense-gas-3d-error-study.v1":
        raise ValueError("unexpected transient error-study schema")
    return payload


def plot(input_path: Path, output_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    payload = _load(input_path)
    names = ["coarse", "medium", "fine"]
    levels = payload["levels"]
    cells = [levels[name]["cells"] for name in names]
    grid_size = [int(nx) * int(ny) * int(nz) for nx, ny, nz in cells]
    rmse = [float(levels[name]["aggregate_rmse_kg_m3"]) for name in names]
    maximum = [float(levels[name]["aggregate_max_absolute_error_kg_m3"]) for name in names]
    residual = [float(levels[name]["maximum_mass_residual_kg"]) for name in names]
    order = float(payload["apparent_order_from_aggregate_rmse"]["medium_to_fine"])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)
    ax = axes[0]
    ax.loglog(grid_size, rmse, "o-", label="aggregate RMSE")
    ax.loglog(grid_size, maximum, "s--", label="aggregate max abs")
    for x, y, label in zip(grid_size, rmse, names):
        ax.annotate(label, (x, y), xytext=(5, 5), textcoords="offset points", fontsize=9)
    ax.set_xlabel("Grid cells (Nx × Ny × Nz)")
    ax.set_ylabel("Error vs. high-resolution reference (kg H₂ m⁻³)")
    ax.set_title("Resolution error")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    ax.text(
        0.03,
        0.04,
        f"medium→fine apparent order = {order:.2f}",
        transform=ax.transAxes,
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.9},
    )

    ax = axes[1]
    ax.loglog(grid_size, residual, "o-", color="#7f1d1d")
    for x, y, label in zip(grid_size, residual, names):
        ax.annotate(label, (x, y), xytext=(5, 5), textcoords="offset points", fontsize=9)
    ax.set_xlabel("Grid cells (Nx × Ny × Nz)")
    ax.set_ylabel("Maximum mass residual (kg)")
    ax.set_title("Mass conservation residual")
    ax.grid(True, which="both", alpha=0.25)
    ax.text(
        0.03,
        0.04,
        "numerical conservation check",
        transform=ax.transAxes,
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.9},
    )

    fig.suptitle("Transient 3-D solver: numerical verification", fontsize=14, fontweight="bold")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="svg", metadata={"Title": "Transient 3-D numerical verification"})
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="error_study.json from the benchmark")
    parser.add_argument("--output", type=Path, required=True, help="SVG output path")
    args = parser.parse_args(argv)
    plot(args.input, args.output)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
