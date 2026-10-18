"""Verify the TI021 resolved-inner to outer-plane conservative remap operator."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.addons.conservative_nested_handoff_3d import (
    remap_conservative_nested_handoff_3d,
)


ROOT = Path(__file__).resolve().parents[1]


def _integral(field: np.ndarray, z_edges: np.ndarray, y_edges: np.ndarray) -> np.ndarray:
    area = np.diff(z_edges)[:, None] * np.diff(y_edges)[None, :]
    return np.sum(field * area, axis=(-2, -1))


def build_audit() -> tuple[dict[str, object], dict[str, np.ndarray]]:
    source_dx = 1.25
    target_dx = 5.0
    source_y = np.arange(-100.0, 100.0 + source_dx, source_dx)
    source_z = np.arange(0.0, 40.0 + source_dx, source_dx)
    target_y = np.arange(-320.0, 320.0 + target_dx, target_dx)
    target_z = np.arange(0.0, 40.0 + target_dx, target_dx)
    time_s = np.asarray((0.0, 60.0, 120.0, 240.0, 600.0))
    yc = 0.5 * (source_y[:-1] + source_y[1:])
    zc = 0.5 * (source_z[:-1] + source_z[1:])
    spatial = np.exp(
        -(yc[None, :] / 28.0) ** 2
        - ((zc[:, None] - 2.5) / 4.0) ** 2
    )
    amplitude = np.asarray((0.0, 0.08, 0.12, 0.04, 0.0))[:, None, None]
    species = amplitude * spatial[None, :, :]
    mass_anomaly = 0.65 * species
    enthalpy = -8.0e5 * species
    lateral = (
        np.asarray((0.0, 1.0, -0.6, 0.2, 0.0))[:, None, None]
        * 0.15 * spatial[None, :, :]
    )
    vertical = (
        np.asarray((0.0, 0.4, 0.8, -0.2, 0.0))[:, None, None]
        * 0.08 * spatial[None, :, :]
    )
    result = remap_conservative_nested_handoff_3d(
        time_s=time_s,
        source_z_edges_m=source_z, source_y_edges_m=source_y,
        target_z_edges_m=target_z, target_y_edges_m=target_y,
        total_mass_flux_anomaly_kg_m2_s=mass_anomaly,
        species_mass_flux_kg_m2_s=species,
        sensible_enthalpy_flux_w_m2=enthalpy,
        lateral_momentum_flux_anomaly_pa=lateral,
        vertical_momentum_flux_anomaly_pa=vertical,
        source_id="TI021 nested-plane analytic concentration-blind pulse",
    )
    source_fields = {
        "total_mass_flux_anomaly": mass_anomaly,
        "species_mass_flux": species,
        "sensible_enthalpy_flux": enthalpy,
        "lateral_momentum_flux_anomaly": lateral,
        "vertical_momentum_flux_anomaly": vertical,
    }
    target_fields = {
        "total_mass_flux_anomaly": result.total_mass_flux_anomaly_kg_m2_s,
        "species_mass_flux": result.species_mass_flux_kg_m2_s,
        "sensible_enthalpy_flux": result.sensible_enthalpy_flux_w_m2,
        "lateral_momentum_flux_anomaly": result.lateral_momentum_flux_anomaly_pa,
        "vertical_momentum_flux_anomaly": result.vertical_momentum_flux_anomaly_pa,
    }
    integrals = {
        name: {
            "source": _integral(field, source_z, source_y).tolist(),
            "target": _integral(target_fields[name], target_z, target_y).tolist(),
        }
        for name, field in source_fields.items()
    }
    residuals = result.diagnostics["relative_integral_residuals"]
    gates = {
        "species_mass_flux_relative_residual_at_most_1e_12": (
            residuals["species_mass_flux"] <= 1.0e-12
        ),
        "total_mass_anomaly_relative_residual_at_most_1e_12": (
            residuals["total_mass_flux_anomaly"] <= 1.0e-12
        ),
        "enthalpy_flux_relative_residual_at_most_1e_12": (
            residuals["sensible_enthalpy_flux"] <= 1.0e-12
        ),
        "lateral_momentum_flux_relative_residual_at_most_1e_12": (
            residuals["lateral_momentum_flux_anomaly"] <= 1.0e-12
        ),
        "vertical_momentum_flux_relative_residual_at_most_1e_12": (
            residuals["vertical_momentum_flux_anomaly"] <= 1.0e-12
        ),
        "no_clipping_or_posthoc_renormalisation": not result.diagnostics[
            "clipping_or_posthoc_renormalisation"
        ],
        "outer_plane_covers_inner_plane": bool(
            target_y[0] <= source_y[0] and target_y[-1] >= source_y[-1]
            and target_z[0] <= source_z[0] and target_z[-1] >= source_z[-1]
        ),
    }
    record = {
        "schema": "degali.thorney-ti021-nested-handoff-audit.v1",
        "status": (
            "accepted_conservative_handoff_remap_operator_only"
            if all(gates.values()) else "withheld_conservative_handoff_remap"
        ),
        "scope": (
            "analytic y-z downstream-plane flux remap only; no TI021 dose or "
            "concentration is read, and outer transient boundary coupling is pending"
        ),
        "configuration": {
            "source_plane_crosswind_bounds_m": [source_y[0], source_y[-1]],
            "source_plane_height_bounds_m": [source_z[0], source_z[-1]],
            "source_cell_size_m": source_dx,
            "target_plane_crosswind_bounds_m": [target_y[0], target_y[-1]],
            "target_plane_height_bounds_m": [target_z[0], target_z[-1]],
            "target_cell_size_m": target_dx,
            "time_s": time_s.tolist(),
            "observation_value_fields_accessed": [],
        },
        "diagnostics": result.diagnostics,
        "integrated_flux_histories": integrals,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
        "next_required": (
            "connect time-interpolated mapped anomalies to the outer finite-volume "
            "inlet, add the declared ambient flux, and verify a coupled overlap case"
        ),
    }
    arrays = {
        "time": time_s, "source_y": source_y, "source_z": source_z,
        "target_y": target_y, "target_z": target_z,
        "source_species": species, "target_species": result.species_mass_flux_kg_m2_s,
    }
    return record, arrays


def _plot(record: dict[str, object], arrays: dict[str, np.ndarray], output: Path) -> None:
    index = 2
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    for axis, y_edges, z_edges, field, title in (
        (axes[0], arrays["source_y"], arrays["source_z"],
         arrays["source_species"][index], "resolved-inner plane (dx=1.25 m)"),
        (axes[1], arrays["target_y"], arrays["target_z"],
         arrays["target_species"][index], "outer plane remap (dx=5 m)"),
    ):
        image = axis.pcolormesh(y_edges, z_edges, field, shading="flat", cmap="magma")
        axis.set_xlabel("crosswind y (m)")
        axis.set_ylabel("height z (m)")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.78, label="species mass flux (kg/m2/s)")
    histories = record["integrated_flux_histories"]
    for name, values in histories.items():
        source = np.asarray(values["source"], dtype=float)
        scale = max(float(np.max(np.abs(source))), 1.0e-30)
        axes[2].plot(arrays["time"], source / scale, "o-", label=f"{name}: inner")
        axes[2].plot(
            arrays["time"], np.asarray(values["target"], dtype=float) / scale,
            "x--", linewidth=1.2, label=f"{name}: outer",
        )
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("normalised integrated plane flux")
    axes[2].set_title("Integral histories overlap after remap")
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=6, ncol=2)
    fig.suptitle("TI021 conservative nested-handoff operator audit")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-directory", type=Path,
        default=(ROOT / "outputs" / "thorney-ti021-nested-handoff-audit-2026-10-11-v1"),
    )
    args = parser.parse_args()
    record, arrays = build_audit()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-nested-handoff-audit.json"
    plot_path = args.output_directory / "ti021-nested-handoff-audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(record, arrays, plot_path)
    print(json.dumps({
        "status": record["status"], "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()), "gates": record["gates"],
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
