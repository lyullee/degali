"""Audit diagnostic-only dense-layer stratification in a completed TI021 run.

Only model density, horizontal velocity, geometry, and coordinate-only
receptor locations are read.  Concentration, dose, arrival, and departure
observations are not accessed, and this audit never modifies transport fluxes.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

# The audit is batch/reproducibility tooling; never require a desktop Tk
# session merely to render its diagnostic figures.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from degali.addons.atmospheric_surface_layer_mixing import (
    NeutralSurfaceLayerMixing,
)
from degali.addons.dense_gas_stratification_3d import (
    diagnose_gradient_richardson_3d,
    summarize_gradient_richardson_regions,
)

try:
    from tools.compare_thorney_ti021_inner_grid_convergence import (
        obstacle_region_masks,
    )
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from compare_thorney_ti021_inner_grid_convergence import (
        obstacle_region_masks,
    )
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-stratification-diagnostic.v1"
CHECKPOINT_SCHEMA = "degali.thorney-ti021-inner-checkpoint.v1"
POLICY_SCHEMA = "degali.thorney-ti021-turbulent-mixing-policy.v1"
PLAN_SCHEMA = "degali.thorney-ti021-stratification-diagnostic-plan.v1"


def _nearest_cell(
    x_m: np.ndarray,
    y_m: np.ndarray,
    z_m: np.ndarray,
    position_m: tuple[float, float, float],
) -> tuple[int, int, int]:
    x, y, z = (float(value) for value in position_m)
    return (
        int(np.argmin(np.abs(z_m - z))),
        int(np.argmin(np.abs(y_m - y))),
        int(np.argmin(np.abs(x_m - x))),
    )


def _load_policy(path: Path) -> tuple[dict[str, object], dict[str, float]]:
    policy = json.loads(path.read_text(encoding="utf-8"))
    if policy.get("schema") != POLICY_SCHEMA:
        raise ValueError("unsupported TI021 turbulent-mixing policy schema")
    if not policy.get("all_gates_pass"):
        raise ValueError("TI021 turbulent-mixing policy gates did not pass")
    if policy["evidence_provenance"].get(
        "concentration_observation_fields_accessed"
    ) != []:
        raise ValueError("mixing policy is not concentration-blind")
    closure = policy["closure"]
    frozen = policy["frozen_policy"]
    values = {
        "friction_velocity_m_s": float(closure["friction_velocity_m_s"]),
        "roughness_length_m": float(closure["roughness_length_m"]),
        "von_karman_constant": float(closure["von_karman_constant"]),
        "effective_schmidt_number": float(
            frozen["primary_turbulent_schmidt_number"]
        ),
    }
    return policy, values


def _validate_plan(
    path: Path,
    *,
    mixing_policy_path: Path,
    shear_regularization_s_inv: float,
) -> dict[str, object]:
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema") != PLAN_SCHEMA or not plan.get("all_gates_pass"):
        raise ValueError("accepted TI021 stratification diagnostic plan is required")
    if plan["provenance"].get("observation_value_fields_accessed") != []:
        raise ValueError("stratification diagnostic plan is not observation-blind")
    if plan["provenance"].get("mixing_policy_sha256") != _sha256(
        mixing_policy_path
    ):
        raise ValueError("diagnostic plan mixing-policy SHA-256 mismatch")
    allowed = np.asarray([
        float(item["shear_regularization_s_inv"])
        for item in plan["members"]
    ])
    if not np.any(np.isclose(
        float(shear_regularization_s_inv), allowed, rtol=0.0, atol=1.0e-15,
    )):
        raise ValueError("shear regularization is not a frozen plan member")
    return plan


def audit_stratification(
    *,
    run_directory: Path,
    mixing_policy_path: Path,
    diagnostic_plan_path: Path,
    shear_regularization_s_inv: float,
    output_directory: Path,
) -> tuple[dict[str, object], dict[str, Path]]:
    run = run_directory.resolve()
    policy_path = mixing_policy_path.resolve()
    plan_path = diagnostic_plan_path.resolve()
    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=False)
    journal_path = run / "ti021-inner-checkpoint.json"
    run_record_path = run / "ti021-inner-transient.json"
    journal = json.loads(journal_path.read_text(encoding="utf-8"))
    run_record = json.loads(run_record_path.read_text(encoding="utf-8"))
    if journal.get("schema") != CHECKPOINT_SCHEMA:
        raise ValueError("unsupported TI021 checkpoint schema")
    if not run_record.get("all_gates_pass"):
        raise ValueError("TI021 run numerical gates did not pass")
    if run_record.get("configuration", {}).get(
        "observation_values_used_to_configure_run"
    ) != []:
        raise ValueError("TI021 run is not observation-blind")
    duration = float(journal["configuration"]["duration_s"])
    entries = list(journal.get("segments", []))
    if not entries or not np.isclose(float(entries[-1]["stop_time_s"]), duration):
        raise ValueError("TI021 checkpoint journal is not complete")

    _policy, policy_values = _load_policy(policy_path)
    _validate_plan(
        plan_path,
        mixing_policy_path=policy_path,
        shear_regularization_s_inv=shear_regularization_s_inv,
    )
    mixing = NeutralSurfaceLayerMixing(
        friction_velocity_m_s=policy_values["friction_velocity_m_s"],
        roughness_length_m=policy_values["roughness_length_m"],
        von_karman_constant=policy_values["von_karman_constant"],
        source_id="hash-pinned TI021 turbulent-mixing policy",
    )
    schmidt = policy_values["effective_schmidt_number"]
    dx = float(journal["configuration"]["cell_size_m"])
    receptor_configuration = run_record["configuration"]["receptors"]
    segment_records: list[dict[str, object]] = []
    receptor_cells: list[dict[str, object]] | None = None
    first_grid: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None

    for entry in entries:
        segment_path = run / str(entry["segment_file"])
        if _sha256(segment_path) != entry["segment_sha256"]:
            raise ValueError(f"segment SHA-256 mismatch: {segment_path}")
        with np.load(segment_path, allow_pickle=False) as payload:
            x_m = np.asarray(payload["x_m"], dtype=float)
            y_m = np.asarray(payload["y_m"], dtype=float)
            z_m = np.asarray(payload["z_m"], dtype=float)
            solid = np.asarray(payload["solid_mask"], dtype=bool)
            density = np.asarray(
                payload["final_mixture_density_kg_m3"], dtype=float,
            )
            velocity_u = np.asarray(payload["final_velocity_u_m_s"], dtype=float)
            velocity_v = np.asarray(payload["final_velocity_v_m_s"], dtype=float)
        if first_grid is None:
            first_grid = (x_m.copy(), y_m.copy(), z_m.copy(), solid.copy())
        else:
            first_x, first_y, first_z, first_solid = first_grid
            if not (
                np.array_equal(x_m, first_x)
                and np.array_equal(y_m, first_y)
                and np.array_equal(z_m, first_z)
                and np.array_equal(solid, first_solid)
            ):
                raise ValueError("grid or solid mask changed between checkpoints")

        diagnostic = diagnose_gradient_richardson_3d(
            density,
            velocity_u,
            velocity_v,
            solid,
            grid_spacing_m=dx,
            shear_regularization_s_inv=shear_regularization_s_inv,
        )
        neutral_viscosity = mixing.eddy_viscosity_m2_s(z_m, solid)
        neutral_scalar_diffusivity = neutral_viscosity / schmidt
        regions = summarize_gradient_richardson_regions(
            diagnostic,
            obstacle_region_masks(x_m, y_m, z_m, solid),
            neutral_scalar_eddy_diffusivity_m2_s=neutral_scalar_diffusivity,
            cell_volume_m3=dx**3,
        )["regions"]

        if receptor_cells is None:
            receptor_cells = []
            for receptor in receptor_configuration:
                position = tuple(float(value) for value in receptor["solver_xyz_m"])
                iz, iy, ix = _nearest_cell(x_m, y_m, z_m, position)
                if solid[iz, iy, ix]:
                    raise ValueError(
                        f"receptor {receptor['label']} maps to a solid cell"
                    )
                center = (float(x_m[ix]), float(y_m[iy]), float(z_m[iz]))
                receptor_cells.append({
                    "label": receptor["label"],
                    "solver_position_m": list(position),
                    "cell_index_zyx": [iz, iy, ix],
                    "cell_center_xyz_m": list(center),
                    "cell_center_distance_m": float(
                        np.linalg.norm(np.asarray(center) - np.asarray(position))
                    ),
                })
        receptor_values = []
        for receptor in receptor_cells:
            iz, iy, ix = receptor["cell_index_zyx"]
            receptor_values.append({
                "label": receptor["label"],
                "gradient_richardson_number": float(
                    diagnostic.gradient_richardson_number[iz, iy, ix]
                ),
                "stable_buoyancy_frequency_squared_s_inv2": float(
                    diagnostic.stable_buoyancy_frequency_squared_s_inv2[
                        iz, iy, ix
                    ]
                ),
                "horizontal_shear_squared_s_inv2": float(
                    diagnostic.horizontal_shear_squared_s_inv2[iz, iy, ix]
                ),
                "neutral_scalar_eddy_diffusivity_m2_s": float(
                    neutral_scalar_diffusivity[iz, iy, ix]
                ),
            })
        segment_records.append({
            "segment_index": int(entry["segment_index"]),
            "stop_time_s": float(entry["stop_time_s"]),
            "segment_path": str(segment_path),
            "segment_sha256": entry["segment_sha256"],
            "regions": regions,
            "receptors": receptor_values,
        })

    assert receptor_cells is not None
    gates = {
        "completed_run_numerical_gates_pass": bool(run_record["all_gates_pass"]),
        "checkpoint_journal_covers_full_declared_duration": True,
        "all_checkpoint_segment_hashes_match": True,
        "grid_and_solid_mask_constant_between_checkpoints": True,
        "all_receptor_sampling_cells_are_fluid": True,
        "run_and_policy_use_no_concentration_observation_for_configuration": True,
        "diagnostic_changes_no_transport_flux": True,
    }
    record = {
        "schema": SCHEMA,
        "status": "accepted_diagnostic_only_stratification_audit",
        "scope": (
            "model-only gradient Richardson diagnosis; not a damping-law "
            "validation and not hydrogen accuracy evidence"
        ),
        "provenance": {
            "run_directory": str(run),
            "run_record_sha256": _sha256(run_record_path),
            "checkpoint_journal_sha256": _sha256(journal_path),
            "mixing_policy_path": str(policy_path),
            "mixing_policy_sha256": _sha256(policy_path),
            "diagnostic_plan_path": str(plan_path),
            "diagnostic_plan_sha256": _sha256(plan_path),
            "observation_value_fields_accessed": [],
        },
        "configuration": {
            "cell_size_m": dx,
            "duration_s": duration,
            "shear_regularization_s_inv": float(
                shear_regularization_s_inv
            ),
            **policy_values,
            "neutral_scalar_diffusivity_formula": "kappa*u_star*z/Sc_t",
            "damping_applied": False,
        },
        "receptor_sampling_cells": receptor_cells,
        "segments": segment_records,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
        "promotion_policy": {
            "closure_promotion_allowed": False,
            "reason": (
                "diagnostic exposure must be evaluated after numerical "
                "grid/time checks and against independent evidence"
            ),
        },
    }
    json_path = output / "ti021-stratification-diagnostic.json"
    figure_path = output / "ti021-stratification-diagnostic.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    times = np.asarray([item["stop_time_s"] for item in segment_records])
    region_names = list(segment_records[0]["regions"])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for name in region_names:
        axes[0, 0].plot(times, [
            item["regions"][name]["gradient_richardson_p50"]
            for item in segment_records
        ], marker="o", label=name)
        axes[0, 1].plot(times, [
            item["regions"][name]["gradient_richardson_p95"]
            for item in segment_records
        ], marker="o", label=name)
        axes[1, 0].plot(times, [
            item["regions"][name]["fraction_at_or_above_0p25"]
            for item in segment_records
        ], marker="o", label=name)
    receptor_matrix = np.asarray([
        [value["gradient_richardson_number"] for value in item["receptors"]]
        for item in segment_records
    ])
    axes[1, 1].plot(times, np.quantile(receptor_matrix, 0.50, axis=1),
                    marker="o", label="sensor-cell p50")
    axes[1, 1].plot(times, np.quantile(receptor_matrix, 0.95, axis=1),
                    marker="o", label="sensor-cell p95")
    axes[0, 0].set(title="Regional Ri_g median", ylabel="Ri_g")
    axes[0, 1].set(title="Regional Ri_g p95", ylabel="Ri_g")
    axes[1, 0].set(title="Regional fraction Ri_g >= 0.25",
                   xlabel="Model time (s)", ylabel="Fluid-cell fraction")
    axes[1, 1].set(title="Receptor sampling-cell Ri_g",
                   xlabel="Model time (s)", ylabel="Ri_g")
    for axis in axes.flat:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    fig.suptitle("TI021 diagnostic-only dense-layer stratification")
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    return record, {"json": json_path, "figure": figure_path}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-directory", type=Path, required=True)
    parser.add_argument("--mixing-policy", type=Path, required=True)
    parser.add_argument("--diagnostic-plan", type=Path, required=True)
    parser.add_argument(
        "--shear-regularization-s-inv", type=float, required=True,
    )
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = audit_stratification(
        run_directory=args.run_directory,
        mixing_policy_path=args.mixing_policy,
        diagnostic_plan_path=args.diagnostic_plan,
        shear_regularization_s_inv=args.shear_regularization_s_inv,
        output_directory=args.output_directory,
    )
    print(json.dumps({
        "status": record["status"],
        "all_gates_pass": record["all_gates_pass"],
        "outputs": {name: str(path) for name, path in paths.items()},
    }, indent=2))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
