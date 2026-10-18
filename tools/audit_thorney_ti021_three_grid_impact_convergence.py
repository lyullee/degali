"""Audit the observation-blind TI021 5/2.5/1.25 m impact window.

The three runs must start from the same report-priority initial condition and
cover the frozen 0--20 s interval.  This audit characterises numerical
sensitivity only; it cannot promote the uncomputed 1.25 m, 600 s event.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from tools.compare_thorney_ti021_inner_grid_convergence import (
        coarsen_fluid_mask_two_to_one,
        coarsen_two_to_one,
        obstacle_region_metrics,
        trace_metrics,
    )
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from compare_thorney_ti021_inner_grid_convergence import (
        coarsen_fluid_mask_two_to_one,
        coarsen_two_to_one,
        obstacle_region_metrics,
        trace_metrics,
    )
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-three-grid-impact-convergence.v1"
GRID_RATIO = 2.0
GCI_SAFETY_FACTOR = 1.25
TI021_WALL_CENTER_M = (60.0, 100.0)
TI021_WALL_RADIUS_M = 50.0


def _observed_order_for_ratios(
    difference_ratio: float, first_ratio: float, second_ratio: float,
) -> float | None:
    """Solve the Richardson order for two possibly unequal refinement ratios."""
    if not all(
        math.isfinite(value) and value > 1.0
        for value in (first_ratio, second_ratio)
    ) or not math.isfinite(difference_ratio) or difference_ratio <= 0.0:
        return None
    limiting_ratio = math.log(first_ratio) / math.log(second_ratio)
    if difference_ratio <= limiting_ratio:
        return None

    def residual(order: float) -> float:
        return (first_ratio**order - 1.0) / (1.0 - second_ratio**(-order)) - difference_ratio

    lower = 1.0e-12
    upper = 1.0
    while residual(upper) < 0.0 and upper < 256.0:
        upper *= 2.0
    if residual(upper) < 0.0:
        return None
    for _ in range(100):
        midpoint = 0.5 * (lower + upper)
        if residual(midpoint) < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    order = 0.5 * (lower + upper)
    return float(order) if math.isfinite(order) and order > 0.0 else None


def scalar_triplet_convergence(
    coarse: float, medium: float, fine: float, *, grid_ratio: float = GRID_RATIO,
    next_grid_ratio: float | None = None,
) -> dict[str, float | str | bool | None]:
    """Return an observed-order/GCI record only for a defensible triplet."""

    values = np.asarray((coarse, medium, fine), dtype=float)
    if not np.all(np.isfinite(values)):
        return {
            "coarse": float(coarse), "medium": float(medium), "fine": float(fine),
            "classification": "non_finite", "monotonic": False,
            "observed_order": None, "fine_gci_absolute": None,
            "fine_gci_relative": None, "richardson_extrapolated": None,
        }
    delta_32 = float(coarse - medium)
    delta_21 = float(medium - fine)
    if delta_32 == 0.0 and delta_21 == 0.0:
        return {
            "coarse": float(coarse), "medium": float(medium), "fine": float(fine),
            "classification": "exact_agreement", "monotonic": True,
            "observed_order": None, "fine_gci_absolute": 0.0,
            "fine_gci_relative": 0.0, "richardson_extrapolated": float(fine),
        }
    monotonic = bool(delta_32 * delta_21 > 0.0)
    if not monotonic or delta_32 == 0.0 or delta_21 == 0.0:
        return {
            "coarse": float(coarse), "medium": float(medium), "fine": float(fine),
            "classification": "oscillatory_or_stationary", "monotonic": False,
            "observed_order": None, "fine_gci_absolute": None,
            "fine_gci_relative": None, "richardson_extrapolated": None,
        }
    second_ratio = grid_ratio if next_grid_ratio is None else next_grid_ratio
    order = _observed_order_for_ratios(
        abs(delta_32 / delta_21), grid_ratio, second_ratio,
    )
    if order is None:
        return {
            "coarse": float(coarse), "medium": float(medium), "fine": float(fine),
            "classification": "monotonic_nonconvergent", "monotonic": True,
            "observed_order": None,
            "fine_gci_absolute": None, "fine_gci_relative": None,
            "richardson_extrapolated": None,
        }
    denominator = second_ratio**order - 1.0
    fine_error = abs(float(fine - medium)) / denominator
    gci_absolute = GCI_SAFETY_FACTOR * fine_error
    return {
        "coarse": float(coarse), "medium": float(medium), "fine": float(fine),
        "classification": "monotonic_convergent", "monotonic": True,
        "observed_order": float(order),
        "fine_gci_absolute": float(gci_absolute),
        "fine_gci_relative": (
            None if abs(float(fine)) <= 1.0e-30
            else float(gci_absolute / abs(float(fine)))
        ),
        "richardson_extrapolated": float(
            fine + (fine - medium) / denominator
        ),
    }


def _error_norms(left: np.ndarray, right: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    difference = np.asarray(left, dtype=float)[mask] - np.asarray(right, dtype=float)[mask]
    return {
        "l1_mean_vol_percent": float(np.mean(np.abs(difference))),
        "rmse_vol_percent": float(np.sqrt(np.mean(difference**2))),
        "linf_vol_percent": float(np.max(np.abs(difference))),
    }


def _successive_field_convergence(
    coarse: np.ndarray, medium_on_coarse: np.ndarray,
    fine_on_coarse: np.ndarray, common_fluid: np.ndarray,
    *, refinement_ratio: float = GRID_RATIO,
) -> dict[str, object]:
    error_32 = _error_norms(coarse, medium_on_coarse, common_fluid)
    error_21 = _error_norms(medium_on_coarse, fine_on_coarse, common_fluid)
    result: dict[str, object] = {
        "common_fully_fluid_coarse_cell_count": int(np.count_nonzero(common_fluid)),
        "coarse_to_medium_error": error_32,
        "medium_to_fine_error_on_coarse_grid": error_21,
    }
    effective = {}
    for key in error_32:
        older = error_32[key]
        newer = error_21[key]
        effective[key] = {
            "error_contraction_ratio": (
                None if older <= 0.0 else float(newer / older)
            ),
            "effective_order_from_norms": (
                None if older <= 0.0 or newer <= 0.0
                else float(math.log(older / newer) / math.log(refinement_ratio))
            ),
            "contracts_with_refinement": bool(newer < older),
        }
    result["effective_convergence"] = effective
    return result


def hazard_geometry(field: dict[str, np.ndarray], *, cell_size_m: float) -> dict[str, object]:
    maximum = np.asarray(field["maximum_concentration_vol_percent"], dtype=float)
    mask = (~np.asarray(field["solid_mask"], dtype=bool)) & (maximum >= 4.0)
    count = int(np.count_nonzero(mask))
    result: dict[str, object] = {
        "cell_count": count,
        "volume_m3": float(count * cell_size_m**3),
    }
    for axis, coordinates, reduction in (
        ("x", np.asarray(field["x_m"], dtype=float), (0, 1)),
        ("y", np.asarray(field["y_m"], dtype=float), (0, 2)),
        ("z", np.asarray(field["z_m"], dtype=float), (1, 2)),
    ):
        occupied = np.any(mask, axis=reduction)
        if not np.any(occupied):
            result[f"{axis}_extent_m"] = None
            result[f"{axis}_minimum_face_m"] = None
            result[f"{axis}_maximum_face_m"] = None
        else:
            selected = coordinates[occupied]
            lower = float(np.min(selected) - 0.5 * cell_size_m)
            upper = float(np.max(selected) + 0.5 * cell_size_m)
            result[f"{axis}_extent_m"] = float(upper - lower)
            result[f"{axis}_minimum_face_m"] = lower
            result[f"{axis}_maximum_face_m"] = upper
    return result


def solid_geometry_metrics(
    field: dict[str, np.ndarray], *, cell_size_m: float,
) -> dict[str, object]:
    """Summarise voxelized solid geometry so grid-dependent blockage is visible."""

    solid = np.asarray(field["solid_mask"], dtype=bool)
    indices = np.nonzero(solid)
    result: dict[str, object] = {
        "solid_cell_count": int(indices[0].size),
        "represented_solid_voxel_volume_m3": float(
            indices[0].size * cell_size_m**3
        ),
    }
    horizontal = np.any(solid, axis=0)
    horizontal_indices = np.argwhere(horizontal)
    if horizontal_indices.size:
        y_indices, x_indices = horizontal_indices[:, 0], horizontal_indices[:, 1]
        radii = np.hypot(
            np.asarray(field["x_m"], dtype=float)[x_indices]
            - TI021_WALL_CENTER_M[0],
            np.asarray(field["y_m"], dtype=float)[y_indices]
            - TI021_WALL_CENTER_M[1],
        )
        arc_length = float(horizontal_indices.shape[0] * cell_size_m)
        target_arc_length = math.pi * TI021_WALL_RADIUS_M
        result.update({
            "horizontal_solid_cell_count": int(horizontal_indices.shape[0]),
            "maximum_solid_center_radial_error_m": float(np.max(
                np.abs(radii - TI021_WALL_RADIUS_M)
            )),
            "maximum_solid_center_radial_error_cells": float(np.max(
                np.abs(radii - TI021_WALL_RADIUS_M)
            ) / cell_size_m),
            "estimated_arc_length_m": arc_length,
            "target_semicircle_arc_length_m": target_arc_length,
            "relative_arc_length_error": abs(
                arc_length - target_arc_length
            ) / target_arc_length,
        })
    else:
        result.update({
            "horizontal_solid_cell_count": 0,
            "maximum_solid_center_radial_error_m": None,
            "maximum_solid_center_radial_error_cells": None,
            "estimated_arc_length_m": 0.0,
            "target_semicircle_arc_length_m": math.pi * TI021_WALL_RADIUS_M,
            "relative_arc_length_error": 1.0,
        })
    for axis, coordinates, index_axis in (
        ("x", np.asarray(field["x_m"], dtype=float), 2),
        ("y", np.asarray(field["y_m"], dtype=float), 1),
        ("z", np.asarray(field["z_m"], dtype=float), 0),
    ):
        selected = indices[index_axis]
        if selected.size:
            lower = float(np.min(coordinates[selected]) - 0.5 * cell_size_m)
            upper = float(np.max(coordinates[selected]) + 0.5 * cell_size_m)
            result[f"{axis}_bounds_m"] = [lower, upper]
            result[f"{axis}_extent_m"] = upper - lower
        else:
            result[f"{axis}_bounds_m"] = None
            result[f"{axis}_extent_m"] = None
    return result


def _load_run(directory: Path) -> dict[str, object]:
    resolved = directory.resolve()
    record_path = resolved / "ti021-inner-transient.json"
    field_path = resolved / "ti021-inner-final-field.npz"
    trace_path = resolved / "ti021-inner-receptor-traces.npz"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    with np.load(field_path, allow_pickle=False) as source:
        field = {key: np.asarray(source[key]).copy() for key in source.files}
    with np.load(trace_path, allow_pickle=False) as source:
        trace = {key: np.asarray(source[key]).copy() for key in source.files}
    return {
        "directory": resolved, "record_path": record_path, "field_path": field_path,
        "trace_path": trace_path, "record": record, "field": field, "trace": trace,
    }


def _receptor_convergence(
    coarse_trace: dict[str, np.ndarray], medium_trace: dict[str, np.ndarray],
    fine_trace: dict[str, np.ndarray], *, refinement_ratios: tuple[float, float] = (GRID_RATIO, GRID_RATIO),
) -> dict[str, object]:
    for candidate in (medium_trace, fine_trace):
        if not np.array_equal(coarse_trace["labels"], candidate["labels"]):
            raise ValueError("all receptor labels must match exactly")
        if not np.allclose(
            coarse_trace["solver_positions_m"], candidate["solver_positions_m"],
            rtol=0.0, atol=1.0e-12,
        ):
            raise ValueError("all receptor coordinates must match")
    metrics = [
        trace_metrics(item["time_s"], item["concentration_vol_percent"], threshold=0.1)
        for item in (coarse_trace, medium_trace, fine_trace)
    ]
    outputs: dict[str, object] = {}
    eligibility = {
        "peak_vol_percent": metrics[2]["peak_vol_percent"] >= 0.1,
        "dose_vol_percent_s": metrics[2]["dose_vol_percent_s"] >= 0.1,
        "arrival_s": np.logical_and.reduce([
            np.isfinite(item["arrival_s"]) for item in metrics
        ]),
    }
    labels = [str(value) for value in coarse_trace["labels"]]
    for metric_name, eligible in eligibility.items():
        per_receptor = []
        relative_gci = []
        monotonic_convergent = 0
        for index in np.flatnonzero(eligible):
            audit = scalar_triplet_convergence(
                metrics[0][metric_name][index], metrics[1][metric_name][index],
                metrics[2][metric_name][index],
                grid_ratio=refinement_ratios[0], next_grid_ratio=refinement_ratios[1],
            )
            audit["label"] = labels[int(index)]
            per_receptor.append(audit)
            if audit["classification"] in {"monotonic_convergent", "exact_agreement"}:
                monotonic_convergent += 1
            if audit["fine_gci_relative"] is not None:
                relative_gci.append(100.0 * float(audit["fine_gci_relative"]))
        count = len(per_receptor)
        outputs[metric_name] = {
            "eligibility_floor": 0.1 if metric_name != "arrival_s" else "finite_on_all_grids",
            "eligible_receptor_count": count,
            "monotonic_convergent_count": monotonic_convergent,
            "monotonic_convergent_fraction": (
                None if count == 0 else float(monotonic_convergent / count)
            ),
            "fine_gci_relative_percent_median": (
                None if not relative_gci else float(np.median(relative_gci))
            ),
            "fine_gci_relative_percent_p90": (
                None if not relative_gci else float(np.percentile(relative_gci, 90.0))
            ),
            "per_receptor": per_receptor,
        }
    return outputs


def audit_three_grid_runs(
    *, coarse_run: Path, medium_run: Path, fine_run: Path, output_directory: Path,
    allow_short_prefix_diagnostic: bool = False,
) -> tuple[dict[str, object], dict[str, Path]]:
    runs = [_load_run(path) for path in (coarse_run, medium_run, fine_run)]
    records = [item["record"] for item in runs]
    if not all(record.get("all_gates_pass") for record in records):
        raise ValueError("all three grid runs must pass numerical gates")
    cell_sizes = [float(record["configuration"]["cell_size_m"]) for record in records]
    time_steps = [float(record["configuration"]["time_step_s"]) for record in records]
    durations = [float(record["configuration"]["duration_s"]) for record in records]
    representations = [
        str(record["configuration"].get(
            "circular_arc_representation", "finite_thickness_voxel",
        ))
        for record in records
    ]
    zero_surface_geometry_held = all(
        representation == "zero_thickness_surface"
        for representation in representations
    )
    space_ratios = np.asarray(cell_sizes[:-1]) / np.asarray(cell_sizes[1:])
    time_ratios = np.asarray(time_steps[:-1]) / np.asarray(time_steps[1:])
    if not np.allclose(space_ratios, GRID_RATIO):
        raise ValueError("spatial refinement ratios must be exactly 2:1")
    if not np.all(np.isfinite(time_ratios)) or np.any(time_ratios < 1.0 - 1.0e-12):
        raise ValueError(
            "time steps must be matched or refined with each spatial refinement"
        )
    full_impact_window = bool(np.allclose(durations, durations[0]) and np.isclose(durations[0], 20.0))
    short_prefix = bool(
        np.allclose(durations, durations[0]) and 0.0 < durations[0] < 20.0
    )
    if not full_impact_window and not (
        allow_short_prefix_diagnostic and short_prefix
    ):
        raise ValueError(
            "runs must cover the frozen 0--20 s window; shorter matched prefixes "
            "require allow_short_prefix_diagnostic=True"
        )

    fields = [item["field"] for item in runs]
    if tuple(2 * size for size in fields[0]["solid_mask"].shape) != fields[1]["solid_mask"].shape:
        raise ValueError("medium grid dimensions must be twice coarse dimensions")
    if tuple(2 * size for size in fields[1]["solid_mask"].shape) != fields[2]["solid_mask"].shape:
        raise ValueError("fine grid dimensions must be twice medium dimensions")
    medium_fluid_coarse = coarsen_fluid_mask_two_to_one(fields[1]["solid_mask"])
    fine_fluid_medium = coarsen_fluid_mask_two_to_one(fields[2]["solid_mask"])
    fine_fluid_coarse = coarsen_fluid_mask_two_to_one(~fine_fluid_medium)
    common_fluid = (~fields[0]["solid_mask"]) & medium_fluid_coarse & fine_fluid_coarse

    field_convergence = {}
    for source_key, label in (
        ("concentration_vol_percent", "final_concentration"),
        ("maximum_concentration_vol_percent", "maximum_concentration"),
    ):
        medium_on_coarse = coarsen_two_to_one(fields[1][source_key])
        fine_on_medium = coarsen_two_to_one(fields[2][source_key])
        fine_on_coarse = coarsen_two_to_one(fine_on_medium)
        field_convergence[label] = _successive_field_convergence(
            fields[0][source_key], medium_on_coarse, fine_on_coarse, common_fluid,
            refinement_ratio=GRID_RATIO,
        )

    hazard = [hazard_geometry(field, cell_size_m=dx) for field, dx in zip(fields, cell_sizes)]
    solid_geometry = [
        solid_geometry_metrics(field, cell_size_m=dx)
        for field, dx in zip(fields, cell_sizes)
    ]
    hazard_convergence = {}
    for key in ("volume_m3", "x_extent_m", "y_extent_m", "z_extent_m"):
        hazard_convergence[key] = scalar_triplet_convergence(
            hazard[0][key], hazard[1][key], hazard[2][key],
        ) if all(item[key] is not None for item in hazard) else {
            "classification": "hazard_region_absent", "monotonic": False,
            "observed_order": None, "fine_gci_absolute": None,
            "fine_gci_relative": None, "richardson_extrapolated": None,
        }

    regions = [
        obstacle_region_metrics(field, cell_size_m=dx)
        for field, dx in zip(fields, cell_sizes)
    ]
    regional_convergence = {}
    for region_name in regions[0]:
        regional_convergence[region_name] = {}
        for metric_name in (
            "maximum_ever_mean_vol_percent", "maximum_ever_peak_vol_percent",
            "maximum_ever_four_vol_percent_volume_m3",
        ):
            regional_convergence[region_name][metric_name] = scalar_triplet_convergence(
                regions[0][region_name][metric_name],
                regions[1][region_name][metric_name],
                regions[2][region_name][metric_name],
            )

    receptor_convergence = _receptor_convergence(
        runs[0]["trace"], runs[1]["trace"], runs[2]["trace"],
        refinement_ratios=(GRID_RATIO, GRID_RATIO),
    )
    field_contracts = all(
        field_convergence[field_name]["effective_convergence"][norm]["contracts_with_refinement"]
        for field_name in field_convergence for norm in ("l1_mean_vol_percent", "rmse_vol_percent")
    )
    hazard_volume = hazard_convergence["volume_m3"]
    hazard_gci_ok = bool(
        hazard_volume.get("fine_gci_relative") is not None
        and float(hazard_volume["fine_gci_relative"]) <= 0.10
    )
    peak_summary = receptor_convergence["peak_vol_percent"]
    dose_summary = receptor_convergence["dose_vol_percent_s"]
    receptor_gate = all(
        summary["monotonic_convergent_fraction"] is not None
        and float(summary["monotonic_convergent_fraction"]) >= 0.70
        and summary["fine_gci_relative_percent_median"] is not None
        and float(summary["fine_gci_relative_percent_median"]) <= 20.0
        and summary["fine_gci_relative_percent_p90"] is not None
        and float(summary["fine_gci_relative_percent_p90"]) <= 50.0
        for summary in (peak_summary, dose_summary)
    )
    gates = {
        "all_run_numerical_gates_pass": True,
        "frozen_observation_blind_20s_window": full_impact_window,
        "exact_two_to_one_spatial_refinement": True,
        "time_steps_matched_or_refined_with_grid": True,
        # An ideal embedded surface permits numerical grid convergence of the
        # stated model, but the report does not identify the real fence's
        # physical thickness. Keep those two claims separate.
        "numerical_obstacle_representation_held_constant": zero_surface_geometry_held,
        "physical_fence_thickness_known": False,
        "physical_obstacle_geometry_identified_and_held_constant": False,
        "field_l1_and_rmse_contract": field_contracts,
        "hazard_volume_monotonic_gci_le_10_percent": hazard_gci_ok,
        "receptor_peak_and_dose_monotonic_fraction_ge_70_percent": receptor_gate,
    }
    impact_supported = bool(all(
        value for key, value in gates.items()
        if key != "physical_fence_thickness_known"
    ))
    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=False)
    record = {
        "schema": SCHEMA,
        "status": (
            "impact_window_grid_sensitivity_supported_full_event_withheld"
            if impact_supported else
            "short_prefix_diagnostic_only" if short_prefix else
            "impact_window_grid_sensitivity_not_converged_full_event_withheld"
        ),
        "scope": (
            (
                f"observation-blind 0--{durations[0]:g} s matched prefix diagnostic; "
                "not a frozen 0--20 s impact-window qualification"
                if short_prefix else
                "observation-blind 0--20 s initial-cloud/fence-impact grid audit; "
            )
            +
            "the 600 s finest-grid event and quantitative facility use remain withheld"
        ),
        "predeclared_engineering_diagnostic_policy": {
            "field": "L1 and RMSE must contract for final and maximum fields",
            "hazard_volume": "monotonic convergent and fine-grid GCI <= 10%",
            "receptors": (
                "peak and dose each require >=70% monotonic-convergent eligible "
                "triplets, median fine GCI <=20%, and p90 fine GCI <=50%"
            ),
            "gci_safety_factor": GCI_SAFETY_FACTOR,
            "note": "diagnostic acceptance policy, not a universal ASME certification claim",
        },
        "provenance": {
            label: {
                "record_path": str(run["record_path"]),
                "record_sha256": _sha256(run["record_path"]),
                "field_sha256": _sha256(run["field_path"]),
                "trace_sha256": _sha256(run["trace_path"]),
            }
            for label, run in zip(("coarse", "medium", "fine"), runs)
        } | {"observation_value_fields_accessed": []},
        "configuration": {
            "cell_sizes_m": cell_sizes, "time_steps_s": time_steps,
            "duration_s": durations[0], "spatial_grid_ratio": GRID_RATIO,
            "window_classification": (
                "frozen_0_20s_impact_window" if full_impact_window else
                "matched_short_prefix_diagnostic_only"
            ),
            "temporal_refinement_ratios": time_ratios.tolist(),
            "temporal_comparison_mode": (
                "matched_time_step_isolates_spatial_grid_effects"
                if np.allclose(time_ratios, 1.0)
                else "time_and_space_refined_together"
            ),
            "receptor_gci_refinement_ratios": [GRID_RATIO, GRID_RATIO],
        },
        "voxelized_obstacle_geometry": {
            "reference_geometry": {
                "shape": "semi-circular impermeable tarpaulin sheets hung over a scaffold frame",
                "center_m_solver_xy": list(TI021_WALL_CENTER_M),
                "radius_m": TI021_WALL_RADIUS_M,
                "height_m": 5.0,
                "arc_extent_deg": 180.0,
                "published_construction": "tarpaulin sheets hung over a scaffold frame",
                "physical_thickness": "tarpaulin gauge not reported; scaffold member sections not reported",
                "numerical_idealization": "impermeable zero-thickness surface; scaffold frame omitted",
                "construction_source": "Davies and Singh (1985), The Phase II Trials: A Data Set on the Effect of Obstructions, doi:10.1016/0304-3894(85)85043-3",
            },
            "by_grid": {
                label: value for label, value in zip(
                    ("coarse", "medium", "fine"), solid_geometry,
                )
            },
            "boundary_representation_by_grid": dict(zip(
                ("coarse", "medium", "fine"), representations,
            )),
            "declared_numerical_wall_thickness_by_grid_m": [
                0.0 if item == "zero_thickness_surface" else dx
                for item, dx in zip(representations, cell_sizes)
            ],
            "declared_wall_thickness_by_grid_m": [
                0.0 if item == "zero_thickness_surface" else dx
                for item, dx in zip(representations, cell_sizes)
            ],
            "wall_thickness_to_cell_size_ratio_by_grid": [
                0.0 if item == "zero_thickness_surface" else 1.0
                for item in representations
            ],
            "numerical_geometry_held_constant": zero_surface_geometry_held,
            "physical_geometry_held_constant": False,
            "physical_fence_thickness_known": False,
            "geometry_comparison_classification": (
                "grid_independent_zero_thickness_embedded_surface"
                if zero_surface_geometry_held else
                "one_cell_surface_regularisation_refinement_not_fixed_physical_geometry"
            ),
            "geometry_builder_source": "tools/audit_thorney_ti021_handoff_streaming_preflight.py",
            "geometry_builder_source_sha256": _sha256(
                Path(__file__).resolve().parents[1]
                / "tools/audit_thorney_ti021_handoff_streaming_preflight.py"
            ),
            "interpretation": (
                "All runs use the same ideal zero-thickness circular surface; "
                "numerical grid sensitivity is interpretable for that model, "
                "but the real fence thickness remains unknown."
                if zero_surface_geometry_held else
                "The voxel wall thickness follows dx, so the comparison mixes "
                "spatial discretisation and geometry regularisation."
            ),
        },
        "grid_comparison_interpretation": (
            "Matched time steps and the grid-independent ideal surface isolate "
            "numerical spatial sensitivity for the zero-thickness model; this "
            "does not establish the unknown physical fence thickness."
            if zero_surface_geometry_held else
            "Matched time steps isolate temporal discretisation, but the wall "
            "thickness also changes with dx. Field, hazard, and receptor "
            "differences combine spatial discretisation with the one-cell "
            "wall-regularisation path."
        ),
        "field_convergence": field_convergence,
        "four_vol_percent_geometry_by_grid": {
            label: value for label, value in zip(("coarse", "medium", "fine"), hazard)
        },
        "four_vol_percent_geometry_convergence": hazard_convergence,
        "obstacle_region_metrics_by_grid": {
            label: value for label, value in zip(("coarse", "medium", "fine"), regions)
        },
        "obstacle_region_convergence": regional_convergence,
        "receptor_convergence": receptor_convergence,
        "gates": gates,
        "impact_window_grid_sensitivity_supported": impact_supported,
        "full_event_grid_convergence_promoted": False,
        "next_required": (
            "the matched-dt 0--20 s spatial audit is complete for the ideal zero-thickness surface, but hazard-volume GCI and receptor peak/dose convergence still fail; improve those numerical convergence results and retain the 600-s finest-grid holdback. Separately register scaffold member geometry and tarp gauge, or justify their omission, before facility-level physical interpretation"
            if zero_surface_geometry_held else
            "register a physically justified obstacle thickness or implement a "
            "grid-independent embedded-surface boundary; rerun matched-geometry "
            "spatial convergence before interpreting GCI, then retain the full-event "
            "finest-grid holdback and combine with the 600 s nested outer transient"
        ),
    }
    json_path = output / "ti021-three-grid-impact-convergence.json"
    plot_path = output / "ti021-three-grid-impact-convergence.png"
    json_path.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    figure, axes = plt.subplots(2, 2, figsize=(12.5, 9.0), constrained_layout=True)
    x = np.arange(3)
    for field_name, color in (("final_concentration", "#4C78A8"), ("maximum_concentration", "#E45756")):
        values = [
            field_convergence[field_name]["coarse_to_medium_error"]["rmse_vol_percent"],
            field_convergence[field_name]["medium_to_fine_error_on_coarse_grid"]["rmse_vol_percent"],
        ]
        axes[0, 0].plot((0, 1), values, marker="o", label=field_name.replace("_", " "), color=color)
    axes[0, 0].set(xticks=(0, 1), xticklabels=("5→2.5 m", "2.5→1.25 m"), ylabel="RMSE (vol%)", title="Common-grid field error contraction")
    axes[0, 0].grid(alpha=0.25); axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(cell_sizes, [item["volume_m3"] for item in hazard], marker="o")
    axes[0, 1].invert_xaxis(); axes[0, 1].set(xlabel="Cell size (m)", ylabel="C ≥ 4 vol% volume (m³)", title="Maximum-ever flammable volume")
    axes[0, 1].grid(alpha=0.25)
    traces = [trace_metrics(run["trace"]["time_s"], run["trace"]["concentration_vol_percent"], threshold=0.1) for run in runs]
    axes[1, 0].scatter(traces[2]["peak_vol_percent"], traces[0]["peak_vol_percent"], s=18, label="5 m", alpha=0.75)
    axes[1, 0].scatter(traces[2]["peak_vol_percent"], traces[1]["peak_vol_percent"], s=18, label="2.5 m", alpha=0.75)
    peak_limit = max(float(np.max(item["peak_vol_percent"])) for item in traces)
    axes[1, 0].plot((0, peak_limit), (0, peak_limit), "k--", linewidth=1)
    axes[1, 0].set(xlabel="1.25 m peak (vol%)", ylabel="Coarser-grid peak (vol%)", title="Receptor peak convergence")
    axes[1, 0].grid(alpha=0.25); axes[1, 0].legend()
    for region_name in regions[0]:
        axes[1, 1].plot(cell_sizes, [item[region_name]["maximum_ever_mean_vol_percent"] for item in regions], marker="o", label=region_name)
    axes[1, 1].invert_xaxis(); axes[1, 1].set(xlabel="Cell size (m)", ylabel="Maximum-ever mean (vol%)", title="Fixed obstacle-region sensitivity")
    axes[1, 1].grid(alpha=0.25); axes[1, 1].legend(fontsize=7)
    figure.suptitle("TI021 three-grid obstacle-impact audit (full-event convergence withheld)")
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    return record, {"json": json_path, "plot": plot_path}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse-run", type=Path, required=True)
    parser.add_argument("--medium-run", type=Path, required=True)
    parser.add_argument("--fine-run", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument(
        "--allow-short-prefix-diagnostic", action="store_true",
        help="allow a matched <20 s prefix while keeping the 20 s gate false",
    )
    args = parser.parse_args()
    record, paths = audit_three_grid_runs(
        coarse_run=args.coarse_run, medium_run=args.medium_run,
        fine_run=args.fine_run, output_directory=args.output_directory,
        allow_short_prefix_diagnostic=args.allow_short_prefix_diagnostic,
    )
    print(json.dumps({
        "status": record["status"], "gates": record["gates"],
        "impact_window_grid_sensitivity_supported": record["impact_window_grid_sensitivity_supported"],
        "full_event_grid_convergence_promoted": record["full_event_grid_convergence_promoted"],
        "outputs": {key: str(value) for key, value in paths.items()},
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
