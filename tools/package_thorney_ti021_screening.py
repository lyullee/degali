"""Package a completed TI021 full-horizon screening run for review.

This reads an existing run only; it never changes the solver state or accesses
observation values.  The output deliberately labels the result as a 2.5 m
screening tier while exposing the same front/roof/side/wake regions used by the
three-grid convergence audit.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

try:
    from tools.compare_thorney_ti021_inner_grid_convergence import (
        obstacle_region_masks,
        obstacle_region_metrics,
        trace_metrics,
    )
except ModuleNotFoundError:
    from compare_thorney_ti021_inner_grid_convergence import (
        obstacle_region_masks,
        obstacle_region_metrics,
        trace_metrics,
    )


def package(
    run: Path,
    output: Path,
    *,
    threshold: float = 4.0,
    score: Path | None = None,
) -> dict[str, object]:
    run = run.resolve()
    record_path = run / "ti021-inner-transient.json"
    field_path = run / "ti021-inner-final-field.npz"
    trace_path = run / "ti021-inner-receptor-traces.npz"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    with np.load(field_path, allow_pickle=False) as source:
        field = {key: np.asarray(source[key]).copy() for key in source.files}
    with np.load(trace_path, allow_pickle=False) as source:
        trace = {key: np.asarray(source[key]).copy() for key in source.files}
    dx = float(record["configuration"]["cell_size_m"])
    maximum = np.asarray(field["maximum_concentration_vol_percent"], dtype=float)
    fluid = ~np.asarray(field["solid_mask"], dtype=bool)
    hazard = fluid & (maximum >= threshold)
    traces = trace_metrics(
        trace["time_s"], trace["concentration_vol_percent"], threshold=0.1,
    )
    regions = obstacle_region_metrics(field, cell_size_m=dx)
    numerics = record.get("numerics", {})
    def _finite_range(name: str) -> dict[str, float] | None:
        values = np.asarray(field[name], dtype=float)[fluid]
        values = values[np.isfinite(values)]
        if values.size == 0:
            return None
        return {"min": float(np.min(values)), "max": float(np.max(values))}

    physics_summary = {
        "state_fields_present": [
            name for name in (
                "pressure_pa_gauge", "mixture_density_kg_m3", "temperature_k",
                "velocity_u_m_s", "velocity_v_m_s", "velocity_w_m_s",
                "concentration_kg_m3", "concentration_vol_percent",
            ) if name in field
        ],
        "pressure_pa_gauge_range": _finite_range("pressure_pa_gauge"),
        "mixture_density_kg_m3_range": _finite_range("mixture_density_kg_m3"),
        "temperature_k_range": _finite_range("temperature_k"),
        "velocity_u_m_s_range": _finite_range("velocity_u_m_s"),
        "velocity_v_m_s_range": _finite_range("velocity_v_m_s"),
        "velocity_w_m_s_range": _finite_range("velocity_w_m_s"),
        "conservation_residuals": {
            "species_mass_relative": numerics.get("maximum_relative_species_mass_residual"),
            "sensible_enthalpy_relative": numerics.get("maximum_relative_sensible_enthalpy_residual"),
            "continuity_kg_m3_s": numerics.get("maximum_continuity_residual_kg_m3_s"),
            "momentum_kg_m_s": numerics.get("maximum_momentum_balance_residual_kg_m_s"),
        },
    }
    density_range = physics_summary["mixture_density_kg_m3_range"]
    temperature_range = physics_summary["temperature_k_range"]
    density_relative_span = (
        None if density_range is None or density_range["max"] == 0.0 else
        float((density_range["max"] - density_range["min"]) / density_range["max"])
    )
    temperature_span = (
        None if temperature_range is None else
        float(temperature_range["max"] - temperature_range["min"])
    )
    physics_summary["thermal_buoyancy_transition"] = {
        "density_relative_span": density_relative_span,
        "temperature_span_k": temperature_span,
        "detectable_in_this_run": bool(
            (density_relative_span is not None and density_relative_span > 1.0e-4)
            or (temperature_span is not None and temperature_span > 0.01)
        ),
        "status": "exercised" if (
            (density_relative_span is not None and density_relative_span > 1.0e-4)
            or (temperature_span is not None and temperature_span > 0.01)
        ) else "not_exercised_in_this_screening_run",
    }
    score_summary: dict[str, object] | None = None
    if score is not None:
        score_path = score.resolve()
        candidates = sorted(score_path.glob("*.json")) if score_path.is_dir() else [score_path]
        if not candidates:
            raise FileNotFoundError(f"no score JSON found under {score_path}")
        score_record = json.loads(candidates[0].read_text(encoding="utf-8"))
        dose = score_record.get("dose_sensor_score", {})
        arrival = score_record.get("arrival_score", {})
        departure = score_record.get("departure_score", {})
        arc = score_record.get("six_arc_score", {})
        score_summary = {
            "source": str(candidates[0]),
            "status": score_record.get("status"),
            "horizon_limited": score_record.get("horizon_limited"),
            "concentration_unit": score_record.get("protocol", {}).get("concentration_unit", "vol%"),
            "dose_unit": score_record.get("protocol", {}).get("dose_unit", "vol% s"),
            "dose_mae_vol_percent_s": dose.get("mean_absolute_error_percent_s"),
            "dose_mape": dose.get("mean_absolute_percentage_error"),
            "dose_normalized_mean_bias": dose.get("normalized_mean_bias"),
            "dose_pearson_r": dose.get("pearson_r"),
            "dose_spearman_rho": dose.get("spearman_rho"),
            "arrival_mae_s": arrival.get("mean_absolute_error_s"),
            "departure_mae_s": departure.get("mean_absolute_error_s"),
            "six_arc_peak_mape": arc.get("peak_mean_absolute_percentage_error"),
            "six_arc_dose_mape": arc.get("dose_mean_absolute_percentage_error"),
            "promotion_allowed": False,
            "interpretation": "post-run comparison only; observation values were not accessed during configuration",
        }
    output.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    top = np.max(maximum, axis=0)
    front = np.max(maximum, axis=1)
    side = np.max(maximum, axis=2)
    region_masks = obstacle_region_masks(
        field["x_m"], field["y_m"], field["z_m"], field["solid_mask"],
    )
    wake_cells = region_masks["exterior_lee_wake_shell"]
    wake = np.max(np.where(wake_cells, maximum, -np.inf), axis=0)
    wake[~np.isfinite(wake)] = np.nan
    top_solid = np.any(field["solid_mask"], axis=0)
    front_solid = np.any(field["solid_mask"], axis=1)
    side_solid = np.any(field["solid_mask"], axis=2)
    levels = np.asarray([0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32])
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_under("#f8fafc")
    norm = LogNorm(vmin=0.1, vmax=max(32.0, float(np.nanmax(maximum))))
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0), constrained_layout=True)
    axes = axes.ravel()
    panels = (
        (top, top_solid, field["x_m"], field["y_m"], "Plan / top (max over z)", "x (m)", "y (m)"),
        (front, front_solid, field["x_m"], field["z_m"], "Front (max over y)", "x (m)", "z (m)"),
        (side, side_solid, field["y_m"], field["z_m"], "Side (max over x)", "y (m)", "z (m)"),
        (wake, np.zeros_like(wake, dtype=bool), field["x_m"], field["y_m"], "Lee wake shell (max over z)", "x (m)", "y (m)"),
    )
    image = None
    for axis, (values, solid, xx, yy, title, xlabel, ylabel) in zip(axes, panels):
        image = axis.contourf(xx, yy, values, levels=levels, cmap=cmap, norm=norm, extend="max")
        if float(np.max(values)) >= threshold:
            axis.contour(xx, yy, values, levels=[threshold], colors="white", linewidths=3.0)
            axis.contour(xx, yy, values, levels=[threshold], colors="#111827", linewidths=1.1)
        if np.any(solid):
            axis.contour(xx, yy, solid.astype(float), levels=[0.5], colors="#111827", linewidths=1.2)
        axis.set(title=title, xlabel=xlabel, ylabel=ylabel)
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.15)
    fig.suptitle("DEGALI TI021 full-horizon obstacle screening field\nTime-maximum H₂ concentration (vol%)", weight="bold", fontsize=15)
    cbar = fig.colorbar(image, ax=axes, shrink=0.86, pad=0.02)
    cbar.set_label("H₂ concentration (vol%)")
    fig.text(0.5, 0.01, f"White/black contour = {threshold:g} vol%; hatched/outlined cells = resolved obstacle", ha="center", fontsize=9)
    figure_path = output / "ti021-600s-screening-front-top-side.png"
    fig.savefig(figure_path, dpi=220)
    fig.savefig(figure_path.with_suffix(".svg"))
    plt.close(fig)

    result = {
        "schema": "degali.thorney-ti021-screening-package.v1",
        "classification": "full-horizon 2.5 m obstacle-resolved screening tier; finest-grid validation withheld",
        "run": {
            "path": str(run),
            "duration_s": float(record["configuration"]["duration_s"]),
            "cell_size_m": dx,
            "time_step_s": float(record["configuration"]["time_step_s"]),
            "wall_clock_s": float(record["runtime"]["wall_clock_s"]),
            "solver_reported_total_s": float(record["runtime"]["solver_reported_total_s"]),
            "all_gates_pass": bool(record["all_gates_pass"]),
        },
        "concentration_unit": "vol%",
        "physics": physics_summary,
        "hazard": {
            "threshold_vol_percent": threshold,
            "maximum_ever_volume_m3": float(np.count_nonzero(hazard) * dx**3),
            "maximum_ever_peak_vol_percent": float(np.max(maximum[fluid])),
            "final_field_volume_m3": float(np.count_nonzero(fluid & (field["concentration_vol_percent"] >= threshold)) * dx**3),
        },
        "obstacle_regions": regions,
        "receptors": {
            "count": int(trace["concentration_vol_percent"].shape[0]),
            "time_s_start": float(trace["time_s"][0]),
            "time_s_end": float(trace["time_s"][-1]),
            "peak_vol_percent_max": float(np.max(traces["peak_vol_percent"])),
            "dose_vol_percent_s_max": float(np.max(traces["dose_vol_percent_s"])),
            "detected_count_at_0p1_vol_percent": int(np.count_nonzero(np.isfinite(traces["arrival_s"]))),
        },
        "independent_error_metrics": score_summary,
        "files": [figure_path.name, figure_path.with_suffix(".svg").name],
        "promotion_allowed": False,
        "promotion_boundary": "1.25 m full-horizon convergence and independent calibrated obstacle-H2 observations remain required",
    }
    (output / "ti021-600s-screening-package.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--score", type=Path, default=None)
    args = parser.parse_args()
    print(json.dumps(package(args.run, args.output, score=args.score), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
