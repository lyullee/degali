"""Conditionally replay a public PRESLHY E3.4 H2 obstacle trial.

The source window is frozen from corrected liquid mass only and the wind speed
is read from the synchronized sonic channel. Hydrogen observations are loaded
only after the solver configuration is fixed, then scored against the three
usable sensor channels. The report geometry is still approximate, so the
result is explicitly withheld from facility promotion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from degali.addons.fast_obstacle_transport_3d import (
    FastObstacleReceptor3D,
    FastObstacleTransport3DConfig,
    solve_fast_obstacle_transport_3d,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.addons.site_geometry import AxisAlignedCuboid
from degali.validation.preslhy_e34 import (
    read_e34_dynamic_wind_hydrogen_history,
    read_e34_pool_history,
    read_e34_wind_history,
    source_only_evaporation_candidates,
)
try:
    from tools.audit_preslhy_e34_dynamic_evidence import audit_dynamic_trial
except ModuleNotFoundError:  # direct ``python tools/<script>.py`` execution
    from audit_preslhy_e34_dynamic_evidence import audit_dynamic_trial


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKBOOK = Path(
    r"C:\Users\lyul\Desktop\DocuData\SLABx_LH2\slabx-lh2\액화수소 검증\10.35097-1319\10.35097-1319\data\dataset\20200416-Sand03Wind-Final.xlsx"
)
DEFAULT_REPORT = Path(
    r"C:\Users\lyul\Desktop\DocuData\SLABx_LH2\slabx-lh2\액화수소 검증\10.35097-1319\10.35097-1319\data\dataset\Draft-20210727-PRESLHY_D3_5_E3_4-UNIgnitedPool.docx"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - observed
    above_one = np.abs(observed) >= 1.0
    if np.any(above_one):
        threshold_mae = float(np.mean(np.abs(residual[above_one])))
        threshold_rmse = float(np.sqrt(np.mean(residual[above_one] ** 2)))
        threshold_mape = float(
            np.mean(np.abs(residual[above_one]) / np.abs(observed[above_one]))
            * 100.0
        )
    else:
        threshold_mae = threshold_rmse = threshold_mape = float("nan")
    return {
        "pair_count": int(observed.size),
        "mae_vol_percent_point": float(np.mean(np.abs(residual))),
        "rmse_vol_percent_point": float(np.sqrt(np.mean(residual**2))),
        "bias_vol_percent_point": float(np.mean(residual)),
        "mape_percent_floor_0p1_vol_percent": float(
            np.mean(np.abs(residual) / np.maximum(np.abs(observed), 0.1)) * 100.0
        ),
        "mape_percent_floor_1_vol_percent": float(
            np.mean(np.abs(residual) / np.maximum(np.abs(observed), 1.0)) * 100.0
        ),
        "above_1_vol_percent_pair_count": int(np.count_nonzero(above_one)),
        "above_1_vol_percent_mae": threshold_mae,
        "above_1_vol_percent_rmse": threshold_rmse,
        "above_1_vol_percent_mape": threshold_mape,
    }


def run_replay(
    *, workbook: Path, report: Path, output: Path, window_ordinal: int = 2,
) -> dict[str, object]:
    workbook = workbook.resolve(); report = report.resolve()
    if not workbook.is_file() or not report.is_file():
        raise FileNotFoundError("PRESLHY workbook/report is not available")
    pool = read_e34_pool_history(workbook)
    candidates = source_only_evaporation_candidates(pool)
    if window_ordinal < 1 or window_ordinal > len(candidates):
        raise ValueError("window_ordinal is outside the frozen source-window set")
    window = candidates[window_ordinal - 1]
    # Wind is configured from the synchronized sonic channel only. H2 values
    # are intentionally not loaded until after the solve configuration exists.
    wind = read_e34_wind_history(workbook)
    wind_mask = (
        (wind.synchronized_time_s >= window.start_s)
        & (wind.synchronized_time_s <= window.end_s)
    )
    if not np.any(wind_mask):
        raise ValueError("frozen source window has no synchronized wind samples")
    mean_u = float(np.mean(wind.wind_u_m_s[wind_mask]))
    mean_v = float(np.mean(wind.wind_v_m_s[wind_mask]))
    mean_w = float(np.mean(wind.wind_w_m_s[wind_mask]))
    # The fast flow currently advances along its x-axis.  Use the measured
    # horizontal sonic magnitude for that axis; feeding the vertical sonic
    # component into an x-only speed would spuriously accelerate advection.
    wind_speed = max(float(np.sqrt(mean_u**2 + mean_v**2)), 0.1)
    pool_width_m = 0.50
    turbulent_diffusivity_m2_s = max(
        0.000061, 0.01 * wind_speed * pool_width_m,
    )
    duration_s = float(window.duration_s)
    # The report's obstacle layout is not machine-registered. This declared
    # surrogate is used only to exercise the obstacle solver path and is never
    # treated as the measured geometry.
    obstacle = AxisAlignedCuboid(
        1.00, 1.40, 0.35, 0.95, 0.00, 0.50,
        label="PRESLHY report-figure surrogate block (not registered)",
    )
    source_schedule = SourceRateSchedule(
        time_s=(0.0, duration_s),
        rate_kg_s=(window.fitted_evaporation_rate_kg_s, 0.0),
        source_id=f"{window.window_id}:mass-only source-rate proxy",
    )
    # The default is the fast 20 ms step used by the screening path.  Around
    # a stronger surrogate obstacle the potential-flow acceleration can make
    # that step violate the scalar CFL bound; halve it deterministically rather
    # than silently clipping the solution or failing the whole replay.
    time_step_s = 0.02
    receptors = tuple(
        FastObstacleReceptor3D(label, (0.75, 0.80, height))
        for label, height in (("CH35", 0.35), ("CH45", 0.45), ("CH55", 0.55))
    )
    started = time.perf_counter()
    for _attempt in range(4):
        config = FastObstacleTransport3DConfig(
            domain_m=(2.0, 1.6, 1.0), cells=(20, 16, 10),
            duration_s=duration_s, time_step_s=time_step_s,
            wind_speed_m_s=wind_speed,
            diffusivity_m2_s=turbulent_diffusivity_m2_s,
            source_schedule=source_schedule,
            # The report describes a nominal 0.50 m x 0.50 m pool.  Use its
            # centre and a geometry-derived, anisotropic Gaussian footprint
            # rather than collapsing the entire liquid source into a 5 cm
            # sphere.
            source_position_m=(0.75, 0.80, 0.05), source_sigma_m=0.05,
            # The measured pool footprint fixes the horizontal sigma.  A
            # 0.20 m near-field release/mixing layer is retained as an
            # explicit conditional surrogate because the workbook has no
            # resolved plume height; it is not a calibrated facility value.
            source_sigma_xyz_m=(
                0.50 / np.sqrt(12.0), 0.50 / np.sqrt(12.0), 0.10,
            ),
            obstacles=(obstacle,), store_fields=False,
        )
        try:
            result = solve_fast_obstacle_transport_3d(config, receptors=receptors)
            break
        except ValueError as exc:
            if "stability limit exceeded" not in str(exc) or _attempt == 3:
                raise
            time_step_s *= 0.5
    runtime_s = time.perf_counter() - started

    # Only now load the observed H2 values and apply the frozen measurement
    # operator to the model receptor traces.
    dynamic = read_e34_dynamic_wind_hydrogen_history(workbook)
    evidence_audit = audit_dynamic_trial(workbook, report)
    score_time = result.time_s + window.start_s
    paired: dict[str, dict[str, float]] = {}
    all_observed: list[float] = []; all_predicted: list[float] = []
    for receptor in receptors:
        label = receptor.label
        usable = dynamic.sensor_quality.get(label) == "usable"
        if not usable:
            continue
        observed = np.interp(
            score_time, dynamic.hydrogen_sample_time_s,
            dynamic.hydrogen_vol_percent[label],
        )
        predicted = result.receptor_traces_kg_m3[label] / 0.00201588 / (
            101325.0 / (8.31446261815324 * 293.15)
        ) * 100.0
        metrics = _metrics(observed, predicted)
        paired[label] = metrics
        all_observed.extend(observed.tolist()); all_predicted.extend(predicted.tolist())
    overall = _metrics(np.asarray(all_observed), np.asarray(all_predicted)) if all_observed else None
    vol = result.maximum_concentration_vol_percent
    fluid = np.isfinite(vol) & (~result.solid_mask)
    threshold = fluid & (vol >= 4.0)
    cell_volume = float(np.prod(result.diagnostics["cell_size_m"]))
    trial_id = workbook.stem.replace("-Final", "")
    manifest_path = (
        ROOT / "outputs" / "preslhy-e34-dynamic-evidence-2026-10-11-v2"
        / f"{trial_id}-field-evidence-manifest.json"
    )
    output_record = {
        "schema": "degali.preslhy-e34-h2-fast-obstacle-replay.v1",
        "status": "withheld_conditional_h2_obstacle_replay",
        "quantitative_prediction_allowed": False,
        "promotion_allowed": False,
        "validation_gate": {
            "status": "withheld",
            "reasons": [
                "report obstacle geometry is approximate and not machine-registered",
                "liquid-mass loss is only a source-rate proxy, not atmospheric vapour flux",
                "sensor calibration certificate is absent (specification_only)",
            ],
            "field_evidence_decision": evidence_audit["evidence_decision"],
        },
        "input": {
            "workbook": str(workbook), "workbook_sha256": _sha256(workbook),
            "report": str(report), "report_sha256": _sha256(report),
            "manifest_candidate": str(manifest_path),
            "manifest_id": evidence_audit["field_evidence_manifest"]["manifest_id"],
            "event_id": workbook.stem.replace("-Final", ""),
            "source_window": window.__dict__,
            "wind_mean_m_s": [mean_u, mean_v, mean_w],
            "wind_speed_m_s": wind_speed,
            "wind_speed_definition": "horizontal sqrt(mean_u^2 + mean_v^2); x-axis fast-flow projection",
            "turbulent_diffusivity_m2_s": turbulent_diffusivity_m2_s,
            "turbulent_diffusivity_rule": "max(molecular H2 diffusivity, 0.01 * horizontal wind speed * 0.50 m pool width)",
            "source_near_field_mixing_layer_m": 0.20,
            "time_step_selection": "start at 0.02 s and halve only on the explicit scalar stability gate",
            "accepted_time_step_s": time_step_s,
        },
        "runtime_s": runtime_s,
        "duration_s": duration_s,
        "cells": list(config.cells),
        "cell_size_m": list(result.diagnostics["cell_size_m"]),
        "source_rate_kg_s": window.fitted_evaporation_rate_kg_s,
        "maximum_h2_vol_percent": float(np.nanmax(vol[fluid])),
        "four_vol_percent_volume_m3": float(np.count_nonzero(threshold) * cell_volume),
        "mass_conservation": {
            "maximum_absolute_residual_kg": result.maximum_mass_residual_kg,
            "passed": result.maximum_mass_residual_kg <= 1.0e-9,
        },
        "sensor_scores": paired,
        "overall_score": overall,
        "sensor_quality": dynamic.sensor_quality,
        "physics": result.diagnostics,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(output_record, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    np.savez_compressed(
        output.with_name(output.stem + "-maximum-field.npz"),
        x_m=result.x_m, y_m=result.y_m, z_m=result.z_m,
        maximum_h2_vol_percent=vol, solid_mask=result.solid_mask,
        threshold_mask=threshold,
    )
    return output_record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--window-ordinal", type=int, default=2)
    args = parser.parse_args()
    print(json.dumps(run_replay(
        workbook=args.workbook, report=args.report, output=args.output,
        window_ordinal=args.window_ordinal,
    ), indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
