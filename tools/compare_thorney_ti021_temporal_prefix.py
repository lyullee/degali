"""Compare two materialized TI021 checkpoint prefixes on a common time axis."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finite_difference_metrics(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(a) & np.isfinite(b)
    if not np.any(mask):
        raise ValueError("comparison arrays have no common finite values")
    d = np.asarray(a[mask] - b[mask], dtype=float)
    return {
        "finite_count": int(d.size),
        "bias": float(np.mean(d)),
        "mae": float(np.mean(np.abs(d))),
        "rmse": float(np.sqrt(np.mean(d * d))),
        "maximum_absolute_difference": float(np.max(np.abs(d))),
    }


def _common_fluid_component_metrics(
    reference: dict[str, np.ndarray], comparison: dict[str, np.ndarray],
) -> dict[str, object]:
    """Compare final scalar/flow fields on cells fluid in both prefixes."""
    coordinate_names = ("x_m", "y_m", "z_m")
    for name in coordinate_names:
        if name in reference and name in comparison and not np.array_equal(
            reference[name], comparison[name],
        ):
            raise ValueError(f"grid coordinate {name} does not match")
    sample = next(
        (value for key, value in reference.items()
         if key not in coordinate_names and np.asarray(value).ndim == 3),
        None,
    )
    if sample is None:
        return {"compared_components": {}, "common_fluid_cell_count": 0}
    shape = np.asarray(sample).shape
    reference_solid = np.asarray(
        reference.get("solid_mask", np.zeros(shape, dtype=bool)), dtype=bool,
    )
    comparison_solid = np.asarray(
        comparison.get("solid_mask", np.zeros(shape, dtype=bool)), dtype=bool,
    )
    if reference_solid.shape != shape or comparison_solid.shape != shape:
        raise ValueError("solid masks must match the 3-D field shape")
    common_fluid = ~(reference_solid | comparison_solid)
    names = (
        "concentration_vol_percent", "maximum_concentration_vol_percent",
        "temperature_k", "mixture_density_kg_m3", "velocity_u_m_s",
        "velocity_v_m_s", "velocity_w_m_s", "pressure_pa_gauge",
    )
    metrics: dict[str, object] = {}
    for name in names:
        if name not in reference or name not in comparison:
            continue
        left = np.asarray(reference[name], dtype=float)
        right = np.asarray(comparison[name], dtype=float)
        if left.shape != shape or right.shape != shape:
            raise ValueError(f"physical field {name} does not match the grid shape")
        finite = common_fluid & np.isfinite(left) & np.isfinite(right)
        if not np.any(finite):
            continue
        difference = left[finite] - right[finite]
        metrics[name] = {
            "finite_common_fluid_cell_count": int(difference.size),
            "bias": float(np.mean(difference)),
            "mae": float(np.mean(np.abs(difference))),
            "rmse": float(np.sqrt(np.mean(difference * difference))),
            "maximum_absolute_difference": float(np.max(np.abs(difference))),
        }
    region_metrics: dict[str, object] = {}
    if all(name in reference and name in comparison for name in coordinate_names):
        x = np.asarray(reference["x_m"], dtype=float)
        y = np.asarray(reference["y_m"], dtype=float)
        z = np.asarray(reference["z_m"], dtype=float)
        xx = x[None, None, :]
        yy = y[None, :, None]
        zz = z[:, None, None]
        if (len(z), len(y), len(x)) != shape:
            raise ValueError("coordinate lengths do not match 3-D field shape")
        radius = np.hypot(xx - 60.0, yy - 100.0)
        tip_distance = np.minimum(
            np.hypot(xx - 60.0, yy - 50.0),
            np.hypot(xx - 60.0, yy - 150.0),
        )
        regions = {
            "windward_inner_face_shell": (
                (xx >= 60.0) & (radius >= 40.0) & (radius < 50.0) & (zz < 5.0)
            ),
            "roof_shear_band": (
                (xx >= 60.0) & (radius >= 45.0) & (radius <= 55.0)
                & (zz >= 5.0) & (zz < 10.0)
            ),
            "side_tip_bypass": (tip_distance <= 10.0) & (zz < 5.0),
            "exterior_lee_wake_shell": (
                (xx >= 60.0) & (radius > 50.0) & (radius <= 70.0)
                & (zz < 7.5)
            ),
        }
        for region_name, region in regions.items():
            eligible = region & common_fluid
            by_component: dict[str, object] = {}
            for name in names:
                if name not in reference or name not in comparison:
                    continue
                left = np.asarray(reference[name], dtype=float)
                right = np.asarray(comparison[name], dtype=float)
                finite = eligible & np.isfinite(left) & np.isfinite(right)
                if not np.any(finite):
                    continue
                difference = left[finite] - right[finite]
                by_component[name] = {
                    "finite_cell_count": int(difference.size),
                    "bias": float(np.mean(difference)),
                    "mae": float(np.mean(np.abs(difference))),
                    "rmse": float(np.sqrt(np.mean(difference * difference))),
                    "maximum_absolute_difference": float(np.max(np.abs(difference))),
                }
            region_metrics[region_name] = {
                "common_fluid_cell_count": int(np.count_nonzero(eligible)),
                "components": by_component,
            }
    return {
        "common_fluid_cell_count": int(np.count_nonzero(common_fluid)),
        "solid_mask_disagreement_cell_count": int(
            np.count_nonzero(reference_solid != comparison_solid)
        ),
        "compared_components": metrics,
        "obstacle_region_components": region_metrics,
    }


def compare_temporal_prefixes(
    reference_run: Path,
    comparison_run: Path,
    output_json: Path,
) -> dict[str, object]:
    reference_run = reference_run.resolve()
    comparison_run = comparison_run.resolve()
    with np.load(
        reference_run / "ti021-inner-final-field.npz", allow_pickle=False,
    ) as reference_field, np.load(
        comparison_run / "ti021-inner-final-field.npz", allow_pickle=False,
    ) as comparison_field:
        reference_components = {
            key: np.asarray(reference_field[key]).copy()
            for key in reference_field.files
        }
        comparison_components = {
            key: np.asarray(comparison_field[key]).copy()
            for key in comparison_field.files
        }
        reference_maximum = np.asarray(
            reference_field["maximum_concentration_vol_percent"], dtype=float,
        )
        comparison_maximum = np.asarray(
            comparison_field["maximum_concentration_vol_percent"], dtype=float,
        )
        if reference_maximum.shape != comparison_maximum.shape:
            raise ValueError("field shapes do not match")
        field_metrics = _finite_difference_metrics(
            reference_maximum, comparison_maximum,
        )
        cell_volume_m3 = float(
            reference_field["x_m"].size
            and (reference_field["x_m"][1] - reference_field["x_m"][0]) ** 3
        )
        if not np.isfinite(cell_volume_m3) or cell_volume_m3 <= 0.0:
            raise ValueError("could not infer a positive field cell volume")
        reference_hazard_volume = float(
            np.count_nonzero(reference_maximum >= 4.0) * cell_volume_m3
        )
        comparison_hazard_volume = float(
            np.count_nonzero(comparison_maximum >= 4.0) * cell_volume_m3
        )
        field_summary = {
            **field_metrics,
            "reference_maximum_vol_percent": float(np.nanmax(reference_maximum)),
            "comparison_maximum_vol_percent": float(np.nanmax(comparison_maximum)),
            "reference_hazard_volume_m3": reference_hazard_volume,
            "comparison_hazard_volume_m3": comparison_hazard_volume,
            "hazard_volume_difference_m3": (
                reference_hazard_volume - comparison_hazard_volume
            ),
        }

    with np.load(
        reference_run / "ti021-inner-receptor-traces.npz", allow_pickle=False,
    ) as reference_receptors, np.load(
        comparison_run / "ti021-inner-receptor-traces.npz", allow_pickle=False,
    ) as comparison_receptors:
        reference_time = np.asarray(reference_receptors["time_s"], dtype=float)
        comparison_time = np.asarray(comparison_receptors["time_s"], dtype=float)
        if not np.array_equal(
            reference_receptors["labels"], comparison_receptors["labels"],
        ):
            raise ValueError("receptor labels do not match")
        comparison_values = np.asarray(
            comparison_receptors["concentration_vol_percent"], dtype=float,
        )
        interpolated = np.vstack([
            np.interp(reference_time, comparison_time, row)
            for row in comparison_values
        ])
        receptor_values = np.asarray(
            reference_receptors["concentration_vol_percent"], dtype=float,
        )
        receptor_summary = {
            **_finite_difference_metrics(receptor_values, interpolated),
            "common_frame_count": int(reference_time.size),
            "common_stop_time_s": float(reference_time[-1]),
            "reference_peak_vol_percent": float(np.nanmax(receptor_values)),
            "comparison_peak_vol_percent_interpolated": float(
                np.nanmax(interpolated)
            ),
        }

    record = {
        "schema": "degali.thorney-ti021-temporal-prefix-comparison.v1",
        "status": "computed_window_temporal_comparison",
        "reference_run": str(reference_run),
        "comparison_run": str(comparison_run),
        "provenance": {
            "reference_record_sha256": _sha256(
                reference_run / "ti021-inner-transient.json"
            ),
            "reference_field_sha256": _sha256(
                reference_run / "ti021-inner-final-field.npz"
            ),
            "reference_trace_sha256": _sha256(
                reference_run / "ti021-inner-receptor-traces.npz"
            ),
            "comparison_record_sha256": _sha256(
                comparison_run / "ti021-inner-transient.json"
            ),
            "comparison_field_sha256": _sha256(
                comparison_run / "ti021-inner-final-field.npz"
            ),
            "comparison_trace_sha256": _sha256(
                comparison_run / "ti021-inner-receptor-traces.npz"
            ),
        },
        "field": field_summary,
        "physical_field_components": _common_fluid_component_metrics(
            reference_components, comparison_components,
        ),
        "receptors": receptor_summary,
        "interpretation": {
            "observation_values_used": False,
            "full_event_validation": False,
            "hazard_volume_convergence_claim": False,
            "receptor_time_axis_interpolation": True,
        },
    }
    output_json = output_json.resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-run", type=Path, required=True)
    parser.add_argument("--comparison-run", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    record = compare_temporal_prefixes(
        args.reference_run, args.comparison_run, args.output_json,
    )
    print(json.dumps(record, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
