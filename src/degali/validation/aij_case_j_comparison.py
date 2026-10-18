"""Pointwise AIJ Case J scoring with explicit coordinate and tracer mapping.

The comparison routine does not make a solver run valid by itself. The caller
must provide an independently justified experiment-to-model translation and
the model concentration must represent ethylene mass concentration (kg/m³).
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TYPE_CHECKING

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .aij_case_j import AijCaseJDataset

if TYPE_CHECKING:
    from degali.addons.resolved_obstacle_flow_3d import ResolvedObstacleFlow3DResult


ETHYLENE_MOLAR_MASS_KG_MOL = 28.05316e-3
_R_GAS_J_MOL_K = 8.31446261815324


@dataclass(frozen=True)
class AijCaseJPointComparison:
    point_id: int
    x_experiment_m: float
    y_experiment_m: float
    z_experiment_m: float
    measured_u_m_s: float
    predicted_u_m_s: float
    measured_v_m_s: float
    predicted_v_m_s: float
    measured_w_m_s: float
    predicted_w_m_s: float
    measured_temperature_k: float
    predicted_temperature_k: float
    measured_ethylene_vol_percent: float
    predicted_ethylene_vol_percent: float


@dataclass(frozen=True)
class AijCaseJComparison:
    points: tuple[AijCaseJPointComparison, ...]
    skipped_point_ids: tuple[int, ...]
    skip_reasons: dict[str, int]
    metrics: dict[str, float | int | None]
    coordinate_translation_model_minus_experiment_m: tuple[float, float, float]
    tracer_molar_mass_kg_mol: float
    concentration_mapping: str
    model_statistic: str
    status: str


def _metric(measured: np.ndarray, predicted: np.ndarray, prefix: str) -> dict[str, float | None]:
    difference = predicted - measured
    correlation: float | None = None
    if len(measured) >= 2 and np.std(measured) > 0.0 and np.std(predicted) > 0.0:
        correlation = float(np.corrcoef(measured, predicted)[0, 1])
    return {
        f"{prefix}_bias": float(np.mean(difference)),
        f"{prefix}_mae": float(np.mean(np.abs(difference))),
        f"{prefix}_rmse": float(np.sqrt(np.mean(difference**2))),
        f"{prefix}_correlation": correlation,
    }


def compare_aij_case_j(
    dataset: AijCaseJDataset,
    result: ResolvedObstacleFlow3DResult,
    *,
    coordinate_translation_model_minus_experiment_m: tuple[float, float, float],
    tracer_molar_mass_kg_mol: float = ETHYLENE_MOLAR_MASS_KG_MOL,
    use_transient_average: bool = True,
) -> AijCaseJComparison:
    """Interpolate model fields at RS points and score in physical units.

    ``result.concentration_kg_m3`` is deliberately used instead of the
    hydrogen-specific vol% property. The caller must ensure that this field
    was transported as passive ethylene mass concentration; hydrogen EOS
    feedback or hydrogen-based vol% conversion is not accepted as ethylene
    validation evidence.
    """
    translation = tuple(float(value) for value in coordinate_translation_model_minus_experiment_m)
    if len(translation) != 3 or not all(math.isfinite(value) for value in translation):
        raise ValueError("coordinate translation must contain three finite values")
    molar_mass = float(tracer_molar_mass_kg_mol)
    if not math.isfinite(molar_mass) or molar_mass <= 0.0:
        raise ValueError("tracer molar mass must be finite and positive")
    axes = tuple(np.asarray(axis, dtype=float) for axis in (result.z_m, result.y_m, result.x_m))
    for axis in axes:
        if axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis)) or np.any(np.diff(axis) <= 0.0):
            raise ValueError("model coordinate axes must be finite, increasing vectors")
    shape = (len(axes[0]), len(axes[1]), len(axes[2]))

    def interpolator(field: np.ndarray, *, method: str = "linear"):
        values = np.asarray(field, dtype=float)
        if values.shape != shape:
            raise ValueError(f"model field shape {values.shape} does not match grid {shape}")
        return RegularGridInterpolator(
            axes, values, method=method, bounds_error=False, fill_value=None,
        )

    average_count = int(getattr(result, "time_average_sample_count", 0))
    use_average = use_transient_average and average_count > 0
    velocity_u = (
        result.time_averaged_velocity_u_m_s
        if use_average else result.velocity_u_m_s
    )
    velocity_v = (
        result.time_averaged_velocity_v_m_s
        if use_average else result.velocity_v_m_s
    )
    velocity_w = (
        result.time_averaged_velocity_w_m_s
        if use_average else result.velocity_w_m_s
    )
    temperature_field = (
        result.time_averaged_temperature_k
        if use_average else result.temperature_k
    )
    concentration_field = (
        result.time_averaged_concentration_kg_m3
        if use_average else result.concentration_kg_m3
    )
    field_map = {
        "u": interpolator(velocity_u),
        "v": interpolator(velocity_v),
        "w": interpolator(velocity_w),
        "temperature": interpolator(
            np.full(shape, result.ambient_temperature_k)
            if temperature_field is None else temperature_field
        ),
        "concentration": interpolator(concentration_field),
        "solid": interpolator(np.asarray(result.solid_mask, dtype=float)),
    }
    samples: list[AijCaseJPointComparison] = []
    skipped: list[int] = []
    skip_reasons = {"outside_physical_domain": 0, "solid_interpolation_stencil": 0, "nonfinite_model_field": 0}
    half_cell = np.asarray(tuple((axis[1] - axis[0]) * 0.5 for axis in axes))
    for point in dataset.observations:
        model_xyz = np.asarray((point.x_m, point.y_m, point.z_m)) + np.asarray(translation)
        sample_zyx = np.asarray((model_xyz[2], model_xyz[1], model_xyz[0]))
        if any(
            sample < axis[0] - margin - 1.0e-12
            or sample > axis[-1] + margin + 1.0e-12
            for sample, axis, margin in zip(sample_zyx, axes, half_cell)
        ):
            skipped.append(point.point_id)
            skip_reasons["outside_physical_domain"] += 1
            continue
        values = {
            name: float(np.asarray(field(sample_zyx[None, :])).reshape(-1)[0])
            for name, field in field_map.items()
        }
        # Do not score a point when any part of its interpolation stencil
        # overlaps a solid or lies beyond the model domain.
        if not all(math.isfinite(value) for value in values.values()):
            skipped.append(point.point_id)
            skip_reasons["nonfinite_model_field"] += 1
            continue
        if values["solid"] > 1.0e-12:
            skipped.append(point.point_id)
            skip_reasons["solid_interpolation_stencil"] += 1
            continue
        local_temperature = values["temperature"]
        total_molar_density = result.ambient_pressure_pa / (_R_GAS_J_MOL_K * local_temperature)
        predicted_vol_percent = (
            100.0 * values["concentration"] / molar_mass / total_molar_density
        )
        samples.append(AijCaseJPointComparison(
            point_id=point.point_id,
            x_experiment_m=point.x_m,
            y_experiment_m=point.y_m,
            z_experiment_m=point.z_m,
            measured_u_m_s=point.u_m_s,
            predicted_u_m_s=values["u"],
            measured_v_m_s=point.v_m_s,
            predicted_v_m_s=values["v"],
            measured_w_m_s=point.w_m_s,
            predicted_w_m_s=values["w"],
            measured_temperature_k=point.temperature_k,
            predicted_temperature_k=local_temperature,
            measured_ethylene_vol_percent=point.concentration_vol_percent,
            predicted_ethylene_vol_percent=predicted_vol_percent,
        ))

    def arrays(measured_name: str, predicted_name: str) -> tuple[np.ndarray, np.ndarray]:
        return (
            np.asarray([getattr(item, measured_name) for item in samples], dtype=float),
            np.asarray([getattr(item, predicted_name) for item in samples], dtype=float),
        )

    metrics: dict[str, float | int | None] = {"point_count": len(samples)}
    for prefix, measured_name, predicted_name in (
        ("u_m_s", "measured_u_m_s", "predicted_u_m_s"),
        ("v_m_s", "measured_v_m_s", "predicted_v_m_s"),
        ("w_m_s", "measured_w_m_s", "predicted_w_m_s"),
        ("temperature_k", "measured_temperature_k", "predicted_temperature_k"),
        ("ethylene_vol_percent", "measured_ethylene_vol_percent", "predicted_ethylene_vol_percent"),
    ):
        metrics.update(_metric(*arrays(measured_name, predicted_name), prefix))
    return AijCaseJComparison(
        points=tuple(samples),
        skipped_point_ids=tuple(skipped),
        skip_reasons=skip_reasons,
        metrics=metrics,
        coordinate_translation_model_minus_experiment_m=translation,
        tracer_molar_mass_kg_mol=molar_mass,
        concentration_mapping="passive_species_mass_concentration_to_ideal_gas_mole_percent_at_local_temperature",
        model_statistic=(
            "finite_transient_time_mean"
            if use_average else "single_transient_final_snapshot"
        ),
        status="scored_component_comparison_not_validation" if samples else "no_comparable_points",
    )


__all__ = [
    "AijCaseJComparison", "AijCaseJPointComparison",
    "ETHYLENE_MOLAR_MASS_KG_MOL", "compare_aij_case_j",
]
