"""Strict loader and component comparison for the AIJ Case H benchmark.

AIJ Case H is a neutral wind-tunnel building-flow and ethylene-tracer dataset.
It is useful independent evidence for obstacle velocity and wake structure, but
it is not cryogenic-hydrogen validation.  The current two-dimensional Tier-2
model also lacks lateral transport and resolved turbulence, so this module
withholds concentration validation and reports mean center-plane velocity as
conditional component evidence only.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from degali.addons.resolved_obstacle_flow_2d import ResolvedObstacleFlow2DResult
from degali.addons.fast_obstacle_transport_3d import FastObstacleTransport3DResult
from degali.addons.resolved_obstacle_flow_3d import ResolvedObstacleFlow3DResult


AIJ_CASE_H_DATASET_SCHEMA = "degali.aij-case-h-dataset.v1"
AIJ_CASE_H_FLOW_COMPARISON_SCHEMA = "degali.aij-case-h-flow-comparison.v1"
AIJ_CASE_H_RECORD_DOI = "10.5281/zenodo.15429872"
AIJ_CASE_H_DECLARED_README_DOI = "10.5281/zenodo.15429871"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _number(value: str, name: str, *, missing_allowed: bool = True) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"AIJ Case H {name} must be numeric") from error
    if not math.isfinite(result):
        raise ValueError(f"AIJ Case H {name} must be finite")
    if missing_allowed and math.isclose(result, -999.9, rel_tol=0.0, abs_tol=1.0e-12):
        return None
    return result


@dataclass(frozen=True)
class AijCaseHApproachPoint:
    x_m: float
    y_m: float
    z_m: float
    u_m_s: float
    v_m_s: float
    w_m_s: float
    u_rms_m_s: float
    v_rms_m_s: float
    w_rms_m_s: float
    uw_m2_s2: float
    uv_m2_s2: float
    turbulent_kinetic_energy_m2_s2: float


@dataclass(frozen=True)
class AijCaseHObservation:
    point_id: int
    x_m: float
    y_m: float
    z_m: float
    u_m_s: float | None
    v_m_s: float | None
    w_m_s: float | None
    u_rms_m_s: float | None
    v_rms_m_s: float | None
    w_rms_m_s: float | None
    turbulent_kinetic_energy_m2_s2: float | None
    concentration_ppm: float | None
    concentration_rms_ppm: float | None


@dataclass(frozen=True)
class AijCaseHDataset:
    schema: str
    record_doi: str
    declared_readme_doi: str
    approach_flow_sha256: str
    result_sha256: str
    readme_sha256: str
    approach_points: tuple[AijCaseHApproachPoint, ...]
    observations: tuple[AijCaseHObservation, ...]
    building_x_bounds_m: tuple[float, float] = (-0.10, 0.0)
    building_y_bounds_m: tuple[float, float] = (-0.05, 0.05)
    building_height_m: float = 0.20
    source_position_m: tuple[float, float, float] = (0.05, 0.0, 0.0)
    reference_height_m: float = 0.20
    reference_speed_m_s: float = 3.23
    reynolds_number: float = 43_000.0
    tracer: str = "ethylene"
    warnings: tuple[str, ...] = (
        "Zenodo record DOI and DOI declared inside the supplied readme differ",
        "Case H is neutral ethylene dispersion, not cryogenic hydrogen",
    )

    @property
    def centerplane_flow_observations(self) -> tuple[AijCaseHObservation, ...]:
        return tuple(
            item for item in self.observations
            if abs(item.y_m) <= 1.0e-12
            and item.u_m_s is not None and item.w_m_s is not None
        )


@dataclass(frozen=True)
class AijCaseHFlowComparisonPoint:
    point_id: int
    x_over_h: float
    z_over_h: float
    measured_u_over_uh: float
    predicted_u_over_uh: float
    measured_w_over_uh: float
    predicted_w_over_uh: float


@dataclass(frozen=True)
class AijCaseHFlowComparison:
    schema: str
    status: str
    points: tuple[AijCaseHFlowComparisonPoint, ...]
    skipped_point_ids: tuple[int, ...]
    u_mae_normalized: float
    u_rmse_normalized: float
    u_correlation: float
    w_mae_normalized: float
    w_rmse_normalized: float
    w_correlation: float
    u_recirculation_sign_accuracy: float
    diagnostics: dict[str, object]
    warnings: tuple[str, ...]


def read_aij_case_h(directory: str | Path) -> AijCaseHDataset:
    """Read the official Zenodo CSV package with immutable file hashes."""
    root = Path(directory).resolve()
    paths = {
        "approach": root / "AF_caseH.csv",
        "result": root / "RS_caseH.csv",
        "readme": root / "readme_caseH.md",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing AIJ Case H files: " + ", ".join(missing))

    approach_points = []
    with paths["approach"].open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {
            "x (m)", " y (m)", " z (m)", "U (m/s)", "V (m/s)", "W (m/s)",
            "u_rms (m/s)", "v_rms (m/s)", "w_rms (m/s)", "uw (m2/s2)",
            "uv (m2/s2)", "k (m2/s2)",
        }
        if set(reader.fieldnames or ()) != expected:
            raise ValueError("AIJ Case H approach-flow CSV header is unsupported")
        for row in reader:
            values = [
                _number(row[key], key, missing_allowed=False)
                for key in reader.fieldnames or ()
            ]
            approach_points.append(AijCaseHApproachPoint(*values))

    observations = []
    with paths["result"].open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = (
            "No.", "x (m)", "y (m)", "z (m)", "U (m/s)", "V (m/s)",
            "W (m/s)", "u_rms (m/s)", "v_rms (m/s)", "w_rsm (m/s)",
            "k (m2/s2)", "<c> (ppm)", "c_rms (ppm)",
        )
        if tuple(reader.fieldnames or ()) != expected:
            raise ValueError("AIJ Case H result CSV header is unsupported")
        for row in reader:
            point_id_value = _number(row["No."], "No.", missing_allowed=False)
            assert point_id_value is not None
            observations.append(AijCaseHObservation(
                point_id=int(point_id_value),
                x_m=float(_number(row["x (m)"], "x", missing_allowed=False)),
                y_m=float(_number(row["y (m)"], "y", missing_allowed=False)),
                z_m=float(_number(row["z (m)"], "z", missing_allowed=False)),
                u_m_s=_number(row["U (m/s)"], "U"),
                v_m_s=_number(row["V (m/s)"], "V"),
                w_m_s=_number(row["W (m/s)"], "W"),
                u_rms_m_s=_number(row["u_rms (m/s)"], "u_rms"),
                v_rms_m_s=_number(row["v_rms (m/s)"], "v_rms"),
                w_rms_m_s=_number(row["w_rsm (m/s)"], "w_rsm"),
                turbulent_kinetic_energy_m2_s2=_number(row["k (m2/s2)"], "k"),
                concentration_ppm=_number(row["<c> (ppm)"], "concentration"),
                concentration_rms_ppm=_number(row["c_rms (ppm)"], "concentration rms"),
            ))
    if len(approach_points) != 17:
        raise ValueError(f"AIJ Case H expected 17 approach points, found {len(approach_points)}")
    if len(observations) != 876:
        raise ValueError(f"AIJ Case H expected 876 observations, found {len(observations)}")
    point_ids = tuple(item.point_id for item in observations)
    if len(set(point_ids)) != len(point_ids):
        raise ValueError("AIJ Case H observation point IDs must be unique")
    return AijCaseHDataset(
        schema=AIJ_CASE_H_DATASET_SCHEMA,
        record_doi=AIJ_CASE_H_RECORD_DOI,
        declared_readme_doi=AIJ_CASE_H_DECLARED_README_DOI,
        approach_flow_sha256=_sha256(paths["approach"]),
        result_sha256=_sha256(paths["result"]),
        readme_sha256=_sha256(paths["readme"]),
        approach_points=tuple(approach_points),
        observations=tuple(observations),
    )


def _bilinear_sample(
    x: np.ndarray,
    z: np.ndarray,
    field: np.ndarray,
    solid: np.ndarray,
    x_value: float,
    z_value: float,
) -> float | None:
    if x_value < x[0] or x_value > x[-1] or z_value < z[0] or z_value > z[-1]:
        return None
    i1 = int(np.searchsorted(x, x_value, side="right"))
    j1 = int(np.searchsorted(z, z_value, side="right"))
    i1 = min(max(i1, 1), len(x) - 1)
    j1 = min(max(j1, 1), len(z) - 1)
    i0, j0 = i1 - 1, j1 - 1
    fx = float((x_value - x[i0]) / (x[i1] - x[i0]))
    fz = float((z_value - z[j0]) / (z[j1] - z[j0]))
    stencil = (
        (j0, i0, (1.0 - fx) * (1.0 - fz)),
        (j0, i1, fx * (1.0 - fz)),
        (j1, i0, (1.0 - fx) * fz),
        (j1, i1, fx * fz),
    )
    if any(solid[j, i] and weight > 1.0e-15 for j, i, weight in stencil):
        return None
    return float(sum(weight * field[j, i] for j, i, weight in stencil))


def _correlation(measured: np.ndarray, predicted: np.ndarray) -> float:
    if len(measured) < 2 or np.std(measured) <= 0.0 or np.std(predicted) <= 0.0:
        return 0.0
    return float(np.corrcoef(measured, predicted)[0, 1])


def compare_aij_case_h_centerplane_flow(
    dataset: AijCaseHDataset,
    result: ResolvedObstacleFlow2DResult,
    *,
    model_x_origin_in_experiment_m: float,
) -> AijCaseHFlowComparison:
    """Compare mean U/W only; concentration and turbulence stay withheld."""
    if not isinstance(dataset, AijCaseHDataset):
        raise TypeError("dataset must be AijCaseHDataset")
    if not isinstance(result, ResolvedObstacleFlow2DResult):
        raise TypeError("result must be ResolvedObstacleFlow2DResult")
    if isinstance(model_x_origin_in_experiment_m, bool) or not math.isfinite(
        float(model_x_origin_in_experiment_m)
    ):
        raise ValueError("model_x_origin_in_experiment_m must be finite")
    points = []
    skipped = []
    for observation in dataset.centerplane_flow_observations:
        model_x = observation.x_m - float(model_x_origin_in_experiment_m)
        predicted_u = _bilinear_sample(
            result.x_m, result.z_m, result.velocity_u_m_s, result.solid_mask,
            model_x, observation.z_m,
        )
        predicted_w = _bilinear_sample(
            result.x_m, result.z_m, result.velocity_w_m_s, result.solid_mask,
            model_x, observation.z_m,
        )
        if predicted_u is None or predicted_w is None:
            skipped.append(observation.point_id)
            continue
        if not math.isfinite(predicted_u) or not math.isfinite(predicted_w):
            raise RuntimeError(
                "AIJ Case H comparison encountered a non-finite model velocity"
            )
        assert observation.u_m_s is not None and observation.w_m_s is not None
        points.append(AijCaseHFlowComparisonPoint(
            point_id=observation.point_id,
            x_over_h=observation.x_m / dataset.reference_height_m,
            z_over_h=observation.z_m / dataset.reference_height_m,
            measured_u_over_uh=observation.u_m_s / dataset.reference_speed_m_s,
            predicted_u_over_uh=predicted_u / dataset.reference_speed_m_s,
            measured_w_over_uh=observation.w_m_s / dataset.reference_speed_m_s,
            predicted_w_over_uh=predicted_w / dataset.reference_speed_m_s,
        ))
    if not points:
        raise ValueError("no AIJ Case H center-plane points overlap the model domain")
    measured_u = np.asarray([item.measured_u_over_uh for item in points])
    predicted_u = np.asarray([item.predicted_u_over_uh for item in points])
    measured_w = np.asarray([item.measured_w_over_uh for item in points])
    predicted_w = np.asarray([item.predicted_w_over_uh for item in points])
    recirculation_mask = np.abs(measured_u) > 0.02
    sign_accuracy = float(np.mean(
        np.sign(measured_u[recirculation_mask])
        == np.sign(predicted_u[recirculation_mask])
    )) if np.any(recirculation_mask) else 0.0
    return AijCaseHFlowComparison(
        schema=AIJ_CASE_H_FLOW_COMPARISON_SCHEMA,
        status="component_evidence_conditional",
        points=tuple(points),
        skipped_point_ids=tuple(skipped),
        u_mae_normalized=float(np.mean(np.abs(predicted_u - measured_u))),
        u_rmse_normalized=float(np.sqrt(np.mean((predicted_u - measured_u) ** 2))),
        u_correlation=_correlation(measured_u, predicted_u),
        w_mae_normalized=float(np.mean(np.abs(predicted_w - measured_w))),
        w_rmse_normalized=float(np.sqrt(np.mean((predicted_w - measured_w) ** 2))),
        w_correlation=_correlation(measured_w, predicted_w),
        u_recirculation_sign_accuracy=sign_accuracy,
        diagnostics={
            "comparison_scope": "mean center-plane U and W only",
            "evaluation_only": True,
            "calibration_used": False,
            "quantitative_hydrogen_validation": False,
            "concentration_validation_withheld": True,
            "turbulence_validation_withheld": True,
            "evaluated_point_count": len(points),
            "skipped_point_count": len(skipped),
        },
        warnings=(
            "the Tier-2 model is two-dimensional while Case H is a finite-width building",
            "effective viscosity is prescribed and does not reproduce the measured turbulent Reynolds number",
            "TKE and Reynolds stresses are not solved",
            "ethylene point-source concentration cannot be compared to a per-span 2-D scalar source",
            "this comparison cannot qualify cryogenic hydrogen predictions",
        ),
    )


def _trilinear_sample_3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    field: np.ndarray,
    solid: np.ndarray,
    x_value: float,
    y_value: float,
    z_value: float,
) -> float | None:
    coordinates_and_values = ((x, x_value), (y, y_value), (z, z_value))
    brackets = []
    for coordinates, value in coordinates_and_values:
        if value < coordinates[0] or value > coordinates[-1]:
            return None
        upper = int(np.searchsorted(coordinates, value, side="right"))
        upper = min(max(upper, 1), len(coordinates) - 1)
        lower = upper - 1
        fraction = float(
            (value - coordinates[lower])
            / (coordinates[upper] - coordinates[lower])
        )
        brackets.append((lower, upper, fraction))
    bx, by, bz = brackets
    stencil = []
    for k, wz in ((bz[0], 1.0 - bz[2]), (bz[1], bz[2])):
        for j, wy in ((by[0], 1.0 - by[2]), (by[1], by[2])):
            for i, wx in ((bx[0], 1.0 - bx[2]), (bx[1], bx[2])):
                weight = wx * wy * wz
                if weight > 1.0e-15:
                    stencil.append((k, j, i, weight))
    if any(solid[k, j, i] for k, j, i, _weight in stencil):
        return None
    return float(sum(weight * field[k, j, i] for k, j, i, weight in stencil))


def compare_aij_case_h_centerplane_flow_3d(
    dataset: AijCaseHDataset,
    result: FastObstacleTransport3DResult | ResolvedObstacleFlow3DResult,
    *,
    model_x_origin_in_experiment_m: float,
    model_y_origin_in_experiment_m: float,
) -> AijCaseHFlowComparison:
    """Compare AIJ mean U/W to the finite-width 3-D potential baseline."""
    if not isinstance(dataset, AijCaseHDataset):
        raise TypeError("dataset must be AijCaseHDataset")
    if not isinstance(result, (FastObstacleTransport3DResult, ResolvedObstacleFlow3DResult)):
        raise TypeError(
            "result must be FastObstacleTransport3DResult or "
            "ResolvedObstacleFlow3DResult"
        )
    resolved_viscous = isinstance(result, ResolvedObstacleFlow3DResult)
    for name, value in (
        ("model_x_origin_in_experiment_m", model_x_origin_in_experiment_m),
        ("model_y_origin_in_experiment_m", model_y_origin_in_experiment_m),
    ):
        if isinstance(value, bool) or not math.isfinite(float(value)):
            raise ValueError(f"{name} must be finite")
    points = []
    skipped = []
    for observation in dataset.centerplane_flow_observations:
        model_x = observation.x_m - float(model_x_origin_in_experiment_m)
        model_y = observation.y_m - float(model_y_origin_in_experiment_m)
        predicted_u = _trilinear_sample_3d(
            result.x_m, result.y_m, result.z_m,
            result.velocity_u_m_s, result.solid_mask,
            model_x, model_y, observation.z_m,
        )
        predicted_w = _trilinear_sample_3d(
            result.x_m, result.y_m, result.z_m,
            result.velocity_w_m_s, result.solid_mask,
            model_x, model_y, observation.z_m,
        )
        if predicted_u is None or predicted_w is None:
            skipped.append(observation.point_id)
            continue
        if not math.isfinite(predicted_u) or not math.isfinite(predicted_w):
            raise RuntimeError(
                "AIJ Case H comparison encountered a non-finite 3-D model velocity"
            )
        assert observation.u_m_s is not None and observation.w_m_s is not None
        points.append(AijCaseHFlowComparisonPoint(
            point_id=observation.point_id,
            x_over_h=observation.x_m / dataset.reference_height_m,
            z_over_h=observation.z_m / dataset.reference_height_m,
            measured_u_over_uh=observation.u_m_s / dataset.reference_speed_m_s,
            predicted_u_over_uh=predicted_u / dataset.reference_speed_m_s,
            measured_w_over_uh=observation.w_m_s / dataset.reference_speed_m_s,
            predicted_w_over_uh=predicted_w / dataset.reference_speed_m_s,
        ))
    if not points:
        raise ValueError("no AIJ Case H center-plane points overlap the 3-D model domain")
    measured_u = np.asarray([item.measured_u_over_uh for item in points])
    predicted_u = np.asarray([item.predicted_u_over_uh for item in points])
    measured_w = np.asarray([item.measured_w_over_uh for item in points])
    predicted_w = np.asarray([item.predicted_w_over_uh for item in points])
    recirculation_mask = np.abs(measured_u) > 0.02
    sign_accuracy = float(np.mean(
        np.sign(measured_u[recirculation_mask])
        == np.sign(predicted_u[recirculation_mask])
    )) if np.any(recirculation_mask) else 0.0
    return AijCaseHFlowComparison(
        schema=AIJ_CASE_H_FLOW_COMPARISON_SCHEMA,
        status="component_evidence_conditional",
        points=tuple(points),
        skipped_point_ids=tuple(skipped),
        u_mae_normalized=float(np.mean(np.abs(predicted_u - measured_u))),
        u_rmse_normalized=float(np.sqrt(np.mean((predicted_u - measured_u) ** 2))),
        u_correlation=_correlation(measured_u, predicted_u),
        w_mae_normalized=float(np.mean(np.abs(predicted_w - measured_w))),
        w_rmse_normalized=float(np.sqrt(np.mean((predicted_w - measured_w) ** 2))),
        w_correlation=_correlation(measured_w, predicted_w),
        u_recirculation_sign_accuracy=sign_accuracy,
        diagnostics={
            "comparison_scope": "mean center-plane U and W only",
            "model_dimension": "finite-width 3-D",
            "model_flow_class": (
                "viscous D3Q19"
                if resolved_viscous else "inviscid potential projection"
            ),
            "evaluation_only": True,
            "calibration_used": False,
            "quantitative_hydrogen_validation": False,
            "concentration_validation_withheld": True,
            "turbulence_validation_withheld": True,
            "evaluated_point_count": len(points),
            "skipped_point_count": len(skipped),
        },
        warnings=(
            (
                "D3Q19 prototype has unverified inlet/outlet boundaries and SGS wall treatment"
                if resolved_viscous
                else "three-dimensional potential flow is inviscid and cannot form the measured wake"
            ),
            (
                "isothermal one-way flow omits cryogenic energy and buoyancy coupling"
                if resolved_viscous
                else "slip walls omit the building and ground boundary layers"
            ),
            "TKE and Reynolds stresses are not independently validated",
            "ethylene concentration is outside this zero-source flow comparison",
            "this comparison cannot qualify cryogenic hydrogen predictions",
        ),
    )


def aij_case_h_flow_comparison_record(
    comparison: AijCaseHFlowComparison,
) -> dict[str, object]:
    if not isinstance(comparison, AijCaseHFlowComparison):
        raise TypeError("comparison must be AijCaseHFlowComparison")
    return {
        "schema": comparison.schema,
        "status": comparison.status,
        "metrics": {
            "u_mae_normalized": comparison.u_mae_normalized,
            "u_rmse_normalized": comparison.u_rmse_normalized,
            "u_correlation": comparison.u_correlation,
            "w_mae_normalized": comparison.w_mae_normalized,
            "w_rmse_normalized": comparison.w_rmse_normalized,
            "w_correlation": comparison.w_correlation,
            "u_recirculation_sign_accuracy": comparison.u_recirculation_sign_accuracy,
        },
        "diagnostics": comparison.diagnostics,
        "warnings": list(comparison.warnings),
        "skipped_point_ids": list(comparison.skipped_point_ids),
        "points": [item.__dict__ for item in comparison.points],
    }


__all__ = [
    "AIJ_CASE_H_DATASET_SCHEMA", "AIJ_CASE_H_FLOW_COMPARISON_SCHEMA",
    "AIJ_CASE_H_RECORD_DOI", "AIJ_CASE_H_DECLARED_README_DOI",
    "AijCaseHApproachPoint", "AijCaseHObservation", "AijCaseHDataset",
    "AijCaseHFlowComparisonPoint", "AijCaseHFlowComparison",
    "read_aij_case_h", "compare_aij_case_h_centerplane_flow",
    "compare_aij_case_h_centerplane_flow_3d",
    "aij_case_h_flow_comparison_record",
]
