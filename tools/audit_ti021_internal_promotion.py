"""Audit an internal TI021 engineering-benchmark promotion.

This audit deliberately separates an internal numerical benchmark promotion
from external validation qualification.  A candidate may replace the previous
engineering benchmark only when it is independently reproduced, improves the
same 145-sensor observable, does not degrade any preregistered height or
distance group, satisfies the conservation gates, and runs within the declared
wall-clock limit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


ROOT = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _arrays(record: dict) -> tuple[list[str], np.ndarray, np.ndarray]:
    rows = record["sensor_rows"]
    return (
        [str(row["label"]) for row in rows],
        np.asarray([row["observed_dose_vol_percent_s"] for row in rows], dtype=float),
        np.asarray([row["predicted_dose_vol_percent_s"] for row in rows], dtype=float),
    )


def _mape(predicted: np.ndarray, observed: np.ndarray, mask: np.ndarray) -> tuple[float, int]:
    selected = mask & (observed > 0.0)
    if not np.any(selected):
        raise ValueError("group has no positive observations")
    return (
        float(np.mean(np.abs(predicted[selected] / observed[selected] - 1.0)) * 100.0),
        int(np.count_nonzero(selected)),
    )


def _groups(predicted: np.ndarray, observed: np.ndarray) -> dict[str, dict]:
    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
    if len(dataset.dose_sensors) != predicted.size:
        raise ValueError("TI021 workbook/result sensor count mismatch")
    heights = np.asarray([sensor.height_m for sensor in dataset.dose_sensors], dtype=float)
    downwind = np.asarray([sensor.local_y_m for sensor in dataset.dose_sensors], dtype=float)
    result: dict[str, dict] = {}
    for height in sorted(set(heights)):
        value, count = _mape(predicted, observed, np.isclose(heights, height))
        result[f"height_{height:g}_m"] = {"mape_percent": value, "positive_count": count}
    for name, lower, upper in (
        ("downwind_lt_150_m", -np.inf, 150.0),
        ("downwind_150_to_300_m", 150.0, 300.0),
        ("downwind_ge_300_m", 300.0, np.inf),
    ):
        value, count = _mape(
            predicted, observed, (downwind >= lower) & (downwind < upper)
        )
        result[name] = {"mape_percent": value, "positive_count": count}
    return result


def _runtime(record: dict) -> float:
    for key in ("runtime_including_output_s", "runtime_with_audit_s", "prediction_runtime_s"):
        if key in record:
            return float(record[key])
    raise ValueError("candidate has no recognized runtime field")


def _wall_diffusion_evidence(record: dict) -> dict:
    required = (float(record.get("wall_height_m", 6.)) > 0.
                and float(record.get("horizontal_diffusivity_m2_s", 0.)) > 0.)
    respected = record.get("horizontal_diffusion_respect_walls", False) is True
    return {"finite_wall_diffusion_required": required,
            "finite_wall_diffusion_enabled": respected,
            "passed": bool(not required or respected),
            "scope": "known horizontal-diffusion boundary defect only; not full physical validation"}


def _vertical_observation_evidence(record: dict) -> dict:
    diagnostics = record.get('diagnostics', {})
    enabled = (record.get('vertical_observation_mode') == 'conservative_volume'
               and diagnostics.get('vertical_observation_mode') == 'conservative_volume')
    residual = diagnostics.get('vertical_observation_maximum_relative_volume_residual')
    passed = (enabled and residual is not None and np.isfinite(residual)
              and 0. <= residual <= 1.e-9)
    return {'conservative_species_volume_enabled': enabled,
            'maximum_relative_species_volume_residual': residual,
            'passed': bool(passed),
            'scope': 'represented-column species volume only; not vertical momentum or field validation'}


def _trace_evidence(path: Path, record: dict) -> dict:
    """Reintegrate saved traces; a duration label alone cannot prove 600 s."""
    trace_path = path.with_suffix(".npz")
    if not trace_path.is_file():
        return {"passed": False, "reason": "missing time-history NPZ"}
    _, observed, predicted = _arrays(record)
    with np.load(trace_path, allow_pickle=False) as data:
        time_s = np.asarray(data["time_s"], dtype=float)
        concentration = np.asarray(data["concentration_vol_percent"], dtype=float)
        saved_dose = np.asarray(data["predicted_dose_vol_percent_s"], dtype=float)
        saved_observed = np.asarray(data["observed_dose_vol_percent_s"], dtype=float)
    shape_ok = bool(time_s.ndim == 1 and time_s.size >= 2
                    and concentration.shape == (predicted.size, time_s.size))
    finite = bool(np.all(np.isfinite(time_s))
                  and np.all(np.isfinite(concentration)))
    complete = bool(shape_ok and finite and np.all(np.diff(time_s) > 0.)
                    and abs(float(time_s[0])) <= 1.e-9
                    and abs(float(time_s[-1]) - 600.) <= 1.e-9
                    and record.get("duration_s") == 600.)
    difference = None
    integral_matches = False
    if shape_ok and finite:
        integrated = np.trapezoid(concentration, time_s, axis=1)
        difference = float(np.max(np.abs(integrated - predicted)))
        integral_matches = bool(np.allclose(integrated, predicted, rtol=1.e-12, atol=1.e-10)
                                and np.array_equal(saved_dose, predicted)
                                and np.array_equal(saved_observed, observed))
    bounded = bool(finite and np.all(concentration >= -1.e-10)
                   and np.all(concentration <= 100. + 1.e-10))
    return {"path": str(trace_path.resolve()), "sha256": _sha256(trace_path),
            "complete_0_to_600_s": complete,
            "concentration_bounded_0_to_100_vol_percent": bounded,
            "reintegration_matches_report": integral_matches,
            "maximum_dose_difference_vol_percent_s": difference,
            "passed": bool(complete and bounded and integral_matches)}


def _implementation_evidence(record: dict) -> dict:
    declared = record.get("implementation_sha256")
    if not isinstance(declared, dict) or not declared:
        return {"passed": False, "reason": "missing implementation SHA-256 map"}
    rows = []
    for raw_path, expected in declared.items():
        path = (ROOT / str(raw_path)).resolve()
        try:
            path.relative_to(ROOT.resolve())
        except ValueError:
            rows.append({"path": str(path), "passed": False,
                         "reason": "implementation path leaves repository"})
            continue
        actual = _sha256(path) if path.is_file() else None
        rows.append({"path": str(path), "expected_sha256": expected,
                     "actual_sha256": actual,
                     "passed": bool(actual is not None and actual == expected)})
    return {"files": rows, "passed": bool(rows and all(row["passed"] for row in rows))}


def _trace_identity(candidate_path: Path, reproduction_path: Path) -> dict:
    with np.load(candidate_path.with_suffix('.npz'), allow_pickle=False) as left:
        left_time=np.asarray(left['time_s'])
        left_trace=np.asarray(left['concentration_vol_percent'])
    with np.load(reproduction_path.with_suffix('.npz'), allow_pickle=False) as right:
        right_time=np.asarray(right['time_s'])
        right_trace=np.asarray(right['concentration_vol_percent'])
    identical=bool(np.array_equal(left_time,right_time)
                   and np.array_equal(left_trace,right_trace))
    maximum=(float(np.max(np.abs(left_trace-right_trace)))
             if left_trace.shape==right_trace.shape else None)
    return {"bitwise_identical_time_and_concentration":identical,
            "maximum_absolute_concentration_difference_vol_percent":maximum,
            "passed":identical}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--reproduction", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runtime-limit-s", type=float, default=30.0)
    parser.add_argument("--target-mape-percent", type=float, default=10.0)
    args = parser.parse_args()

    baseline = _load(args.baseline)
    candidate = _load(args.candidate)
    reproduction = _load(args.reproduction)
    baseline_labels, baseline_observed, baseline_predicted = _arrays(baseline)
    candidate_labels, candidate_observed, candidate_predicted = _arrays(candidate)
    reproduction_labels, reproduction_observed, reproduction_predicted = _arrays(reproduction)

    dataset = read_thorney_ti021_workbook(ROOT / "ti21_input_ascii.xls")
    authoritative_observed = np.array([s.dose_percent_s for s in dataset.dose_sensors])
    expected_labels = [f"TI021-S{i+1:03d}" for i in range(len(dataset.dose_sensors))]
    authoritative_identity = bool(
        baseline_labels == expected_labels
        and np.array_equal(baseline_observed, authoritative_observed)
        and all(record.get("workbook_sha256") == dataset.workbook_sha256
                for record in (baseline, candidate, reproduction)))
    trace_evidence = {name: _trace_evidence(path, record) for name, path, record in (
        ("baseline", args.baseline, baseline), ("candidate", args.candidate, candidate),
        ("reproduction", args.reproduction, reproduction))}
    full_event_evidence = all(item["passed"] for item in trace_evidence.values())

    sensor_identity = bool(
        baseline_labels == candidate_labels == reproduction_labels
        and np.array_equal(baseline_observed, candidate_observed)
        and np.array_equal(candidate_observed, reproduction_observed)
    )
    observation_species_mode_identity = bool(
        baseline.get("observation_species_mode", "r12_component")
        == candidate.get("observation_species_mode", "r12_component")
        == reproduction.get("observation_species_mode", "r12_component")
    )
    vertical_observation_mode_identity = bool(
        baseline.get("vertical_observation_mode")
        == candidate.get("vertical_observation_mode")
        == reproduction.get("vertical_observation_mode")
        == "conservative_volume"
    )
    exact_reproduction = bool(np.array_equal(candidate_predicted, reproduction_predicted))
    max_reproduction_difference = float(
        np.max(np.abs(candidate_predicted - reproduction_predicted))
    )

    positive = baseline_observed > 0.0
    baseline_mape, positive_count = _mape(
        baseline_predicted, baseline_observed, np.ones_like(positive, dtype=bool)
    )
    candidate_mape, _ = _mape(
        candidate_predicted, candidate_observed, np.ones_like(positive, dtype=bool)
    )
    baseline_groups = _groups(baseline_predicted, baseline_observed)
    candidate_groups = _groups(candidate_predicted, candidate_observed)
    # Freeze route membership to the baseline; moving a sensor to a different
    # execution path must not remove it from that path's non-regression gate.
    routes = np.array([row["route"] for row in baseline["sensor_rows"]])
    for route in sorted(set(routes)):
        mask = routes == route
        if not np.any(mask & positive):
            continue
        for groups, prediction in ((baseline_groups, baseline_predicted),
                                   (candidate_groups, candidate_predicted)):
            value, count = _mape(prediction, baseline_observed, mask)
            groups[f"baseline_route_{route}"] = {"mape_percent": value, "positive_count": count}
    group_comparison = {
        name: {
            "baseline_mape_percent": baseline_groups[name]["mape_percent"],
            "candidate_mape_percent": candidate_groups[name]["mape_percent"],
            "delta_percentage_points": (
                candidate_groups[name]["mape_percent"]
                - baseline_groups[name]["mape_percent"]
            ),
            "positive_count": candidate_groups[name]["positive_count"],
            "nondegraded": bool(
                candidate_groups[name]["mape_percent"]
                <= baseline_groups[name]["mape_percent"] + 1.0e-12
            ),
        }
        for name in baseline_groups
    }
    all_groups_nondegraded = bool(
        all(item["nondegraded"] for item in group_comparison.values())
    )

    diagnostics = reproduction["diagnostics"]
    conservation = {
        "maximum_relative_species_mass_residual": float(
            diagnostics["maximum_relative_species_mass_residual"]
        ),
        "mixture_mass_residual_kg": float(diagnostics["mixture_mass_residual_kg"]),
        "momentum_conservation_gate": diagnostics["momentum_conservation_gate"] is True,
    }
    conservation["passed"] = bool(
        abs(conservation["maximum_relative_species_mass_residual"]) <= 1.0e-8
        and abs(conservation["mixture_mass_residual_kg"]) <= 1.0e-6
        and conservation["momentum_conservation_gate"]
    )
    runtime_s = _runtime(reproduction)
    runtime_passed = bool(runtime_s <= args.runtime_limit_s)
    wall_diffusion_evidence = {
        name: _wall_diffusion_evidence(record)
        for name, record in (("candidate", candidate), ("reproduction", reproduction))}
    wall_diffusion_passed = all(item["passed"] for item in wall_diffusion_evidence.values())
    vertical_observation_evidence = {
        name: _vertical_observation_evidence(record)
        for name, record in (("candidate", candidate), ("reproduction", reproduction))}
    implementation_evidence = {
        name: _implementation_evidence(record)
        for name, record in (("baseline", baseline), ("candidate", candidate),
                             ("reproduction", reproduction))}
    implementation_passed = all(
        item["passed"] for item in implementation_evidence.values())
    trace_identity = _trace_identity(args.candidate, args.reproduction)
    improves_baseline = bool(candidate_mape < baseline_mape)
    internal_promotion_allowed = bool(
        sensor_identity
        and authoritative_identity
        and full_event_evidence
        and observation_species_mode_identity
        and vertical_observation_mode_identity
        and candidate.get("observation_species_mode") == "released_mixture"
        and args.candidate.resolve() != args.reproduction.resolve()
        and exact_reproduction
        and trace_identity["passed"]
        and implementation_passed
        and improves_baseline
        and all_groups_nondegraded
        and conservation["passed"]
        and wall_diffusion_passed
        and all(item['passed'] for item in vertical_observation_evidence.values())
        and runtime_passed
    )

    report = {
        "schema": "degali.ti021-internal-promotion-audit.v4",
        "status": (
            "internal_engineering_benchmark_promoted"
            if internal_promotion_allowed
            else "internal_engineering_benchmark_rejected"
        ),
        "scope_guard": (
            "Internal engineering benchmark only; this audit does not grant external "
            "field-validation or predictive-use qualification."
        ),
        "files": {
            "baseline": str(args.baseline.resolve()),
            "baseline_sha256": _sha256(args.baseline),
            "candidate": str(args.candidate.resolve()),
            "candidate_sha256": _sha256(args.candidate),
            "independent_reproduction": str(args.reproduction.resolve()),
            "independent_reproduction_sha256": _sha256(args.reproduction),
        },
        "observable": {
            "metric": "600 s integrated concentration dose MAPE in vol%-s",
            "sensor_count": len(candidate_labels),
            "positive_observed_sensor_count": positive_count,
            "sensor_order_and_observations_identical": sensor_identity,
            "observation_species_mode_identical": observation_species_mode_identity,
            "vertical_observation_mode_identical_and_conservative": (
                vertical_observation_mode_identity),
            "matches_authoritative_workbook": authoritative_identity,
        },
        "full_event_trace_evidence": trace_evidence,
        "implementation_evidence": implementation_evidence,
        "score": {
            "baseline_mape_percent": baseline_mape,
            "candidate_mape_percent": candidate_mape,
            "improvement_percentage_points": baseline_mape - candidate_mape,
            "improves_baseline": improves_baseline,
            "target_mape_percent": args.target_mape_percent,
            "target_reached": bool(candidate_mape <= args.target_mape_percent),
        },
        "independent_reproduction": {
            "bitwise_identical_predictions": exact_reproduction,
            "maximum_absolute_prediction_difference_vol_percent_s": max_reproduction_difference,
            **trace_identity,
        },
        "group_comparison": group_comparison,
        "all_preregistered_groups_nondegraded": all_groups_nondegraded,
        "conservation": conservation,
        "known_wall_diffusion_boundary_gate": wall_diffusion_evidence,
        "vertical_observation_species_volume_gate": vertical_observation_evidence,
        "runtime": {
            "runtime_including_output_s": runtime_s,
            "limit_s": args.runtime_limit_s,
            "passed": runtime_passed,
        },
        "internal_engineering_benchmark_promotion_allowed": internal_promotion_allowed,
        "external_validation_promotion_allowed": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
