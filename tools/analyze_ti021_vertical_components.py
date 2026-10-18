"""Attribute TI021 dose error to saved physical observation components.

The input must be a component-diagnostic dense-layer report and its adjacent
NPZ bundle.  This tool only integrates already saved physical branches; it
does not fit sensor gains, offsets, or correction factors.
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


def _summary(
    mask: np.ndarray,
    observed: np.ndarray,
    predicted: np.ndarray,
    component_dose: np.ndarray,
    labels: list[str],
) -> dict:
    selected = mask & (observed > 0.0)
    if not np.any(selected):
        return {"positive_sensor_count": 0}
    obs = observed[selected]
    pred = predicted[selected]
    component = component_dose[:, selected]
    total_predicted = float(np.sum(pred))
    return {
        "positive_sensor_count": int(np.count_nonzero(selected)),
        "mape_percent": float(np.mean(np.abs(pred / obs - 1.0)) * 100.0),
        "normalized_mean_bias": float(np.sum(pred - obs) / np.sum(obs)),
        "observed_dose_sum_vol_percent_s": float(np.sum(obs)),
        "predicted_dose_sum_vol_percent_s": total_predicted,
        "component_dose_sum_vol_percent_s": {
            label: float(np.sum(component[index]))
            for index, label in enumerate(labels)
        },
        "component_share_of_predicted_dose": {
            label: (
                float(np.sum(component[index]) / total_predicted)
                if total_predicted > 0.0
                else 0.0
            )
            for index, label in enumerate(labels)
        },
    }


def analyze(report_path: Path, workbook_path: Path) -> dict:
    report_path = report_path.resolve()
    bundle_path = report_path.with_suffix(".npz")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    dataset = read_thorney_ti021_workbook(workbook_path)
    rows = report["sensor_rows"]
    if len(rows) != len(dataset.dose_sensors):
        raise ValueError("report/workbook sensor count differs")

    observed = np.asarray(
        [row["observed_dose_vol_percent_s"] for row in rows], dtype=float
    )
    predicted = np.asarray(
        [row["predicted_dose_vol_percent_s"] for row in rows], dtype=float
    )
    authoritative = np.asarray(
        [sensor.dose_percent_s for sensor in dataset.dose_sensors], dtype=float
    )
    if not np.array_equal(observed, authoritative):
        raise ValueError("report observations are not the authoritative workbook doses")

    with np.load(bundle_path, allow_pickle=False) as bundle:
        time_s = np.asarray(bundle["time_s"], dtype=float)
        concentration = np.asarray(
            bundle["concentration_vol_percent"], dtype=float
        )
        component_trace = np.asarray(
            bundle["observation_component_concentration_vol_percent"], dtype=float
        )
        labels = [str(item) for item in bundle["observation_component_labels"]]

    expected_shape = (len(labels), observed.size, time_s.size)
    if component_trace.shape != expected_shape:
        raise ValueError(
            f"component trace shape {component_trace.shape} != {expected_shape}"
        )
    reconstructed = np.sum(component_trace, axis=0)
    trace_error = float(np.max(np.abs(reconstructed - concentration)))
    if not np.allclose(reconstructed, concentration, rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("component traces do not reconstruct total concentration")
    integrated = np.trapezoid(concentration, time_s, axis=1)
    component_dose = np.trapezoid(component_trace, time_s, axis=2)
    dose_error = float(np.max(np.abs(np.sum(component_dose, axis=0) - predicted)))
    if not np.allclose(integrated, predicted, rtol=1.0e-12, atol=1.0e-10):
        raise ValueError("total trace does not reconstruct reported dose")
    if not np.allclose(
        np.sum(component_dose, axis=0), predicted, rtol=1.0e-12, atol=1.0e-10
    ):
        raise ValueError("component doses do not reconstruct reported dose")

    heights = np.asarray(
        [sensor.height_m for sensor in dataset.dose_sensors], dtype=float
    )
    downwind = np.asarray(
        [sensor.local_y_m for sensor in dataset.dose_sensors], dtype=float
    )
    groups: dict[str, np.ndarray] = {
        "all": np.ones(observed.size, dtype=bool),
        **{
            f"height_{height:g}_m": np.isclose(heights, height)
            for height in sorted(set(heights))
        },
        "downwind_lt_150_m": downwind < 150.0,
        "downwind_150_to_300_m": (downwind >= 150.0) & (downwind < 300.0),
        "downwind_ge_300_m": downwind >= 300.0,
    }
    positive = observed > 0.0
    relative_error = np.full(observed.shape, np.nan)
    relative_error[positive] = predicted[positive] / observed[positive] - 1.0
    ranked = np.flatnonzero(positive)[
        np.argsort(np.abs(relative_error[positive]))[::-1]
    ]
    worst_sensors = []
    for index in ranked[:20]:
        sensor = dataset.dose_sensors[int(index)]
        worst_sensors.append(
            {
                "label": rows[int(index)]["label"],
                "height_m": sensor.height_m,
                "downwind_m": sensor.local_y_m,
                "observed_dose_vol_percent_s": float(observed[index]),
                "predicted_dose_vol_percent_s": float(predicted[index]),
                "relative_error_percent": float(relative_error[index] * 100.0),
                "component_dose_vol_percent_s": {
                    label: float(component_dose[position, index])
                    for position, label in enumerate(labels)
                },
            }
        )

    return {
        "schema": "degali.ti021-vertical-component-attribution.v1",
        "method": (
            "direct integration of saved physical-component concentration traces; "
            "no sensor-specific fitting"
        ),
        "report": str(report_path),
        "report_sha256": _sha256(report_path),
        "bundle": str(bundle_path),
        "bundle_sha256": _sha256(bundle_path),
        "component_labels": labels,
        "maximum_component_trace_reconstruction_error_vol_percent": trace_error,
        "maximum_component_dose_reconstruction_error_vol_percent_s": dose_error,
        "groups": {
            name: _summary(mask, observed, predicted, component_dose, labels)
            for name, mask in groups.items()
        },
        "worst_sensors": worst_sensors,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--workbook", type=Path, default=ROOT / "ti21_input_ascii.xls")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.report, args.workbook)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
