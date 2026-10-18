"""Score short TI021 resolved runs for initial-condition diagnostics only.

This utility compares the first 10 s threshold crossing against the workbook's
arrival labels.  The workbook does not contain 0--10 s concentration traces,
so the result is a timing diagnostic, not concentration MAPE or validation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook
from degali.addons.field_evidence_manifest import (
    read_field_evidence_manifest_json,
)


def _manifest_clock_alignment(path: Path) -> dict[str, object]:
    """Return a verified clock alignment or a fail-closed withholding record."""
    try:
        manifest = read_field_evidence_manifest_json(path)
    except (OSError, ValueError, TypeError, KeyError) as error:
        return {
            "status": "withheld",
            "reason": f"FieldEvidenceManifest could not be verified: {error}",
        }
    readiness = manifest.quantitative_readiness_record()
    alignment = readiness.get("clock_alignment")
    if not isinstance(alignment, dict):
        return {
            "status": "withheld",
            "reason": "manifest did not emit clock_alignment metadata",
            "readiness": readiness,
        }
    if readiness.get("status") != "accepted" or alignment.get("status") != "accepted":
        return {
            "status": "withheld",
            "reason": "quantitative FieldEvidenceManifest readiness is not accepted",
            "readiness": readiness,
            "clock_alignment": alignment,
        }
    try:
        offset = float(alignment["offset_s"])
        uncertainty = float(alignment["uncertainty_s"])
    except (KeyError, TypeError, ValueError) as error:
        return {
            "status": "withheld",
            "reason": f"accepted clock alignment is not numeric: {error}",
            "readiness": readiness,
            "clock_alignment": alignment,
        }
    return {
        "status": "accepted",
        "offset_s": offset,
        "uncertainty_s": uncertainty,
        "reference_time_utc": alignment.get("reference_time_utc"),
        "synchronization_method": alignment.get("synchronization_method"),
        "covered_channels": alignment.get("covered_channels"),
        "readiness": readiness,
    }


def score_run(
    workbook: Path,
    run: Path,
    *,
    threshold: float = 0.1,
    diagnostic_time_origin_offset_s: float | None = None,
    field_evidence_manifest: Path | None = None,
    relative_peak_fraction: float | None = None,
) -> dict[str, object]:
    if diagnostic_time_origin_offset_s is not None and field_evidence_manifest is not None:
        raise ValueError(
            "diagnostic-time-origin-offset-s and field-evidence-manifest are mutually exclusive"
        )
    record = json.loads((run / "ti021-inner-transient.json").read_text(encoding="utf-8"))
    with np.load(run / "ti021-inner-receptor-traces.npz") as bundle:
        time_s = np.asarray(bundle["time_s"], dtype=float)
        traces = np.asarray(bundle["concentration_vol_percent"], dtype=float)
        labels = np.asarray(bundle["labels"]).astype(str)
    dataset = read_thorney_ti021_workbook(workbook)
    indices = np.asarray([
        item["workbook_row_index_zero_based"]
        for item in record["configuration"]["receptors"]
    ], dtype=int)
    observed = np.asarray([
        np.nan if dataset.dose_sensors[index].arrival_time_s is None
        else dataset.dose_sensors[index].arrival_time_s
        for index in indices
    ])
    predicted = np.asarray([
        time_s[hit[0]] if (hit := np.flatnonzero(trace >= threshold)).size else np.nan
        for trace in traces
    ])
    valid = np.isfinite(observed) & np.isfinite(predicted)
    result = {
        "run": str(run.resolve()),
        "initial_cloud_flat_to_flat_m": record["configuration"].get(
            "initial_cloud_flat_to_flat_m", 14.0,
        ),
        "initial_cloud_z_max_m": record["configuration"].get(
            "initial_cloud_z_max_m", 13.0,
        ),
        "threshold_vol_percent": float(threshold),
        "sensor_count": int(labels.size),
        "observed_arrival_available_count": int(np.count_nonzero(np.isfinite(observed))),
        "model_detected_among_available_count": int(np.count_nonzero(valid)),
        "mean_signed_error_s": None if not np.any(valid) else float(np.mean(predicted[valid] - observed[valid])),
        "mean_absolute_error_s": None if not np.any(valid) else float(np.mean(np.abs(predicted[valid] - observed[valid]))),
        "early_detection_count": int(np.count_nonzero(valid & (predicted < observed))),
        "maximum_receptor_vol_percent": float(np.max(traces)),
        "diagnostic_only": True,
        "limitation": (
            "arrival labels are an observation operator, not measured 0--10 s "
            "concentration traces; do not report these values as MAPE"
        ),
    }
    if relative_peak_fraction is not None:
        fraction = float(relative_peak_fraction)
        if not 0.0 < fraction < 1.0:
            raise ValueError("relative_peak_fraction must lie in (0, 1)")
        peaks = np.max(traces, axis=1)
        relative_predicted = np.asarray([
            time_s[hit[0]] if (hit := np.flatnonzero(
                trace >= fraction * peak,
            )).size and peak > 0.0 else np.nan
            for trace, peak in zip(traces, peaks)
        ])
        relative_valid = np.isfinite(observed) & np.isfinite(relative_predicted)
        result.update({
            "relative_peak_fraction": fraction,
            "relative_model_detected_among_available_count": int(
                np.count_nonzero(relative_valid)
            ),
            "relative_mean_signed_error_s": (
                None if not np.any(relative_valid) else float(
                    np.mean(relative_predicted[relative_valid] - observed[relative_valid])
                )
            ),
            "relative_mean_absolute_error_s": (
                None if not np.any(relative_valid) else float(
                    np.mean(np.abs(relative_predicted[relative_valid] - observed[relative_valid]))
                )
            ),
            "relative_early_detection_count": int(
                np.count_nonzero(relative_valid & (relative_predicted < observed))
            ),
            "relative_operator_status": "diagnostic_only_truncated_peak_window",
        })
    if diagnostic_time_origin_offset_s is not None:
        offset = float(diagnostic_time_origin_offset_s)
        result.update({
            "diagnostic_time_origin_offset_s": offset,
            "time_origin_corrected_mean_absolute_error_s": (
                None if not np.any(valid) else float(
                    np.mean(np.abs(predicted[valid] + offset - observed[valid]))
                )
            ),
            "time_origin_correction_status": "observation_operator_sensitivity_only",
        })
    if field_evidence_manifest is not None:
        alignment = _manifest_clock_alignment(field_evidence_manifest)
        result["field_evidence_manifest"] = str(field_evidence_manifest.resolve())
        result["manifest_clock_alignment"] = alignment
        if alignment["status"] == "accepted":
            offset = float(alignment["offset_s"])
            result.update({
                "manifest_clock_offset_s": offset,
                "manifest_clock_corrected_mean_absolute_error_s": (
                    None if not np.any(valid) else float(
                        np.mean(np.abs(predicted[valid] + offset - observed[valid]))
                    )
                ),
                "manifest_clock_correction_status": (
                    "observation_operator_correction_applied"
                ),
            })
        else:
            result["manifest_clock_correction_status"] = "withheld"
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--run-directory", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold-vol-percent", type=float, default=0.1)
    parser.add_argument(
        "--diagnostic-time-origin-offset-s", type=float,
        help=(
            "optional observation-clock sensitivity; never changes the solver "
            "field or constitutes validation"
        ),
    )
    parser.add_argument(
        "--field-evidence-manifest", type=Path,
        help=(
            "optional verified FieldEvidenceManifest; applies its accepted "
            "common-clock offset only when quantitative readiness is accepted"
        ),
    )
    parser.add_argument(
        "--relative-peak-fraction", type=float,
        help=(
            "optional diagnostic model operator: first crossing of this "
            "fraction of each model trace's local peak"
        ),
    )
    args = parser.parse_args()
    report = {
        "schema": "degali.thorney-ti021-initial-cloud-sweep.v1",
        "status": "diagnostic_only",
        "runs": [score_run(
            args.workbook, run, threshold=args.threshold_vol_percent,
            diagnostic_time_origin_offset_s=args.diagnostic_time_origin_offset_s,
            field_evidence_manifest=args.field_evidence_manifest,
            relative_peak_fraction=args.relative_peak_fraction,
        )
                 for run in args.run_directory],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
