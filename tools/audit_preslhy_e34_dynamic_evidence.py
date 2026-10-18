"""Build fail-closed E3.4 dynamic wind/H2 evidence manifests.

The public workbooks contain synchronized 10 Hz sonic velocity, three local
hydrogen channels and corrected pool mass.  The accompanying report supplies
layout and nominal instrument specifications.  This tool fingerprints those
files, creates one ``FieldEvidenceManifest`` per wind trial, and records what
still prevents obstacle-transport model qualification.  It exports summaries,
not raw measurement arrays.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from degali.addons.field_evidence_manifest import (
    FieldEvidenceManifest,
    FieldEvidenceManifestArtifact,
)
from degali.validation.preslhy_e34 import (
    read_e34_dynamic_wind_hydrogen_history,
    read_e34_pool_history,
)


SCHEMA = "degali.preslhy-e34-dynamic-evidence-audit.v1"
SENSOR_POSITIONS_POOL_CENTRE_M = {
    "CH35": (0.0, 0.0, 0.35),
    "CH45": (0.0, 0.0, 0.45),
    "CH55": (0.0, 0.0, 0.55),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _statistics(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0 or np.any(~np.isfinite(array)):
        raise ValueError("dynamic evidence statistic requires a finite vector")
    return {
        "minimum": float(np.min(array)),
        "mean": float(np.mean(array)),
        "standard_deviation": float(np.std(array)),
        "maximum": float(np.max(array)),
    }


def audit_dynamic_trial(
    workbook_path: str | Path,
    report_path: str | Path,
    *,
    audit_root: str | Path | None = None,
) -> dict[str, Any]:
    workbook = Path(workbook_path).resolve()
    report = Path(report_path).resolve()
    root = Path(audit_root).resolve() if audit_root is not None else workbook.parent
    if not workbook.is_file():
        raise FileNotFoundError(workbook)
    if not report.is_file():
        raise FileNotFoundError(report)
    dynamic = read_e34_dynamic_wind_hydrogen_history(workbook)
    pool = read_e34_pool_history(workbook)
    trial_id = workbook.stem.removesuffix("-Final")
    workbook_sha = _sha256(workbook)
    report_sha = _sha256(report)
    artifact_paths = {
        "source_boundary": (workbook, workbook_sha),
        "weather": (workbook, workbook_sha),
        "obstacle_geometry": (report, report_sha),
        "receptor_observations": (workbook, workbook_sha),
        "common_clock": (workbook, workbook_sha),
    }
    manifest = FieldEvidenceManifest(
        manifest_id=f"preslhy-e34-{trial_id}-dynamic-v1",
        event_id=trial_id,
        audit_root=str(root),
        channel_artifacts=tuple(
            FieldEvidenceManifestArtifact(
                channel=channel, path=str(path), sha256=digest,
            )
            for channel, (path, digest) in artifact_paths.items()
        ),
        dataset_id=f"preslhy-e34:{trial_id}:dynamic-wind-h2",
        observed_dataset_path=str(workbook),
        observed_dataset_sha256=workbook_sha,
        observed_row_count=int(dynamic.original_time_s.size),
        source_boundary_id=(
            f"{trial_id}:corrected-pool-mass-history:vapour-flux-unresolved"
        ),
        weather_id=(
            f"{trial_id}:young-81000-synchronized-10hz:{dynamic.wind_quality}"
        ),
        obstacle_geometry_id=(
            "preslhy-e34-report-figure2:site-obstacles-approximate-not-"
            "machine-registered"
        ),
        receptor_geometry_id=(
            "preslhy-e34-report-figure4:pool-centre-z0p35-z0p45-z0p55"
        ),
        temporal_operator_id=(
            "preslhy-e34-h2-sample-tube-delay-minus-35s-v1"
        ),
        common_clock_id=(
            f"{trial_id}:workbook-synchronized-to-thermocouple-clock-v1"
        ),
        scope="lh2_obstacle_transport",
        sensor_set_id="preslhy-e34-ftc300-ch35-ch45-ch55-v1",
        operator_id="degali-local-evidence-audit-2026-10-11",
        sensor_registry_artifact=FieldEvidenceManifestArtifact(
            channel="sensor_registry", path=str(report), sha256=report_sha,
        ),
        sensor_calibration_status="specification_only",
    )
    manifest.verify_files()

    reasons = [
        "pool_mass_is_measured_but_atmospheric_h2_vapour_flux_is_not_directly_measured",
        "site_obstacles_are_available_only_in_a_report_figure_not_machine_registered_3d_geometry",
        "sensor_accuracy_is_a_nominal_specification_not_an_event_calibration_certificate",
        "hydrogen_sample_transport_delay_is_approximately_not_exactly_35s",
    ]
    defective_sensors = sorted(
        label for label, status in dynamic.sensor_quality.items()
        if status != "usable"
    )
    if defective_sensors:
        reasons.append(
            "workbook_flagged_defective_hydrogen_channels:"
            + ",".join(defective_sensors)
        )
    if dynamic.wind_quality != "usable":
        reasons.append("workbook_flagged_defective_sonic_wind_channels")

    speed = np.sqrt(
        dynamic.wind_u_m_s**2
        + dynamic.wind_v_m_s**2
        + dynamic.wind_w_m_s**2
    )
    usable_hydrogen = tuple(
        label for label, status in dynamic.sensor_quality.items()
        if status == "usable"
    )
    dynamic_replay_allowed = bool(
        dynamic.wind_quality == "usable" and usable_hydrogen
    )
    return {
        "trial_id": trial_id,
        "field_evidence_manifest": manifest.as_record(),
        "evidence_decision": {
            "status": "conditional" if dynamic_replay_allowed else "withheld",
            "dynamic_component_replay_allowed": dynamic_replay_allowed,
            "obstacle_transport_model_qualification_allowed": False,
            "promotion_allowed": False,
            "reasons": reasons,
        },
        "source_history": {
            "sample_count": int(pool.time_s.size),
            "time_range_s": [float(pool.time_s[0]), float(pool.time_s[-1])],
            "corrected_pool_mass_range_kg": [
                float(np.min(pool.liquid_mass_kg)),
                float(np.max(pool.liquid_mass_kg)),
            ],
            "atmospheric_vapour_mass_flow_available": False,
        },
        "dynamic_history": {
            "sample_count": int(dynamic.original_time_s.size),
            "median_sample_interval_s": dynamic.median_sample_interval_s,
            "sampling_frequency_hz": 1.0 / dynamic.median_sample_interval_s,
            "synchronized_time_range_s": [
                float(dynamic.synchronized_time_s[0]),
                float(dynamic.synchronized_time_s[-1]),
            ],
            "hydrogen_sample_time_range_s": [
                float(dynamic.hydrogen_sample_time_s[0]),
                float(dynamic.hydrogen_sample_time_s[-1]),
            ],
            "declared_sample_transport_delay_s": (
                dynamic.declared_sample_transport_delay_s
            ),
            "sensor_positions_pool_centre_m": {
                label: list(position)
                for label, position in SENSOR_POSITIONS_POOL_CENTRE_M.items()
            },
            "sensor_quality": dynamic.sensor_quality,
            "usable_hydrogen_channels": list(usable_hydrogen),
            "hydrogen_vol_percent_statistics_full_record": {
                label: _statistics(values)
                for label, values in dynamic.hydrogen_vol_percent.items()
            },
            "wind_quality": dynamic.wind_quality,
            "wind_statistics_full_record": {
                "u_m_s": _statistics(dynamic.wind_u_m_s),
                "v_m_s": _statistics(dynamic.wind_v_m_s),
                "w_m_s": _statistics(dynamic.wind_w_m_s),
                "speed_m_s": _statistics(speed),
            },
            "full_record_statistics_are_validation_scores": False,
        },
        "instrument_evidence": {
            "hydrogen_sensor": "Messkonzept FTC300 0-100 vol% H2",
            "hydrogen_sensor_nominal_non_linearity": "<=1% full scale at 290 K",
            "cryogenic_accuracy_available": False,
            "sonic_anemometer": "Young Model 81000 0-40 m/s",
            "sonic_velocity_nominal_accuracy_m_s": 0.05,
            "sampling_frequency_hz": 10.0,
            "layout_source": "PRESLHY D3.5 E3.4 report Figures 2 and 4",
        },
    }


def build_dynamic_evidence_package(
    data_directory: str | Path,
    report_path: str | Path,
    output_directory: str | Path,
) -> dict[str, Any]:
    data = Path(data_directory).resolve()
    report = Path(report_path).resolve()
    output = Path(output_directory).resolve()
    workbooks = sorted(data.glob("*Wind-Final.xlsx"))
    if not workbooks:
        raise ValueError("no E3.4 wind workbooks found")
    trials = [
        audit_dynamic_trial(path, report, audit_root=data)
        for path in workbooks
    ]
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    manifest_files: list[str] = []
    for trial in trials:
        destination = output / f"{trial['trial_id']}-field-evidence-manifest.json"
        destination.write_text(
            json.dumps(
                trial["field_evidence_manifest"], indent=2, allow_nan=False,
            ) + "\n",
            encoding="utf-8",
        )
        manifest_files.append(destination.name)
    record = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": "PRESLHY Experiment series 3.4",
        "source_doi": "10.35097/1319",
        "raw_measurements_copied": False,
        "trial_count": len(trials),
        "dynamic_component_replay_allowed_count": sum(
            trial["evidence_decision"]["dynamic_component_replay_allowed"]
            for trial in trials
        ),
        "obstacle_transport_model_qualification_allowed": False,
        "manifest_files": manifest_files,
        "overview_plot": "dynamic-evidence-overview.png",
        "trials": trials,
    }
    plot_dynamic_evidence(
        trials, output / record["overview_plot"],
    )
    (output / "audit.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return record


def plot_dynamic_evidence(
    trials: list[dict[str, Any]], output_path: str | Path,
) -> None:
    """Plot pool mass, measured wind speed and usable H2 channels by trial."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    count = len(trials)
    figure, axes = plt.subplots(
        3, count, figsize=(4.2 * count, 8.4), squeeze=False,
        constrained_layout=True, sharex="col",
    )
    for column, trial in enumerate(trials):
        workbook = Path(
            trial["field_evidence_manifest"]["observed_dataset_path"]
        )
        dynamic = read_e34_dynamic_wind_hydrogen_history(workbook)
        pool = read_e34_pool_history(workbook)
        axes[0, column].plot(pool.time_s, pool.liquid_mass_kg, color="#2563eb")
        axes[0, column].set_title(
            trial["trial_id"] + "\n" + trial["evidence_decision"]["status"]
        )
        speed = np.sqrt(
            dynamic.wind_u_m_s**2
            + dynamic.wind_v_m_s**2
            + dynamic.wind_w_m_s**2
        )
        colour = "#dc2626" if dynamic.wind_quality == "usable" else "#9ca3af"
        axes[1, column].plot(
            dynamic.synchronized_time_s, speed, color=colour, linewidth=0.7,
        )
        for label, sensor_colour in (("CH45", "#059669"), ("CH55", "#7c3aed")):
            if dynamic.sensor_quality[label] == "usable":
                axes[2, column].plot(
                    dynamic.hydrogen_sample_time_s,
                    dynamic.hydrogen_vol_percent[label],
                    label=label, color=sensor_colour, linewidth=0.7,
                )
        axes[2, column].axhline(
            4.0, color="black", linestyle="--", linewidth=0.8,
        )
        axes[2, column].legend(frameon=False, fontsize=8)
        axes[2, column].set_xlabel("synchronized event time (s)")
        for row in range(3):
            axes[row, column].grid(alpha=0.2)
    axes[0, 0].set_ylabel("corrected pool mass (kg)")
    axes[1, 0].set_ylabel("sonic wind speed (m/s)")
    axes[2, 0].set_ylabel("H2 (vol%)")
    figure.suptitle(
        "PRESLHY E3.4 dynamic evidence — raw clocks preserved, no model fitting",
        fontsize=13,
    )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_directory", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record = build_dynamic_evidence_package(
        args.data_directory, args.report, args.output_directory,
    )
    print(json.dumps({
        "schema": record["schema"],
        "trial_count": record["trial_count"],
        "dynamic_component_replay_allowed_count": record[
            "dynamic_component_replay_allowed_count"
        ],
        "obstacle_transport_model_qualification_allowed": record[
            "obstacle_transport_model_qualification_allowed"
        ],
        "output_directory": str(args.output_directory.resolve()),
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
