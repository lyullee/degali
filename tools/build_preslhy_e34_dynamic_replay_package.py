"""Build a concentration-sealed PRESLHY E3.4 dynamic replay input package.

The source interval must already be frozen by
``freeze_preslhy_e34_source_windows.py``.  This command verifies the source and
dynamic evidence hashes, rotates Young U/V into a window-mean instrument wind
frame using wind only, and reduces the point measurement to predeclared 10 s
block means suitable for a whole-inlet mean-flow boundary.  H2 values are not
exported or summarized; only usable channel names, positions, clocks and row
counts are exposed until the numerical configuration is frozen.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from degali.addons.field_evidence_manifest import (
    field_evidence_manifest_from_mapping,
)
from degali.addons.resolved_obstacle_flow_3d import (
    ResolvedInletVelocityHistory3D,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule
from degali.validation.preslhy_e34 import (
    read_e34_dynamic_wind_hydrogen_history,
    read_e34_pool_history,
)


SCHEMA = "degali.preslhy-e34-dynamic-replay-input.v1"
WIND_REDUCTION_ID = "young-uvw-window-mean-horizontal-frame-10s-block-v1"
YOUNG_MANUAL_URL = (
    "https://www.youngusa.com/wp-content/uploads/2008/01/81000-9028I29.pdf"
)
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


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")
    return value


def _trial(record: dict[str, Any], trial_id: str) -> dict[str, Any]:
    matches = [item for item in record["trials"] if item["trial_id"] == trial_id]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one record for trial {trial_id!r}")
    return matches[0]


def _block_mean_history(
    time_s: np.ndarray,
    u_m_s: np.ndarray,
    v_m_s: np.ndarray,
    w_m_s: np.ndarray,
    *,
    start_s: float,
    end_s: float,
    block_duration_s: float,
) -> tuple[ResolvedInletVelocityHistory3D, dict[str, Any]]:
    if end_s <= start_s:
        raise ValueError("dynamic replay window must have positive duration")
    if block_duration_s <= 0.0 or not math.isfinite(block_duration_s):
        raise ValueError("block_duration_s must be positive and finite")
    selected = (time_s >= start_s) & (time_s <= end_s)
    if np.count_nonzero(selected) < 2:
        raise ValueError("dynamic replay window contains fewer than two wind samples")
    uv = np.column_stack((u_m_s[selected], v_m_s[selected]))
    horizontal_mean = np.mean(uv, axis=0)
    horizontal_speed = float(np.linalg.norm(horizontal_mean))
    if horizontal_speed <= 0.0:
        raise ValueError("window-mean horizontal wind vector has zero magnitude")
    stream_axis = horizontal_mean / horizontal_speed
    cross_axis = np.asarray((-stream_axis[1], stream_axis[0]))

    relative = time_s[selected] - start_s
    stream = uv @ stream_axis
    cross = uv @ cross_axis
    vertical = w_m_s[selected]
    duration = end_s - start_s
    edges = np.arange(0.0, duration, block_duration_s)
    if edges.size == 0 or edges[0] != 0.0:
        edges = np.insert(edges, 0, 0.0)
    if not math.isclose(float(edges[-1]), duration, abs_tol=1.0e-12):
        edges = np.append(edges, duration)

    block_u: list[float] = []
    block_v: list[float] = []
    block_w: list[float] = []
    block_counts: list[int] = []
    for index, (left, right) in enumerate(zip(edges, edges[1:])):
        mask = (relative >= left) & (
            relative <= right if index == len(edges) - 2 else relative < right
        )
        count = int(np.count_nonzero(mask))
        if count == 0:
            raise ValueError("one or more wind-reduction blocks contain no samples")
        block_counts.append(count)
        block_u.append(float(np.mean(stream[mask])))
        block_v.append(float(np.mean(cross[mask])))
        block_w.append(float(np.mean(vertical[mask])))

    # Piecewise-constant schedules retain the last value as an explicit end
    # node, mirroring SourceRateSchedule's finite endpoint convention.
    history = ResolvedInletVelocityHistory3D(
        time_s=tuple(float(value) for value in edges),
        u_m_s=tuple(block_u + [block_u[-1]]),
        v_m_s=tuple(block_v + [block_v[-1]]),
        w_m_s=tuple(block_w + [block_w[-1]]),
        history_id=WIND_REDUCTION_ID,
        interpolation="piecewise_constant",
    )
    diagnostics = {
        "instrument_horizontal_stream_axis_uv": [
            float(stream_axis[0]), float(stream_axis[1]),
        ],
        "instrument_horizontal_cross_axis_uv": [
            float(cross_axis[0]), float(cross_axis[1]),
        ],
        "site_azimuth_available": False,
        "vertical_axis": "Young_W",
        "raw_streamwise_minimum_m_s": float(np.min(stream)),
        "raw_streamwise_nonpositive_fraction": float(np.mean(stream <= 0.0)),
        "block_streamwise_minimum_m_s": float(min(block_u)),
        "block_streamwise_nonpositive_count": int(np.count_nonzero(
            np.asarray(block_u) <= 0.0
        )),
        "block_sample_counts": block_counts,
        "block_duration_s": block_duration_s,
    }
    return history, diagnostics


def build_replay_input_package(
    source_freeze_path: str | Path,
    dynamic_audit_path: str | Path,
    *,
    trial_id: str,
    window_ordinal: int,
    output_directory: str | Path,
    wind_block_duration_s: float = 10.0,
) -> dict[str, Any]:
    """Build one hash-bound, observation-sealed dynamic replay package."""

    source_path = Path(source_freeze_path).resolve()
    audit_path = Path(dynamic_audit_path).resolve()
    output = Path(output_directory).resolve()
    source_record = _load_json(source_path)
    dynamic_record = _load_json(audit_path)
    if source_record.get("schema") != "degali.preslhy-e34-source-window-freeze.v1":
        raise ValueError("unsupported source-window freeze schema")
    if dynamic_record.get("schema") != "degali.preslhy-e34-dynamic-evidence-audit.v1":
        raise ValueError("unsupported dynamic evidence audit schema")
    source_trial = _trial(source_record, trial_id)
    dynamic_trial = _trial(dynamic_record, trial_id)
    if window_ordinal < 1 or window_ordinal > source_trial["window_count"]:
        raise ValueError("window_ordinal is outside the frozen source-window set")
    window = source_trial["windows"][window_ordinal - 1]
    if not window.get("source_only_selection", False):
        raise ValueError("source window was not selected by the source-only policy")
    if not dynamic_trial["evidence_decision"]["dynamic_component_replay_allowed"]:
        raise ValueError("dynamic evidence decision withholds this trial")

    manifest = field_evidence_manifest_from_mapping(
        dynamic_trial["field_evidence_manifest"]
    )
    manifest.verify_files()
    workbook = Path(source_trial["workbook_path"]).resolve()
    workbook_sha = _sha256(workbook)
    if workbook_sha != source_trial["workbook_sha256"]:
        raise ValueError("source-window workbook fingerprint changed")
    if workbook_sha != manifest.observed_dataset_sha256:
        raise ValueError("source and dynamic packages do not bind the same workbook")

    dynamic = read_e34_dynamic_wind_hydrogen_history(workbook)
    if dynamic.wind_quality != "usable":
        raise ValueError("workbook flags the wind history defective")
    start_s = float(window["start_s"])
    end_s = float(window["end_s"])
    inlet_history, wind_diagnostics = _block_mean_history(
        dynamic.synchronized_time_s,
        dynamic.wind_u_m_s,
        dynamic.wind_v_m_s,
        dynamic.wind_w_m_s,
        start_s=start_s,
        end_s=end_s,
        block_duration_s=wind_block_duration_s,
    )
    usable_sensors = tuple(
        label for label, status in dynamic.sensor_quality.items()
        if status == "usable"
    )
    observed_mask = (
        (dynamic.hydrogen_sample_time_s >= start_s)
        & (dynamic.hydrogen_sample_time_s <= end_s)
    )
    if np.count_nonzero(observed_mask) < 2:
        raise ValueError("frozen physical observation window has too few H2 rows")

    duration = end_s - start_s
    source_schedule = SourceRateSchedule(
        time_s=(0.0, duration),
        rate_kg_s=(float(window["fitted_evaporation_rate_kg_s"]), 0.0),
        source_id=window["window_id"] + ":linear-mass-loss-proxy",
        rate_operator="piecewise_constant",
    )
    frozen_mass_loss = float(window["smoothed_mass_lost_kg"])
    source_mass_relative_difference = abs(
        source_schedule.released_mass_kg - frozen_mass_loss
    ) / max(frozen_mass_loss, 1.0e-30)
    source_proxy_mass_closure_passed = bool(
        source_mass_relative_difference <= 0.10
    )
    if output.exists():
        raise FileExistsError(output)
    component_replay_allowed = bool(
        wind_diagnostics["block_streamwise_nonpositive_count"] == 0
        and usable_sensors
        and source_proxy_mass_closure_passed
    )
    record = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "trial_id": trial_id,
        "window_ordinal": window_ordinal,
        "window_id": window["window_id"],
        "source_freeze": {
            "path": str(source_path),
            "sha256": _sha256(source_path),
            "policy_id": source_record["policy_id"],
            "selection_frozen_before_model_scoring": source_record[
                "selection_frozen_before_model_scoring"
            ],
        },
        "dynamic_evidence_audit": {
            "path": str(audit_path),
            "sha256": _sha256(audit_path),
        },
        "field_evidence_manifest_id": manifest.manifest_id,
        "workbook_path": str(workbook),
        "workbook_sha256": workbook_sha,
        "source_window_synchronized_clock_s": [start_s, end_s],
        "duration_s": duration,
        "source_schedule": asdict(source_schedule),
        "source_schedule_status": "conditional_liquid_mass_loss_proxy",
        "source_schedule_mass_kg": source_schedule.released_mass_kg,
        "frozen_smoothed_mass_loss_kg": frozen_mass_loss,
        "source_schedule_mass_relative_difference": (
            source_mass_relative_difference
        ),
        "source_proxy_mass_closure_tolerance": 0.10,
        "source_proxy_mass_closure_passed": source_proxy_mass_closure_passed,
        "inlet_velocity_history": asdict(inlet_history),
        "wind_reduction": {
            "id": WIND_REDUCTION_ID,
            "selection_inputs": ["Young_U", "Young_V", "Young_W"],
            "hydrogen_concentration_used": False,
            "model_prediction_used": False,
            "coordinate_frame": (
                "instrument_window_mean_horizontal_wind_aligned_not_site_registered"
            ),
            "manufacturer_coordinate_convention_source": YOUNG_MANUAL_URL,
            "manufacturer_convention": (
                "U/V horizontal orthogonal components and W vertical when "
                "the instrument is installed vertically"
            ),
            "experiment_installation_azimuth_documented": False,
            "whole_inlet_assumption": (
                "one point-history block mean times the declared spatial inlet shape"
            ),
            **wind_diagnostics,
        },
        "sealed_observation_operator": {
            "values_exported": False,
            "values_summarized": False,
            "sealed_until_model_configuration_is_frozen": True,
            "physical_observation_clock": (
                "workbook synchronized time minus declared sample-tube delay"
            ),
            "declared_sample_transport_delay_s": (
                dynamic.declared_sample_transport_delay_s
            ),
            "selected_row_count": int(np.count_nonzero(observed_mask)),
            "usable_channels": list(usable_sensors),
            "excluded_channels": [
                label for label, status in dynamic.sensor_quality.items()
                if status != "usable"
            ],
            "positions_pool_centre_m": {
                label: list(SENSOR_POSITIONS_POOL_CENTRE_M[label])
                for label in usable_sensors
            },
        },
        "evidence_decision": {
            "status": "conditional" if component_replay_allowed else "withheld",
            "dynamic_component_replay_input_allowed": component_replay_allowed,
            "obstacle_transport_model_qualification_allowed": False,
            "reasons": [
                "site_azimuth_and_machine_registered_obstacle_geometry_are_unavailable",
                "liquid_mass_loss_is_not_a_direct_atmospheric_vapour_flux_measurement",
                "point_wind_block_mean_is_applied_over_the_complete_inlet_plane",
                "sensor_evidence_is_nominal_specification_not_event_calibration",
            ] + (
                ["one_or_more_10s_blocks_reverse_the_xmin_normal_flow"]
                if wind_diagnostics["block_streamwise_nonpositive_count"]
                else []
            ) + (
                ["linear_source_proxy_misses_frozen_mass_loss_by_more_than_10pct"]
                if not source_proxy_mass_closure_passed else []
            ),
        },
        "review_plot": "dynamic-replay-input-review.png",
    }
    output.mkdir(parents=True)
    plot_replay_input(
        dynamic, window, inlet_history, wind_diagnostics,
        output / record["review_plot"],
    )
    (output / "dynamic-replay-input.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return record


def plot_replay_input(
    dynamic,
    window: dict[str, Any],
    history: ResolvedInletVelocityHistory3D,
    wind_diagnostics: dict[str, Any],
    output_path: str | Path,
) -> None:
    """Plot raw wind and the exported mean-flow boundary without H2 values."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    start_s, end_s = float(window["start_s"]), float(window["end_s"])
    selected = (
        (dynamic.synchronized_time_s >= start_s)
        & (dynamic.synchronized_time_s <= end_s)
    )
    relative = dynamic.synchronized_time_s[selected] - start_s
    uv = np.column_stack((
        dynamic.wind_u_m_s[selected], dynamic.wind_v_m_s[selected],
    ))
    stream_axis = np.asarray(
        wind_diagnostics["instrument_horizontal_stream_axis_uv"]
    )
    cross_axis = np.asarray(
        wind_diagnostics["instrument_horizontal_cross_axis_uv"]
    )
    stream = uv @ stream_axis
    cross = uv @ cross_axis
    figure, axes = plt.subplots(
        3, 1, figsize=(11.5, 8.0), sharex=True, constrained_layout=True,
    )
    channels = (
        (stream, history.u_m_s, "streamwise", "#2563eb"),
        (cross, history.v_m_s, "cross-stream", "#7c3aed"),
        (dynamic.wind_w_m_s[selected], history.w_m_s, "vertical W", "#059669"),
    )
    for axis, (raw, reduced, label, colour) in zip(axes, channels):
        axis.plot(relative, raw, color="#94a3b8", linewidth=0.55, label="raw 10 Hz")
        axis.step(
            history.time_s, reduced, where="post", color=colour,
            linewidth=1.6, label="10 s whole-inlet mean boundary",
        )
        axis.axhline(0.0, color="black", linewidth=0.7)
        axis.set_ylabel(label + " (m/s)")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, loc="upper right")
    axes[-1].set_xlabel("time from frozen source-window start (s)")
    figure.suptitle(
        "PRESLHY E3.4 dynamic replay input — H2 observations sealed",
        fontsize=13,
    )
    figure.savefig(Path(output_path), dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_freeze", type=Path)
    parser.add_argument("dynamic_audit", type=Path)
    parser.add_argument("--trial-id", required=True)
    parser.add_argument("--window-ordinal", type=int, required=True)
    parser.add_argument("--wind-block-duration-s", type=float, default=10.0)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record = build_replay_input_package(
        args.source_freeze,
        args.dynamic_audit,
        trial_id=args.trial_id,
        window_ordinal=args.window_ordinal,
        output_directory=args.output_directory,
        wind_block_duration_s=args.wind_block_duration_s,
    )
    print(json.dumps({
        "schema": record["schema"],
        "trial_id": record["trial_id"],
        "window_id": record["window_id"],
        "duration_s": record["duration_s"],
        "dynamic_component_replay_input_allowed": record[
            "evidence_decision"
        ]["dynamic_component_replay_input_allowed"],
        "observation_values_exported": record[
            "sealed_observation_operator"
        ]["values_exported"],
        "output_directory": str(args.output_directory.resolve()),
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
