"""Build a mass/wind-only fill-prehistory diagnostic after a failed primary.

The primary score is never replaced.  A fixed 10% smoothed-mass excursion
locates fill onset before the selected post-fill decay.  The already frozen
post-fill evaporation rate is extended backward as an explicit diagnostic
assumption; positive liquid filling is never interpreted as measured vapour
flow.  Only corrected pool mass and Young U/V/W are opened by this command.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from degali.addons.resolved_obstacle_flow_3d import (  # noqa: E402
    ResolvedInletVelocityHistory3D,
)
from degali.addons.semi_fv_obstacle import SourceRateSchedule  # noqa: E402
from degali.validation.preslhy_e34 import (  # noqa: E402
    E34FillPrehistoryPolicy,
    E34SourceWindowPolicy,
    read_e34_pool_history,
    read_e34_wind_history,
    source_only_fill_prehistory_candidate,
)


SCHEMA = "degali.preslhy-e34-fill-prehistory-diagnostic-input.v1"
POLICY_ID = "preslhy-e34-mass-only-10pct-fill-onset-v1"


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


def _fixed_frame_block_history(
    time_s: np.ndarray,
    wind_u_m_s: np.ndarray,
    wind_v_m_s: np.ndarray,
    wind_w_m_s: np.ndarray,
    *,
    start_s: float,
    end_s: float,
    stream_axis_uv: tuple[float, float],
    cross_axis_uv: tuple[float, float],
    block_duration_s: float,
) -> tuple[ResolvedInletVelocityHistory3D, dict[str, Any]]:
    selected = (time_s >= start_s) & (time_s <= end_s)
    if np.count_nonzero(selected) < 2:
        raise ValueError("prehistory interval contains fewer than two wind samples")
    relative = time_s[selected] - start_s
    uv = np.column_stack((wind_u_m_s[selected], wind_v_m_s[selected]))
    stream = uv @ np.asarray(stream_axis_uv)
    cross = uv @ np.asarray(cross_axis_uv)
    vertical = wind_w_m_s[selected]
    duration = end_s - start_s
    edges = np.arange(0.0, duration, block_duration_s)
    if not edges.size or edges[0] != 0.0:
        edges = np.insert(edges, 0, 0.0)
    if not math.isclose(float(edges[-1]), duration, abs_tol=1.0e-12):
        edges = np.append(edges, duration)
    block_u: list[float] = []
    block_v: list[float] = []
    block_w: list[float] = []
    counts: list[int] = []
    for index, (left, right) in enumerate(zip(edges, edges[1:])):
        mask = (relative >= left) & (
            relative <= right if index == edges.size - 2 else relative < right
        )
        count = int(np.count_nonzero(mask))
        if count == 0:
            raise ValueError("prehistory wind block contains no samples")
        counts.append(count)
        block_u.append(float(np.mean(stream[mask])))
        block_v.append(float(np.mean(cross[mask])))
        block_w.append(float(np.mean(vertical[mask])))
    history = ResolvedInletVelocityHistory3D(
        time_s=tuple(float(value) for value in edges),
        u_m_s=tuple(block_u + [block_u[-1]]),
        v_m_s=tuple(block_v + [block_v[-1]]),
        w_m_s=tuple(block_w + [block_w[-1]]),
        history_id="young-uvw-primary-frame-fill-prehistory-10s-block-v1",
        interpolation="piecewise_constant",
    )
    return history, {
        "raw_streamwise_minimum_m_s": float(np.min(stream)),
        "raw_streamwise_nonpositive_fraction": float(np.mean(stream <= 0.0)),
        "block_streamwise_minimum_m_s": float(min(block_u)),
        "block_streamwise_nonpositive_count": int(np.count_nonzero(
            np.asarray(block_u) <= 0.0
        )),
        "block_sample_counts": counts,
    }


def build_fill_prehistory_diagnostic(
    primary_replay_input_path: str | Path,
    primary_score_path: str | Path,
    output_directory: str | Path,
    *,
    relative_fill_onset_fraction: float = 0.10,
    wind_block_duration_s: float = 10.0,
) -> dict[str, Any]:
    replay_path = Path(primary_replay_input_path).resolve()
    score_path = Path(primary_score_path).resolve()
    output = Path(output_directory).resolve()
    if output.exists():
        raise FileExistsError(output)
    replay = _load_json(replay_path)
    score = _load_json(score_path)
    if replay.get("schema") != "degali.preslhy-e34-dynamic-replay-input.v1":
        raise ValueError("unsupported primary replay schema")
    if score.get("schema") != "degali.preslhy-e34-dynamic-component-score.v1":
        raise ValueError("unsupported primary score schema")
    if score.get("all_predeclared_dynamic_accuracy_targets_passed") is not False:
        raise ValueError("fill-prehistory diagnosis requires a failed primary score")
    if score.get("observation_control", {}).get(
        "observation_used_for_parameter_selection"
    ) is not False:
        raise ValueError("primary score indicates observation-based model selection")
    if score["configuration_freeze"]["sha256"] != _sha256(
        Path(score["configuration_freeze"]["path"])
    ):
        raise ValueError("primary score configuration fingerprint changed")

    source_freeze_path = Path(replay["source_freeze"]["path"]).resolve()
    if _sha256(source_freeze_path) != replay["source_freeze"]["sha256"]:
        raise ValueError("source-window freeze fingerprint changed")
    source_freeze = _load_json(source_freeze_path)
    source_policy = E34SourceWindowPolicy(**source_freeze["policy"])
    workbook = Path(replay["workbook_path"]).resolve()
    if _sha256(workbook) != replay["workbook_sha256"]:
        raise ValueError("raw workbook fingerprint changed")
    pool = read_e34_pool_history(workbook)
    wind = read_e34_wind_history(workbook)
    if wind.wind_quality != "usable":
        raise ValueError("workbook flags prehistory wind defective")
    prehistory_policy = E34FillPrehistoryPolicy(
        smoothing_sigma_s=source_policy.smoothing_sigma_s,
        relative_fill_onset_fraction=relative_fill_onset_fraction,
        minimum_prehistory_s=10.0,
    )
    candidate = source_only_fill_prehistory_candidate(
        pool,
        window_ordinal=int(replay["window_ordinal"]),
        source_window_policy=source_policy,
        prehistory_policy=prehistory_policy,
    )
    primary_start, primary_end = replay["source_window_synchronized_clock_s"]
    if not math.isclose(candidate.score_start_s, primary_start, abs_tol=1.0e-8):
        raise ValueError("mass-only candidate score start differs from primary")
    if not math.isclose(candidate.score_end_s, primary_end, abs_tol=1.0e-8):
        raise ValueError("mass-only candidate score end differs from primary")
    stream_axis = tuple(
        float(value) for value in replay["wind_reduction"][
            "instrument_horizontal_stream_axis_uv"
        ]
    )
    cross_axis = tuple(
        float(value) for value in replay["wind_reduction"][
            "instrument_horizontal_cross_axis_uv"
        ]
    )
    history, wind_diagnostics = _fixed_frame_block_history(
        wind.synchronized_time_s,
        wind.wind_u_m_s,
        wind.wind_v_m_s,
        wind.wind_w_m_s,
        start_s=candidate.prehistory_start_s,
        end_s=candidate.score_end_s,
        stream_axis_uv=stream_axis,
        cross_axis_uv=cross_axis,
        block_duration_s=wind_block_duration_s,
    )
    primary_rate = float(replay["source_schedule"]["rate_kg_s"][0])
    source_schedule = SourceRateSchedule(
        time_s=(0.0, candidate.full_duration_s),
        rate_kg_s=(primary_rate, 0.0),
        source_id=replay["source_schedule"]["source_id"] + ":fill-prehistory-diagnostic",
        rate_operator="piecewise_constant",
    )
    allowed = wind_diagnostics["block_streamwise_nonpositive_count"] == 0
    record = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy_id": POLICY_ID,
        "policy": asdict(prehistory_policy),
        "primary_replay": {"path": str(replay_path), "sha256": _sha256(replay_path)},
        "failed_primary_score": {
            "path": str(score_path),
            "sha256": _sha256(score_path),
            "retained_as_primary": True,
            "status": score["evidence_decision"]["status"],
        },
        "workbook": {"path": str(workbook), "sha256": _sha256(workbook)},
        "data_access": {
            "corrected_pool_mass_opened": True,
            "young_uvw_opened": True,
            "hydrogen_concentration_opened_by_builder": False,
            "hydrogen_concentration_used_to_select_onset": False,
            "model_prediction_used_to_select_onset": False,
            "primary_observations_had_already_been_scored_before_this_diagnostic": True,
        },
        "mass_only_prehistory": asdict(candidate),
        "model_clock": {
            "prehistory_start_synchronized_s": candidate.prehistory_start_s,
            "score_start_model_time_s": candidate.prehistory_duration_s,
            "score_end_model_time_s": candidate.full_duration_s,
            "score_observation_clock_synchronized_s": [primary_start, primary_end],
        },
        "source_schedule": asdict(source_schedule),
        "source_assumption": {
            "prehistory_rate_kg_s": primary_rate,
            "basis": "extend_the_frozen_postfill_rate_backward_without_h2_fitting",
            "prehistory_assumed_released_mass_kg": (
                primary_rate * candidate.prehistory_duration_s
            ),
            "score_window_released_mass_kg": primary_rate * candidate.score_duration_s,
            "full_released_mass_kg": source_schedule.released_mass_kg,
            "prehistory_mass_closure_available": False,
            "atmospheric_vapour_flux_qualified": False,
        },
        "inlet_velocity_history": asdict(history),
        "wind_reduction": {
            "coordinate_frame": replay["wind_reduction"]["coordinate_frame"],
            "axes_reused_from_failed_primary": True,
            "stream_axis_uv": list(stream_axis),
            "cross_axis_uv": list(cross_axis),
            "block_duration_s": wind_block_duration_s,
            **wind_diagnostics,
        },
        "sealed_observation_operator_metadata": replay["sealed_observation_operator"],
        "evidence_decision": {
            "status": (
                "diagnostic_input_allowed"
                if allowed else "withheld_nonpositive_block_inflow"
            ),
            "diagnostic_run_allowed": allowed,
            "eligible_as_independent_validation": False,
            "may_replace_failed_primary_score": False,
            "obstacle_transport_model_qualification_allowed": False,
            "quantitative_facility_prediction_allowed": False,
        },
        "review_plot": "fill-prehistory-diagnostic-input.png",
    }
    output.mkdir(parents=True)
    plot_path = output / record["review_plot"]
    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=False)
    axes[0].plot(pool.time_s, pool.liquid_mass_kg, lw=0.8, label="corrected pool mass")
    axes[0].axvspan(
        candidate.prehistory_start_s, candidate.score_start_s,
        color="tab:blue", alpha=0.18, label="mass-only fill prehistory",
    )
    axes[0].axvspan(
        candidate.score_start_s, candidate.score_end_s,
        color="tab:orange", alpha=0.18, label="retained primary score window",
    )
    axes[0].set_xlim(candidate.prehistory_start_s - 20.0, candidate.score_end_s + 20.0)
    axes[0].set_ylabel("m(LH2) (kg)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].step(history.time_s, history.u_m_s, where="post", label="stream U")
    axes[1].step(history.time_s, history.v_m_s, where="post", label="cross V")
    axes[1].step(history.time_s, history.w_m_s, where="post", label="vertical W")
    axes[1].axvline(candidate.prehistory_duration_s, color="tab:orange", ls="--")
    axes[1].set_xlabel("Diagnostic model time (s)")
    axes[1].set_ylabel("Applied velocity (m/s)")
    axes[1].legend(ncol=3)
    axes[1].grid(alpha=0.25)
    figure.suptitle(
        "PRESLHY E3.4 fill-prehistory diagnostic input — mass/wind only, no H2 fit"
    )
    figure.tight_layout()
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)
    (output / "fill-prehistory-diagnostic-input.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-replay-input", type=Path, required=True)
    parser.add_argument("--primary-score", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    arguments = parser.parse_args(argv)
    record = build_fill_prehistory_diagnostic(
        arguments.primary_replay_input,
        arguments.primary_score,
        arguments.output_directory,
    )
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if record["evidence_decision"]["diagnostic_run_allowed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
