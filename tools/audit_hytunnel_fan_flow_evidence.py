"""Audit what the public HyTunnel fan-flow workbook can actually support.

The workbook supplies spatial one-minute average/minimum/maximum velocity
matrices.  Those extrema are not a velocity time series and cannot uniquely
define RMS turbulence, intensity, integral scale, spectrum or dissipation.
This audit permits a mean-profile candidate while fail-closing any synthetic
turbulence qualification.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from degali.validation.hytunnel_ps import (
    DATASET_DOI,
    read_hytunnel_flow_workbook,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKBOOK = (
    ROOT
    / "reference/datasets/hytunnel-ps-4668554/raw/HTE2440PS005FLMT00300321.xlsx"
)
DEFAULT_OUTPUT = ROOT / "outputs/hytunnel-fan-flow-evidence-2026-10-11-v1.json"
SCHEMA = "degali.hytunnel-fan-flow-evidence.v1"


def audit_fan_flow_evidence(
    workbook_path: str | Path = DEFAULT_WORKBOOK,
) -> dict[str, Any]:
    source = Path(workbook_path)
    dataset = read_hytunnel_flow_workbook(source)
    selected_id = dataset.configuration_ids[0]
    selected_points = tuple(
        point for point in dataset.points
        if point.configuration_id == selected_id
    )
    first_x = min(
        point.x_m for point in selected_points if point.statistic == "average"
    )
    first_average_plane = sorted(
        (
            point for point in selected_points
            if point.statistic == "average"
            and abs(point.x_m - first_x) <= 1.0e-12
        ),
        key=lambda point: point.y_m,
    )
    velocities = np.asarray(
        [point.velocity_x_m_s for point in first_average_plane], dtype=float,
    )
    statistics = sorted({point.statistic for point in selected_points})
    return {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "conditional_mean_profile_only_turbulence_withheld",
        "promotion_allowed": False,
        "mean_profile_candidate_allowed": True,
        "turbulence_inlet_qualification_allowed": False,
        "dataset": {
            "doi": DATASET_DOI,
            "workbook": source.name,
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "configuration_count": len(dataset.configuration_ids),
            "point_count": len(dataset.points),
            "selected_configuration_id": selected_id,
            "available_spatial_statistics": statistics,
        },
        "selected_first_average_plane": {
            "x_m": first_x,
            "relative_y_m": [float(point.y_m) for point in first_average_plane],
            "velocity_x_m_s": velocities.tolist(),
            "minimum_velocity_m_s": float(np.min(velocities)),
            "mean_velocity_m_s": float(np.mean(velocities)),
            "maximum_velocity_m_s": float(np.max(velocities)),
        },
        "missing_required_turbulence_evidence": [
            "velocity_time_series",
            "rms_velocity_or_reynolds_stress",
            "turbulence_intensity",
            "integral_length_scale",
            "energy_spectrum",
            "dissipation_rate",
        ],
        "decision_basis": (
            "spatial one-minute minima and maxima do not determine temporal RMS "
            "or turbulence scales; only the independently measured mean profile "
            "may be transferred without inventing closure inputs"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    record = audit_fan_flow_evidence(args.workbook)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": record["status"],
        "selected_first_average_plane": record[
            "selected_first_average_plane"
        ],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
