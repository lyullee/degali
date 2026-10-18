"""Audit the local Thorney Island Trial-021 curved-fence evidence package."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degali.validation.thorney_ti021 import (
    read_thorney_ti021_workbook,
    thorney_ti021_evidence_record,
)


ROOT = Path(__file__).resolve().parents[1]


def _plot(record: dict[str, object], output: Path) -> None:
    observations = record["observations"]
    sensors = observations["dose_sensors"]
    x = np.asarray([item["local_x_m"] for item in sensors], dtype=float)
    y = np.asarray([item["local_y_m"] for item in sensors], dtype=float)
    dose = np.asarray([item["dose_percent_s"] for item in sensors], dtype=float)
    positive = dose > 0.0

    fig, axis = plt.subplots(figsize=(8.2, 7.0), constrained_layout=True)
    axis.scatter(
        x[~positive], y[~positive], s=13, facecolors="none", edgecolors="0.65",
        linewidths=0.7, label="zero reported dose",
    )
    scatter = axis.scatter(
        x[positive], y[positive], c=np.log10(dose[positive]), s=24,
        cmap="viridis", edgecolors="black", linewidths=0.25,
        label="positive dose sensor",
    )
    theta = np.linspace(0.0, np.pi, 361)
    axis.plot(50.0 * np.cos(theta), 50.0 * np.sin(theta), color="crimson",
              linewidth=3.0, label="report: 50 m radius impermeable fence")
    axis.scatter([0.0], [0.0], marker="*", s=150, color="orange",
                 edgecolor="black", label="release centre")
    wind_angle_rad = math.radians(-6.1)
    wind_dx = 60.0 * math.sin(wind_angle_rad)
    wind_dy = 60.0 * math.cos(wind_angle_rad)
    axis.arrow(0.0, 0.0, wind_dx, wind_dy, width=0.8, head_width=5.0,
               length_includes_head=True, color="tab:blue", alpha=0.8)
    axis.text(wind_dx - 5.0, wind_dy + 3.0,
              "report mean wind\n-6.1 deg from +site-y",
              color="tab:blue", ha="center", va="bottom")
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("site x relative to release (m)")
    axis.set_ylabel("site y relative to release (m)")
    axis.set_title("Thorney Island TI021 evidence geometry (not a model result)")
    axis.grid(alpha=0.2)
    axis.legend(loc="upper right", fontsize=8)
    bar = fig.colorbar(scatter, ax=axis, shrink=0.82)
    bar.set_label("log10 reported dose (%·s)")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--curved-wall-geometry-audit", type=Path)
    parser.add_argument("--initial-cloud-audit", type=Path)
    parser.add_argument("--r12-property-audit", type=Path)
    parser.add_argument("--full-scale-preflight", type=Path)
    parser.add_argument("--turbulent-mixing-policy", type=Path)
    parser.add_argument("--full-transient-plan", type=Path)
    parser.add_argument("--nested-handoff-audit", type=Path)
    parser.add_argument("--outer-inlet-coupling-audit", type=Path)
    parser.add_argument("--outer-momentum-overlap-audit", type=Path)
    parser.add_argument("--handoff-streaming-preflight", type=Path)
    parser.add_argument("--coarse-inner-transient", type=Path)
    parser.add_argument("--coarse-inner-score", type=Path)
    parser.add_argument("--completed-inner-transient", type=Path)
    parser.add_argument("--inner-transient-score", type=Path)
    parser.add_argument("--corrector-tolerance-audit", type=Path)
    parser.add_argument("--strict-peak-window-audit", type=Path)
    parser.add_argument("--full-event-grid-comparison", type=Path)
    parser.add_argument(
        "--output-directory", type=Path,
        default=ROOT / "outputs" / "thorney-ti021-evidence-audit-2026-10-11-v1",
    )
    args = parser.parse_args()

    dataset = read_thorney_ti021_workbook(args.workbook)
    record = thorney_ti021_evidence_record(
        dataset,
        report_path=args.report,
        curved_wall_geometry_audit_path=args.curved_wall_geometry_audit,
        initial_cloud_audit_path=args.initial_cloud_audit,
        r12_property_audit_path=args.r12_property_audit,
        full_scale_preflight_path=args.full_scale_preflight,
        turbulent_mixing_policy_path=args.turbulent_mixing_policy,
        full_transient_plan_path=args.full_transient_plan,
        nested_handoff_audit_path=args.nested_handoff_audit,
        outer_inlet_coupling_audit_path=args.outer_inlet_coupling_audit,
        outer_momentum_overlap_audit_path=args.outer_momentum_overlap_audit,
        handoff_streaming_preflight_path=args.handoff_streaming_preflight,
        coarse_inner_transient_path=args.coarse_inner_transient,
        coarse_inner_score_path=args.coarse_inner_score,
        completed_inner_transient_path=args.completed_inner_transient,
        inner_transient_score_path=args.inner_transient_score,
        corrector_tolerance_audit_path=args.corrector_tolerance_audit,
        strict_peak_window_audit_path=args.strict_peak_window_audit,
        full_event_grid_comparison_path=args.full_event_grid_comparison,
    )
    args.output_directory.mkdir(parents=True, exist_ok=True)
    json_path = args.output_directory / "ti021-evidence-audit.json"
    plot_path = args.output_directory / "ti021-evidence-geometry.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _plot(record, plot_path)
    print(json.dumps({
        "status": record["status"],
        "json": str(json_path.resolve()),
        "plot": str(plot_path.resolve()),
        "dose_sensor_count": record["observations"]["dose_sensor_count"],
        "positive_dose_sensor_count": record["observations"][
            "positive_dose_sensor_count"
        ],
        "conflicts": record["conflicts_and_missing_evidence"],
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
