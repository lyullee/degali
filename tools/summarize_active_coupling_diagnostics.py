"""Build the reproducible decision matrix for active pressure coupling."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "obstacle-physics-rebuild-2026-10-10"

CASES = (
    (
        "baseline equilibrium replacement",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85.json",
        True,
        "rejected: final H2/T non-monotonic",
    ),
    (
        "Strang source split",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85_strang.json",
        True,
        "rejected: source ordering is not causal defect",
    ),
    (
        "zero gravity",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85_g0.json",
        True,
        "rejected: buoyancy force is not causal defect",
    ),
    (
        "variable-density projection only",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85_variable_density_only.json",
        True,
        "rejected: removing pre-projection does not cure feedback",
    ),
    (
        "deferred Guo pressure force",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85_deferred_guo_force.json",
        True,
        "rejected: final H2/T non-monotonic",
    ),
    (
        "unified FV momentum/projection",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85_finite_volume_projection.json",
        True,
        "numerical component accepted; external H2 validation still required",
    ),
    (
        "transport-only pressure correction",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85_transport_only.json",
        False,
        "diagnostic pass; rejected for missing pressure-to-momentum feedback",
    ),
    (
        "no mass-flux correction",
        "active_thermal_corrector_time_audit_sc0p7_pr0p85_no_mass_flux.json",
        False,
        "diagnostic pass; rejected for missing variable-density continuity",
    ),
)


def main() -> int:
    records = []
    for label, filename, momentum_feedback, decision in CASES:
        path = OUTPUT / filename
        source = json.loads(path.read_text(encoding="utf-8"))
        ratios = {
            name: float(metric["medium_fine_rms"])
            / float(metric["coarse_medium_rms"])
            for name, metric in source["field_difference_metrics"].items()
        }
        numerical_pass = bool(
            source["numerical_gates"]["successive_field_differences_decrease"]
        )
        records.append({
            "label": label,
            "source_artifact": filename,
            "source_status": source["status"],
            "pressure_to_momentum_feedback": momentum_feedback,
            "successive_difference_ratios_medium_fine_over_coarse_medium": ratios,
            "numerical_time_gate_passed": numerical_pass,
            "physical_coupling_complete": bool(momentum_feedback),
            "promoted": False,
            "decision": decision,
        })

    fixed_path = OUTPUT / "active_time_fixed_carrier_diagnostic.json"
    fixed = json.loads(fixed_path.read_text(encoding="utf-8"))
    fixed_ratios = {
        name: float(metric["medium_fine_rms"])
        / float(metric["coarse_medium_rms"])
        for name, metric in fixed["field_metrics"].items()
    }
    grid_path = OUTPUT / (
        "active_thermal_corrector_grid_audit_cell_integrated_source_"
        "sc0p7_pr0p85_finite_volume_projection.json"
    )
    grid = json.loads(grid_path.read_text(encoding="utf-8"))
    record = {
        "schema": "degali.active-pressure-coupling-decision-matrix.v1",
        "status": "accepted_numerical_candidate_external_validation_required",
        "acceptance_rule": (
            "Every field ratio must be below 1 and pressure-to-momentum feedback "
            "must remain physically complete. Numerical convergence alone cannot "
            "promote a deliberately decoupled diagnostic."
        ),
        "fixed_carrier_corrector_control": {
            "source_artifact": fixed_path.name,
            "status": fixed["status"],
            "successive_difference_ratios_medium_fine_over_coarse_medium": (
                fixed_ratios
            ),
        },
        "unified_fv_grid_gate": {
            "source_artifact": grid_path.name,
            "status": grid["status"],
            "field_difference_metrics": grid["field_difference_metrics"],
            "finest_pair_changes": grid["finest_pair_changes"],
            "gate_reasons": grid["gate_reasons"],
        },
        "coupling_candidates": records,
        "decision": (
            "The scalar/EOS/mass-flux corrector and carrier-flow components "
            "converge separately. Feedback into isothermal BGK fails, whereas "
            "the unified finite-volume variable-density momentum/projection "
            "candidate passes the active component time and grid gates. It is "
            "not promoted until external hydrogen validation is passed."
        ),
    }
    json_path = OUTPUT / "active_pressure_coupling_decision_matrix.json"
    figure_path = OUTPUT / "active_pressure_coupling_decision_matrix.png"
    json_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    metric_order = (
        "maximum_h2_vol_percent", "final_h2_vol_percent", "final_temperature_k",
    )
    labels = [item["label"] for item in records]
    values = np.asarray([
        [
            item["successive_difference_ratios_medium_fine_over_coarse_medium"][name]
            for name in metric_order
        ]
        for item in records
    ])
    figure, axis = plt.subplots(figsize=(12.5, 6.5), constrained_layout=True)
    image = axis.imshow(values, cmap="RdYlGn_r", vmin=0.0, vmax=1.5, aspect="auto")
    axis.set_xticks(np.arange(3), ("maximum H2", "final H2", "final temperature"))
    axis.set_yticks(np.arange(len(labels)), labels)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            axis.text(
                column, row, f"{value:.2f}", ha="center", va="center",
                color="white" if value > 0.95 else "black", fontweight="bold",
            )
    axis.axhline(4.5, color="black", linewidth=2.0)
    axis.axhline(5.5, color="black", linewidth=2.0)
    axis.text(
        2.48, 4.35, "stacked BGK feedback above / unified FV candidate below",
        ha="right", va="bottom", fontsize=9,
    )
    axis.text(
        2.48, 5.35, "unified candidate above / diagnostic controls below",
        ha="right", va="bottom", fontsize=9,
    )
    axis.set_title(
        "Active pressure-coupling time refinement\n"
        "unified FV passes numerical gates; external H2 validation pending"
    )
    colorbar = figure.colorbar(image, ax=axis, label="medium–fine / coarse–medium RMS")
    colorbar.ax.axhline(1.0, color="black", linewidth=1.5)
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(json.dumps({
        "json": str(json_path.resolve()),
        "figure": str(figure_path.resolve()),
        "status": record["status"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
