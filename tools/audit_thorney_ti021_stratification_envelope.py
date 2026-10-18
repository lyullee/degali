"""Run every frozen TI021 Richardson regularisation diagnostic member."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

# Keep envelope audits runnable in CI and headless Windows sessions.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from tools.audit_thorney_ti021_stratification_diagnostic import (
        PLAN_SCHEMA,
        audit_stratification,
    )
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from audit_thorney_ti021_stratification_diagnostic import (
        PLAN_SCHEMA,
        audit_stratification,
    )
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-stratification-diagnostic-envelope.v1"


def audit_stratification_envelope(
    *,
    run_directory: Path,
    mixing_policy_path: Path,
    diagnostic_plan_path: Path,
    output_directory: Path,
) -> tuple[dict[str, object], dict[str, Path]]:
    output = output_directory.resolve()
    output.mkdir(parents=True, exist_ok=False)
    plan_path = diagnostic_plan_path.resolve()
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema") != PLAN_SCHEMA or not plan.get("all_gates_pass"):
        raise ValueError("accepted TI021 stratification diagnostic plan is required")

    member_records = []
    member_artifacts = []
    for member in plan["members"]:
        member_id = str(member["member_id"])
        member_output = output / member_id
        record, paths = audit_stratification(
            run_directory=run_directory,
            mixing_policy_path=mixing_policy_path,
            diagnostic_plan_path=plan_path,
            shear_regularization_s_inv=float(
                member["shear_regularization_s_inv"]
            ),
            output_directory=member_output,
        )
        member_records.append(record)
        member_artifacts.append({
            "member_id": member_id,
            "role": member["role"],
            "shear_regularization_s_inv": member[
                "shear_regularization_s_inv"
            ],
            "record_path": str(paths["json"]),
            "record_sha256": _sha256(paths["json"]),
            "figure_path": str(paths["figure"]),
            "figure_sha256": _sha256(paths["figure"]),
        })

    segment_count = len(member_records[0]["segments"])
    if any(len(record["segments"]) != segment_count for record in member_records):
        raise ValueError("diagnostic members do not share a segment count")
    region_names = list(member_records[0]["segments"][0]["regions"])
    receptor_labels = [
        item["label"] for item in member_records[0]["segments"][0]["receptors"]
    ]
    aggregate_segments = []
    for segment_index in range(segment_count):
        member_segments = [
            record["segments"][segment_index] for record in member_records
        ]
        times = [float(item["stop_time_s"]) for item in member_segments]
        if not np.allclose(times, times[0], rtol=0.0, atol=1.0e-12):
            raise ValueError("diagnostic members do not share a model time")
        region_envelopes = {}
        for region_name in region_names:
            region_envelopes[region_name] = {}
            for metric in (
                "gradient_richardson_mean",
                "gradient_richardson_p50",
                "gradient_richardson_p90",
                "gradient_richardson_p95",
                "gradient_richardson_maximum",
                "fraction_at_or_above_0p25",
            ):
                values = np.asarray([
                    item["regions"][region_name][metric]
                    for item in member_segments
                ], dtype=float)
                region_envelopes[region_name][metric] = {
                    "minimum": float(np.min(values)),
                    "maximum": float(np.max(values)),
                    "member_values": {
                        artifact["member_id"]: float(value)
                        for artifact, value in zip(member_artifacts, values)
                    },
                }
        receptor_envelopes = []
        for receptor_index, label in enumerate(receptor_labels):
            values = np.asarray([
                item["receptors"][receptor_index][
                    "gradient_richardson_number"
                ]
                for item in member_segments
            ], dtype=float)
            receptor_envelopes.append({
                "label": label,
                "gradient_richardson_minimum": float(np.min(values)),
                "gradient_richardson_maximum": float(np.max(values)),
                "member_values": {
                    artifact["member_id"]: float(value)
                    for artifact, value in zip(member_artifacts, values)
                },
            })
        aggregate_segments.append({
            "stop_time_s": times[0],
            "regions": region_envelopes,
            "receptors": receptor_envelopes,
        })

    gates = {
        "all_frozen_regularization_members_executed": (
            len(member_records) == len(plan["members"])
        ),
        "all_member_diagnostic_gates_pass": all(
            record["all_gates_pass"] for record in member_records
        ),
        "member_model_times_match": True,
        "post_execution_member_selection_prohibited": (
            plan["execution_policy"]["selection_after_execution"]
            == "prohibited"
        ),
        "no_damping_or_flux_modification": all(
            not record["configuration"]["damping_applied"]
            for record in member_records
        ),
    }
    record = {
        "schema": SCHEMA,
        "status": "accepted_diagnostic_only_regularization_envelope",
        "scope": (
            "model-only Ri_g sensitivity to frozen regularisation; not a "
            "damping closure or physical-validation result"
        ),
        "provenance": {
            "run_directory": str(run_directory.resolve()),
            "mixing_policy_sha256": _sha256(mixing_policy_path.resolve()),
            "diagnostic_plan_sha256": _sha256(plan_path),
            "observation_value_fields_accessed": [],
        },
        "member_artifacts": member_artifacts,
        "segments": aggregate_segments,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
        "promotion_policy": {
            "closure_promotion_allowed": False,
            "member_selection_allowed": False,
        },
    }
    json_path = output / "ti021-stratification-envelope.json"
    figure_path = output / "ti021-stratification-envelope.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    times = np.asarray([item["stop_time_s"] for item in aggregate_segments])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for region_name in region_names:
        for axis, metric, title in (
            (axes[0, 0], "gradient_richardson_p50", "Regional Ri_g median"),
            (axes[0, 1], "gradient_richardson_p95", "Regional Ri_g p95"),
            (axes[1, 0], "fraction_at_or_above_0p25", "Fraction Ri_g >= 0.25"),
        ):
            low = np.asarray([
                item["regions"][region_name][metric]["minimum"]
                for item in aggregate_segments
            ])
            high = np.asarray([
                item["regions"][region_name][metric]["maximum"]
                for item in aggregate_segments
            ])
            axis.plot(times, 0.5 * (low + high), label=region_name)
            axis.fill_between(times, low, high, alpha=0.15)
            axis.set_title(title)
    receptor_low = np.asarray([
        np.quantile([
            receptor["gradient_richardson_minimum"]
            for receptor in item["receptors"]
        ], 0.50)
        for item in aggregate_segments
    ])
    receptor_high = np.asarray([
        np.quantile([
            receptor["gradient_richardson_maximum"]
            for receptor in item["receptors"]
        ], 0.95)
        for item in aggregate_segments
    ])
    axes[1, 1].fill_between(
        times, receptor_low, receptor_high, alpha=0.25,
        label="sensor-cell regularisation envelope",
    )
    axes[1, 1].set_title("Receptor-cell Ri_g envelope")
    for axis in axes.flat:
        axis.set_xlabel("Model time (s)")
        axis.set_ylabel("Ri_g or fraction")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    fig.suptitle("TI021 frozen gradient-Richardson regularisation envelope")
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    return record, {"json": json_path, "figure": figure_path}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-directory", type=Path, required=True)
    parser.add_argument("--mixing-policy", type=Path, required=True)
    parser.add_argument("--diagnostic-plan", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record, paths = audit_stratification_envelope(
        run_directory=args.run_directory,
        mixing_policy_path=args.mixing_policy,
        diagnostic_plan_path=args.diagnostic_plan,
        output_directory=args.output_directory,
    )
    print(json.dumps({
        "status": record["status"],
        "all_gates_pass": record["all_gates_pass"],
        "outputs": {name: str(path) for name, path in paths.items()},
    }, indent=2))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
