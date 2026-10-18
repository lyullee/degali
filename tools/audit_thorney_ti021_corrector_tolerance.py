"""Audit TI021 nonlinear-corrector tolerance at the 2.5 m grid."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.thorney-ti021-corrector-tolerance-audit.v1"

FIELD_NAMES = (
    "concentration_vol_percent",
    "velocity_u_m_s",
    "velocity_v_m_s",
    "velocity_w_m_s",
    "pressure_pa_gauge",
)
HANDOFF_NAMES = (
    "total_mass_flux_anomaly_kg_m2_s",
    "species_mass_flux_kg_m2_s",
    "sensible_enthalpy_flux_w_m2",
    "lateral_momentum_flux_anomaly_pa",
    "vertical_momentum_flux_anomaly_pa",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _difference(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float]:
    left = np.asarray(reference, dtype=float)
    right = np.asarray(candidate, dtype=float)
    if left.shape != right.shape:
        raise ValueError("comparison arrays must have identical shapes")
    valid = np.isfinite(left) & np.isfinite(right)
    if not np.any(valid):
        raise ValueError("comparison arrays have no jointly finite values")
    absolute = np.abs(right[valid] - left[valid])
    return {
        "maximum_absolute_difference": float(np.max(absolute)),
        "root_mean_square_difference": float(np.sqrt(np.mean(absolute**2))),
        "relative_l1_difference": float(
            np.sum(absolute) / max(np.sum(np.abs(left[valid])), 1.0e-30)
        ),
    }


def build_audit(
    reference_directory: Path,
    candidate_directory: Path,
    *,
    rejected_reference_directory: Path | None = None,
    rejected_candidate_directory: Path | None = None,
) -> dict[str, object]:
    reference = reference_directory.resolve()
    candidate = candidate_directory.resolve()
    reference_record_path = reference / "ti021-inner-transient.json"
    candidate_record_path = candidate / "ti021-inner-transient.json"
    reference_record = json.loads(reference_record_path.read_text(encoding="utf-8"))
    candidate_record = json.loads(candidate_record_path.read_text(encoding="utf-8"))
    for record in (reference_record, candidate_record):
        if record.get("status") != "accepted_inner_transient_qualification":
            raise ValueError("tolerance audit inputs must be accepted qualifications")
        config = record["configuration"]
        if (
            config["cell_size_m"] != 2.5
            or config["duration_s"] != 10.0
            or config["time_step_s"] != 0.0625
            or config["observation_values_used_to_configure_run"] != []
        ):
            raise ValueError("tolerance audit requires matched no-fit 2.5 m/10 s runs")
    reference_tolerance = reference_record["configuration"]["corrector_tolerance"]
    candidate_tolerance = candidate_record["configuration"]["corrector_tolerance"]
    if reference_tolerance != 1.0e-6 or candidate_tolerance != 1.0e-4:
        raise ValueError("expected 1e-6 reference and 1e-4 candidate")

    with np.load(reference / "ti021-inner-final-field.npz") as left, np.load(
        candidate / "ti021-inner-final-field.npz"
    ) as right:
        field_differences = {
            name: _difference(left[name], right[name]) for name in FIELD_NAMES
        }
    with np.load(reference / "ti021-inner-receptor-traces.npz") as left, np.load(
        candidate / "ti021-inner-receptor-traces.npz"
    ) as right:
        if not np.array_equal(left["time_s"], right["time_s"]):
            raise ValueError("tolerance runs do not share the same event clock")
        receptor_difference = _difference(
            left["concentration_vol_percent"],
            right["concentration_vol_percent"],
        )
    with h5py.File(reference / "ti021-inner-handoff.h5", "r") as left, h5py.File(
        candidate / "ti021-inner-handoff.h5", "r"
    ) as right:
        if not np.array_equal(left["time_s"][:], right["time_s"][:]):
            raise ValueError("handoff histories do not share the same event clock")
        handoff_differences = {
            name: _difference(left[name][:], right[name][:])
            for name in HANDOFF_NAMES
        }

    runtime_reference = reference_record["runtime"]["wall_clock_s"]
    runtime_candidate = candidate_record["runtime"]["wall_clock_s"]
    rejected = None
    if (
        rejected_reference_directory is not None
        and rejected_candidate_directory is not None
    ):
        rejected_reference = rejected_reference_directory.resolve()
        rejected_candidate = rejected_candidate_directory.resolve()
        with np.load(
            rejected_reference / "ti021-inner-final-field.npz"
        ) as left, np.load(
            rejected_candidate / "ti021-inner-final-field.npz"
        ) as right:
            rejected_concentration = _difference(
                left["concentration_vol_percent"],
                right["concentration_vol_percent"],
            )
            rejected_velocity = max(
                _difference(left[name], right[name])[
                    "maximum_absolute_difference"
                ]
                for name in (
                    "velocity_u_m_s", "velocity_v_m_s", "velocity_w_m_s",
                )
            )
        rejected_record_path = rejected_candidate / "ti021-inner-transient.json"
        rejected_record = json.loads(
            rejected_record_path.read_text(encoding="utf-8")
        )
        rejected = {
            "tolerance": rejected_record["configuration"]["corrector_tolerance"],
            "horizon_s": rejected_record["configuration"]["duration_s"],
            "concentration": rejected_concentration,
            "maximum_velocity_component_difference_m_s": rejected_velocity,
            "maximum_concentration_vol_percent": rejected_record["numerics"][
                "maximum_vol_percent"
            ],
            "decision": "rejected_excess_initial_maximum_principle_error",
        }
    gates = {
        "candidate_concentration_relative_l1_at_most_5e_4": (
            field_differences["concentration_vol_percent"][
                "relative_l1_difference"
            ] <= 5.0e-4
        ),
        "candidate_concentration_linf_at_most_0p05vol_percent": (
            field_differences["concentration_vol_percent"][
                "maximum_absolute_difference"
            ] <= 0.05
        ),
        "candidate_receptor_linf_at_most_0p01vol_percent": (
            receptor_difference["maximum_absolute_difference"] <= 0.01
        ),
        "candidate_velocity_component_linf_at_most_0p002m_s": (
            max(
                field_differences[name]["maximum_absolute_difference"]
                for name in (
                    "velocity_u_m_s", "velocity_v_m_s", "velocity_w_m_s",
                )
            ) <= 0.002
        ),
        "candidate_handoff_mass_anomaly_relative_l1_at_most_5e_5": (
            handoff_differences[
                "total_mass_flux_anomaly_kg_m2_s"
            ]["relative_l1_difference"] <= 5.0e-5
        ),
        "candidate_continuity_gate_retained": (
            candidate_record["numerics"][
                "maximum_continuity_residual_kg_m3_s"
            ] <= 1.0e-7
        ),
        "candidate_mass_and_energy_gates_retained": (
            candidate_record["gates"]["species_mass_conservative"]
            and candidate_record["gates"]["sensible_enthalpy_conservative"]
        ),
        "candidate_corrector_converged_every_step": (
            candidate_record["numerics"]["corrector_nonconverged_steps"] == 0
        ),
        "candidate_at_least_1p5x_faster": (
            runtime_reference / runtime_candidate >= 1.5
        ),
        "observation_values_not_used": (
            reference_record["configuration"][
                "observation_values_used_to_configure_run"
            ] == []
            and candidate_record["configuration"][
                "observation_values_used_to_configure_run"
            ] == []
        ),
    }
    accepted = all(gates.values())
    return {
        "schema": SCHEMA,
        "status": (
            "accepted_1e4_corrector_tolerance_for_dx2p5_production_with_postcheck"
            if accepted else "withheld_corrector_tolerance_relaxation"
        ),
        "scope": (
            "observation-blind 10 s numerical equivalence audit; full 600 s "
            "conservation and solution-difference postcheck remains mandatory"
        ),
        "configuration": {
            "cell_size_m": 2.5,
            "time_step_s": 0.0625,
            "comparison_horizon_s": 10.0,
            "reference_tolerance": reference_tolerance,
            "candidate_tolerance": candidate_tolerance,
            "observation_value_fields_accessed": [],
        },
        "provenance": {
            "reference_record": {
                "path": str(reference_record_path),
                "sha256": _sha256(reference_record_path),
            },
            "candidate_record": {
                "path": str(candidate_record_path),
                "sha256": _sha256(candidate_record_path),
            },
        },
        "runtime": {
            "reference_wall_clock_s": runtime_reference,
            "candidate_wall_clock_s": runtime_candidate,
            "speedup": runtime_reference / runtime_candidate,
            "reference_corrector_mean_iterations_per_step": (
                reference_record["numerics"][
                    "corrector_mean_iterations_per_step"
                ]
            ),
            "candidate_corrector_mean_iterations_per_step": (
                candidate_record["numerics"][
                    "corrector_mean_iterations_per_step"
                ]
            ),
            "candidate_linear_600s_estimate_hours": candidate_record[
                "runtime"
            ]["linear_extrapolation_to_600s_hours"],
        },
        "final_field_differences": field_differences,
        "receptor_trace_difference": receptor_difference,
        "handoff_history_differences": handoff_differences,
        "rejected_1e3_diagnostic": rejected,
        "gates": gates,
        "all_gates_pass": accepted,
        "decision": {
            "production_tolerance": 1.0e-4 if accepted else 1.0e-6,
            "mandatory_postcheck": (
                "repeat or checkpoint a strict 1e-6 reference window near the "
                "plume peak and compare final conservation before promotion"
            ),
        },
    }


def _plot(record: dict[str, object], output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    metrics = record["final_field_differences"]
    labels = ["C vol%", "u", "v", "w", "pressure"]
    keys = list(FIELD_NAMES)
    values = [metrics[key]["relative_l1_difference"] for key in keys]
    axes[0].bar(labels, values, color="tab:blue")
    axes[0].set_yscale("log")
    axes[0].axhline(5.0e-4, color="crimson", linestyle="--",
                    label="concentration gate")
    axes[0].set_ylabel("relative L1 difference")
    axes[0].set_title("1e-4 candidate vs 1e-6 reference")
    axes[0].legend(fontsize=8)
    runtime = record["runtime"]
    axes[1].bar(
        ["1e-6 reference", "1e-4 candidate"],
        [runtime["reference_wall_clock_s"], runtime["candidate_wall_clock_s"]],
        color=["0.35", "tab:green"],
    )
    axes[1].set_ylabel("10 s run wall clock (s)")
    axes[1].set_title(
        f"speedup {runtime['speedup']:.2f}x; mean corrector iterations "
        f"{runtime['reference_corrector_mean_iterations_per_step']:.2f} -> "
        f"{runtime['candidate_corrector_mean_iterations_per_step']:.2f}"
    )
    fig.suptitle("TI021 dx=2.5 m nonlinear-corrector tolerance audit")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-directory", type=Path, required=True)
    parser.add_argument("--candidate-directory", type=Path, required=True)
    parser.add_argument("--rejected-reference-directory", type=Path)
    parser.add_argument("--rejected-candidate-directory", type=Path)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    record = build_audit(
        args.reference_directory,
        args.candidate_directory,
        rejected_reference_directory=args.rejected_reference_directory,
        rejected_candidate_directory=args.rejected_candidate_directory,
    )
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "ti021-corrector-tolerance-audit.json"
    plot_path = output / "ti021-corrector-tolerance-audit.png"
    json_path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    _plot(record, plot_path)
    print(json.dumps({
        "status": record["status"],
        "all_gates_pass": record["all_gates_pass"],
        "runtime": record["runtime"],
        "concentration_difference": record["final_field_differences"][
            "concentration_vol_percent"
        ],
        "receptor_difference": record["receptor_trace_difference"],
        "json": str(json_path),
        "plot": str(plot_path),
    }, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
