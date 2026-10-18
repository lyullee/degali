"""Audit segmented TI021 restart output against a matched one-shot run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

try:
    from tools.run_thorney_ti021_inner_transient import _sha256
except ModuleNotFoundError:
    from run_thorney_ti021_inner_transient import _sha256


SCHEMA = "degali.thorney-ti021-checkpoint-equivalence.v1"

RECEPTOR_FIELDS = (
    "time_s", "concentration_kg_m3", "concentration_vol_percent",
    "temperature_k",
)
FINAL_FIELDS = (
    "concentration_kg_m3", "concentration_vol_percent",
    "maximum_concentration_vol_percent", "temperature_k",
    "mixture_density_kg_m3", "velocity_u_m_s", "velocity_v_m_s",
    "velocity_w_m_s", "pressure_pa_gauge",
)
HANDOFF_FIELDS = (
    "time_s", "total_mass_flux_anomaly_kg_m2_s",
    "species_mass_flux_kg_m2_s", "sensible_enthalpy_flux_w_m2",
    "lateral_momentum_flux_anomaly_pa",
    "vertical_momentum_flux_anomaly_pa",
)


def _npz_differences(left: Path, right: Path, fields: tuple[str, ...]):
    result = {}
    with np.load(left, allow_pickle=False) as a, np.load(
        right, allow_pickle=False,
    ) as b:
        for field in fields:
            if a[field].shape != b[field].shape:
                raise ValueError(f"shape mismatch for {field}")
            delta = np.asarray(a[field], dtype=float) - np.asarray(
                b[field], dtype=float,
            )
            result[field] = {
                "maximum_absolute_difference": float(np.nanmax(np.abs(delta))),
                "root_mean_square_difference": float(np.sqrt(np.nanmean(delta**2))),
            }
    return result


def _hdf5_differences(left: Path, right: Path):
    result = {}
    with h5py.File(left, "r") as a, h5py.File(right, "r") as b:
        for field in HANDOFF_FIELDS:
            if a[field].shape != b[field].shape:
                raise ValueError(f"handoff shape mismatch for {field}")
            delta = np.asarray(a[field][:], dtype=float) - np.asarray(
                b[field][:], dtype=float,
            )
            result[field] = {
                "maximum_absolute_difference": float(np.max(np.abs(delta))),
                "root_mean_square_difference": float(np.sqrt(np.mean(delta**2))),
            }
    return result


def audit_checkpoint_equivalence(
    *, checkpointed_run: Path, one_shot_run: Path, output_json: Path,
) -> dict[str, object]:
    checkpointed = checkpointed_run.resolve()
    one_shot = one_shot_run.resolve()
    receptor = _npz_differences(
        checkpointed / "ti021-inner-receptor-traces.npz",
        one_shot / "ti021-inner-receptor-traces.npz",
        RECEPTOR_FIELDS,
    )
    final = _npz_differences(
        checkpointed / "ti021-inner-final-field.npz",
        one_shot / "ti021-inner-final-field.npz",
        FINAL_FIELDS,
    )
    handoff = _hdf5_differences(
        checkpointed / "ti021-inner-handoff.h5",
        one_shot / "ti021-inner-handoff.h5",
    )
    gates = {
        "time_axes_identical": (
            receptor["time_s"]["maximum_absolute_difference"] == 0.0
            and handoff["time_s"]["maximum_absolute_difference"] == 0.0
        ),
        "receptor_vol_percent_linf_at_most_1e_10": (
            receptor["concentration_vol_percent"][
                "maximum_absolute_difference"
            ] <= 1.0e-10
        ),
        "field_vol_percent_linf_at_most_1e_10": (
            final["concentration_vol_percent"][
                "maximum_absolute_difference"
            ] <= 1.0e-10
        ),
        "velocity_component_linf_at_most_1e_10_m_s": max(
            final[field]["maximum_absolute_difference"]
            for field in (
                "velocity_u_m_s", "velocity_v_m_s", "velocity_w_m_s",
            )
        ) <= 1.0e-10,
        "handoff_field_linf_at_most_1e_10": max(
            values["maximum_absolute_difference"]
            for field, values in handoff.items() if field != "time_s"
        ) <= 1.0e-10,
    }
    gates = {key: bool(value) for key, value in gates.items()}
    record = {
        "schema": SCHEMA,
        "status": (
            "accepted_checkpoint_equivalent_to_one_shot"
            if all(gates.values()) else "withheld_checkpoint_not_equivalent"
        ),
        "scope": (
            "matched configuration segmented-restart versus one-shot numerical "
            "equivalence; this is a restart integrity test, not physical validation"
        ),
        "provenance": {
            "checkpointed_record": {
                "path": str(checkpointed / "ti021-inner-transient.json"),
                "sha256": _sha256(checkpointed / "ti021-inner-transient.json"),
            },
            "one_shot_record": {
                "path": str(one_shot / "ti021-inner-transient.json"),
                "sha256": _sha256(one_shot / "ti021-inner-transient.json"),
            },
        },
        "receptor_differences": receptor,
        "final_field_differences": final,
        "handoff_differences": handoff,
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }
    destination = output_json.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpointed-run", type=Path, required=True)
    parser.add_argument("--one-shot-run", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    record = audit_checkpoint_equivalence(
        checkpointed_run=args.checkpointed_run,
        one_shot_run=args.one_shot_run,
        output_json=args.output_json,
    )
    print(json.dumps(record, indent=2, ensure_ascii=False))
    return 0 if record["all_gates_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
