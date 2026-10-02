"""Compare a declared E3.4 mass-loss window with the conduction-only pool source.

The manifest is deliberately external to the repository.  It records a
predeclared observation window and independently selected substrate values;
the tool neither searches for the most favourable interval nor fits a heat
transfer coefficient.  E3.4 workbooks and the resulting JSON are user-held
evidence, never package data.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from degali.addons.pool_evaporation import (
    SolidSubstrate,
    substrate_conduction_evaporation,
)
from degali.validation.preslhy_e34 import evaporation_window, read_e34_pool_history


def _cumulative_at(result, time_s: float) -> float:
    time = np.asarray([step.elapsed_s for step in result.steps], dtype=float)
    mass = np.asarray([step.evaporated_mass_kg for step in result.steps], dtype=float)
    return float(np.interp(time_s, np.r_[0.0, time], np.r_[0.0, mass]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="external JSON declaration")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite source-validation evidence")
    declaration = json.loads(args.manifest.read_text(encoding="utf-8"))
    workbook = Path(declaration["workbook"]).expanduser().resolve()
    if not workbook.is_file():
        raise FileNotFoundError(workbook)
    trial_id = str(declaration["trial_id"]).strip()
    if not trial_id or not workbook.stem.startswith(trial_id):
        raise ValueError(
            "trial_id must exactly identify the supplied E3.4 workbook prefix; "
            "do not map a report example to a similarly named trial"
        )
    report_trial = declaration.get("report_trial_id")
    if report_trial is not None:
        report_trial = str(report_trial).strip()
        if report_trial != trial_id:
            raise ValueError(
                "report_trial_id differs from workbook trial_id; "
                "report-to-workbook numerical comparison is not established"
            )
    contact = float(declaration["pool_contact_start_s"])
    start = float(declaration["window_start_s"])
    end = float(declaration["window_end_s"])
    if not (contact <= start < end):
        raise ValueError("require pool_contact_start_s <= window_start_s < window_end_s")
    properties = declaration["substrate"]
    substrate = SolidSubstrate(
        conductivity_w_m_k=float(properties["conductivity_w_m_k"]),
        density_kg_m3=float(properties["density_kg_m3"]),
        heat_capacity_j_kg_k=float(properties["heat_capacity_j_kg_k"]),
        initial_temperature_k=float(properties["initial_temperature_k"]),
        depth_m=float(properties["depth_m"]),
        cells=int(properties.get("cells", 96)),
    )
    history = read_e34_pool_history(workbook)
    observed = evaporation_window(
        history, start_s=start, end_s=end,
        thermocouples=tuple(declaration.get("thermocouples", ())),
    )
    model = substrate_conduction_evaporation(
        substrate,
        area_m2=float(declaration["area_m2"]), duration_s=end - contact,
        time_step_s=float(declaration.get("time_step_s", 1.0)),
        saturation_temperature_k=float(declaration.get("saturation_temperature_k", 20.27)),
        latent_heat_j_kg=float(declaration.get("latent_heat_j_kg", 4.46e5)),
        initial_liquid_mass_kg=declaration.get("initial_liquid_mass_kg"),
    )
    elapsed_start, elapsed_end = start - contact, end - contact
    predicted_loss = _cumulative_at(model, elapsed_end) - _cumulative_at(model, elapsed_start)
    predicted_rate = predicted_loss / (end - start)
    result = {
        "completed": True,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": "PRESLHY E3.4",
        "source_doi": "10.35097/1319",
        "source_files_distributed": False,
        "coefficient_fitted": False,
        "default_model_changed": False,
        "workbook": workbook.name,
        "workbook_trial_id": trial_id,
        "report_trial_id": report_trial,
        "report_to_workbook_identity_confirmed": report_trial == trial_id,
        "declared_window_s": [start, end],
        "pool_contact_start_s": contact,
        "observed": {
            "mass_loss_kg": observed.mass_lost_kg,
            "evaporation_rate_kg_s": observed.evaporation_rate_kg_s,
            "mass_slope_r_squared": observed.mass_slope_r_squared,
            "temperature_span_k": observed.temperature_span_k,
            "notes": observed.notes,
        },
        "conduction_model": {
            "mass_loss_kg": predicted_loss,
            "mean_evaporation_rate_kg_s": predicted_rate,
            "finite_depth_ratio": model.finite_depth_ratio,
            "effectively_semi_infinite": model.depth_is_effectively_semi_infinite,
            "declared_substrate": properties,
        },
        "interpretation": (
            "A mismatch is a source-boundary finding, not permission to fit an "
            "interfacial boiling or atmospheric-transfer coefficient."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
