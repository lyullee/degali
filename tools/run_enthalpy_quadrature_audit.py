"""Check all five section fluxes at 2x and 4x the saved quadrature order."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from run_preslhy_ambient_profile_audit import ROOT, replay


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("field", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite an existing quadrature audit")
    field = json.loads(args.field.read_text(encoding="utf-8"))
    data = json.loads((ROOT / "reference/preslhy/e35_reduced.json").read_text(encoding="utf-8"))
    rows = []
    for trial in data["trials"]:
        n = trial["trial"]
        if n not in field["selected_trials"]:
            continue
        trajectory, _ = replay(field, trial)
        if trajectory is None:
            raise ValueError("quadrature audit requires complete stored trajectories")
        model = trajectory.model
        q = model.quadrature_points
        indices = np.unique(np.linspace(0, len(trajectory.states)-1, 9, dtype=int))
        for index in indices:
            state = trajectory.states[index]
            fluxes = [model._as_array(model.integral_fluxes(state, quadrature_points=points))
                      for points in (q, 2*q, 4*q)]
            buoyancy = [float(9.81 * state[2] * np.sum(
                (model.rhoa-model.profiles(state, quadrature_points=points)[1])
                * model._quadrature(points)[1])) for points in (q, 2*q, 4*q)]
            scales = np.maximum(np.abs(fluxes[-1]), [1e-12, 1e-12, 1., 1., 1.])
            rows.append({"trial": n, "x_m": float(state[5]), "orders": [q, 2*q, 4*q],
                         "fluxes": [f.tolist() for f in fluxes],
                         "buoyancy_force_N_per_m": buoyancy,
                         "buoyancy_scaled_1x_to_4x": abs(buoyancy[0]-buoyancy[2])/max(abs(buoyancy[2]), 1.),
                         "relative_1x_to_4x": (np.abs(fluxes[0]-fluxes[2])/scales).tolist(),
                         "relative_2x_to_4x": (np.abs(fluxes[1]-fluxes[2])/scales).tolist()})
        print(f"quadrature trial {n}: {len(indices)} sections checked", flush=True)
    maximum = max(max(row["relative_1x_to_4x"]) for row in rows)
    maximum_buoyancy = max(row["buoyancy_scaled_1x_to_4x"] for row in rows)
    payload = {"flux_order": ["mass", "H2", "momentum_x", "momentum_z", "total_energy"],
               "field_file": args.field.name, "field_sha256": hashlib.sha256(args.field.read_bytes()).hexdigest(),
               "threshold": 1e-5, "maximum_relative_change": maximum, "passed": maximum <= 1e-5,
               "maximum_buoyancy_scaled_change": maximum_buoyancy,
               "buoyancy_note": "Diagnostic force convergence; scale is max(abs(fine force), 1 N/m). Not a fitted physical coefficient.",
               "sections": rows}
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"maximum_relative_change": maximum, "passed": payload["passed"]}))


if __name__ == "__main__":
    main()
