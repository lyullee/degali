"""Re-run a small, physically declared TI021 600 s physics sweep.

The sweep starts from a previously recorded candidate configuration, changes
only global obstacle/wake residence parameters, and reports dose MAPE plus
conservation gates. It never reads observed dose values to create a sensor
specific correction. Every member is a complete 600 s run.
"""

from __future__ import annotations

import argparse
import contextlib
import inspect
import io
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.score_ti021_dense_layer_hybrid import run


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base",
        type=Path,
        default=Path("outputs/ti021-core-meander-s005-l150-w50-600s-2026-10-18.json"),
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("outputs/ti021-current-physics-sweep"),
    )
    args = parser.parse_args()
    base_record = json.loads(args.base.read_text(encoding="utf-8"))
    accepted = set(inspect.signature(run).parameters)
    base = {key: value for key, value in base_record.items() if key in accepted}
    args.output_directory.mkdir(parents=True, exist_ok=True)

    trials = {
        "reattach040-wake025": {
            "wall_impulse_reattached_fraction": 0.40,
            "wall_impulse_wake_capture_rate_s_inv": 0.025,
        },
        "reattach050-wake020": {
            "wall_impulse_reattached_fraction": 0.50,
            "wall_impulse_wake_capture_rate_s_inv": 0.020,
        },
        "reattach050-wake025": {
            "wall_impulse_reattached_fraction": 0.50,
            "wall_impulse_wake_capture_rate_s_inv": 0.025,
        },
        "reattach040-wake025-elevated020": {
            "wall_impulse_reattached_fraction": 0.40,
            "wall_impulse_wake_capture_rate_s_inv": 0.025,
            "wall_impulse_elevated_transport_fraction": 0.20,
            "wall_impulse_elevated_advection_duration_s": 60.0,
        },
    }
    summary = []
    for name, changes in trials.items():
        parameters = dict(base)
        parameters.update(changes)
        parameters["output"] = args.output_directory / f"{name}.json"
        # The solver prints its complete report; keep the sweep output compact.
        with contextlib.redirect_stdout(io.StringIO()):
            result = run(**parameters)
        observed = np.asarray(
            [row["observed_dose_vol_percent_s"] for row in result["sensor_rows"]],
            dtype=float,
        )
        predicted = np.asarray(
            [row["predicted_dose_vol_percent_s"] for row in result["sensor_rows"]],
            dtype=float,
        )
        positive = observed > 0.0
        summary.append(
            {
                "name": name,
                "dose_mape_percent": float(
                    np.mean(np.abs(predicted[positive] / observed[positive] - 1.0))
                    * 100.0
                ),
                "runtime_s": float(result["runtime_including_output_s"]),
                "maximum_relative_species_mass_residual": float(
                    result["diagnostics"]["maximum_relative_species_mass_residual"]
                ),
                "momentum_conservation_gate": bool(
                    result["diagnostics"]["momentum_conservation_gate"]
                ),
                "configuration": changes,
            }
        )
    report = {
        "schema": "degali.ti021-current-physics-sweep.v1",
        "base_record": str(args.base.resolve()),
        "full_event_horizon_s": 600.0,
        "sensor_specific_correction": False,
        "results": summary,
    }
    report_path = args.output_directory / "sweep-summary.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
