"""Summarize checkpoint-prefix evolution in vol% and receptor space."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


SCHEMA = "degali.thorney-ti021-prefix-evolution-summary.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _prefix_record(prefix: Path) -> dict[str, object]:
    field_path = prefix / "ti021-inner-final-field.npz"
    trace_path = prefix / "ti021-inner-receptor-traces.npz"
    record_path = prefix / "ti021-inner-transient.json"
    with np.load(field_path, allow_pickle=False) as field:
        solid = np.asarray(field["solid_mask"], dtype=bool)
        peak_field = np.asarray(field["maximum_concentration_vol_percent"], dtype=float)
        final_field = np.asarray(field["concentration_vol_percent"], dtype=float)
    with np.load(trace_path, allow_pickle=False) as trace:
        time_s = np.asarray(trace["time_s"], dtype=float)
        labels = np.asarray(trace["labels"])
        receptors = np.asarray(trace["concentration_vol_percent"], dtype=float)
    record = json.loads(record_path.read_text(encoding="utf-8"))
    fluid = ~solid
    return {
        "path": str(prefix.resolve()),
        "record_sha256": _sha256(record_path),
        "stop_time_s": float(time_s[-1]),
        "frame_count": int(time_s.size),
        "field_maximum_ever_vol_percent": float(np.nanmax(peak_field[fluid])),
        "field_final_maximum_vol_percent": float(np.nanmax(final_field[fluid])),
        "field_final_p99_vol_percent": float(np.nanpercentile(final_field[fluid], 99.0)),
        "hazard_volume_final_m3": float(np.count_nonzero((final_field >= 4.0) & fluid) * 1.25**3),
        "hazard_volume_ever_m3": float(np.count_nonzero((peak_field >= 4.0) & fluid) * 1.25**3),
        "receptor_labels_sha256": hashlib.sha256(labels.tobytes()).hexdigest(),
        "receptor_global_peak_vol_percent": float(np.nanmax(receptors)),
        "receptor_endpoint_max_vol_percent": float(np.nanmax(receptors[:, -1])),
        "receptor_endpoint_mean_vol_percent": float(np.nanmean(receptors[:, -1])),
        "numerics": record.get("numerics", {}),
    }


def summarize(earlier: Path, later: Path) -> dict[str, object]:
    first = _prefix_record(earlier.resolve())
    second = _prefix_record(later.resolve())
    if second["stop_time_s"] <= first["stop_time_s"]:
        raise ValueError("later prefix must have a strictly larger stop time")
    return {
        "schema": SCHEMA,
        "status": "accepted_computed_window_evolution",
        "earlier": first,
        "later": second,
        "delta_later_minus_earlier": {
            key: second[key] - first[key]
            for key in (
                "stop_time_s", "field_maximum_ever_vol_percent",
                "field_final_maximum_vol_percent", "field_final_p99_vol_percent",
                "hazard_volume_final_m3", "hazard_volume_ever_m3",
                "receptor_global_peak_vol_percent",
                "receptor_endpoint_max_vol_percent",
                "receptor_endpoint_mean_vol_percent",
            )
        },
        "promotion_allowed": False,
        "interpretation": (
            "Computed-window evolution only; no full-event or independent-field "
            "validation claim."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--earlier-prefix", type=Path, required=True)
    parser.add_argument("--later-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.earlier_prefix, args.later_prefix)
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=False)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
