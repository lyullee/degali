"""Audit vertical arrival ordering in TI021 sensor groups.

This is a data-structure diagnostic, not a fit.  A single hydrostatic layer
or a monotone vertical concentration profile cannot represent groups whose
arrival order reverses with height.  Such groups are reported explicitly so
that a calibration cannot hide a missing vertical/recirculation state.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from degali.validation.thorney_ti021 import read_thorney_ti021_workbook


def _inversion_count(values: list[float]) -> int:
    return sum(
        1 for i, left in enumerate(values)
        for right in values[i + 1:]
        if left > right
    )


def audit(workbook: Path) -> dict[str, object]:
    dataset = read_thorney_ti021_workbook(workbook)
    groups: dict[tuple[float, float], list[dict[str, float]]] = defaultdict(list)
    for index, sensor in enumerate(dataset.dose_sensors):
        if sensor.arrival_time_s is None or sensor.dose_percent_s <= 0.0:
            continue
        groups[(float(sensor.local_x_m), float(sensor.local_y_m))].append({
            "workbook_row_index_zero_based": float(index),
            "height_m": float(sensor.height_m),
            "arrival_s": float(sensor.arrival_time_s),
        })

    records: list[dict[str, object]] = []
    for (local_x, local_y), entries in sorted(groups.items()):
        entries.sort(key=lambda item: item["height_m"])
        heights = [float(item["height_m"]) for item in entries]
        arrivals = [float(item["arrival_s"]) for item in entries]
        nondecreasing = all(left <= right for left, right in zip(arrivals, arrivals[1:]))
        nonincreasing = all(left >= right for left, right in zip(arrivals, arrivals[1:]))
        records.append({
            "local_x_m": local_x,
            "local_y_m": local_y,
            "sensor_count": len(entries),
            "height_m": heights,
            "arrival_s": arrivals,
            "arrival_range_s": float(max(arrivals) - min(arrivals)),
            "nondecreasing_with_height": nondecreasing,
            "nonincreasing_with_height": nonincreasing,
            "arrival_inversion_count": _inversion_count(arrivals),
            "entries": entries,
        })

    multi = [record for record in records if record["sensor_count"] >= 3]
    contradictory = [record for record in multi if not (
        record["nondecreasing_with_height"] or record["nonincreasing_with_height"]
    )]
    ranges = np.asarray([record["arrival_range_s"] for record in multi], dtype=float)
    return {
        "schema": "degali.ti021-vertical-ordering-audit.v1",
        "status": "diagnostic_only",
        "workbook": str(workbook.resolve()),
        "positive_dose_arrival_group_count": len(records),
        "groups_with_at_least_three_heights": len(multi),
        "non_monotone_group_count": len(contradictory),
        "non_monotone_group_fraction": (
            None if not multi else float(len(contradictory) / len(multi))
        ),
        "vertical_arrival_range_summary_s": {
            "median": None if not ranges.size else float(np.median(ranges)),
            "maximum": None if not ranges.size else float(np.max(ranges)),
        },
        "groups": records,
        "interpretation": (
            "Non-monotone groups cannot be represented by a single uniform layer "
            "or a monotone vertical concentration profile. They require a resolved "
            "vertical/recirculation state or an explicitly withheld observation "
            "operator; no vertical coefficient was fitted from these labels."
        ),
        "diagnostic_only": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.workbook)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        key: result[key] for key in (
            "groups_with_at_least_three_heights",
            "non_monotone_group_count",
            "non_monotone_group_fraction",
            "vertical_arrival_range_summary_s",
        )
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
