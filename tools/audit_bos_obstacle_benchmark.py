"""Fail-closed audit for the published BOS obstacle benchmark descriptor."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


REQUIRED_SIGNATURES = {
    "flat_plate": "gaussian_like_radial_profile",
    "flat_plate_with_sidewalls": "W_radial_profile",
    "dome": "W_like_radial_profile",
}


def audit(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    physical = data.get("physical_scope", {})
    machine = data.get("machine_data", {})
    promotion = data.get("promotion", {})
    signatures = physical.get("reported_shape_signatures", {})
    errors: list[str] = []
    if data.get("schema") != "degali.public-obstacle-bos-benchmark.v1":
        errors.append("schema_mismatch")
    if data.get("citation", {}).get("doi") != "10.1016/j.est.2024.110764":
        errors.append("doi_mismatch")
    if physical.get("measurement") != "background-oriented schlieren with concentration inversion":
        errors.append("measurement_definition_missing")
    for key, expected in REQUIRED_SIGNATURES.items():
        if signatures.get(key) != expected:
            errors.append(f"shape_signature_missing:{key}")
    if physical.get("reported_inversion_error_percent") != 3.81:
        errors.append("reported_error_missing")
    for key in ("raw_concentration_coordinates_available", "obstacle_cad_available", "common_clock_available", "sensor_calibration_available"):
        if machine.get(key) is not False:
            errors.append(f"machine_data_gate_not_explicit:{key}")
    if promotion.get("promotion_allowed") is not False:
        errors.append("promotion_not_fail_closed")
    return {
        "schema": "degali.public-obstacle-bos-audit.v1",
        "benchmark_id": data.get("benchmark_id"),
        "status": "withheld_literature_aggregate_only" if not errors else "invalid_descriptor",
        "promotion_allowed": False,
        "errors": errors,
        "reported_inversion_error_percent": physical.get("reported_inversion_error_percent"),
        "shape_signatures": signatures,
        "machine_data": machine,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("descriptor", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.descriptor)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {args.output} status={result['status']} promotion_allowed={result['promotion_allowed']}")
    return 0 if result["status"] == "withheld_literature_aggregate_only" else 1


if __name__ == "__main__":
    raise SystemExit(main())
