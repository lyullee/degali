"""Reconcile existing PRESLHY source flows with D4.8 Appendix 1."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "reference/preslhy/measured_pipe_source_2026-09-05.json"
REPORT = ROOT / "reference/lh2/preslhy-d4.8-condensed-phases-2020.pdf"
PROTOCOL = ROOT / "docs/prereg-d48-source-boundary-reconciliation.md"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_screen(value: float, reference: float, tolerance: float = 0.10) -> dict:
    relative = (value - reference) / reference
    return {
        "status": "pass" if abs(relative) <= tolerance else "fail",
        "reference_g_s": reference,
        "relative_difference": relative,
    }


def range_screen(value: float, lower: float, upper: float, tolerance: float = 0.10) -> dict:
    if lower <= value <= upper:
        status = "pass"
        nearest = value
    else:
        nearest = lower if value < lower else upper
        status = "near_range" if abs(value - nearest) / nearest <= tolerance else "fail"
    return {
        "status": status,
        "reference_range_g_s": [lower, upper],
        "relative_distance_from_nearest_bound": (value - nearest) / nearest,
    }


def audit(source: dict) -> dict:
    rows = {int(row["trial"]): row for row in source["trials"]}
    if not {10, 11, 12}.issubset(rows):
        raise ValueError("source table is missing a D4.8 comparison trial")
    values = {n: float(rows[n]["pressure_loss_mass_flow_g_s"]) for n in (10, 11, 12)}
    if not all(math.isfinite(v) and v > 0 for v in values.values()):
        raise ValueError("source flows must be finite and positive")
    screens = {
        "10": {"pressure_loss_g_s": values[10], **exact_screen(values[10], 298.0)},
        "11": {"pressure_loss_g_s": values[11], **exact_screen(values[11], 265.0)},
        "12": {"pressure_loss_g_s": values[12], **range_screen(values[12], 90.0, 100.0)},
    }
    stop = any(screens[str(n)]["status"] == "fail" for n in (10, 11))
    return {
        "screens": screens,
        "stop_before_downstream": stop,
        "decision": "retain_trial_specific_pressure_loss" if not stop else "stop",
        "trials_22_to_25_inferred_from_summary": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    result = audit(source)
    result.update({
        "completed": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": str(PROTOCOL.relative_to(ROOT)).replace("\\", "/"),
        "inputs_sha256": {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
            for path in (SOURCE, REPORT, PROTOCOL, Path(__file__))
        },
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

