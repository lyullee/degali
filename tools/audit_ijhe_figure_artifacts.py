"""Audit IJHE graphical and derived SVG artifacts without opening raw data."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "degali.ijhe-figure-artifact-audit.v1"

FIGURES = (
    ("IJHE_GRAPHICAL_ABSTRACT.svg", True),
    ("outputs/ijhe-operator-comparison-2026-10-08.svg", False),
    ("outputs/transient-3d-error-study-2026-10-08/transient-3d-convergence.svg", False),
    ("tmp/ijhe-lh2-paper-figures/preslhy-concentration-residual-map.svg", False),
    ("tmp/ijhe-lh2-paper-figures/ffi-arc-residual-map.svg", False),
    ("tmp/ijhe-time-aligned/e35_lfl_uncertainty_envelope.svg", False),
    ("tmp/ijhe-decision-scenarios/decision_lfl_distances.svg", False),
)
GRAPHICAL_ABSTRACT_MIN_WIDTH = 1328
GRAPHICAL_ABSTRACT_MIN_HEIGHT = 531
GRAPHICAL_ABSTRACT_PREFERRED_ASPECT = 2.5


def _dimension_px(value: str) -> float:
    match = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]+)?)(?:px)?\s*", value)
    if not match:
        raise ValueError(f"unsupported SVG dimension: {value!r}")
    return float(match.group(1))


def _record(relative: str, required: bool) -> dict[str, object]:
    path = ROOT / relative
    record: dict[str, object] = {"path": relative, "required": required, "present": path.is_file()}
    if not path.is_file():
        record["status"] = "missing_required" if required else "not_regenerated_on_this_host"
        return record
    try:
        root = ET.parse(path).getroot()
        width = root.attrib.get("width")
        height = root.attrib.get("height")
        view_box = root.attrib.get("viewBox")
        if not width or not height or not view_box:
            raise ValueError("width, height and viewBox are required")
    except (OSError, ET.ParseError, ValueError) as error:
        record["status"] = "invalid_xml"
        record["error"] = str(error)
        return record
    record.update(
        {
            "status": "valid",
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "width": width,
            "height": height,
            "viewBox": view_box,
        }
    )
    if relative == "IJHE_GRAPHICAL_ABSTRACT.svg":
        try:
            pixel_width = _dimension_px(width)
            pixel_height = _dimension_px(height)
        except ValueError as error:
            record["status"] = "invalid_dimensions"
            record["error"] = str(error)
            return record
        record.update(
            {
                "pixel_width": pixel_width,
                "pixel_height": pixel_height,
                "aspect_ratio": round(pixel_width / pixel_height, 4),
                "preferred_aspect_ratio": GRAPHICAL_ABSTRACT_PREFERRED_ASPECT,
                "minimum_pixels": [GRAPHICAL_ABSTRACT_MIN_WIDTH, GRAPHICAL_ABSTRACT_MIN_HEIGHT],
            }
        )
        if pixel_width < GRAPHICAL_ABSTRACT_MIN_WIDTH or pixel_height < GRAPHICAL_ABSTRACT_MIN_HEIGHT:
            record["status"] = "invalid_dimensions"
            record["error"] = (
                f"graphical abstract is {pixel_width:g}x{pixel_height:g}px; "
                f"minimum is {GRAPHICAL_ABSTRACT_MIN_WIDTH}x{GRAPHICAL_ABSTRACT_MIN_HEIGHT}px"
            )
    return record


def audit() -> dict[str, object]:
    records = [_record(relative, required) for relative, required in FIGURES]
    required_failures = [item for item in records if item["required"] and item.get("status") != "valid"]
    invalid_optional = [item for item in records if not item["required"] and item.get("status") == "invalid_xml"]
    if required_failures or invalid_optional:
        status = "fail"
    elif any(item.get("status") != "valid" for item in records):
        status = "warning"
    else:
        status = "pass"
    return {"schema": SCHEMA, "status": status, "figures": records}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args(argv)
    payload = audit()
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"schema": payload["schema"], "status": payload["status"], "figures": len(payload["figures"])}))
    return 0 if payload["status"] != "fail" else 1


if __name__ == "__main__":
    raise SystemExit(main())
